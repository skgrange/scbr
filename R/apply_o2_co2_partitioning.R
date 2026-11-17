#' Function to apply the oxygen-carbon dioxide (O2-CO2) partitioning method to
#' atmospheric observations to calculate a carbon budget where carbon emissions
#' are allocated to the atmospheric reservoir, the land sink, and the ocean sink
#' using stoichiometry. 
#' 
#' @author Stuart K. Grange
#' 
#' @param df A data frame or tibble containing atmospheric observations and 
#' emissions data, with columns for `id`, `date`, `co2`, `o2`, 
#' `carbon_emissions_fossil_fuels`, and `carbon_emissions_land_use_change`.
#' 
#' @param conversion_factor_gigatonnes_carbon_to_ppm Conversion factor for 
#' transforming gigatonnes of carbon (GtC) to atmospheric mole fraction 
#' (ppm).
#' 
#' @param conversion_factor_gigatonnes_oxygen_to_ppm Conversion factor for 
#' transforming gigatonnes of oxygen to atmospheric mole fraction (ppm).
#' 
#' @param oxidation_ratio_fossil_fuels The stoichiometric oxygen-to-carbon 
#' ratio for fossil fuels ((\eqn{\alpha_{\text{FF}}})).
#' 
#' @param oxidation_ratio_biosphere The stoichiometric oxygen-to-carbon ratio 
#' for terrestrial biosphere exchange ((\eqn{\alpha_{\text{B}}})).
#' 
#' @param o2_ocean_outgasing_year Annual rate of ocean oxygen outgassing.
#' 
#' @return A tibble containing the partitioned carbon budget variables.
#' 
#' @export
apply_o2_co2_partitioning <- function(
    df, conversion_factor_gigatonnes_carbon_to_ppm, 
    conversion_factor_gigatonnes_oxygen_to_ppm, oxidation_ratio_fossil_fuels,
    oxidation_ratio_biosphere, o2_ocean_outgasing_year) {
  
  # Check input tibble
  check_o2_co2_partitioning_input(df)
  
  # Ensure table is ordered correctly
  df <- arrange(df, id, date)
  
  # Apply function to each id
  df |> 
    dplyr::group_split(id) |> 
    purrr::map(
      \(x) apply_o2_co2_partitioning_worker(
        df = x,
        o2_ocean_outgasing_year = o2_ocean_outgasing_year,
        conversion_factor_gigatonnes_carbon_to_ppm = 
          conversion_factor_gigatonnes_carbon_to_ppm,
        conversion_factor_gigatonnes_oxygen_to_ppm = 
          conversion_factor_gigatonnes_oxygen_to_ppm,
        oxidation_ratio_fossil_fuels = oxidation_ratio_fossil_fuels,
        oxidation_ratio_biosphere = oxidation_ratio_biosphere
      )
    ) |> 
    purrr::list_rbind() |> 
    relocate(inputs_perturbed)
  
}


apply_o2_co2_partitioning_worker <- function(
    df, o2_ocean_outgasing_year, conversion_factor_gigatonnes_carbon_to_ppm, 
    conversion_factor_gigatonnes_oxygen_to_ppm, oxidation_ratio_fossil_fuels,
    oxidation_ratio_biosphere) {
  
  # Split into two-row tibbles for year by year partitioning
  list_df <- purrr::map(seq_len(nrow(df) - 1), \(x) df[x:(x + 1), ])
  
  # Apply o2-co2 partitioning logic to add new variables to the input table
  list_df |> 
    purrr::map(apply_o2_co2_partitioning_paired_worker,
               o2_ocean_outgasing_year,
               conversion_factor_for_oxygen = 32.00876 / 1e15,
               conversion_factor_gigatonnes_carbon_to_ppm,
               conversion_factor_gigatonnes_oxygen_to_ppm,
               oxidation_ratio_fossil_fuels,
               oxidation_ratio_biosphere) |> 
    purrr::list_rbind() |> 
    mutate(inputs_perturbed = FALSE)
  
}


apply_o2_co2_partitioning_paired_worker <- function(
    df, o2_ocean_outgasing_year, conversion_factor_for_oxygen, 
    conversion_factor_gigatonnes_carbon_to_ppm, 
    conversion_factor_gigatonnes_oxygen_to_ppm, oxidation_ratio_fossil_fuels,
    oxidation_ratio_biosphere) {
  
  # Check if two rows have passed
  stopifnot(nrow(df) == 2L)
  
  # Get first and last observations, useful for the many calculations to avoid 
  # many head and tail calls
  df_first <- df[1, ]
  df_last <- df[2, ]
  
  # Calculate expected values based on stoichiometric ratios
  co2_expected <- df_first$co2 + df_first$carbon_emissions_fossil_fuels /
    conversion_factor_gigatonnes_carbon_to_ppm
  
  o2_expected <- df_first$o2 + df_first$carbon_emissions_fossil_fuels /
    conversion_factor_gigatonnes_carbon_to_ppm * oxidation_ratio_fossil_fuels
  
  # Change unit system for oxygen
  o2_ocean_outgasing_year_ppm <- o2_ocean_outgasing_year * 
    conversion_factor_for_oxygen / conversion_factor_gigatonnes_oxygen_to_ppm
  
  # Subtract ocean outgassing from final o2 value
  o2_last_without_outgassing <- df_last$o2 - o2_ocean_outgasing_year_ppm
  
  # Land use change calculations
  # co2
  co2_last_land_use_change <- co2_expected + 
    df_first$carbon_emissions_land_use_change /
    conversion_factor_gigatonnes_carbon_to_ppm
  
  # o2
  o2_last_land_use_change <- o2_expected +
    df_first$carbon_emissions_land_use_change /
    conversion_factor_gigatonnes_carbon_to_ppm * oxidation_ratio_biosphere
  
  # Calculate the ocean uptake
  co2_after_land_partition <- (o2_last_land_use_change - o2_last_without_outgassing) /
    oxidation_ratio_biosphere + df_last$co2
  
  co2_ocean_uptake <- (co2_last_land_use_change - co2_after_land_partition) * 
    conversion_factor_gigatonnes_carbon_to_ppm
  
  # Land uptake
  co2_land_uptake <- (o2_last_land_use_change - o2_last_without_outgassing) *
    conversion_factor_gigatonnes_carbon_to_ppm / oxidation_ratio_biosphere
  
  # Sum the two sinks
  carbon_both_sinks <- co2_ocean_uptake + co2_land_uptake
  
  # How much carbon goes into the atmosphere? 
  carbon_atmosphere <- (df_last$co2 - df_first$co2) * 
    conversion_factor_gigatonnes_carbon_to_ppm
  
  # The anthropogenic emissions
  fossil_fuel_emissions <- df_first$carbon_emissions_fossil_fuels
  land_use_change_emissions <- df_first$carbon_emissions_land_use_change
  sum_of_emissions <- fossil_fuel_emissions + land_use_change_emissions
  
  # Calculate the budget's imbalance, should be very very small
  imbalance <- fossil_fuel_emissions + land_use_change_emissions - 
    carbon_atmosphere - co2_land_uptake - co2_ocean_uptake
  
  # Add variables to input table
  df <- df_first |> 
    rename(fossil_fuel_emissions = carbon_emissions_fossil_fuels,
           land_use_change_emissions = carbon_emissions_land_use_change) |> 
    mutate(co2_expected = !!co2_expected,
           o2_expected = !!o2_expected,
           sum_of_emissions,
           atmosphere = carbon_atmosphere,
           ocean_sink = co2_ocean_uptake,
           land_sink = co2_land_uptake,
           both_sinks = carbon_both_sinks,
           imbalance)
  
  # # Build tables for vector plots
  # if (include_vectors) {
  #   # Coordinate tables for the vector plot
  #   df_expected <- tibble::tribble(
  #     ~plotting_group,   ~co2,                  ~o2,
  #     "expected_values", df_first$co2, df_first$o2,
  #     "expected_values", co2_expected,          o2_expected
  #   )
  #   
  #   df_land_use_change <- tibble::tribble(
  #     ~plotting_group,    ~co2,                     ~o2,
  #     "land_use_change",  co2_expected,             o2_expected,
  #     "land_use_change",  co2_last_land_use_change, o2_last_land_use_change
  #   )
  #   
  #   df_land_sink <- tibble::tribble(
  #     ~plotting_group, ~co2,                      ~o2,
  #     "land_sink",      df_last$co2,     o2_last_without_outgassing,
  #     "land_sink",      co2_after_land_partition, o2_last_land_use_change
  #   )
  #   
  #   df_ocean_sink <- tibble::tribble(
  #     ~plotting_group, ~co2,                      ~o2,
  #     "ocean_sink",     co2_after_land_partition, o2_last_land_use_change,
  #     "ocean_sink",     co2_last_land_use_change, o2_last_land_use_change
  #   )
  #   
  #   df_ocean_o2_outgassing <- tibble::tribble(
  #     ~plotting_group,       ~co2,                 ~o2,
  #     "ocean_o2_outgassing", df_last$co2, o2_last_without_outgassing,
  #     "ocean_o2_outgassing", df_last$co2, df_last$o2
  #   )
  #   
  #   df_atmospheric_co2 <- tibble::tribble(
  #     ~plotting_group,   ~co2,                  ~o2,
  #     "atmospheric_co2", df_first$co2, df_last$o2,
  #     "atmospheric_co2", df_last$co2,  df_last$o2
  #   )
  #   
  #   df_atmospheric_o2 <- tibble::tribble(
  #     ~plotting_group,   ~co2,                  ~o2,
  #     "atmospheric_o2",  df_first$co2, df_first$o2,
  #     "atmospheric_o2",  df_first$co2, df_last$o2
  #   )
  #   
  #   # Bind
  #   df_for_plot <- dplyr::bind_rows(
  #     df_expected,
  #     df_land_use_change,
  #     df_ocean_o2_outgassing,
  #     df_land_sink,
  #     df_ocean_sink,
  #     df_atmospheric_co2,
  #     df_atmospheric_o2
  #   )
  #   
  # }
  
  return(df)
  
}


check_o2_co2_partitioning_input <- function(df) {
  
  # Define required columns for the partitioning pipeline
  variables <- c(
    "id", "date", "co2", "o2", "carbon_emissions_fossil_fuels", 
    "carbon_emissions_land_use_change"
  )
  
  # Check for missing variables
  missing_variables <- setdiff(variables, names(df))
  
  if (length(missing_variables) > 0) {
    missing_variables_text <- stringr::str_c("`", missing_variables, "`")
    cli::cli_abort("There are missing variables in `df`: {missing_variables_text}...")
  }
  
  return(invisible(NULL))
  
}
