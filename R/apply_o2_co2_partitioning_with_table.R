#' Function to apply the oxygen-carbon dioxide partitioning method using a table
#' of Monte Carlo simulation parameters. 
#' 
#' Wrapper function to execute the oxygen-carbon dioxide 
#' (\eqn{\text{O}_2}-\eqn{\text{CO}_2}) stoichiometric partitioning framework 
#' across an ensemble of perturbed input parameters and observational 
#' uncertainties. Each row of `df_sim` defines a Monte Carlo realization. 
#' 
#' \code{apply_o2_co2_partitioning_with_table} will conduct budgeting in 
#' parallel if \strong{mirai} daemons have been registered. 
#' 
#' @author Stuart K. Grange
#' 
#' @param df_sim A data frame or tibble defining the Monte Carlo simulation 
#' ensemble, where each row represents a single realization and contains 
#' columns for `simulation_id`, `co2_uncertainty`, `o2_uncertainty`, `span`, 
#' `carbon_emissions_fossil_fuels_uncertainty`, 
#' `carbon_land_use_change_emissions_uncertainty`, `o2_ocean_outgasing_year`, 
#' `conversion_factor_gigatonnes_carbon_to_ppm`, 
#' `conversion_factor_gigatonnes_oxygen_to_ppm`, `oxidation_ratio_fossil_fuels`, 
#' and `oxidation_ratio_biosphere`. 
#' 
#' @param df_obs A data frame or tibble containing monthly atmospheric 
#' observations (`id`, `date`, `co2`, `o2`). 
#' 
#' @param df_emissions A data frame or tibble containing annual global carbon 
#' emissions (`date`, `carbon_emissions_fossil_fuels`, 
#' `carbon_emissions_land_use_change`). 
#' 
#' @param progress Should a progress bar be displayed? 
#' 
#' @return Tibble. 
#' 
#' @seealso [apply_o2_co2_partitioning()]
#' 
#' @export
apply_o2_co2_partitioning_with_table <- function(df_sim, df_obs, df_emissions,
                                                 progress) {
  
  # Split mapping table into a list
  list_sim <- split(df_sim, seq_len(nrow(df_sim)))
  
  # Do a budget for each id within the observations with the variables passed 
  # in list_map, possibly in parallel, using map is easier than pmap here with
  # the parallel call
  list_sim |> 
    purrr::map(
      purrr::in_parallel(
        \(x) apply_o2_co2_partitioning_with_input_perturbation(
          df_obs = df_obs,
          df_emissions = df_emissions,
          simulation_id = x$simulation_id,
          co2_uncertainty = x$co2_uncertainty,
          o2_uncertainty = x$o2_uncertainty,
          span = x$span,
          carbon_emissions_fossil_fuels_uncertainty = 
            x$carbon_emissions_fossil_fuels_uncertainty,
          carbon_land_use_change_emissions_uncertainty = 
            x$carbon_land_use_change_emissions_uncertainty,
          o2_ocean_outgasing_year = x$o2_ocean_outgasing_year,
          conversion_factor_gigatonnes_carbon_to_ppm = 
            x$conversion_factor_gigatonnes_carbon_to_ppm,
          conversion_factor_gigatonnes_oxygen_to_ppm = 
            x$conversion_factor_gigatonnes_oxygen_to_ppm,
          oxidation_ratio_fossil_fuels = x$oxidation_ratio_fossil_fuels,
          oxidation_ratio_biosphere = x$oxidation_ratio_biosphere
        ),
        df_obs = df_obs,
        df_emissions = df_emissions,
        apply_o2_co2_partitioning_with_input_perturbation = 
          apply_o2_co2_partitioning_with_input_perturbation,
        apply_o2_co2_partitioning = apply_o2_co2_partitioning,
        apply_o2_co2_partitioning_worker = apply_o2_co2_partitioning_worker,
        apply_o2_co2_partitioning_paired_worker =
          apply_o2_co2_partitioning_paired_worker,
        decompose_time_series_and_prepare_for_budget = 
          decompose_time_series_and_prepare_for_budget
      ),
      .progress = progress
    ) |> 
    purrr::list_rbind() |> 
    relocate(inputs_perturbed,
             simulation_id)
  
}


apply_o2_co2_partitioning_with_input_perturbation <- function(
    df_obs, df_emissions, simulation_id, co2_uncertainty, o2_uncertainty, span,
    carbon_emissions_fossil_fuels_uncertainty, 
    carbon_land_use_change_emissions_uncertainty, o2_ocean_outgasing_year,
    conversion_factor_gigatonnes_carbon_to_ppm, 
    conversion_factor_gigatonnes_oxygen_to_ppm, oxidation_ratio_fossil_fuels,
    oxidation_ratio_biosphere) {
  
  # For parallel processing
  requireNamespace("dplyr")
  
  # Add uncertainty and sample the two carbon emissions variables, these
  # uncertainties are one sigma and will overwrite the values, different logic
  # used for the different variables
  df_emissions <- df_emissions |> 
    mutate(
      carbon_emissions_fossil_fuels = rnorm(
        n = dplyr::n(),
        mean = carbon_emissions_fossil_fuels,
        sd = carbon_emissions_fossil_fuels *
          carbon_emissions_fossil_fuels_uncertainty
      ),
      carbon_emissions_land_use_change = rnorm(
        n = dplyr::n(),
        mean = carbon_emissions_land_use_change,
        sd = carbon_land_use_change_emissions_uncertainty
      )
    )
  
  # Prepare atmospheric observations. Monthly means are decomposed, the trend
  # component is smoothed and aggregated to annual means, and also join 
  # emissions data to observations
  df_obs_annual <- decompose_time_series_and_prepare_for_budget(
    df_obs,
    df_emissions = df_emissions,
    co2_uncertainty = co2_uncertainty,
    o2_uncertainty = o2_uncertainty,
    span = span
  )
  
  # Do the partitioning
  df <- apply_o2_co2_partitioning(
    df_obs_annual,
    conversion_factor_gigatonnes_carbon_to_ppm = 
      conversion_factor_gigatonnes_carbon_to_ppm, 
    conversion_factor_gigatonnes_oxygen_to_ppm = 
      conversion_factor_gigatonnes_oxygen_to_ppm, 
    oxidation_ratio_fossil_fuels = oxidation_ratio_fossil_fuels,
    oxidation_ratio_biosphere = oxidation_ratio_biosphere,
    o2_ocean_outgasing_year = o2_ocean_outgasing_year
  ) |> 
    mutate(inputs_perturbed = TRUE,
           simulation_id = simulation_id)
  
  return(df)
  
}
