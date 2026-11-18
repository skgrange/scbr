#' Function to squash R check's global variable notes. 
#' 
#' @name zzz
if (getRversion() >= "2.15.1") {
  
  # What variables are causing issues?
  variables <- c(
    "id", "inputs_perturbed", "carbon_emissions_fossil_fuels", 
    "carbon_emissions_land_use_change", "co2", "date_end", "o2", "observations",
    "remainder", "seasonal", "sub.labels", "trend", "trend_smooth", "trends",
    "value", "value_type", "variable", "year_min", "year", "weights", "n",
    "simulation_id"
  )
  
  # Squash the note
  utils::globalVariables(variables)
  
}
