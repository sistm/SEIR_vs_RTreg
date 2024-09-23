sim_SIR_Simulx_ME_2params <- function(tmax = 121, 
                                      path_to_ind_params, 
                                      regressor_df_path, 
                                      ME_I_a = 0.0408, ME_I_b = 0.1){
  
  
  # this function requires a Simulx project already loaded!
  defineIndividualElement(name = "mlx_Pop_ind", element = path_to_ind_params)
  defineRegressorElement(name = "mlx_Reg", element = regressor_df_path)
  
  
  defineOutputElement(name = "I_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "I"))
  defineOutputElement(name = "S_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "S"))
  defineOutputElement(name = "R_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "R"))
  defineOutputElement(name = "transmission", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "transmission"))
  defineOutputElement(name = "predI", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predI"))
   
  
  setGroupElement(group = "simulationGroup1", 
                  elements = c("mlx_Pop_ind"))
  
  setGroupElement(group = "simulationGroup1", 
                  elements = c("S_mlx", "R_mlx", "I_mlx", "transmission", "predI"))
  
  setGroupRemaining(group = "simulationGroup1",
                    remaining = list(a = 0, b = 0))
  
  runSimulation()
  sim <- getSimulationResults()
  
  
  # bring results into required format
  res_df <- data.frame()
  for(k in 1:length(sim$res)){
    if(k == 1){
      res_df <- sim$res[[k]]
    }else{
      res_df <- left_join(res_df, sim$res[[k]], by = c("id", "time"))
    }
  }
  
  res_df <- res_df  %>%
    rename(IncI = predI) %>%
    rowwise() %>%
    mutate(IncI_ME = max(rnorm(1, IncI, (ME_I_a + ME_I_b*IncI)/1.96), 0)) %>%
    ungroup()
  
}


SIR_monolix_data_creation_ME <- function(simulation_results, 
                                     popsize_df, 
                                     start_date = as.Date("2020-03-02"), 
                                     end_date = as.Date("2020-06-30"), 
                                     ld1_start = 16, ld1_end = 70, 
                                     BG1_start = 71){
  
  popsize_variable <- popsize_df %>% pull(popsize)
  
  n_days <- difftime(end_date, start_date) + 1
  
  data_for_sim <- data.frame(dept_id = rep(1:94, each = n_days), 
                             day = rep(1:n_days, 94), 
                             popsize = rep(popsize_variable, each = n_days)) %>%
    mutate(lockdown1 = ifelse(between(day, ld1_start, ld1_end), 1, 0),
           BG1 = ifelse(day >=  BG1_start, 1, 0))
  
  monolix_data <- data_for_sim %>%
    left_join(., simulation_results %>% 
                select(id, time, IncI_ME), 
              by = c("dept_id" = "id", "day" = "time")) %>%
    rename(IncI = IncI_ME) %>%
    mutate(IncI = ifelse(IncI < 0, 0, IncI)) # correct data in case some observations are < 0
  
  return(monolix_data)
}


