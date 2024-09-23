sim_SEIRAHD_Simulx_ME_2params <- function(tmax = 121, 
                                          path_to_ind_params, 
                                          regressor_df_path, 
                                          ME_I_a = 0.0408, ME_I_b = 0.1, 
                                          ME_Hin_a = 3.62e-03, ME_Hin_b = 0.05, 
                                          ME_Hprev_b = 0.139,  
                                          ME_D_a = 0.000191, ME_D_b = 0.0754){
  

  # this function requires a Simulx project already loaded!
  defineIndividualElement(name = "mlx_Pop_ind", element = path_to_ind_params)
  defineRegressorElement(name = "mlx_Reg", element = regressor_df_path)
  
  
  defineOutputElement(name = "A_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "A"))
  defineOutputElement(name = "E_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "E"))
  defineOutputElement(name = "I_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "I"))
  defineOutputElement(name = "S_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "S"))
  defineOutputElement(name = "H_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "H"))
  defineOutputElement(name = "R_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "R"))
  defineOutputElement(name = "D_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "D"))
  defineOutputElement(name = "transmission", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "transmission"))
  defineOutputElement(name = "predI", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predI"))
  defineOutputElement(name = "predHin", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predHin"))
  defineOutputElement(name = "predHprev", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predHprev"))
  defineOutputElement(name = "predD", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predD"))
  
  
  setGroupElement(group = "simulationGroup1", 
                  elements = c("mlx_Pop_ind"))
  
  setGroupElement(group = "simulationGroup1", 
                  elements = c("A_mlx", "S_mlx", "H_mlx", "R_mlx", "I_mlx", "E_mlx", "D_mlx", 
                               "transmission", "predI", "predHin", "predHprev", "predD"))
  
  setGroupRemaining(group = "simulationGroup1",
                    remaining = list(apredHin = 0, bpredHin = 0, 
                                     apredI = 0, bpredI = 0, 
                                     apredD = 0, bpredD = 0, 
                                     bpredHprev = 0))
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
    rename(IncI = predI, 
           IncD = predD, 
           IncH = predHin,
           PrevH = predHprev) %>%
    rowwise() %>%
    mutate(IncI_ME = max(rnorm(1, IncI, (ME_I_a + ME_I_b*IncI)/1.96), 0),
           IncH_ME = max(rnorm(1, IncH, (ME_Hin_a + ME_Hin_b*IncH)/1.96), 0),
           PrevH_ME = max(rnorm(1, PrevH, (ME_Hprev_b*PrevH)/1.96), 0),
           IncD_ME = max(rnorm(1, IncD, (ME_D_b*IncD)/1.96), 0)) %>%
    ungroup()

}


monolix_data_creation_ME <- function(simulation_results, 
                                     popsize_df, 
                                     start_date = as.Date("2020-03-02"), 
                                     end_date = as.Date("2020-06-30")){
  
  popsize_variable <- popsize_df %>% pull(popsize)
  
  data_for_sim <- data.frame(dept_id = rep(1:94, each = 121), 
                             day = rep(1:121, 94), 
                             popsize = rep(popsize_variable, each = 121)) %>%
    mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0),
           BG1 = ifelse(day > 70, 1, 0))
  
  monolix_data <- data_for_sim %>%
    left_join(., simulation_results %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("dept_id" = "id", "day" = "time")) %>%
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  return(monolix_data)
}


sim_SEIRAHD_Simulx_ME2_2params <- function(tmax = 121, 
                                           path_to_ind_params, 
                                           regressor_df_path, 
                                           ME_I_a = 0.0408, ME_I_b = 0.1, 
                                           ME_Hin_a = 3.62e-03, ME_Hin_b = 0.05, 
                                           ME_Hprev_b = 0.139,  
                                           ME_D_a = 0.000191, ME_D_b = 0.0754){
  
  # this function requires a Simulx project already loaded!
  defineIndividualElement(name = "mlx_Pop_ind", element = path_to_ind_params)
  defineRegressorElement(name = "mlx_Reg", element = regressor_df_path)
  
  
  defineOutputElement(name = "A_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "A"))
  defineOutputElement(name = "E_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "E"))
  defineOutputElement(name = "I_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "I"))
  defineOutputElement(name = "S_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "S"))
  defineOutputElement(name = "H_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "H"))
  defineOutputElement(name = "R_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "R"))
  defineOutputElement(name = "D_mlx", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "D"))
  defineOutputElement(name = "transmission", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "transmission"))
  defineOutputElement(name = "predI", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predI"))
  defineOutputElement(name = "predHin", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predHin"))
  defineOutputElement(name = "predHprev", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predHprev"))
  defineOutputElement(name = "predD", 
                      element = list(data = data.frame(time = seq(1:tmax)), output = "predD"))
  
  
  setGroupElement(group = "simulationGroup1", 
                  elements = c("mlx_Pop_ind"))
  
  setGroupElement(group = "simulationGroup1", 
                  elements = c("A_mlx", "S_mlx", "H_mlx", "R_mlx", "I_mlx", "E_mlx", "D_mlx", 
                               "transmission", "predI", "predHin", "predHprev", "predD"))
  
  setGroupRemaining(group = "simulationGroup1",
                    remaining = list(apredHin = 0, bpredHin = 0, 
                                     apredI = 0, bpredI = 0, 
                                     apredD = 0, bpredD = 0, 
                                     bpredHprev = 0))

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
    rename(IncI = predI, 
           IncD = predD, 
           IncH = predHin,
           PrevH = predHprev) 
  
}


monolix_data_creation <- function(simulation_results, 
                                  popsize_df, 
                                  start_date = as.Date("2020-03-02"), 
                                  end_date = as.Date("2020-06-30")){
  
  popsize_variable <- popsize_df %>% pull(popsize)
  
  data_for_sim <- data.frame(dept_id = rep(1:94, each = 121), 
                             day = rep(1:121, 94), 
                             popsize = rep(popsize_variable, each = 121)) %>%
    mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0),
           BG1 = ifelse(day > 70, 1, 0))
  
  monolix_data <- data_for_sim %>%
    left_join(., simulation_results %>% 
                select(id, time, IncI, IncH, PrevH, IncD), 
              by = c("dept_id" = "id", "day" = "time")) %>%
    group_by(dept_id) %>%
    mutate(initH = PrevH[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI, IncH, PrevH, IncD), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH" ~ 1, 
                              obs_id == "PrevH" ~ 2,
                              obs_id == "IncI" ~ 3, 
                              obs_id == "IncD" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  return(monolix_data)
}