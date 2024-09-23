EpiEstim_reg_fun <- function(data_for_est, Inc_name, rep_num, lag_NPIs = FALSE, lag_days = 0, 
                             meansi = 7.5, stdsi = 5, meanprior = 2, stdprior = 4, cut_days = 0){

  
  Rt_list <- foreach(i = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
    Inc_series <- data_for_est %>% 
      filter(dept_id == i) 
    
    Inc <- as.numeric(na.omit(Inc_series[[Inc_name]]*Inc_series$popsize/10^4))
    
    Rt_estim <- estimate_R(Inc,
                           method = "parametric_si", 
                           config = make_config(list(
                             mean_si = meansi,
                             std_si = stdsi, 
                             mean_prior = meanprior, 
                             std_prior = stdprior)))$R
    Rt_estim <- Rt_estim %>%
      mutate(t = (t_end+t_start)/2,
             id = i, .before = t_start) %>%
      dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
      rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  if(lag_NPIs){
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = lag_days, default = 0), 
             BG1 = lag(BG1, n = lag_days, default = 0)) %>%
      ungroup() %>%
      unique() %>%
      left_join(., Rt_df, by = c("dept_id", "day")) %>%
      filter(day > cut_days)
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  } else{
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      unique() %>%
      left_join(., Rt_df, by = c("dept_id", "day")) %>%
      filter(day > cut_days)
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  }
  
  
  coefs_Rt_reg <- coefficients(reg_res)$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_res, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  coefs_Rt_reg <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("Lockdown 1", "Barrier gestures"))) %>%
    mutate(rep = rep_num)
  
  return(coefs_Rt_reg)
  
}


popparam_cleaning <- function(df){
  df_clean <- df %>%
    mutate(parameter = str_remove(parameter, "_pop"), 
           value = 0 - value) %>% # invert the sign of the parameter because it's estimated different in Monolix
    filter(grepl("beta", parameter) & !grepl("omega", parameter)) 
  
  return(df_clean)
}

bootstrap_cleaning <- function(df, boot_rep = i, sim_rep = j){
  df_clean <- df %>%
    mutate(parameter = str_remove(parameter, "_pop"), 
          value = 0 - value) %>% # invert the sign of the parameter because it's estimated different in Monolix
    filter(grepl("beta", parameter) & !grepl("omega", parameter)) %>%
    mutate(boot_rep = boot_rep, sim_rep = sim_rep)
  
  return(df_clean)
}

bootstrap_CI_calc <- function(boot_list){
  boot_df <- do.call("rbind.data.frame", boot_list) %>%
    group_by(parameter, sim_rep) %>%
    summarize(mean_est1 = mean(value), 
              sd_est = sd(value), 
              CI_LL2 = quantile(value, probs = 0.025), 
              CI_UL2 = quantile(value, probs = 0.975)) %>%
    mutate(CI_LL1 = mean_est1 - 1.96*sd_est, 
           CI_UL1 = mean_est1 + 1.96*sd_est)
  
  return(boot_df)
}


bootstrap_summary <- function(bootstrap_list, true_val_NPI1 = -1.45, true_val_NPI2 = -0.5){
  boot_df <- do.call("rbind.data.frame", bootstrap_list) %>%
    ungroup() %>% 
    mutate(parameter = factor(parameter, 
                              levels = c("beta_ld1", "beta_BG1"),
                              labels = c("NPI 1", "NPI 2")), 
           true_value = ifelse(parameter == "NPI 1", true_val_NPI1, true_val_NPI2), 
           unique_sims = n()/2, 
           CI_covers = ifelse(between(true_value, CI_LL2, CI_UL2), 1, 0), 
           bias = true_value - mean_est2, 
           rel_bias = abs(true_value - mean_est2)/abs(true_value)*100) %>%
    group_by(parameter) %>%
    mutate(perc_CI_covers = sum(CI_covers)/unique_sims*100, 
           mean_bias = mean(bias, na.rm = TRUE), 
           mean_rel_bias = mean(rel_bias, na.rm = TRUE))
  
  return(boot_df)
}


reg_summary <- function(reg_df, true_val_NPI1 = -1.45, true_val_NPI2 = -0.5){
  reg_df %<>% 
    mutate(parameter = factor(parameter, 
                              levels = c("Lockdown 1", "Barrier gestures"), 
                              labels = c("NPI 1", "NPI 2")), 
           true_value = ifelse(parameter == "NPI 1", true_val_NPI1, true_val_NPI2), 
           unique_sims = length(unique(rep)), 
           CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
           bias = true_value - value, 
           rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
    group_by(parameter, model) %>%
    mutate(perc_CI_covers = sum(CI_covers)/unique_sims*100, 
              mean_bias = mean(bias), 
              mean_rel_bias = mean(rel_bias))
    
  return(reg_df)
}


SEIR_summary_2params <- function(SEIR_list){
  df <- do.call("rbind.data.frame", SEIR_list) %>%
    mutate(parameter = str_remove(parameter, "_pop"))
  
  df_boot <- df %>%
    group_by(parameter, model) %>%
    mutate(value = 0-value) %>%
    summarize(median_est = median(value),
              CI_LL = quantile(value, probs = 0.025),
              CI_UL = quantile(value, probs = 0.975)) %>%
    mutate(parameter = str_remove(parameter, "_pop")) %>%
    filter(grepl("beta", parameter)) %>%
    filter(!grepl("omega", parameter)) %>%
    mutate(parameter = str_remove(parameter, "beta_")) %>%
    mutate(parameter = factor(parameter, 
                              levels = c("ld1", "BG1"), 
                              labels = c("Lockdown 1", "Barrier gestures"))) %>%
    rename(value = median_est)
}


SEIR_summary_2params2 <- function(SEIR_list){
  df <- do.call("rbind.data.frame", SEIR_list) %>%
    mutate(parameter = str_remove(parameter, "_pop")) %>%
    filter(grepl("beta", parameter) & !grepl("omega", parameter)) %>%
    mutate(parameter = str_remove(parameter, "beta_")) %>%
    mutate(parameter = factor(parameter, 
                              levels = c("ld1", "BG1"), 
                              labels = c("Lockdown 1", "Barrier gestures")), 
           value = 0-value)
}


calc_Rt_SEIR <- function(b1, gamma, S, N){
  R_eff <- b1*S/(gamma*N)
  return(R_eff)
}


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


EpiEstim_only_fun <- function(data_for_est, Inc_name, meansi = 7.5, stdsi = 5, 
                              meanprior = 2, stdprior = 4){
  
  
  Rt_list <- foreach(i = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
    Inc_series <- data_for_est %>% 
      filter(dept_id == i) 
    
    Inc <- as.numeric(na.omit(Inc_series[[Inc_name]]*Inc_series$popsize/10^4))
    
    Rt_estim <- estimate_R(Inc,
                           method = "parametric_si",
                           config = make_config(list(
                             mean_si = meansi,
                             std_si = stdsi, 
                             mean_prior = meanprior, 
                             std_prior = stdprior)))$R
    Rt_estim <- Rt_estim %>%
      mutate(t = (t_end+t_start)/2,
             id = i, .before = t_start) %>%
      dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
      rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  return(Rt_df)
}




Rt_calc_fun <- function(Inc_df, id_col_name, time_col_name, Inc_col_name, model_name, 
                        tstart = NULL, tend = NULL, Rt_ref_df, cut_days = 0,
                        Rt_prior = 1, Rt_sd_prior = 2, meansi = 10.1, stdsi = 8.75){
  # select input data according to desired input for Rt calculation
  Inc_data <- Inc_df[c(id_col_name, time_col_name, Inc_col_name)]
  colnames(Inc_data) <- c("id", "time", "Inc")
  
  # Rt calculation
  Rt_list <- vector(mode = "list")
  
  if(is_null(tstart) | (is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(mean_si = meansi, 
                                                         std_si = stdsi, 
                                                         mean_prior = Rt_prior, 
                                                         std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = (t_end+t_start)/2 + cut_days,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  } 
  
  if(!is_null(tstart) & !(is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(
                                 mean_si = meansi, std_si = stdsi, 
                                 t_start = tstart, t_end = tend, 
                                 mean_prior = Rt_prior, 
                                 std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = (t_end+t_start)/2 + cut_days,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  Rt_comp_df <- Rt_df %>%
    left_join(., Rt_ref_df, by = c("dept_id", "day")) %>%
    mutate(Rt_residual = Rt - Rt_real, 
           model = model_name)
  
  RMSE_df <- Rt_comp_df %>%
    filter(!is.na(Rt)) %>%
    group_by(dept_id, model) %>%
    summarize(RMSE = sqrt(sum(Rt_residual^2)/n())) %>%
    ungroup() %>%
    mutate(RMSE_total = mean(RMSE))
  
  return(list(Rt_comp = Rt_comp_df, RMSE = RMSE_df))
  
}


Rt_reg_only_fun <- function(data_for_est, lag_NPIs = FALSE, lag_days = 0, model_name, 
                            fits = FALSE){
  
  if(lag_NPIs){
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, Rt, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = lag_days, default = 0), 
             BG1 = lag(BG1, n = lag_days, default = 0)) %>%
      ungroup()
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  } else{
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, Rt, lockdown1, BG1) 
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  }
  
  
  coefs_Rt_reg <- coefficients(reg_res)$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_res, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  coefs_Rt_reg <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(lag = lag_days, model = model_name)
  
  
  if(fits){
    fitted_vals <- exp(fitted(reg_res))
    fitted_df <- data_for_est %>%
      filter(!is.na(Rt)) %>%
      mutate(Rt_fitted = fitted_vals, lag = lag_days, model = model_name)
    
    return(list(coefs = coefs_Rt_reg, fits = fitted_df))
  } else{
    
    return(coefs_Rt_reg)
  }
  
}



Rt_calc_fun2 <- function(Inc_df, id_col_name, time_col_name, Inc_col_name, 
                         tstart = NULL, tend = NULL, cut_days = 0,
                         Rt_prior = 1, Rt_sd_prior = 2, meansi = 10.1, stdsi = 8.75){
  # select input data according to desired input for Rt calculation
  Inc_data <- Inc_df[c(id_col_name, time_col_name, Inc_col_name)]
  colnames(Inc_data) <- c("id", "time", "Inc")
  
  # Rt calculation
  Rt_list <- vector(mode = "list")
  
  if(is_null(tstart) | (is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(mean_si = meansi, 
                                                         std_si = stdsi, 
                                                         mean_prior = Rt_prior, 
                                                         std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = (t_end+t_start)/2 + cut_days,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Median(R)`, `Std(R)`,  `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, Rt_median = `Median(R)`, sd = `Std(R)`,  CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  } 
  
  if(!is_null(tstart) & !(is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(
                                 mean_si = meansi, std_si = stdsi, 
                                 t_start = tstart, t_end = tend, 
                                 mean_prior = Rt_prior, 
                                 std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = (t_end+t_start)/2,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(id, t, t_start, t_end, `Mean(R)`, `Median(R)`, `Std(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, Rt_median = `Median(R)`, sd = `Std(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  
  return(Rt_df)
  
}


Rt_reg_only_fun_7d <- function(data_for_est, lag_NPIs = FALSE, lag_days = 0, model_name, 
                               fits = FALSE){
  
  first_day_with_R_est <- min(data_for_est$day[!is.na(data_for_est$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(data_for_est$day), 7)
  
  if(lag_NPIs){
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, Rt, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = lag_days, default = 0), 
             BG1 = lag(BG1, n = lag_days, default = 0)) %>%
      ungroup() %>%
      filter(day %in% seven_day_seq) 
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  } else{
    
    reg_data <- data_for_est %>%
      filter(day %in% seven_day_seq) %>%
      dplyr::select(dept_id, day, Rt, lockdown1, BG1) 
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  }
  
  
  coefs_Rt_reg <- coefficients(reg_res)$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_res, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  coefs_Rt_reg <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(lag = lag_days, model = model_name)
  
  
  if(fits){
    fitted_vals <- exp(fitted(reg_res))
    fitted_df <- data_for_est %>%
      filter(!is.na(Rt) & day %in% seven_day_seq) %>%
      mutate(Rt_fitted = fitted_vals, lag = lag_days, model = model_name)
    
    return(list(coefs = coefs_Rt_reg, fits = fitted_df))
  } else{
    
    return(coefs_Rt_reg)
  }
  
}


my_plot_theme_big <- function(){
  
  theme_bw() %+replace%
    theme(
      plot.title = element_text(family = "serif", size = 20), 
      axis.title = element_text(family = "serif", size = 15), 
      axis.text.x = element_text(family = "serif", size = 13), 
      axis.text.y = element_text(family = "serif", size = 13), 
      legend.text = element_text(family = "serif", size = 14.5), 
      strip.text = element_text(family = "serif", size = 14.5)
    )
}


my_plot_theme_smaller <- function(){
  
  theme_bw() %+replace%
    theme(
      plot.title = element_text(family = "serif", size = 16, hjust = 0), 
      axis.title = element_text(family = "serif", size = 13), 
      axis.text.x = element_text(family = "serif", size = 11), 
      axis.text.y = element_text(family = "serif", size = 11), 
      legend.text = element_text(family = "serif", size = 12), 
      strip.text = element_text(family = "serif", size = 12)
    )
}



two_step_fun <- function(data_list, dataset_number, popsize_df, model_name, smooth = FALSE, 
                         weekly_datapoints = FALSE, NPI1_start, NPI1_stop, Inc_name, cut_days = 0){
  
  reg_data <- data_list[[dataset_number]]  %>% 
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  true_Rt_df <- reg_data %>%
    mutate(Rt_real = transmission*S/(0.2*10000)) %>%
    dplyr::select(dept_id, day, Rt_real)
  
  
  if(smooth){
    
    Rt_df <- Rt_calc_fun(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = Inc_name, model_name = model_name, 
                         Rt_ref_df = true_Rt_df, 
                         meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)
  } else {

    new_tstart <- 2:100
    new_tend <- 2:100
    
    Rt_df <- Rt_calc_fun(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = Inc_name, model_name = model_name, 
                         Rt_ref_df = true_Rt_df, 
                         meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                         tstart = new_tstart, tend = new_tend)
  }
  

  if(weekly_datapoints){
    reg_res <- Rt_reg_only_fun_7d(data_for_est = Rt_df$Rt_comp %>%
                                    mutate(lockdown1 = ifelse(between(day, NPI1_start, NPI1_stop), 1, 0), 
                                           BG1 = ifelse(day > NPI1_stop, 1, 0)) %>%
                                    filter(day > cut_days), 
                                  model_name = model_name, lag_NPIs = FALSE, fits = FALSE)
  } else {
    reg_res <- Rt_reg_only_fun(data_for_est = Rt_df$Rt_comp %>%
                                 mutate(lockdown1 = ifelse(between(day, NPI1_start, NPI1_stop), 1, 0), 
                                        BG1 = ifelse(day > NPI1_stop, 1, 0)) %>%
                                 filter(day > cut_days), 
                               model_name = model_name, lag_NPIs = FALSE, fits = FALSE)
  }
  
  res_df <- reg_res %>%
    mutate(smoothed = ifelse(smooth, 1, 0), 
           weekly = ifelse(weekly_datapoints, 1, 0)) %>%
    ungroup()
  
  return(res_df)
}


calc_Rt <- function(b1, S, Di = 5, alpha = 0.55, Dq, re = 0.844, risk_hosp, 
                    VE_I, VE_H, ne = 10000){
  R_eff <- b1*(1-VE_I)*S/ne*(Di*alpha*(1-re)+(Di*Dq*re)/((1-risk_hosp)*Dq+Di*(1-VE_H)*risk_hosp))
  return(R_eff)
}
