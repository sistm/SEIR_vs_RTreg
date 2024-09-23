library(tidyverse)
library(lme4)
library(EpiEstim)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(progress)


dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SIR_simulx_data_creation_2params"
setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs")
source("gamma.pars.from.quantiles.fun.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

popsize_df <- read.csv("popsize_df.csv") %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id-1))


load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI_neg3_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI_neg1.2_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI_neg0.4_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI0.6_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI1.6_list.RData"))



cl <- makeCluster(12)
registerDoParallel(cl)

new_tstart <- 2:100
new_tend <- 2:100
seeds <- seq(101, 600)


# Late NPI 1 datasets -----------------------------------------------------

#### initI = -3 ####
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)

res_boot_SIR_late_initI_neg3_list <- list()
for(i in 1:100){

  # load data
  reg_data <- sim_res_late_initI_neg3_list[[i]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  Rt_df <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                       Inc_col_name = "IncI_unscaled", 
                       meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                       tstart = new_tstart, tend = new_tend)
  
  
  Rt_df %<>%
    mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
           BG1 = ifelse(day > 55, 1, 0)) 
  
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = Rt_df, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  reg_res_simple_7d <- Rt_reg_only_fun(data_for_est = Rt_df %>% filter(day %in% seven_day_seq), 
                                       model_name = "simple regression 7d", 
                                       lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)

    Rt_df %<>% 
      filter(!is.na(Rt)) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    # regression with sampled values
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
    reg_res_qsampled_7d <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), 
                                data = Rt_df %>% filter(day %in% seven_day_seq))
    
    # get reg coefficients
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_Rt_qsampled_7d <- coefficients(reg_res_qsampled_7d)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_qsampled, coefs_Rt_qsampled_7d)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    mutate(median = median(value), 
           CI_LL = quantile(value, probs = 0.025), 
           CI_UL = quantile(value, probs = 0.975), 
           parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    bind_rows(reg_res_simple_7d) %>%
    mutate(rep = i)
  
  res_boot_SIR_late_initI_neg3_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_SIR_late_initI_neg3_df <- do.call("rbind.data.frame", res_boot_SIR_late_initI_neg3_list)  
save(res_boot_SIR_late_initI_neg3_df, file = "res_boot_SIR_late_initI_neg3_df.RData")



#### initI = -1.2 ####
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)

res_boot_SIR_late_initI_neg1.2_list <- list()
for(i in 1:100){

  # load data
  reg_data <- sim_res_late_initI_neg1.2_list[[i]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  Rt_df <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI_unscaled", 
                        meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                        tstart = new_tstart, tend = new_tend)
  
  
  Rt_df %<>%
    mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
           BG1 = ifelse(day > 55, 1, 0)) 
  
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = Rt_df, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  reg_res_simple_7d <- Rt_reg_only_fun(data_for_est = Rt_df %>% filter(day %in% seven_day_seq), 
                                       model_name = "simple regression 7d", 
                                       lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    Rt_df %<>% 
      filter(!is.na(Rt)) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    # regression with sampled values
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
    reg_res_qsampled_7d <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), 
                                data = Rt_df %>% filter(day %in% seven_day_seq))
    
    # get reg coefficients
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_Rt_qsampled_7d <- coefficients(reg_res_qsampled_7d)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_qsampled, coefs_Rt_qsampled_7d)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    mutate(median = median(value), 
           CI_LL = quantile(value, probs = 0.025), 
           CI_UL = quantile(value, probs = 0.975), 
           parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    bind_rows(reg_res_simple_7d) %>%
    mutate(rep = i)
  
  res_boot_SIR_late_initI_neg1.2_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_SIR_late_initI_neg1.2_df <- do.call("rbind.data.frame", res_boot_SIR_late_initI_neg1.2_list)  
save(res_boot_SIR_late_initI_neg1.2_df, file = "res_boot_SIR_late_initI_neg1.2_df.RData")


#### initI = -0.4 ####
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)

res_boot_SIR_late_initI_neg0.4_list <- list()
for(i in 1:100){

  # load data
  reg_data <- sim_res_late_initI_neg0.4_list[[i]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  Rt_df <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI_unscaled", 
                        meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                        tstart = new_tstart, tend = new_tend)
  
  
  Rt_df %<>%
    mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
           BG1 = ifelse(day > 55, 1, 0)) 
  
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = Rt_df, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  reg_res_simple_7d <- Rt_reg_only_fun(data_for_est = Rt_df %>% filter(day %in% seven_day_seq), 
                                       model_name = "simple regression 7d", 
                                       lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    Rt_df %<>% 
      filter(!is.na(Rt)) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    # regression with sampled values
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
    reg_res_qsampled_7d <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), 
                                data = Rt_df %>% filter(day %in% seven_day_seq))
    
    # get reg coefficients
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_Rt_qsampled_7d <- coefficients(reg_res_qsampled_7d)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_qsampled, coefs_Rt_qsampled_7d)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    mutate(median = median(value), 
           CI_LL = quantile(value, probs = 0.025), 
           CI_UL = quantile(value, probs = 0.975), 
           parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    bind_rows(reg_res_simple_7d) %>%
    mutate(rep = i)
  
  res_boot_SIR_late_initI_neg0.4_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_SIR_late_initI_neg0.4_df <- do.call("rbind.data.frame", res_boot_SIR_late_initI_neg0.4_list)  
save(res_boot_SIR_late_initI_neg0.4_df, file = "res_boot_SIR_late_initI_neg0.4_df.RData")



#### initI = 0.6 ####
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)

res_boot_SIR_late_initI0.6_list <- list()
for(i in 1:100){
  
  # load data
  reg_data <- sim_res_late_initI0.6_list[[i]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  Rt_df <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI", 
                        meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                        tstart = new_tstart, tend = new_tend)
  
  
  Rt_df %<>%
    mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
           BG1 = ifelse(day > 55, 1, 0)) 
  
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = Rt_df, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  reg_res_simple_7d <- Rt_reg_only_fun(data_for_est = Rt_df %>% filter(day %in% seven_day_seq), 
                                       model_name = "simple regression 7d", 
                                       lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    Rt_df %<>% 
      filter(!is.na(Rt)) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    # regression with sampled values
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
    reg_res_qsampled_7d <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), 
                                data = Rt_df %>% filter(day %in% seven_day_seq))
    
    # get reg coefficients
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_Rt_qsampled_7d <- coefficients(reg_res_qsampled_7d)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_qsampled, coefs_Rt_qsampled_7d)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    mutate(median = median(value), 
           CI_LL = quantile(value, probs = 0.025), 
           CI_UL = quantile(value, probs = 0.975), 
           parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    bind_rows(reg_res_simple_7d) %>%
    mutate(rep = i)
  
  res_boot_SIR_late_initI0.6_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_SIR_late_initI0.6_df <- do.call("rbind.data.frame", res_boot_SIR_late_initI0.6_list)  
save(res_boot_SIR_late_initI0.6_df, file = "res_boot_SIR_late_initI0.6_df.RData")



#### initI = 1.6 ####
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)

res_boot_SIR_late_initI1.6_list <- list()
for(i in 1:100){
  
  # load data
  reg_data <- sim_res_late_initI1.6_list[[i]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  Rt_df <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI", 
                        meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                        tstart = new_tstart, tend = new_tend)
  
  
  Rt_df %<>%
    mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
           BG1 = ifelse(day > 55, 1, 0)) 
  
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = Rt_df, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  reg_res_simple_7d <- Rt_reg_only_fun(data_for_est = Rt_df %>% filter(day %in% seven_day_seq), 
                                       model_name = "simple regression 7d", 
                                       lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    Rt_df %<>% 
      filter(!is.na(Rt)) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    # regression with sampled values
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
    reg_res_qsampled_7d <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), 
                                data = Rt_df %>% filter(day %in% seven_day_seq))
    
    # get reg coefficients
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_Rt_qsampled_7d <- coefficients(reg_res_qsampled_7d)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_qsampled, coefs_Rt_qsampled_7d)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    mutate(median = median(value), 
           CI_LL = quantile(value, probs = 0.025), 
           CI_UL = quantile(value, probs = 0.975), 
           parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    bind_rows(reg_res_simple_7d) %>%
    mutate(rep = i)
  
  res_boot_SIR_late_initI1.6_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_SIR_late_initI1.6_df <- do.call("rbind.data.frame", res_boot_SIR_late_initI1.6_list)  
save(res_boot_SIR_late_initI1.6_df, file = "res_boot_SIR_late_initI1.6_df.RData")

