library(tidyverse)
library(lme4)
library(EpiEstim)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(progress)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs")
dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once10"
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/gamma.pars.from.quantiles.fun.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id-1))


# ABM10 hybrid 7d --------------------------------------------------------

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)


cl <- makeCluster(16)
registerDoParallel(cl)

seeds <- seq(101, 600)

res_boot_ABM_h_list_7d <- list()
for(i in 1:100){
  
  pb$tick(0)
  # load data
  reg_data_ABM_h <- read.csv(paste0(dir1, "/data_covasim_hybrid10_Rt_", i, ".csv")) %>%
    filter(day > 29) %>%
    mutate(day = day - 29)
  
  
  # Estimate Rt 
  Rt_df <- Rt_calc_fun2(reg_data_ABM_h, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI", meansi = 7.8, stdsi = 4.4, 
                        Rt_prior = 2, Rt_sd_prior = 4)
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  reg_data_ABM_h_new <-  reg_data_ABM_h %>%
    dplyr::select(dept_id, day, lockdown1, BG1) %>%
    group_by(dept_id) %>%
    mutate(lockdown1 = lag(lockdown1, n = 5, default = 0), 
           BG1 = lag(BG1, n = 5, default = 0)) %>%
    ungroup() %>%
    unique() %>%
    left_join(., Rt_df, by = c("dept_id", "day")) %>%
    filter(day %in% seven_day_seq) 
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = reg_data_ABM_h_new, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    # randomly sample values for Rt from its distribution: 
    # 1. by randomly sampling one value each day
    # 2. by sampling a constant quantile of the distribution
    Rt_df %<>% 
      filter(!is.na(Rt))  %>%
      filter(day %in% seven_day_seq) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate), 
        Rt_rsampled = rgamma(1, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    
    reg_data_ABM_h_new <- reg_data_ABM_h  %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = 5, default = 0), 
             BG1 = lag(BG1, n = 5, default = 0)) %>%
      ungroup() %>%
      unique() %>%
      filter(day %in% seven_day_seq) %>%
      left_join(., Rt_df, by = c("dept_id", "day"))
    
    # regression with sampled values
    reg_res_rsampled <- lmer(log(Rt_rsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM_h_new)
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM_h_new)
    
    
    # get reg coefficients
    coefs_Rt_rsampled <- coefficients(reg_res_rsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "random")
    
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_sampled <- rbind(coefs_Rt_rsampled, coefs_Rt_qsampled)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    summarize(median = median(value), 
              CI_LL = quantile(value, probs = 0.025), 
              CI_UL = quantile(value, probs = 0.975)) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    mutate(rep = i)
  
  res_boot_ABM_h_list_7d[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_df_ABM10_h_7d <- do.call("rbind.data.frame", res_boot_ABM_h_list_7d)  
save(res_boot_df_ABM10_h_7d, file = "res_boot_df_ABM10_h_7d.RData")




# ABM10 random mixing 7d --------------------------------------------------------

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)


cl <- makeCluster(16)
registerDoParallel(cl)

seeds <- seq(101, 600)

res_boot_ABM_rm_list_7d <- list()
for(i in 1:100){
  
  pb$tick(0)
  # load data
  reg_data_ABM_rm <- read.csv(paste0(dir1, "/data_covasim_rm10_Rt_", i, ".csv")) %>%
    filter(day > 29) %>%
    mutate(day = day - 29)
  
  
  # Estimate Rt 
  Rt_df <- Rt_calc_fun2(reg_data_ABM_rm, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI", meansi = 8.45, stdsi = 5.45, Rt_prior = 1, 
                        Rt_sd_prior = 2)
  
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  reg_data_ABM_rm_new <-  reg_data_ABM_rm %>%
    dplyr::select(dept_id, day, lockdown1, BG1) %>%
    group_by(dept_id) %>%
    mutate(lockdown1 = lag(lockdown1, n = 5, default = 0), 
           BG1 = lag(BG1, n = 5, default = 0)) %>%
    ungroup() %>%
    unique() %>%
    left_join(., Rt_df, by = c("dept_id", "day")) %>%
    filter(day %in% seven_day_seq) 

  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = reg_data_ABM_rm_new, model_name = "simple regression", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    # randomly sample values for Rt from its distribution: 
    # 1. by randomly sampling one value each day
    # 2. by sampling a constant quantile of the distribution
    Rt_df %<>% 
      filter(!is.na(Rt))  %>% 
      filter(day %in% seven_day_seq) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate), 
        Rt_rsampled = rgamma(1, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    
    reg_data_ABM_rm_new <- reg_data_ABM_rm %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = 5, default = 0), 
             BG1 = lag(BG1, n = 5, default = 0)) %>%
      ungroup() %>%
      unique() %>% 
      filter(day %in% seven_day_seq) %>%
      left_join(., Rt_df, by = c("dept_id", "day"))
    
    # regression with sampled values
    reg_res_rsampled <- lmer(log(Rt_rsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM_rm_new)
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM_rm_new)
    
    
    # get reg coefficients
    coefs_Rt_rsampled <- coefficients(reg_res_rsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "random")
    
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile")
    
    coefs_sampled <- rbind(coefs_Rt_rsampled, coefs_Rt_qsampled)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    summarize(median = median(value), 
              CI_LL = quantile(value, probs = 0.025), 
              CI_UL = quantile(value, probs = 0.975)) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    mutate(rep = i)
  
  res_boot_ABM_rm_list_7d[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_df_ABM10_rm_7d <- do.call("rbind.data.frame", res_boot_ABM_rm_list_7d)  
save(res_boot_df_ABM10_rm_7d, file = "res_boot_df_ABM10_rm_7d.RData")


# SEIRAHD 3 data 7d IncI ------------------------------------------------------------
load(paste0(dir2, "/sim_res_Simulx_2params_new3_list.RData"))

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)


cl <- makeCluster(16)
registerDoParallel(cl)

seeds <- seq(101, 600)

res_boot_Simulx3_IncI_7d_list <- list()
for(i in 1:100){
  
  pb$tick(0)
  # load data
  reg_data_Simulx3 <- sim_res_Simulx_2params_new3_list[[i]] %>%
    dplyr::select(id, time, IncI_ME) %>%
    rename(IncI = IncI_ME, dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI = IncI*popsize/10^4, 
           lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
           BG1 = ifelse(day > 70, 1, 0))
  
  
  # Estimate Rt 
  Rt_df <- Rt_calc_fun2(reg_data_Simulx3, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncI", meansi = 10.1, stdsi = 8.75, Rt_prior = 1, 
                        Rt_sd_prior = 2)
  
  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  
  reg_data_Simulx3_new <- reg_data_Simulx3 %>%
    dplyr::select(dept_id, day, lockdown1, BG1) %>%
    group_by(dept_id) %>%
    mutate(lockdown1 = lag(lockdown1, n = 5, default = 0), 
           BG1 = lag(BG1, n = 5, default = 0)) %>%
    ungroup() %>%
    unique() %>%
    left_join(., Rt_df, by = c("dept_id", "day")) %>%
    filter(day %in% seven_day_seq) 
  
  
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est = reg_data_Simulx3_new, model_name = "simple regression 7d", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    # randomly sample values for Rt from its distribution: 
    # 1. by randomly sampling one value each day
    # 2. by sampling a constant quantile of the distribution
    Rt_df %<>% 
      filter(!is.na(Rt)) %>%
      filter(day %in% seven_day_seq) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate), 
        Rt_rsampled = rgamma(1, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    reg_data_Simulx3_new <- reg_data_Simulx3 %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = 5, default = 0), 
             BG1 = lag(BG1, n = 5, default = 0)) %>%
      ungroup() %>%
      unique()  %>%
      filter(day %in% seven_day_seq) %>%
      left_join(., Rt_df, by = c("dept_id", "day"))
    
    # regression with sampled values
    reg_res_rsampled <- lmer(log(Rt_rsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_Simulx3_new)
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_Simulx3_new)
    
    
    # get reg coefficients
    coefs_Rt_rsampled <- coefficients(reg_res_rsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "random 7d")
    
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_rsampled, coefs_Rt_qsampled)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    summarize(median = median(value), 
              CI_LL = quantile(value, probs = 0.025), 
              CI_UL = quantile(value, probs = 0.975)) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    mutate(rep = i)
  
  res_boot_Simulx3_IncI_7d_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_df_Simulx3_IncI_7d <- do.call("rbind.data.frame", res_boot_Simulx3_IncI_7d_list)  
save(res_boot_df_Simulx3_IncI_7d, file = "res_boot_df_Simulx3_IncI_7d.RData")



# SEIRAHD 3 data 7d IncH --------------------------------------------------


pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = 100,    
                       clear = TRUE,    
                       width = 100)


cl <- makeCluster(16)
registerDoParallel(cl)

seeds <- seq(101, 600)

res_boot_Simulx3_IncH_7d_list <- list()
for(i in 1:100){
  
  pb$tick(0)
  # load data
  reg_data_Simulx3 <- sim_res_Simulx_2params_new3_list[[i]] %>%
    dplyr::select(id, time, IncH_ME) %>%
    rename(IncH = IncH_ME, dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncH = IncH*popsize/10^4, 
           lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
           BG1 = ifelse(day > 70, 1, 0))
  
  
  # Estimate Rt 
  Rt_df <- Rt_calc_fun2(reg_data_Simulx3, id_col_name = "dept_id", time_col_name = "day", 
                        Inc_col_name = "IncH", meansi = 10.1, stdsi = 8.75, Rt_prior = 1, 
                        Rt_sd_prior = 2)

  first_day_with_R_est <- min(Rt_df$day[!is.na(Rt_df$Rt)]) + 3
  seven_day_seq <- seq(first_day_with_R_est, max(Rt_df$day), 7)
  
  
  reg_data_Simulx3_new <- reg_data_Simulx3 %>%
    dplyr::select(dept_id, day, lockdown1, BG1) %>%
    group_by(dept_id) %>%
    mutate(lockdown1 = lag(lockdown1, n = 10, default = 0), 
           BG1 = lag(BG1, n = 10, default = 0)) %>%
    ungroup() %>%
    unique() %>%
    left_join(., Rt_df, by = c("dept_id", "day")) %>%
    filter(day %in% seven_day_seq) 
  
  # run normal regression
  reg_res_simple <- Rt_reg_only_fun(data_for_est =  reg_data_Simulx3_new, model_name = "simple regression 7d", 
                                    lag_NPIs = FALSE, fits = FALSE) %>%
    select(-lag) %>%
    rename(median = value, method = model)
  
  
  # include uncertainty in regression estimation by repeatedly sampling from Rt distribution
  reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
    
    # select quantile for Rt quantile sampling
    set.seed(seeds[k])
    sel_quantile <- round(rbeta(1, shape1 = 2, shape2 = 2), 3)
    
    # randomly sample values for Rt from its distribution: 
    # 1. by randomly sampling one value each day
    # 2. by sampling a constant quantile of the distribution
    Rt_df %<>% 
      filter(!is.na(Rt))  %>%
      filter(day %in% seven_day_seq) %>%
      rowwise() %>%
      mutate(
        # derive underlying gamma distribution
        gshape = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$shape, 
        grate = gamma.parms.from.quantiles(q = c(Rt, CI_UL), p = c(0.5, 0.975))$rate, 
        # sample
        Rt_qsampled = qgamma(sel_quantile, shape = gshape, rate = grate), 
        Rt_rsampled = rgamma(1, shape = gshape, rate = grate)
      ) %>%
      ungroup()
    
    reg_data_Simulx3_new <- reg_data_Simulx3 %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = 10, default = 0), 
             BG1 = lag(BG1, n = 10, default = 0)) %>%
      ungroup() %>%
      unique() %>%
      left_join(., Rt_df, by = c("dept_id", "day")) %>%
      filter(day %in% seven_day_seq) 
    
    # regression with sampled values
    reg_res_rsampled <- lmer(log(Rt_rsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_Simulx3_new)
    reg_res_qsampled <- lmer(log(Rt_qsampled) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_Simulx3_new)
    
    
    # get reg coefficients
    coefs_Rt_rsampled <- coefficients(reg_res_rsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "random 7d")
    
    coefs_Rt_qsampled <- coefficients(reg_res_qsampled)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(boot_rep = k, 
             method = "quantile 7d")
    
    coefs_sampled <- rbind(coefs_Rt_rsampled, coefs_Rt_qsampled)
    
  }
  
  # assemble parameters
  reg_sampled_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
    group_by(parameter, method) %>%
    summarize(median = median(value), 
              CI_LL = quantile(value, probs = 0.025), 
              CI_UL = quantile(value, probs = 0.975)) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    bind_rows(reg_res_simple) %>%
    mutate(rep = i)
  
  res_boot_Simulx3_IncH_7d_list[[i]] <- reg_sampled_boot_df
  
  pb$tick()
}

res_boot_df_Simulx3_IncH_7d <- do.call("rbind.data.frame", res_boot_Simulx3_IncH_7d_list)  
save(res_boot_df_Simulx3_IncH_7d, file = "res_boot_df_Simulx3_IncH_7d.RData")

