library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(RColorBrewer)
library(ggside)
library(kableExtra)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
#source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/deSolve_SIR_model_function.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"



# Load data ---------------------------------------------------------------

popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

indparams_SEIRAHD <- read.table(paste0(dir2, "/ind_params/ind_2params_new3_1.txt"), 
                                  header = TRUE, sep = " ")



data_list_SEIRAHD <- loadRData(paste(dir2, "sim_res_Simulx_2params_new3_list.RData", sep = "/"))

dataset1_SEIRAHD <- data_list_SEIRAHD[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_SEIRAHD, by = "id") %>%
  mutate(Rt_real = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncE = transmission*S*(I+0.55*A)/10^4,
         IncE_unscaled = IncE*popsize/10^4,
         IncI_unscaled = IncI_ME*popsize/10^4, 
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
         BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_SEIRAHD <- dataset1_SEIRAHD %>% dplyr::select(dept_id, day, Rt_real)



# Estimation on 100 datasets  -------------------------------------------------------------

cl <- makeCluster(14)
registerDoParallel(cl)

reg_Rt_IncI_list <- list()
reg_Rt_IncE_list <- list()
for(i in 1:100){
  
  indparams <- read.table(paste0(dir2, "/ind_params/ind_2params_new3_", i, ".txt"), 
                                  header = TRUE, sep = " ")
  
  data_for_est <- data_list_SEIRAHD[[i]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(indparams, by = "id") %>%
    mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = IncI_ME*popsize/10^4, 
           IncE = transmission*S*(I+0.55*A)/10^4,
           IncE_unscaled = IncE*popsize/10^4,
           lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  
  Rt_comp_res_I <- Rt_calc_fun2(Inc_df = data_for_est, id_col_name = "dept_id", time_col_name = "day", 
                               Inc_col_name = "IncI_unscaled", 
                               meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
  
  Rt_comp_res_E <- Rt_calc_fun2(Inc_df = data_for_est, id_col_name = "dept_id", time_col_name = "day", 
                               Inc_col_name = "IncE_unscaled", 
                               meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
  
  first_Rt_day <- min(Rt_comp_res_E$day[!is.na(Rt_comp_res_E$Rt)]) + 3
  seven_day_seq <- seq(first_Rt_day, max(Rt_comp_res_E$day), 7)
  
  
  reg_data_I <- Rt_comp_res_I %>%
    mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
           BG1 = ifelse(day > 70, 1, 0)) %>%
    group_by(dept_id) %>%
    mutate(across(lockdown1:BG1, function(x) lag(x, 5, default = 0))) %>%
    ungroup() %>%
    filter(day %in% seven_day_seq)
  
  reg_res_I <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_I)
 
  reg_data_E <- Rt_comp_res_E %>%
    mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
           BG1 = ifelse(day > 70, 1, 0)) %>%
    filter(day %in% seven_day_seq)
  
  reg_res_E <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_E)

  
  coefs_I <- coefficients(reg_res_I)$dept_id[1, 2:3]
  coefs_E <- coefficients(reg_res_E)$dept_id[1, 2:3]
  
  
  confint_I <- data.frame(confint(reg_res_I, method="Wald"))[-c(1:2), ]
  names(confint_I) <- c("CI_LL", "CI_UL")
  
  confint_E <- data.frame(confint(reg_res_E, method="Wald"))[-c(1:2), ]
  names(confint_E) <- c("CI_LL", "CI_UL")
  
  
  reg_Rt_I <- coefs_I %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_I[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = "IncI")
  
  reg_Rt_E <- coefs_E %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_E[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = "IncE")

  reg_Rt_IncE_list[[i]] <- reg_Rt_E
}

stopCluster(cl)


summary_regs_IncE <- bind_rows(reg_Rt_IncI_list) %>%
  bind_rows(bind_rows(reg_Rt_IncE_list)) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         unique_sims = 100, 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  group_by(parameter, model) %>%
  mutate(perc_CI_covers = sum(CI_covers)/unique_sims*100, 
         mean_bias = mean(bias), 
         mean_rel_bias = mean(rel_bias))


metrics_regs_IncE <- summary_regs_IncE %>%
  dplyr::select(parameter, model, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(coverage = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metrics_regs_IncE %>% 
  ungroup() %>%
  dplyr::select(-parameter) %>%
  kable(digits = 2, format = "html"#, table.attr = "style='width:62%;'"
  ) %>%
  kable_styling(bootstrap_options = "striped") %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)



fitted_vals_I <- exp(fitted(reg_res_I))
fitted_df_I <- reg_data_I %>%
  mutate(Rt_fitted = fitted_vals_I) %>%
  select(day, dept_id, Rt_fitted)

dataset1_SEIRAHD_fitted_I <- dataset1_SEIRAHD %>%
  select(dept_id, day, Rt_real, IncI_unscaled) %>%
  left_join(Rt_comp_res_I %>% select(dept_id, day, Rt), by = c("dept_id", "day")) %>%
  left_join(fitted_df_I, by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


fitted_vals_E <- exp(fitted(reg_res_E))
fitted_df_E <- reg_data_E %>%
  mutate(Rt_fitted = fitted_vals_E) %>%
  select(day, dept_id, Rt_fitted)

dataset1_SEIRAHD_fitted_E <- dataset1_SEIRAHD %>%
  select(dept_id, day, Rt_real, IncE_unscaled) %>%
  left_join(Rt_comp_res_E %>% select(dept_id, day, Rt), by = c("dept_id", "day")) %>%
  left_join(fitted_df_E, by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))

selected_depts <- c(2, 14, 24, 51, 76, 90) 

ggplot(dataset1_SEIRAHD_fitted_E %>% filter(dept_id %in% selected_depts), 
       aes(x = day)) +
  geom_line(aes(y = Rt, col = "Epi Estim"), linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit")) + 
  geom_line(aes(y = Rt_real, col = "True"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncE_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(title = "IncE fits", 
       col = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_color_manual(values = c("blue", "deeppink", "black")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        legend.key.size = unit(1, "cm"),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggplot(dataset1_SEIRAHD_fitted_I %>% filter(dept_id %in% selected_depts), 
       aes(x = day)) +
  geom_line(aes(y = Rt, col = "Epi Estim"), linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit")) + 
  geom_line(aes(y = Rt_real, col = "True"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(title = "IncI fits", 
       col = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_color_manual(values = c("blue", "deeppink", "black")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        legend.key.size = unit(1, "cm"),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))
