library(tidyverse)
library(kableExtra)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/SIR_simulx_data_creation_2params")
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


# Load 2-step regression results --------------------------------------------------------------

dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"

load(paste0(dir3, "/res_boot_SIR_late_InitI_neg3_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI_neg1.2_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI_neg0.4_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI0.6_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI1.6_df.RData"))



# Load bootstrap results from mechanistic model ---------------------------

dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SIR_bootstrap"

load(paste0(dir4, "/SIR_boot_list_initI_neg3.RData"))
load(paste0(dir4, "/SIR_boot_list_initI_neg1.2.RData"))
load(paste0(dir4, "/SIR_boot_list_initI_neg0.4.RData"))
load(paste0(dir4, "/SIR_boot_list_initI0.6.RData"))
load(paste0(dir4, "/SIR_boot_list_initI1.6.RData"))


SIR_boot_df <- do.call("rbind.data.frame", SIR_boot_list_initI_neg3) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI_neg1.2)) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI_neg0.4)) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI0.6)) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI1.6)) %>%
  mutate(parameter = ifelse(parameter == "beta_ld1", "NPI 1", "NPI 2"), 
         true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         CI_covers = ifelse(between(true_value, CI_LL1, CI_UL1), 1, 0), 
         bias = true_value - mean_est1, 
         rel_bias = abs(true_value - mean_est1)/abs(true_value)*100) %>%
  group_by(parameter, initI) %>%
  mutate(perc_CI_covers = sum(CI_covers), 
         mean_bias = mean(bias), 
         mean_rel_bias = mean(rel_bias)) %>%
  ungroup()


metric_df_SIR <- SIR_boot_df %>%
  dplyr::select(parameter, initI, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(perc_bootstrap_CI_covers = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  mutate(coverage = NA) %>%
  pivot_longer(cols = c(perc_bootstrap_CI_covers, coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(model = "Mechanistic") %>%
  pivot_wider(names_from = c(initI, model), values_from = value) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)",
                            metric == "perc_bootstrap_CI_covers" ~ "bootstrap CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))



# Combine results into table ----------------------------------------------

metrics_df_reg_late <- res_boot_SIR_late_initI_neg3_df %>% 
  ungroup() %>%
  mutate(initI = -3) %>%
  bind_rows(res_boot_SIR_late_initI_neg1.2_df %>% ungroup() %>% mutate(initI = -1.2)) %>%
  bind_rows(res_boot_SIR_late_initI_neg0.4_df %>% ungroup() %>% mutate(initI = -0.4)) %>%
  bind_rows(res_boot_SIR_late_initI0.6_df %>% ungroup() %>% mutate(initI = 0.6)) %>%
  bind_rows(res_boot_SIR_late_initI1.6_df %>% ungroup() %>% mutate(initI = 1.6)) %>%
  dplyr::select(parameter, median, CI_LL, CI_UL, initI, method, rep) %>%
  unique() %>%
  rename(value = median)  %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  group_by(parameter, initI, method) %>%
  summarize(perc_CI_covers = sum(CI_covers), 
            mean_bias = mean(bias), 
            mean_rel_bias = mean(rel_bias)) %>%
  group_by(parameter, initI) %>%
  mutate(perc_bootstrap_CI_covers = perc_CI_covers[method == "quantile 7d"]) %>%
  filter(method == "simple regression 7d") %>%
  pivot_longer(cols = c(mean_bias, mean_rel_bias, perc_CI_covers, perc_bootstrap_CI_covers), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = c(initI, method), values_from = value) %>%
  mutate(metric = factor(metric, levels = c("mean_bias", "mean_rel_bias", "perc_CI_covers", 
                                             "perc_bootstrap_CI_covers"))) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "mean_bias" ~ "absolute bias", 
                            metric == "perc_CI_covers" ~ "CI coverage (%)",
                            metric == "perc_bootstrap_CI_covers" ~ "bootstrap CI coverage (%)", 
                            metric == "mean_rel_bias" ~ "relative bias (%)"))

# print metrics table
metrics_df_reg_late %>% 
  ungroup() %>%
  left_join(metric_df_SIR, by = c("parameter", "metric")) %>%
  dplyr::select(-parameter) %>%
  dplyr::select(`-3_simple regression 7d`, `-3_Mechanistic`, `-1.2_simple regression 7d`, `-1.2_Mechanistic`, 
                `-0.4_simple regression 7d`, `-0.4_Mechanistic`, `0.6_simple regression 7d`, `0.6_Mechanistic`,
                `1.6_simple regression 7d`, `1.6_Mechanistic`) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 4) %>%
  pack_rows("NPI 2", 5, 8)

