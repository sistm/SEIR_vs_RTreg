library(tidyverse)
library(kableExtra)


source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")
dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_new3_Simulx_SEIR"
dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_new3_Simulx_SEIRAHD"


# Load data ---------------------------------------------------------------

load(paste0(dir, "/res_boot_df_Simulx3_IncI_7d.RData"))
load(paste0(dir, "/res_boot_df_Simulx3_IncH_7d.RData"))

load(paste0(dir3, "/SEIR_boot_list_2params_new3_Simulx.RData"))
load(paste0(dir4, "/SEIRAHD_boot_list_2params_new3_Simulx.RData"))



# Bring data into right format and assemble -------------------------------

SEIR_boot_df_comp3 <- bootstrap_summary(SEIR_boot_list_all3, true_val_NPI2 = -0.8) 
SEIRAHD_boot_df_comp3 <- bootstrap_summary(SEIRAHD_boot_list_all3, true_val_NPI2 = -0.8) 

metric_df_mech <- SEIR_boot_df_comp3 %>%
  mutate(model = "SEIR model") %>%
  bind_rows(SEIRAHD_boot_df_comp3 %>% mutate(model = "SEIRAHD model")) %>%
  rename(mean_est = mean_est2, CI_LL = CI_LL2, CI_UL = CI_UL2) %>%
  dplyr::select(parameter, model, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(perc_bootstrap_CI_covers = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  mutate(coverage = NA) %>%
  pivot_longer(cols = c(perc_bootstrap_CI_covers, coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)",
                            metric == "perc_bootstrap_CI_covers" ~ "bootstrap CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))


metric_df_reg <- res_boot_df_Simulx3_IncI_7d %>%
  ungroup() %>%
  mutate(model = "regression IncI") %>%
  bind_rows(res_boot_df_Simulx3_IncH_7d %>% mutate(model = "regression IncH") %>% ungroup()) %>%
  filter(method != "random 7d") %>%
  dplyr::select(parameter, median, CI_LL, CI_UL, model, method, rep) %>%
  unique() %>%
  rename(value = median) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  group_by(parameter, model, method) %>%
  summarize(perc_CI_covers = sum(CI_covers), 
            mean_bias = mean(bias), 
            mean_rel_bias = mean(rel_bias)) %>%
  mutate(perc_bootstrap_CI_covers = perc_CI_covers[method == "quantile 7d"]) %>%
  filter(method == "simple regression 7d") %>%
  pivot_longer(cols = c(perc_CI_covers, perc_bootstrap_CI_covers, mean_bias, mean_rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("mean_bias", "mean_rel_bias", "perc_CI_covers", 
                                            "perc_bootstrap_CI_covers"))) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "mean_bias" ~ "absolute bias", 
                            metric == "perc_CI_covers" ~ "CI coverage (%)", 
                            metric == "perc_bootstrap_CI_covers" ~ "bootstrap CI coverage (%)", 
                            metric == "mean_rel_bias" ~ "relative bias (%)")) %>%
  dplyr::select(-method)

# print metrics table
metric_df_reg %>% 
  ungroup() %>%
  left_join(metric_df_mech, by = c("parameter", "metric")) %>%
  dplyr::select(-parameter) %>%
  select(metric, `SEIR model`, `SEIRAHD model`, `regression IncI`, `regression IncH`) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 4) %>%
  pack_rows("NPI 2", 5, 8)

