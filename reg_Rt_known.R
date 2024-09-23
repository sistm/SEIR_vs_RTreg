library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(RColorBrewer)
library(kableExtra)
library(ggside)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir10 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once10"


# load data and calculate Rt ----------------------------------------------
popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

indparams_Simulx3_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_new3_1.txt"), 
                                  header = TRUE, sep = " ")

data_list_Simulx_new3 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new3_list.RData", sep = "/"))


dataset1_Simulx_new3 <- data_list_Simulx_new3[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx3_1, by = "id") %>%
  mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4,
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
         BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new3 <- dataset1_Simulx_new3 %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)

data_ABM_rm_cov10 <- read.csv(paste0(dir10, "/data_covasim_rm10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)
data_ABM_hybrid_cov10 <- read.csv(paste0(dir10, "/data_covasim_hybrid10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)


# Known Rt regressions ----------------------------------------------------

#### fits figures ####
reg_Rt_known_Simulx3 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new3)
reg_Rt_known_ABM_rm10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov10)
reg_Rt_known_ABM_hybrid10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov10)

data_list <- list(dataset1_Simulx_new3, data_ABM_rm_cov10, data_ABM_hybrid_cov10)
reg_models_list <- list(reg_Rt_known_Simulx3, reg_Rt_known_ABM_rm10, reg_Rt_known_ABM_hybrid10)
model_names <- c("SEIRAHD", "random mixing ABM", "multi-layer ABM")

fit_data_list <- list()
for(i in 1:length(reg_models_list)){
  
  fitted <- exp(fitted(reg_models_list[[i]]))
  
  fit_data_list[[i]] <- data_list[[i]] %>%
    mutate(fitted_Rt = fitted, 
           model = model_names[i])
}

reg_fits_all <- do.call("bind_rows", fit_data_list) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         IncI_unscaled = ifelse(is.na(IncI_unscaled), IncI, IncI_unscaled))


selected_depts <- c(14, 31, 42)

ggplot(reg_fits_all %>% filter(dept_id %in% selected_depts & grepl("ABM", model)), 
       aes(x = day, y = Rt)) + 
  geom_line(aes(col = "Real Rt")) + 
  geom_line(aes(y = fitted_Rt, col = "Fitted Rt")) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y") + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_grid(rows = vars(model), cols = vars(dept_id2)) +
  scale_color_brewer(palette = "Set1", 
                     labels = c(expression("Fitted" ~R[t]), expression("Real" ~R[t]))) + 
  theme_bw() + 
  labs(y = expression(R[t]), col = "", x = "Day") +
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 11))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Rt_known_ABMs.jpeg", dpi = 400, width = 10, height = 6)



ggplot(reg_fits_all %>% filter(dept_id %in% selected_depts & !grepl("ABM", model)), 
       aes(x = day, y = Rt)) + 
  geom_line(aes(col = "Real Rt")) + 
  geom_line(aes(y = fitted_Rt, col = "Fitted Rt")) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y") + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  scale_color_brewer(palette = "Set1", 
                     labels = c(expression("Fitted" ~R[t]), expression("Real" ~R[t]))) + 
  theme_bw() + 
  labs(y = expression(R[t]), col = "", x = "Day") +
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 11))

save(plot_Rt_known_SEIRAHD, file = "plot_Rt_known_SEIRAHD.RData")



# Bias table --------------------------------------------------------------

#### regressions ####
reg_Rt_known_Simulx3_list <- list()
for(j in 1:100){
  indparams_Simulx3 <- read.table(paste0(dir2, "/ind_params/ind_2params_new3_", j, ".txt"), 
                                  header = TRUE, sep = " ")
  
  data_Simulx_new3 <- data_list_Simulx_new3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(indparams_Simulx3, by = "id") %>%
    mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = IncI*popsize/10^4, 
           IncH_unscaled = IncH*popsize/10^4,
           lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  reg_Rt_known <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_Simulx_new3)
  
  coefs_Rt_reg <- coefficients(reg_Rt_known)$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_Rt_known, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  reg_Rt_known_comp <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = "Simulx3", rep = j)
  
  reg_Rt_known_Simulx3_list[[j]] <- reg_Rt_known_comp
}


reg_Rt_known_ABM10_rm_list <- list()
reg_Rt_known_ABM10_h_list <- list()

for(j in 1:100){
  data_ABM_rm_cov10 <- read.csv(paste0(dir10, "/data_covasim_rm10_Rt_", j, ".csv")) %>%
    filter(day > 29) %>%
    mutate(day = day - 29) %>%
    filter(Rt > 0)
  data_ABM_hybrid_cov10 <- read.csv(paste0(dir10, "/data_covasim_hybrid10_Rt_", j, ".csv")) %>%
    filter(day > 29) %>%
    mutate(day = day - 29) %>%
    filter(Rt > 0)
  
  reg_Rt_known_rm <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov10)
  reg_Rt_known_h <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov10)
  
  coefs_Rt_reg_rm <- coefficients(reg_Rt_known_rm)$dept_id
  coefs_Rt_reg_h <- coefficients(reg_Rt_known_h)$dept_id
  
  confint_Rt_reg_rm <- data.frame(confint(reg_Rt_known_rm, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg_rm) <- c("CI_LL", "CI_UL")
  
  confint_Rt_reg_h <- data.frame(confint(reg_Rt_known_h, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg_h) <- c("CI_LL", "CI_UL")
  
  reg_Rt_known_comp_rm <- coefs_Rt_reg_rm %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg_rm[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = "ABM10 rm", rep = j) 
  
  reg_Rt_known_comp_h <- coefs_Rt_reg_h %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg_h[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = "ABM10 hybrid", rep = j)
  
  reg_Rt_known_ABM10_rm_list[[j]] <- reg_Rt_known_comp_rm
  reg_Rt_known_ABM10_h_list[[j]] <- reg_Rt_known_comp_h
}


#### results summary ####
summary_regs_Rt_known <- do.call("rbind.data.frame", reg_Rt_known_Simulx3_list) %>%
  bind_rows(do.call("rbind.data.frame", reg_Rt_known_ABM10_rm_list)) %>%
  bind_rows(do.call("rbind.data.frame", reg_Rt_known_ABM10_h_list)) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         unique_sims = length(unique(rep)), 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  group_by(parameter, model) %>%
  mutate(perc_CI_covers = sum(CI_covers)/unique_sims*100, 
         mean_bias = mean(bias), 
         mean_rel_bias = mean(rel_bias))


metrics_regs_Rt_known <- summary_regs_Rt_known %>%
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
metrics_regs_Rt_known %>% 
  ungroup() %>%
  dplyr::select(-parameter) %>%
  kable(digits = 2, format = "html"#, table.attr = "style='width:62%;'"
  ) %>%
  kable_styling(bootstrap_options = "striped") %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)




