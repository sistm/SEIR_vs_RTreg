library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(kableExtra)
library(patchwork)
library(ggside)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
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


true_Rt_df_ABM_rm_10 <- data_ABM_rm_cov10 %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)
true_Rt_df_ABM_hybrid_10 <- data_ABM_hybrid_cov10 %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)

set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)


# known Rt regressions ------------------------------------------------

# run regressions
reg_Rt_known_SEIRAHD <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new3)

reg_Rt_known_ABM_rm10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov10)
reg_Rt_known_ABM_hybrid10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov10)

reg_models_list <- list(reg_Rt_known_SEIRAHD, reg_Rt_known_ABM_rm10, reg_Rt_known_ABM_hybrid10)
model_names <- c("SEIRAHD", "random mixing ABM", "multi-layer ABM")


reg_Rt_known_comp_list <- list()

# summarize results
for(i in 1:length(reg_models_list)){
  coefs_Rt_reg <- coefficients(reg_models_list[[i]])$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_models_list[[i]], method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  reg_Rt_known_comp_list[[i]] <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = model_names[i])
}

reg_Rt_known_comp <- do.call("rbind.data.frame", reg_Rt_known_comp_list) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8))


# fits
data_list <- list(dataset1_Simulx_new3, data_ABM_rm_cov10, data_ABM_hybrid_cov10)
fit_data_list <- list()

for(i in 1:length(reg_models_list)){
  
  fitted <- exp(fitted(reg_models_list[[i]]))
  
  fit_data_list[[i]] <- data_list[[i]] %>%
    mutate(fitted_Rt = fitted, 
           model = model_names[i])
}

reg_fits <- do.call("bind_rows", fit_data_list) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         IncI_unscaled = ifelse(is.na(IncI_unscaled), IncI, IncI_unscaled))


selected_depts <- c(14, 31, 42)

ggplot(reg_fits %>% filter(dept_id %in% selected_depts & grepl("ABM", model)), 
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



plot_Rt_known_SEIRAHD <- ggplot(reg_fits %>% filter(dept_id %in% selected_depts & !grepl("ABM", model)), 
                                aes(x = day, y = Rt)) + 
  geom_line(aes(col = "Real Rt"), linewidth = 0.8) + 
  geom_line(aes(y = fitted_Rt, col = "Fitted Rt"), linewidth = 0.8) + 
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



#### on the 100 datasets ####

reg_Rt_known_SEIRAHD_list <- list()
for(j in 1:100){
  indparams_SEIRAHD <- read.table(paste0(dir2, "/ind_params/ind_2params_new3_", j, ".txt"), 
                                  header = TRUE, sep = " ")
  
  data_SEIRAHD <- data_list_Simulx_new3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(indparams_SEIRAHD, by = "id") %>%
    mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = IncI*popsize/10^4, 
           IncH_unscaled = IncH*popsize/10^4,
           lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  reg_Rt_known <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_SEIRAHD)
  
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
    mutate(model = "SEIRAHD", rep = j)
  
  reg_Rt_known_SEIRAHD_list[[j]] <- reg_Rt_known_comp
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


summary_regs_Rt_known <- do.call("rbind.data.frame", reg_Rt_known_SEIRAHD_list) %>%
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



# Fits EpiEstim + reg -----------------------------------------------------


## Rt estimation -----------------------------------------------------------

cl <- makeCluster(12)
registerDoParallel(cl)

Rt_comp_SEIRAHD <- Rt_calc_fun(Inc_df = dataset1_Simulx_new3, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "SEIRAHD", 
                            Rt_ref_df = true_Rt_df_Simulx_new3,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)

Rt_comp_res_ABM10_h <- Rt_calc_fun(Inc_df = data_ABM_hybrid_cov10, id_col_name = "dept_id", time_col_name = "day", 
                                   Inc_col_name = "IncI", model_name = "multi-layer ABM", 
                                   Rt_ref_df = true_Rt_df_ABM_hybrid_10,
                                   meansi = 7.8, stdsi = 4.4, Rt_prior = 1)

Rt_comp_res_ABM10_rm <- Rt_calc_fun(Inc_df = data_ABM_rm_cov10, id_col_name = "dept_id", time_col_name = "day", 
                                    Inc_col_name = "IncI", model_name = "random mixing ABM", 
                                    Rt_ref_df = true_Rt_df_ABM_rm_10,
                                    meansi = 8.45, stdsi = 5.45, Rt_prior = 1)
stopCluster(cl)

first_Rt_day <- min(Rt_comp_SEIRAHD$Rt_comp$day) + 3
seven_day_seq <- seq(first_Rt_day, max(Rt_comp_SEIRAHD$Rt_comp$day), 7)


## Regressions -------------------------------------------------------------
# all with NPIs lagged by 5 days

# summary SEIRAHD dataset
reg_data_SEIRAHD <- Rt_comp_SEIRAHD$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  group_by(dept_id) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  ungroup() %>%
  filter(day %in% seven_day_seq) 

reg_res_SEIRAHD <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_SEIRAHD)

fitted_vals_SEIRAHD <- exp(fitted(reg_res_SEIRAHD))
fitted_df_SEIRAHD <- reg_data_SEIRAHD %>%
  filter(day %in% seven_day_seq) %>%
  mutate(Rt_fitted = fitted_vals_SEIRAHD) %>%
  select(day, dept_id, Rt_fitted)

dataset1_SEIRAHD_fitted <- dataset1_SEIRAHD %>%
  select(dept_id, day, IncI_unscaled) %>%
  left_join(Rt_comp_SEIRAHD$Rt_comp, by = c("dept_id", "day")) %>%
  left_join(fitted_df_SEIRAHD, by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


# summary ABM 10
reg_data_ABM10_h <- Rt_comp_res_ABM10_h$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  group_by(dept_id) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  ungroup() %>%
  filter(day %in% seven_day_seq) 

reg_res_ABM10_h <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM10_h)

fitted_vals_ABM10_h <- exp(fitted(reg_res_ABM10_h))
fitted_df_ABM10_h <- reg_data_ABM10_h %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals_ABM10_h) %>%
  select(day, dept_id, Rt_fitted)

data_ABM_hybrid_cov10_fitted <- data_ABM_hybrid_cov10 %>%
  select(day, dept_id, IncI) %>%
  left_join(Rt_comp_res_ABM10_h$Rt_comp, by = c("dept_id", "day")) %>%
  left_join(fitted_df_ABM10_h, by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))



reg_data_ABM10_rm <- Rt_comp_res_ABM10_rm$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  group_by(dept_id) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  ungroup() %>%
  filter(day %in% seven_day_seq) 

reg_res_ABM10_rm <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM10_rm)

fitted_vals_ABM10_rm <- exp(fitted(reg_res_ABM10_rm))
fitted_df_ABM10_rm <- reg_data_ABM10_rm %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals_ABM10_rm) %>%
  select(dept_id, day, Rt_fitted)

data_ABM_rm_cov10_fitted <- data_ABM_rm_cov10  %>%
  select(day, dept_id, IncI) %>%
  left_join(Rt_comp_res_ABM10_rm$Rt_comp, by = c("dept_id", "day")) %>%
  left_join(fitted_df_ABM10_rm, by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


## plots fits ---------------------------------------------------

rect_cols <- sequential_hcl(5, palette = "BluYl")

plot_fits_SEIRAHD <- ggplot(dataset1_SEIRAHD_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
       aes(x = day)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(y = Rt_real, col = "True"), linewidth = 0.8) +
  geom_line(aes(y = Rt, col = "EpiEstim"), linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_color_manual(values = c("blue", "deeppink", "black")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(#title = "Simulx infections 3, short NPI-free period, no days cut off, no NPI lag", 
    y = expression(R[t]), x = "Day", col = expression(R[t]~ "type")) +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.title = element_text(family = "serif", size = 14, hjust = 0),
        legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

plot_fits_SEIRAHD
ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx3_9panels.jpeg", 
       dpi = 400, width = 10, height = 7.5)

plot_Rt_known_SEIRAHD + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  plot_fits_SEIRAHD +
  plot_layout(nrow = 2, heights = c(1, 3)) + 
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(family = "serif", face = "bold", size = 18, hjust = 0, vjust = 0))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/Rt_known_and_reg_fits.jpeg", 
       dpi = 400, width = 10, height = 10)


ggplot(data_ABM_hybrid_cov10_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
       aes(x = day)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(y = Rt_real, col = "Real Rt"), linewidth = 0.8) +
  geom_line(aes(y = Rt, col = "EpiEstim Rt"), linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit Rt")) +
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_color_manual(values = c("blue", "black", "deeppink")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(#title = "Simulx infections 3, short NPI-free period, no days cut off, no NPI lag", 
    y = expression(R[t]), x = "Day", col = expression(R[t]~ "type")) +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.title = element_text(family = "serif", size = 14, hjust = 0),
        legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

ggplot(data_ABM_rm_cov10_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
       aes(x = day)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(y = Rt_real, col = "Real Rt"), linewidth = 0.8) +
  geom_line(aes(y = Rt, col = "EpiEstim Rt"), linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit Rt")) +
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_color_manual(values = c("blue", "black", "deeppink")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(#title = "Simulx infections 3, short NPI-free period, no days cut off, no NPI lag", 
    y = expression(R[t]), x = "Day", col = expression(R[t]~ "type")) +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.title = element_text(family = "serif", size = 14, hjust = 0),
        legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


