library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(RColorBrewer)
library(ggridges)
library(kableExtra)
library(ggside)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir10 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once10"


# load data and calculate Rt ----------------------------------------------
popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

# SEIRAHD-created data
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


# ABM-created data
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


# Fits EpiEstim + reg -----------------------------------------------------

## Rt estimation -----------------------------------------------------------

cl <- makeCluster(12)
registerDoParallel(cl)

Rt_comp_res3 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new3, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "SEIRAHD", 
                            Rt_ref_df = true_Rt_df_Simulx_new3,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)

Rt_comp_res_ABM10_h <- Rt_calc_fun(Inc_df = data_ABM_hybrid_cov10, id_col_name = "dept_id", time_col_name = "day", 
                                   Inc_col_name = "IncI", model_name = "ABM hybrid", 
                                   Rt_ref_df = true_Rt_df_ABM_hybrid_10,
                                   meansi = 7.8, stdsi = 4.4, Rt_prior = 1)

Rt_comp_res_ABM10_rm <- Rt_calc_fun(Inc_df = data_ABM_rm_cov10, id_col_name = "dept_id", time_col_name = "day", 
                                    Inc_col_name = "IncI", model_name = "ABM random mixing", 
                                    Rt_ref_df = true_Rt_df_ABM_rm_10,
                                    meansi = 7.8, stdsi = 4.4, Rt_prior = 1)
stopCluster(cl)

first_Rt_day <- min(Rt_comp_SEIRAHD$Rt_comp$day) + 3
seven_day_seq <- seq(first_Rt_day, max(Rt_comp_SEIRAHD$Rt_comp$day), 7)

## Regressions -------------------------------------------------------------
# SEIRAHD dataset
reg_data3 <- Rt_comp_res3$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  ungroup() %>%
  filter(day %in% seven_day_seq) 

reg_res3 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data3)

fitted_vals3 <- exp(fitted(reg_res3))
fitted_df3 <- reg_data3 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals3)

dataset1_Simulx_new3_fitted <- dataset1_Simulx_new3 %>%
  left_join(fitted_df3 %>% rename(Rt_EE = Rt), by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


# ABM multi-layer dataset
reg_data_ABM10_h <- Rt_comp_res_ABM10_h$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  filter(day %in% seven_day_seq)


reg_res_ABM10_h <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM10_h)

fitted_vals_ABM10_h <- exp(fitted(reg_res_ABM10_h))
fitted_df_ABM10_h <- reg_data_ABM10_h %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals_ABM10_h)

data_ABM_hybrid_cov10_fitted <- data_ABM_hybrid_cov10 %>%
  left_join(fitted_df_ABM10_h %>% rename(Rt_EE = Rt), by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


# ABM random mixing dataset
reg_data_ABM10_rm <- Rt_comp_res_ABM10_rm$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  filter(day %in% seven_day_seq)

reg_res_ABM10_rm <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_ABM10_rm)

fitted_vals_ABM10_rm <- exp(fitted(reg_res_ABM10_rm))
fitted_df_ABM10_rm <- reg_data_ABM10_rm %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals_ABM10_rm)

data_ABM_rm_cov10_fitted <- data_ABM_rm_cov10 %>%
  left_join(fitted_df_ABM10_rm %>% rename(Rt_EE = Rt), by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


## plots fits ---------------------------------------------------
rect_cols <- sequential_hcl(5, palette = "BluYl")

ggplot(dataset1_Simulx_new3_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
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
  labs(y = expression(R[t]), x = "Day", col = expression(R[t]~ "type")) +
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

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx3_9panels.jpeg", 
       dpi = 400, width = 10, height = 7.5)

# ABM multi-layer
ggplot(data_ABM_hybrid_cov10_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
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
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_color_manual(values = c("blue", "deeppink", "black")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(y = expression(R[t]), x = "Day", col = expression(R[t]~ "type")) +
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

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_ABM_h_9panels.jpeg", 
       dpi = 400, width = 10, height = 7.5)


# ABM random mixing
ggplot(data_ABM_rm_cov10_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
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
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_color_manual(values = c("blue", "deeppink", "black")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(y = expression(R[t]), x = "Day", col = expression(R[t]~ "type")) +
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

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_ABM_rm_9panels.jpeg", 
       dpi = 400, width = 10, height = 7.5)
