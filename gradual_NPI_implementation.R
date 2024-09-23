library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(ggside)
library(RColorBrewer)
library(kableExtra)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")
dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"


# load necessary data
load("sim_res3_slow_NPI.RData")
load("sim_res_1week_slow_list.RData")
load("sim_res_2week_slow_list.RData")
load("sim_res_1week_slow_late_list.RData")
load("sim_res_2week_slow_late_list.RData")


sim_res_slow_names <- c("1week_slow", "2week_slow", "1week_slow_late", "2week_slow_late")
names2 <- c("1 week gradient early", "2 week gradient early", 
            "1 week gradient late", "2 week gradient late")

popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

ind_params3_lowb1 <- read.table("ind_params3_lowb1_gradual.txt", header = TRUE, sep = " ")
ind_params3 <- read.table("ind_params3_gradual.txt", header = TRUE, sep = " ")

NPI_df_1week <- read.csv("ld1_reg_1week_slow.csv")
NPI_df_2week <- read.csv("ld1_reg_2week_slow.csv")
NPI_df_1week_late <- read.csv("ld1_reg_1week_slow_late.csv")
NPI_df_2week_late <- read.csv("ld1_reg_2week_slow_late.csv")
  
  

# Bias table --------------------------------------------------------------

cl <- makeCluster(12)
registerDoParallel(cl)

#### 1 week slow ####
reg_res_list_1week_slow_7d <- vector(mode = "list")
for(j in 1:100){
  
  reg_data <- sim_res_1week_slow_list[[j]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI_ME*popsize/10^4)
  
  Rt_res <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = "IncI_unscaled", meansi = 10.1, stdsi = 8,75, 
                         Rt_prior = 1, Rt_sd_prior = 2)
  
  reg_res <- Rt_reg_only_fun_7d(data_for_est = Rt_res %>% 
                                  left_join(NPI_df_1week %>% select(id, time, lockdown1, BG1), 
                                              by = c("dept_id" = "id", "day" = "time")), 
                                model_name = "1week slow", lag_NPIs = TRUE, lag_days = 5, fits = FALSE)
  
  
  reg_res_list_1week_slow_7d[[j]] <- reg_res
  
}


#### 2 weeks slow ####
reg_res_list_2week_slow_7d <- vector(mode = "list")
for(j in 1:100){
  
  reg_data <- sim_res_2week_slow_list[[j]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI_ME*popsize/10^4)
  
  Rt_res <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = "IncI_unscaled", meansi = 10.1, stdsi = 8,75, 
                         Rt_prior = 1, Rt_sd_prior = 2)
  
  reg_res <- Rt_reg_only_fun_7d(data_for_est = Rt_res %>% 
                                  left_join(NPI_df_2week %>% select(id, time, lockdown1, BG1), 
                                            by = c("dept_id" = "id", "day" = "time")), 
                                model_name = "2week slow", lag_NPIs = TRUE, lag_days = 5, fits = FALSE)
  
  
  reg_res_list_2week_slow_7d[[j]] <- reg_res
  
}


#### 1 week slow late ####
reg_res_list_1week_slow_late_7d <- vector(mode = "list")
for(j in 1:100){
  
  reg_data <- sim_res_1week_slow_late_list[[j]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI_ME*popsize/10^4)
  
  Rt_res <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = "IncI_unscaled", meansi = 10.1, stdsi = 8,75, 
                         Rt_prior = 1, Rt_sd_prior = 2)
  
  reg_res <- Rt_reg_only_fun_7d(data_for_est = Rt_res %>% 
                                  left_join(NPI_df_1week_late %>% select(id, time, lockdown1, BG1), 
                                            by = c("dept_id" = "id", "day" = "time")), 
                                model_name = "1week slow late", lag_NPIs = TRUE, lag_days = 5, fits = FALSE)
  
  
  reg_res_list_1week_slow_late_7d[[j]] <- reg_res
  
}


#### 2 weeks slow late ####
reg_res_list_2week_slow_late_7d <- vector(mode = "list")
for(j in 1:100){
  
  reg_data <- sim_res_2week_slow_late_list[[j]] %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    mutate(IncI_unscaled = IncI_ME*popsize/10^4)
  
  Rt_res <- Rt_calc_fun2(Inc_df = reg_data, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = "IncI_unscaled", meansi = 10.1, stdsi = 8,75, 
                         Rt_prior = 1, Rt_sd_prior = 2)
  
  reg_res <- Rt_reg_only_fun_7d(data_for_est = Rt_res %>% 
                                  left_join(NPI_df_2week_late %>% select(id, time, lockdown1, BG1), 
                                            by = c("dept_id" = "id", "day" = "time")), 
                                model_name = "2week slow late", lag_NPIs = TRUE, lag_days = 5, fits = FALSE)
  
  
  reg_res_list_2week_slow_late_7d[[j]] <- reg_res
  
}

stopCluster(cl)

reg_res_1week_slow_7d_df <- do.call("rbind.data.frame", reg_res_list_1week_slow_7d)
reg_res_2week_slow_7d_df <- do.call("rbind.data.frame", reg_res_list_2week_slow_7d)
reg_res_1week_slow_7d_late_df <- do.call("rbind.data.frame", reg_res_list_1week_slow_late_7d)
reg_res_2week_slow_7d_late_df <- do.call("rbind.data.frame", reg_res_list_2week_slow_late_7d)

save(reg_res_1week_slow_7d_df, file = "reg_res_1week_slow_7d_df.RData")
save(reg_res_2week_slow_7d_df, file = "reg_res_2week_slow_7d_df.RData")
save(reg_res_1week_slow_7d_late_df, file = "reg_res_1week_slow_7d_late_df.RData")
save(reg_res_2week_slow_7d_late_df, file = "reg_res_2week_slow_7d_late_df.RData")


#### table ####

metric_df_reg <- reg_res_1week_slow_7d_df %>%
  bind_rows(reg_res_2week_slow_7d_df, reg_res_1week_slow_7d_late_df, reg_res_2week_slow_7d_late_df) %>%
  dplyr::select(parameter, value, CI_LL, CI_UL, model) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  group_by(parameter, model) %>%
  summarize(perc_CI_covers = sum(CI_covers), 
            mean_bias = mean(bias), 
            mean_rel_bias = mean(rel_bias)) %>%
  pivot_longer(cols = c(perc_CI_covers, mean_bias, mean_rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("mean_bias", "mean_rel_bias", "perc_CI_covers"))) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "mean_bias" ~ "absolute bias", 
                            metric == "perc_CI_covers" ~ "CI coverage (%)", 
                            metric == "mean_rel_bias" ~ "relative bias (%)")) 

# print metrics table
metric_df_reg %>% 
  ungroup() %>%
  dplyr::select(-parameter) %>%
  kable(digits = 2, format = "latex") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)



# Plots -------------------------------------------------------

cl <- makeCluster(10)
registerDoParallel(cl)

res_EE_slow_NPI_list3 <- list()
res_reg_slow_NPI_list3 <- list()
fits_reg_slow_NPI_list3 <- list()

for(j in 1:4){

  indprms <- ind_params3
  
  # assemble data
  data_SEIR_monolix <- read.table(paste0("data_sim_SEIR_Simulx3_2params_", sim_res_slow_names[j], ".txt"), 
                                  sep = ",", header = TRUE) %>%
    select(dept_id, day, lockdown1, BG1)
  
  data_for_est <- sim_res3_slow_NPI[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(indprms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4)) %>% 
    rename(dept_id = id, day = time) %>%
    left_join(data_SEIR_monolix, by = c("dept_id", "day"))
  
  # EpiEstim only
  res_EE <- Rt_calc_fun2(Inc_df = data_for_est, id_col_name = "dept_id", time_col_name = "day", 
                         Inc_col_name = "IncI_unscaled", 
                         meansi = 10.1, stdsi = 8.75, Rt_prior = 1) %>%
    mutate(model = sim_res_slow_names[j]) %>%
    left_join(data_for_est %>% select(dept_id, day, lockdown1, BG1), by = c("dept_id", "day"))
  
  res_EE_slow_NPI_list3[[j]] <- res_EE %>%
    mutate(model = names2[j])
  
  
  first_Rt_day <- min(res_EE$day[!is.na(res_EE$Rt)]) + 3
  seven_day_seq <- seq(first_Rt_day, max(res_EE$day), 7)
  
  
  reg_data <- res_EE %>%
    group_by(dept_id) %>%
    mutate(across(lockdown1:BG1, function(x) lag(x, 5, default = 0))) %>%
    ungroup() %>%
    filter(day %in% seven_day_seq)
  
  
  # regression
  reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
  
  # coefficients
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
                              labels = c("Lockdown 1", "Barrier gestures")), 
           model = sim_res_slow_names[j])
  
  res_reg_slow_NPI_list3[[j]] <- coefs_Rt_reg %>%
    mutate(model = names2[j])
  
  
  # regression fits
  fitted_vals <- reg_data %>% 
    mutate(Rt_fitted = exp(fitted(reg_res))) %>%
    select(dept_id, day, Rt_fitted)
  
  fitted_df <- res_EE %>%
    left_join(fitted_vals, by = c("dept_id", "day")) %>%
    left_join(data_for_est %>% select(dept_id, day, Rt_SEIRAHD, IncI_unscaled)) %>%
    select(dept_id, day, IncI_unscaled, Rt_SEIRAHD, Rt, CI_LL, CI_UL, Rt_fitted, model, lockdown1, BG1) %>%
    rename(Rt_true = Rt_SEIRAHD)
  
  fits_reg_slow_NPI_list3[[j]] <- fitted_df %>%
    mutate(model = names2[j])
}

stopCluster(cl)


res_EE_slow_NPI_df3 <- do.call("rbind.data.frame", res_EE_slow_NPI_list3)  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94))) 

fits_reg_slow_NPI_df3 <- bind_rows(fits_reg_slow_NPI_list3) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94))) 



selected_depts <- c(2, 51, 76) 

ggplot(fits_reg_slow_NPI_df3 %>% filter(dept_id %in% selected_depts), 
       aes(x = day)) +
  geom_line(aes(y = Rt, col = "Epi Estim"), linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit")) + 
  geom_line(aes(y = Rt_true, col = "True"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_grid(rows = vars(dept_id2), cols = vars(model)) +
  labs(col = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
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

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/fits_gradient_NPIs.jpeg", 
       dpi = 400, width = 12, height = 8.5)
