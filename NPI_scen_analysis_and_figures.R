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
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")
dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"


# load necessary data
load("sim_res_params3_ld1_start.RData")
load("sim_res_params3_ld1_start_high.RData")
load("sim_res_params3_ld1_strength.RData")

popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

ind_params3_lowb1 <- read.table("ind_params3_lowb1.txt", header = TRUE, sep = " ")
ind_params3 <- read.table("ind_params3.txt", header = TRUE, sep = " ")


ld1_start <- c(20, 30, 40, 50, 60)
ld1_strength <- seq(0.5, 2, 0.1)


set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)



# Simulx 3 ----------------------------------------------------------------


### EpiEstim ####
cl <- makeCluster(10)
registerDoParallel(cl)

res_EE_ld1_start_list3 <- list()

for(j in 1:length(sim_res_ld1_start3)){
  data_x <- sim_res_ld1_start3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params3_lowb1, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 50 + ld1_start[j]), 1, 0), 
           BG1 = ifelse(time > 50 + ld1_start[j], 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_EE_ld1_start <- EpiEstim_only_fun(data_for_est = data_x, 
                                        Inc_name = "IncI_unscaled", 
                                        meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_ld1_start_list3[[j]] <- res_EE_ld1_start
}



res_EE_ld1_start_high_list3 <- list()

for(j in 1:length(sim_res_ld1_start3)){
  data_x <- sim_res_ld1_start_high3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params3, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 50 + ld1_start[j]), 1, 0), 
           BG1 = ifelse(time > 50 + ld1_start[j], 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_EE_ld1_start_high <- EpiEstim_only_fun(data_for_est = data_x, 
                                             Inc_name = "IncI_unscaled", 
                                             meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_ld1_start_high_list3[[j]] <- res_EE_ld1_start_high
}


res_EE_ld1_strength_list3 <- list()

for(j in 1:length(sim_res_ld1_strength3)){
  ind_prms <- read.table(paste0("ind_params3_", ld1_strength[j], ".txt"), header = TRUE, sep = " ")
  data_x <- sim_res_ld1_strength3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_prms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, 20, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_EE_ld1_strength <- EpiEstim_only_fun(data_for_est = data_x, 
                                           Inc_name = "IncI_unscaled", 
                                           meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_strength = ld1_strength[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_ld1_strength_list3[[j]] <- res_EE_ld1_strength
}

stopCluster(cl)


### compare and plot results ####
res_EE_ld1_start_df3 <- do.call("rbind.data.frame", res_EE_ld1_start_list3) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))

res_EE_ld1_start_high_df3 <- do.call("rbind.data.frame", res_EE_ld1_start_high_list3) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))

res_EE_ld1_strength_df3 <- do.call("rbind.data.frame", res_EE_ld1_strength_list3) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)), 
         ld1_strength_cat = cut(ld1_strength, 4, 
                                labels = c("very weak", "weak", "moderate", "strong")), 
         ld1_strength = 0 - ld1_strength)


p4 <- ggplot(res_EE_ld1_start_df3 %>% filter(dept_id %in% selected_depts[c(1:3, 5)]),
       aes(x = day, y = Rt, col = as.factor(ld1_start))) + 
  geom_line(aes(linetype = "EpiEstim Rt"), linewidth = 0.8) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_start)), alpha = 0.3, colour = NA) +
  geom_line(aes(y = Rt_SEIRAHD, linetype = "Real Rt"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) + 
  scale_x_continuous(expand = c(0.05, 0.05), breaks = seq(0, 150, 50)) + 
  labs(title = "Low baseline transmission", 
       col = "NPI 1 start day", fill = "NPI 1 start day", x = "Day",
       y = expression(R[t]), linetype = "") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  scale_linetype(labels = c(expression(EpiEstim~R[t]), expression(Real~R[t]))) +
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
         linetype = guide_legend(order = 2)) +
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


p5 <- ggplot(res_EE_ld1_start_high_df3 %>% filter(dept_id %in% selected_depts[c(1:3, 5)]),
       aes(x = day, y = Rt, col = as.factor(ld1_start))) + 
  geom_line(aes(linetype = "EpiEstim Rt"), linewidth = 0.8) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_start)), alpha = 0.3, colour = NA) +
  geom_line(aes(y = Rt_SEIRAHD, linetype = "Real Rt"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) + 
  scale_x_continuous(expand = c(0.05, 0.05), breaks = seq(0, 150, 50)) + 
  labs(title = "High baseline transmission", 
       col = "NPI 1 start day", fill = "NPI 1 start day", x = "Day",
       y = expression(R[t]), linetype = "") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  scale_linetype(labels = c(expression(EpiEstim~R[t]), expression(Real~R[t]))) +
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
         linetype = guide_legend(order = 2)) +
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0),
        legend.title = element_text(family = "serif", size = 13), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


viridis_pal <- sequential_hcl(16, palette = "Viridis")
p6 <- ggplot(res_EE_ld1_strength_df3 %>% filter(dept_id %in% selected_depts[c(1:3, 5)]),
       aes(x = day, y = Rt, col = as.factor(ld1_strength))) + 
  geom_line(aes(y = Rt_SEIRAHD, linetype = "Real Rt"), linewidth = 0.8) + 
  geom_line(aes(linetype = "EpiEstim Rt"), linewidth = 0.8) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_strength)), alpha = 0.3, colour = NA) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  #ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_grid(ld1_strength_cat~dept_id2) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(col = "NPI 1 coefficient", fill = "NPI 1 coefficient", x = "Day",
       y = expression(R[t]), linetype = "") +
  theme_bw() +
  scale_color_manual(values = rev(viridis_pal), breaks = unique(res_EE_ld1_strength_df3$ld1_strength)) +
  scale_fill_manual(values = rev(viridis_pal), breaks = unique(res_EE_ld1_strength_df3$ld1_strength)) + 
  scale_linetype(labels = c(expression(EpiEstim~R[t]), expression(Real~R[t]))) +
  guides(linetype = "none", fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

plot4a <- p4 + p5 + plot_layout(guides = 'collect', tag_level = 'new')
plot4a / p6 + 
  plot_layout(heights = c(1, 1.3)) + 
  plot_annotation(tag_levels = c('A', '1')) & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/Rt_traj_NPI_scen.jpeg", dpi = 300, width = 14, height = 9)


### regressions ####
cl <- makeCluster(10)
registerDoParallel(cl)

res_reg_ld1_start_list3 <- list()

for(j in 1:length(sim_res_ld1_start3)){
  data_x <- sim_res_ld1_start3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params3_lowb1, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 50 + ld1_start[j]), 1, 0), 
           BG1 = ifelse(time > 50 + ld1_start[j], 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_reg_ld1_start <- EpiEstim_reg_fun(data_for_est = data_x, 
                                        Inc_name = "IncI_unscaled", 
                                        rep_num = 1,  lag_NPIs = TRUE, lag_days = 5,
                                        meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) 
  
  res_reg_ld1_start_list3[[j]] <- res_reg_ld1_start
}


res_reg_ld1_start_high_list3 <- list()

for(j in 1:length(sim_res_ld1_start3)){
  data_x <- sim_res_ld1_start_high3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params3, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 50 + ld1_start[j]), 1, 0), 
           BG1 = ifelse(time > 50 + ld1_start[j], 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_reg_ld1_start_high <- EpiEstim_reg_fun(data_for_est = data_x, 
                                             Inc_name = "IncI_unscaled", 
                                             rep_num = 1, lag_NPIs = TRUE, lag_days = 5,
                                             meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) 
  
  res_reg_ld1_start_high_list3[[j]] <- res_reg_ld1_start_high
}


res_reg_ld1_strength_list3 <- list()

for(j in 1:length(sim_res_ld1_strength3)){
  ind_prms <- read.table(paste0("ind_params3_", ld1_strength[j], ".txt"), header = TRUE, sep = " ")
  data_x <- sim_res_ld1_strength3[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_prms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, 20, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_reg_ld1_strength <- EpiEstim_reg_fun(data_for_est = data_x, 
                                           Inc_name = "IncI_unscaled", 
                                           rep_num = 1, lag_NPIs = TRUE, lag_days = 5,
                                           meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_strength = ld1_strength[j]) 
  
  res_reg_ld1_strength_list3[[j]] <- res_reg_ld1_strength
}

stopCluster(cl)


res_reg_ld1_start_df3 <- do.call("rbind.data.frame", res_reg_ld1_start_list3) %>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"), 
         true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100)
res_reg_ld1_start_high_df3 <- do.call("rbind.data.frame", res_reg_ld1_start_high_list3) %>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"), 
         true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100)
res_reg_ld1_strength_df3 <- do.call("rbind.data.frame", res_reg_ld1_strength_list3) %>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"),  
         true_value = ifelse(parameter == "NPI 1", 0-ld1_strength, -0.8), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100, 
         ld1_strength = 0-ld1_strength)


# plots
p1 <- ggplot(res_reg_ld1_start_df3, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_start, 
                                 y = value, col = as.factor(ld1_start))) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 2, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(col = "NPI 1 start day", x = "NPI 1 start day", y = "Coefficient value", 
       title = "Low baseline transmission") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13))

p2 <- ggplot(res_reg_ld1_start_high_df3, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_start, 
                                      y = value, col = as.factor(ld1_start))) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 2, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "High baseline transmission", col = "NPI 1 start day",
       x = "NPI 1 start day", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13))


p3 <- ggplot(res_reg_ld1_strength_df3, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_strength, 
                                      y = value, col = as.factor(ld1_strength))) + 
  geom_pointrange() + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 2, scale = "free_y") +
  geom_segment(aes(y = true_value, yend = true_value, x = ld1_strength-0.05, xend = ld1_strength+0.05), 
               linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(col = "NPI 1 coefficient", x = "True NPI 1 coefficient", y = "NPI 1 coefficient") +
  theme_bw() +
  scale_color_manual(values = rev(viridis_pal), breaks = unique(res_reg_ld1_strength_df3$ld1_strength)) + 
  theme(legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13))

plot1a <- p1 + p2 + plot_layout(guides = 'collect', tag_level = 'new')
plot1a / p3 + 
  plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = c('A', '1')) & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/Rt_traj_NPI_scen_coefs.jpeg", dpi = 300, width = 14, height = 9)



# tables
metric_table_reg_ld1_start3 <- res_reg_ld1_start_df3 %>%
  select(parameter, ld1_start, bias, rel_bias) %>%
  unique() %>%
  pivot_longer(cols = c(bias, rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = ld1_start, values_from = value, names_prefix = "day ") %>%
  mutate(metric = factor(metric, levels = c("bias", "rel_bias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "rel_bias" ~ "relative bias (%)"))

metric_table_reg_ld1_start3 %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 2) %>%
  pack_rows("NPI 2", 3, 4)


metric_table_reg_ld1_start_high3 <- res_reg_ld1_start_high_df3 %>%
  select(parameter, ld1_start, bias, rel_bias) %>%
  unique() %>%
  pivot_longer(cols = c(bias, rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = ld1_start, values_from = value, names_prefix = "day ") %>%
  mutate(metric = factor(metric, levels = c("bias", "rel_bias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "rel_bias" ~ "relative bias (%)"))

metric_table_reg_ld1_start_high3 %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 2) %>%
  pack_rows("NPI 2", 3, 4)



metric_table_reg_ld1_strength3 <- res_reg_ld1_strength_df3 %>%
  select(parameter, ld1_strength, bias, rel_bias) %>%
  unique() %>%
  pivot_longer(cols = c(bias, rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = ld1_strength, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("bias", "rel_bias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "rel_bias" ~ "relative bias (%)"))

metric_table_reg_ld1_strength3 %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 2) %>%
  pack_rows("NPI 2", 3, 4)


# save results
save(res_reg_ld1_start_df3, res_reg_ld1_start_high_df3, res_reg_ld1_strength_df3, 
     res_EE_ld1_start_df3 , res_EE_ld1_start_high_df3, res_EE_ld1_strength_df3, 
     file = "Rt_NPI_scen_res_dfs3.RData")
