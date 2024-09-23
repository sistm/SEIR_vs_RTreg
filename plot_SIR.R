library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(ggside)
library(kableExtra)
library(patchwork)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")
dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SIR_simulx_data_creation_2params"
dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SIR_bootstrap"
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")



# Load data ---------------------------------------------------------------

popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id-1))

# simulated datasets
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI_neg3_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI_neg1.2_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI_neg0.4_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI0.6_list.RData"))
load(paste0(dir, "/Datasets_depl_sus_late/sim_res_late_initI1.6_list.RData"))


# regression bootstraps
load(paste0(dir3, "/res_boot_SIR_late_InitI_neg3_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI_neg1.2_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI_neg0.4_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI0.6_df.RData"))
load(paste0(dir3, "/res_boot_SIR_late_InitI1.6_df.RData"))


# mechanistic model bootstraps
load(paste0(dir1, "/SIR_boot_list_initI_neg3.RData"))
load(paste0(dir1, "/SIR_boot_list_initI_neg1.2.RData"))
load(paste0(dir1, "/SIR_boot_list_initI_neg0.4.RData"))
load(paste0(dir1, "/SIR_boot_list_initI0.6.RData"))
load(paste0(dir1, "/SIR_boot_list_initI1.6.RData"))




# CI coverage figure -----------------------------------------------------------------

df_reg_late <- res_boot_SIR_late_initI_neg3_df %>% ungroup() %>% mutate(initI = -3) %>%
  bind_rows(res_boot_SIR_late_initI_neg1.2_df %>% ungroup() %>% mutate(initI = -1.2)) %>%
  bind_rows(res_boot_SIR_late_initI_neg0.4_df %>% ungroup() %>% mutate(initI = -0.4)) %>%
  bind_rows(res_boot_SIR_late_initI0.6_df %>% ungroup() %>% mutate(initI = 0.6)) %>%
  bind_rows(res_boot_SIR_late_initI1.6_df %>% ungroup() %>% mutate(initI = 1.6)) %>%
  select(parameter, median, CI_LL, CI_UL, initI, method, rep) %>%
  filter(method %in% c("quantile", "simple regression")) %>%
  unique() %>%
  rename(value = median) %>%
  mutate(method = ifelse(method == "quantile", "bootstrap regression", "simple regression"))

SIR_boot_df <- do.call("rbind.data.frame", SIR_boot_list_initI_neg3) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI_neg1.2)) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI_neg0.4)) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI0.6)) %>%
  bind_rows(do.call("rbind.data.frame", SIR_boot_list_initI1.6)) %>%
  select(-contains("1")) %>%
  rename_with(~str_remove(., "2"), everything()) %>%
  rename(value = mean_est, rep = sim_rep) %>%
  mutate(parameter = ifelse(parameter == "beta_ld1", "NPI 1", "NPI 2"), 
         method = "mechanistic model")

df_all_SIR <- bind_rows(df_reg_late, SIR_boot_df) %>%
  mutate(method = factor(method, 
                         levels = c("mechanistic model", "simple regression", "bootstrap regression")), 
         depl_sus = case_when(initI == -3 ~ "2%", 
                              initI == -1.2 ~ "10%",
                              initI == -0.4 ~ "20%",
                              initI == 0.6 ~ "40%",
                              initI == 1.6 ~ "60%"), 
         depl_sus = factor(depl_sus, 
                           levels = c("2%", "10%", "20%", "40%", "60%")))

br_palette <- diverging_hcl("Blue-Red", n = 20)
SIR_cols <- c(br_palette[c(20, 5, 1)])

ggplot(df_all_SIR, aes(x = rep, y = value, ymin = CI_LL, ymax = CI_UL, col = method)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  facet_grid(cols = vars(parameter), rows = vars(depl_sus)) + 
  geom_hline(data = data.frame(parameter = c("NPI 1", "NPI 2"), 
                               true_value = c(-1.45, -0.8)), 
             aes(yintercept = true_value), linetype = "dashed") +
  scale_color_manual(values = SIR_cols) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(#title = "Late NPI 1, 60% depletion of susceptibles", 
       x = "Data set", y = "Regression coefficient", col = "Estimation method")  +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16, hjust = 0), 
        axis.title = element_text(family = "serif", size = 15), 
        axis.text.x = element_text(family = "serif", size = 13.5), 
        axis.text.y = element_text(family = "serif", size = 13.5), 
        legend.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(family = "serif", size = 15), 
        strip.text = element_text(family = "serif", size = 15), 
        strip.background = element_rect(fill = "grey90"),
        legend.position = "bottom")

ggsave("SIR_depl_sus_plot_all.jpeg", width = 10, height = 12)

# Evaluation Metrics -----------------------------------------------------------------

metrics_df <- df_all_SIR %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  group_by(parameter, method, depl_sus) %>%
  summarize(perc_CI_covers = sum(CI_covers), 
            mean_bias = mean(bias), 
            mean_rel_bias = mean(rel_bias)) %>%
  pivot_longer(cols = c(perc_CI_covers, mean_bias, mean_rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = c(method, depl_sus), values_from = value) %>%
  mutate(metric = factor(metric, levels = c("mean_bias", "mean_rel_bias", "perc_CI_covers"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "mean_bias" ~ "absolute bias", 
                            metric == "perc_CI_covers" ~ "CI coverage (%)", 
                            metric == "mean_rel_bias" ~ "relative bias (%)"))


# print metrics table
metrics_df %>% 
  ungroup() %>%
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)



# Regression fits ---------------------------------------------------------

cl <- makeCluster(12)
registerDoParallel(cl)

new_tstart <- 2:100
new_tend <- 2:100

#### 2% depl sus ####
reg_data2 <- sim_res_late_initI_neg3_list[[1]] %>% 
  rename(dept_id = id, day = time) %>%
  left_join(popsize_df, by = "dept_id") %>%
  mutate(IncI_unscaled = IncI*popsize/10^4)

true_Rt_df2 <- reg_data2 %>%
  mutate(Rt_real = transmission*S/(0.2*10000)) %>%
  dplyr::select(dept_id, day, Rt_real)


Rt_comp_res2 <- Rt_calc_fun(Inc_df = reg_data2, id_col_name = "dept_id", time_col_name = "day", 
                           Inc_col_name = "IncI_unscaled", model_name = "2% depl_sus", 
                           Rt_ref_df = true_Rt_df2, 
                           meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                           tstart = new_tstart, tend = new_tend)

first_Rt_day <- min(Rt_comp_res2$Rt_comp$day) + 3 # +3 because this is when Rt estimates are available
seven_day_seq <- seq(first_Rt_day, max(Rt_comp_res2$Rt_comp$day), 7)

reg_data2_7d <- Rt_comp_res2$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
         BG1 = ifelse(day > 55, 1, 0)) %>%
  filter(day %in% seven_day_seq)

reg_res2 <- Rt_reg_only_fun(data_for_est = reg_data2_7d, 
                           model_name = "2% depl_sus", lag_NPIs = FALSE, fits = TRUE)


res_df_fitted2 <-  reg_data2 %>%
  left_join(true_Rt_df2, by = c("dept_id", "day")) %>%
  left_join(Rt_comp_res2$Rt_comp %>% select(dept_id, day, Rt), by = c("dept_id", "day")) %>%
  left_join(reg_res2$fits %>% select(dept_id, day, Rt_fitted), by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))

# set.seed(12345)
# selected_depts <- sample(unique(reg_data$dept_id), 9)



#### 60% depl sus ####
reg_data60 <- sim_res_late_initI1.6_list[[1]] %>% 
  rename(dept_id = id, day = time) %>%
  left_join(popsize_df, by = "dept_id") %>%
  mutate(IncI_unscaled = IncI*popsize/10^4)

true_Rt_df60 <- reg_data60 %>%
  mutate(Rt_real = transmission*S/(0.2*10000)) %>%
  dplyr::select(dept_id, day, Rt_real)


Rt_comp_res60 <- Rt_calc_fun(Inc_df = reg_data60, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "60% depl_sus", 
                            Rt_ref_df = true_Rt_df60, 
                            meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                            tstart = new_tstart, tend = new_tend)

reg_data60_7d <- Rt_comp_res60$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 35, 55), 1, 0), 
         BG1 = ifelse(day > 55, 1, 0)) %>%
  filter(day %in% seven_day_seq)

reg_res60 <- Rt_reg_only_fun(data_for_est = reg_data60_7d,
                            model_name = "60% depl_sus", lag_NPIs = FALSE, fits = TRUE)


res_df_fitted60 <-  reg_data60 %>%
  left_join(true_Rt_df60, by = c("dept_id", "day")) %>%
  left_join(Rt_comp_res60$Rt_comp %>% select(dept_id, day, Rt), by = c("dept_id", "day")) %>%
  left_join(reg_res60$fits %>% select(dept_id, day, Rt_fitted), by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


#### plots ####
selected_depts <- c(2, 24, 51, 90)

rect_cols <- sequential_hcl(5, palette = "BluYl")

p2 <- ggplot(res_df_fitted2 %>% filter(dept_id %in% selected_depts & day <= 75), 
       aes(x = day)) +
  annotate("rect", xmin = 35, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 55, xmax = 75, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 45, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 3, fontface = 2, family = "serif") + 
  annotate("label", x = 65, y = Inf, label = "NPI 2", 
           hjust = 0.5, vjust = 1, size = 3, fontface = 2, family = "serif") + 
  geom_line(aes(y = Rt, col = "Epi Estim"), linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit")) + 
  geom_line(aes(y = Rt_real, col = "True"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(title = "A - 2% depletion of susceptibles", 
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

p60 <- ggplot(res_df_fitted60 %>% filter(dept_id %in% selected_depts & day <= 75), 
             aes(x = day)) +
  annotate("rect", xmin = 35, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 55, xmax = 75, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 45, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 3, fontface = 2, family = "serif") + 
  annotate("label", x = 65, y = Inf, label = "NPI 2", 
           hjust = 0.5, vjust = 1, size = 3, fontface = 2, family = "serif") + 
  geom_line(aes(y = Rt, col = "Epi Estim"), linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(y = Rt_fitted, col = "Regression fit")) + 
  geom_line(aes(y = Rt_real, col = "True"), linewidth = 0.8) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(title = "B - 60% depletion of susceptibles", 
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

p2 / p60 + plot_layout(guides = 'collect')

ggsave("SIR_depl_sus_fits.jpeg", width = 7.5, height = 10)




# Data creation figs ------------------------------------------------------

initI_late <- c(-3, -1.2, -0.4, 0.6, 1.6)
data_list_late <- list(sim_res_late_initI_neg3_list, sim_res_late_initI_neg1.2_list, 
                       sim_res_late_initI_neg0.4_list, sim_res_late_initI0.6_list, 
                       sim_res_late_initI1.6_list)
depl_sus_NPI1 <- c(2, 10, 20, 40, 60) 

plot_data <- sim_res_late_initI_neg3_list[[1]] %>% mutate(depl_sus = "2% depletion") %>%
  bind_rows(sim_res_late_initI_neg1.2_list[[1]] %>% mutate(depl_sus = "10% depletion")) %>%
  bind_rows(sim_res_late_initI_neg0.4_list[[1]] %>% mutate(depl_sus = "20% depletion")) %>%
  bind_rows(sim_res_late_initI0.6_list[[1]] %>% mutate(depl_sus = "40% depletion")) %>%
  bind_rows(sim_res_late_initI1.6_list[[1]] %>% mutate(depl_sus = "60% depletion")) %>%
  mutate(depl_sus = factor(depl_sus, levels = paste0(depl_sus_NPI1, "% depletion")))

#### parameter ranges ####
S_start_list_late <- list()
R0_list_late <- list()
transmission_list_late <- list()

for(i in 1:5){
  S_start <- c()
  R0 <- c()
  b0 <- c()
  for(j in 1:100){
    S_start <- c(S_start, data_list_late[[i]][[j]]$S[1])
    R0 = c(R0, data_list_late[[i]][[j]]$transmission[1]*5)
    b0 = c(b0, data_list_late[[i]][[j]]$transmission[1])
  }
  S_start_list_late[[i]] <- range(S_start)
  R0_list_late[[i]] <- range(R0)
  transmission_list_late[[i]] <- range(b0)
}
S_start_list_late
R0_list_late
transmission_list_late


#### plots ####
rect_cols <- sequential_hcl(5, palette = "BluYl")
pS <- ggplot(plot_data, aes(x = time, y = S, group = id))  + 
  annotate("rect", xmin = 35, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 55, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 45, y = 0, label = "NPI 1", 
           hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 77, y = 0, label = "NPI 2", 
           hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif") + 
  geom_line(linewidth = 0.2) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(1, 100), breaks = seq(25, 100, 25)) + 
  scale_y_continuous(limits = c(0, 10000)) + 
  facet_wrap(~depl_sus, nrow = 1) + 
  labs(x = "Day", y = "Susceptibles") +
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13))

pI <- ggplot(plot_data, aes(x = time, y = IncI, group = id)) + 
  annotate("rect", xmin = 35, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 55, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 45, y = 450, label = "NPI 1", 
           hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 77, y = 450, label = "NPI 2", 
           hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif") +
  geom_line(linewidth = 0.2) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(1, 100), breaks = seq(25, 100, 25)) + 
  scale_y_continuous(limits = c(0, 470)) + 
  facet_wrap(~depl_sus, nrow = 1) + 
  labs(x = "Day", y = "Incident cases") +
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13))

pS / pI + plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(family = "serif", size = 18, hjust = 0, vjust = 0))

ggsave("SIR_data_creation_plot.jpeg", width = 12, height = 7, dpi = 300)
                  