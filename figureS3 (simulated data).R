library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)
library(patchwork)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")

dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir10 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once10"

data_SEIRAHD3 <- read.table(paste0(dir1, "/data_sim_SEIRAHD_Simulx_2params_new3_ME1.txt"), 
                            header = TRUE, sep = ",")
data_ABM_rm_cov10 <- read.csv(paste0(dir10, "/data_covasim_rm10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)
data_ABM_hybrid_cov10 <- read.csv(paste0(dir10, "/data_covasim_hybrid10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)

rect_cols <- sequential_hcl(5, palette = "BluYl")

#### SEIRAHD 3 & ABM 6/7 ####
# plot with higher NPI 2 data

p1 <- ggplot(data_SEIRAHD3 %>% filter(obs_id == 3), 
             aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "SEIRAHD model cases", 
       x = "Day", y = "Cases") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p2 <- ggplot(data_SEIRAHD3 %>% filter(obs_id == 1), 
             aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "SEIRAHD model hospitalizations", 
       x = "Day", y = "Hospital admissions") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

#### ABM 10 ####
p3 <- ggplot(data_ABM_rm_cov10, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Random mixing ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p4 <- ggplot(data_ABM_hybrid_cov10, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line()  + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Multi-layer ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("Data generation BG0.8 models.jpeg", dpi = 400, width = 16, height = 9.5)

