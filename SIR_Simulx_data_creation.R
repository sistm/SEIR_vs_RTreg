library(tidyverse)
library(colorspace)
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx")

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/SIR_simulx_data_creation_2params")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/SIR_Simulx_function_2params.R")

dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/"
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"

popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id-1))


# Depletion of susceptibles scenarios -------------------------------

#### data generation #### 

project.file <- "sim_SIR_2params.mlxtran"
importProject(project.file)

initI_vals <- c(-3, -1.2, -0.4, 0.6, 1.6)
seeds <- 101:200

sim_res_late_initI_neg3_list <- list()
sim_res_late_initI_neg1.2_list <- list()
sim_res_late_initI_neg0.4_list <- list()
sim_res_late_initI0.6_list <- list()
sim_res_late_initI1.6_list <- list()

for(j in 1:100){
  for(i in 1:5){
    set.seed(seeds[j])
    inits <- data.frame(id = 1:94, 
                        b1 = rnorm(94, 0.425, 0.01),
                        log_initI = rnorm(94, initI_vals[i], 0.04), 
                        beta_ld1 = 1.45, 
                        beta_BG1 = 0.8)
    
    write.table(inits, file = paste0(getwd(), "/inits/inits_late_initI", initI_vals[i], "_", j, ".txt"), 
                row.names = FALSE)
    
    
    
    sim_res <- sim_SIR_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/inits/inits_late_initI", initI_vals[i], "_", j, ".txt"), 
                                         regressor_df_path = "ld1_reg_df_2params_shortLD_long_0.8.csv",
                                         tmax = 100)
    
    if(i == 1){
      sim_res_late_initI_neg3_list[[j]] <- sim_res
    } else if(i == 2){
      sim_res_late_initI_neg1.2_list[[j]] <- sim_res
    } else if(i == 3){
      sim_res_late_initI_neg0.4_list[[j]] <- sim_res
    } else if(i == 4){
      sim_res_late_initI0.6_list[[j]] <- sim_res
    } else {
      sim_res_late_initI1.6_list[[j]] <- sim_res
    }


    monolix_SIR <- SIR_monolix_data_creation_ME(simulation_results = sim_res,
                                                popsize_df = popsize_df,
                                                start_date = as.Date("2020-03-02"),
                                                end_date = as.Date("2020-06-09"),
                                                ld1_start = 35, ld1_end = 55,
                                                BG1_start = 56)
    write.table(monolix_SIR, file = paste0(getwd(), "/Datasets_depl_sus_late/data_sim_late_initI", initI_vals[i], "_", j, ".txt"),
                row.names = FALSE, sep = ",")
  }
}

save(sim_res_late_initI_neg3_list, file = paste0(getwd(), "/Datasets_depl_sus_late/sim_res_late_initI_neg3_list.RData"))
save(sim_res_late_initI_neg1.2_list, file = paste0(getwd(), "/Datasets_depl_sus_late/sim_res_late_initI_neg1.2_list.RData"))
save(sim_res_late_initI_neg0.4_list, file = paste0(getwd(), "/Datasets_depl_sus_late/sim_res_late_initI_neg0.4_list.RData"))
save(sim_res_late_initI0.6_list, file = paste0(getwd(), "/Datasets_depl_sus_late/sim_res_late_initI0.6_list.RData"))
save(sim_res_late_initI1.6_list, file = paste0(getwd(), "/Datasets_depl_sus_late/sim_res_late_initI1.6_list.RData"))


#### plots ####
initI_late <- c(-3, -1.2, -0.4, 0.6, 1.6)
data_list_late <- list(sim_res_late_initI_neg3_list, sim_res_late_initI_neg1.2_list, 
                  sim_res_late_initI_neg0.4_list, sim_res_late_initI0.6_list, 
                  sim_res_late_initI1.6_list)

S_start_list_late <- list()
R0_list_late <- list()

for(i in 1:5){
  S_start <- c()
  R0 <- c()
  for(j in 1:100){
    S_start <- c(S_start, data_list_late[[i]][[j]]$S[1])
    R0 = c(R0, data_list_late[[i]][[j]]$transmission[1]*5)
  }
  S_start_list_late[[i]] <- range(S_start)
  R0_list_late[[i]] <- range(R0)
}
S_start_list_late
R0_list_late


plot_S <- function(df){
  ggplot(df, aes(x = time, y = S, group = id)) + 
    geom_line() +
    scale_x_continuous(expand = c(0.01, 0.01)) + 
    scale_y_continuous(limits = c(0, 10000)) + 
    labs(x = "Day", y = "Susceptibles") +
    theme_bw() + 
    theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
          legend.title = element_text(family = "serif", size = 13),
          plot.title = element_text(family = "serif", size = 16), 
          axis.title = element_text(family = "serif", size = 13), 
          axis.text.x = element_text(family = "serif", size = 12), 
          axis.text.y = element_text(family = "serif", size = 12))
  
}

plot_IncI <- function(df, popsize_df = popsize_df, model_name){
  plot_df <- df %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    mutate(IncI_unscaled = IncI*popsize/10^4)
  
  ggplot(plot_df, aes(x = time, y = IncI, group = id)) + 
    geom_line() +
    scale_x_continuous(expand = c(0.01, 0.01)) + 
    labs(x = "Day", y = "Incident cases") +
    #scale_y_continuous(limits = c(0, 800)) + 
    theme_bw() + 
    theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
          legend.title = element_text(family = "serif", size = 13),
          plot.title = element_text(family = "serif", size = 16), 
          axis.title = element_text(family = "serif", size = 13), 
          axis.text.x = element_text(family = "serif", size = 12), 
          axis.text.y = element_text(family = "serif", size = 12))
  
}

# depletion of susceptibles before NPI 1 implementation
depl_sus_NPI1 <- c(2, 10, 20, 40, 60) 
# depletion of susceptibles at the end of simulation
depl_sus_end <- c(7, 20, 30, 50, 65) 

rect_cols <- sequential_hcl(5, palette = "BluYl")
for(i in 1:5){
  plot <- plot_S(df = data_list_late[[i]][[1]]) + 
    labs(title = paste0("Late NPI 1, S depletion before NPI1: ~", depl_sus_NPI1[i], "%")) + 
    annotate("rect", xmin = 35, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
    annotate("rect", xmin = 55, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
    annotate("label", x = 45, y = 0, label = "NPI 1", 
             hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif") + 
    annotate("label", x = 77, y = 0, label = "NPI 2", 
             hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif")
  print(plot)
}

for(i in 1:5){
  plot <- plot_IncI(df = data_list_late[[i]][[1]], popsize_df = popsize_df) + 
    scale_y_continuous(limits = c(0, 460)) +
    labs(title = paste0("Late NPI 1, S depletion before NPI1: ~", depl_sus_NPI1[i], "%")) + 
    annotate("rect", xmin = 35, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
    annotate("rect", xmin = 55, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
    annotate("label", x = 45, y = 460, label = "NPI 1", 
             hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif") + 
    annotate("label", x = 77, y = 460, label = "NPI 2", 
             hjust = 0.5, vjust = 0, size = 4, fontface = 2, family = "serif")
  print(plot)
}
