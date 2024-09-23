library(dplyr)
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "/cm/shared/dev/modules/generic/apps/tools/monolix/2021R2/")


dir <- "~/MonolixFiles/2023-10-14_boot_sim_2params_all_new3_Simulx_SEIRAHD"


for(j in 1:2){
  # resample data
  data <- read.table(paste0("data_sim_SEIRAHD_Simulx_2params_new3_ME", j, ".txt"), header = TRUE, sep = ",") 
  
  for(i in 1:100){
    sampled_depts <- sample(data$dept_id, size = length(unique(data$dept_id)), replace = TRUE)
    resampled_data <- lapply(sampled_depts, function(x) filter(data, dept_id == x))
    resampled_data <- do.call("rbind", resampled_data)
    write.table(resampled_data, paste0(dir, "/resampled_data/data_2params_all_", j, "boot", i, ".txt"), 
                row.names = FALSE, sep = ",")
    
    ProjectName = paste0("SEIRAHD_Simulx_new3_2params_", j, "boot", i)
    data_file = paste0(dir, "/resampled_data/data_2params_all_", j, "boot", i, ".txt")
    model_file = "code_SEIRAHD_BG1.txt"
    
    newProject(data = list(dataFile = data_file, 
                           headerTypes = c("id", "time", "observation", "obsid", rep("regressor", 4)), 
                           observationTypes = list(predHin = "continuous", predHprev = "continuous", predI = "continuous", predD = "continuous"), 
                           mapping = list("1" = 'predHin', "2" = 'predHprev', "3" = 'predI', "4" = 'predD')), 
               modelFile = model_file)
    
    setErrorModel(ypredHin = "combined1", ypredHprev = "proportional", ypredI = "combined1", ypredD = "combined1")
    
    setPopulationParameterInformation(b1_pop = list(initialValue = 0.6),
                                      omega_b1 = list(initialValue = 0.1),
                                      beta_ld1_pop = list(initialValue = 1.45),
                                      beta_BG1_pop = list(initialValue = 0.8), 
                                      log_initE_pop = list(initialValue = 3), 
                                      omega_log_initE = list(initialValue = 0.3))
    
    setPopulationParameterInformation(apredHin = list(initialValue = 0.04),
                                      bpredHin = list(initialValue = 1.22),
                                      bpredHprev = list(initialValue = 1.65), 
                                      apredI = list(initialValue = 0.174), 
                                      bpredI = list(initialValue = 0.782), 
                                      apredD = list(initialValue = 0.01), 
                                      bpredD = list(initialValue = 0.8))
    
    
    setIndividualParameterDistribution(list(beta_ld1 = "logNormal", beta_BG1 = "logNormal", 
                                            b1 = "logitNormal", log_initE = "logNormal"))
    setIndividualParameterVariability(list(beta_ld1 = FALSE, beta_BG1 = FALSE, 
                                           b1 = TRUE, log_initE = TRUE))
    
    setIndividualLogitLimits(b1 = c(0, 1))
    
    
    
    popparams <- getPopulationParameterInformation()  
    
    popini <- sapply(1:(nrow(popparams)-7), function(j){
      if(popparams$initialValue[j] >= 0){
        runif(n=1, min=popparams$initialValue[j]/2, max=popparams$initialValue[j]*2)
      }else if(popparams$initialValue[j] < 0){
        runif(n=1, min=popparams$initialValue[j]*2, max=popparams$initialValue[j]/2)
      }
    })
    
    # set sampled values as new initial estimates (except for weather on position 9)
    newpopparams <- popparams
    newpopparams$initialValue[1:(nrow(popparams)-7)] <- popini
    
    setIndividualLogitLimits(b1 = c(min(0, newpopparams$initialValue[newpopparams$name == "b1_pop"]-0.01), max(1, newpopparams$initialValue[newpopparams$name == "b1_pop"]+0.01)))
    
    setPopulationParameterInformation(newpopparams)
    
    
    
    setPopulationParameterEstimationSettings(variability = "firstStage")
    
    setPopulationParameterEstimationSettings(nbburningiterations = 20, 
                                             nbexploratoryiterations = 500, 
                                             nbsmoothingiterations = 200,
                                             exploratoryautostop = TRUE, 
                                             smoothingautostop = TRUE) 
    
    saveProject(projectFile = paste0(ProjectName, ".mlxtran"))
    
    scenario = getScenario()
    scenario$tasks[4:5] <- FALSE
    scenario$tasks[6] <- FALSE
    setScenario(scenario)
    
    lixoftConnectors::runScenario()
    saveProject()
  }

}  



