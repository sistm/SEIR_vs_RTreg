library(dplyr)
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "/cm/shared/dev/modules/generic/apps/tools/monolix/2021R2/")



ProjectName = "sim_SEIRAHD_2params_init_est"
data_file = "data_SEIRAHD_2params_init_est.txt"
model_file = "code_SEIRAHD_BG1.txt"
  
newProject(data = list(dataFile = data_file, 
                         headerTypes = c("id", "time", "observation", "obsid", rep("regressor", 4)), 
                         observationTypes = list(predHin = "continuous", predHprev = "continuous", predI = "continuous", predD = "continuous"), 
                         mapping = list("1" = 'predHin', "2" = 'predHprev', "3" = 'predI', "4" = 'predD')), 
             modelFile = model_file)
  
setErrorModel(ypredHin = "combined1", ypredHprev = "proportional", ypredI = "combined1", ypredD = "combined1")
  
setPopulationParameterInformation(b1_pop = list(initialValue = 0.6, method = "MAP", priorValue = 0.6, priorSD = 0.05),
                                  omega_b1 = list(initialValue = 0.1), 
                                  beta_ld1_pop = list(initialValue = 1.45),
                                  beta_BG1_pop = list(initialValue = 0.5), 
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
                                        b1 = "logNormal", log_initE = "logNormal"))
setIndividualParameterVariability(list(beta_ld1 = FALSE, beta_BG1 = FALSE, 
                                        b1 = TRUE, log_initE = TRUE))
  
  
setPopulationParameterEstimationSettings(variability = "firstStage")
  
setPopulationParameterEstimationSettings(nbburningiterations = 25, 
                                         nbexploratoryiterations = 500, 
                                         nbsmoothingiterations = 250,
                                         exploratoryautostop = TRUE, 
                                         smoothingautostop = TRUE) 
  
saveProject(projectFile = paste0(ProjectName, ".mlxtran"))
  
scenario = getScenario()
scenario$tasks[4:5] <- TRUE
scenario$tasks[6] <- FALSE
setScenario(scenario)
  
lixoftConnectors::runScenario()
saveProject()
 

