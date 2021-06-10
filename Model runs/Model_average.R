library(StateSpaceSIR)

# Load all the models
file_names <- c("Base/Base",
                "Sensitivity_1/Sensitivity_1",
                "Sensitivity_2/Sensitivity_2",
                "Sensitivity_3/Sensitivity_3",
                "Sensitivity_4/Sensitivity_4",
                "Sensitivity_5/Sensitivity_5",
                "Sensitivity_6/Sensitivity_6",
                "Sensitivity_7/Sensitivity_7",
                "Sensitivity_8/Sensitivity_8",
                "Sensitivity_9/Sensitivity_9",
                "Sensitivity_10/Sensitivity_10",
                "Model_average/Model_average")

for(i in 1:length(file_names)){
  load(file = paste0("Model runs/",file_names[i], ".Rdata"))
}

# Densities
plot_density(SIR = list(sir_base[[1]]),  file_name = paste0("Model runs/",file_names[1]),   priors = list(sir_base[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_1[[1]]),  file_name = paste0("Model runs/",file_names[2]),   priors = list(sensitivity_1[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_2[[1]]),  file_name = paste0("Model runs/",file_names[3]),   priors = list(sensitivity_2[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_3[[1]]),  file_name = paste0("Model runs/",file_names[4]),   priors = list(sensitivity_3[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_4[[1]]),  file_name = paste0("Model runs/",file_names[5]),   priors = list(sensitivity_4[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_5[[1]]),  file_name = paste0("Model runs/",file_names[6]),   priors = list(sensitivity_5[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_6[[1]]),  file_name = paste0("Model runs/",file_names[7]),   priors = list(sensitivity_6[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_7[[1]]),  file_name = paste0("Model runs/",file_names[8]),   priors = list(sensitivity_7[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_8[[1]]),  file_name = paste0("Model runs/",file_names[9]),   priors = list(sensitivity_8[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_9[[1]]),  file_name = paste0("Model runs/",file_names[10]),   priors = list(sensitivity_9[[2]]), inc_reference = FALSE)
plot_density(SIR = list(sensitivity_10[[1]]),  file_name = paste0("Model runs/",file_names[11]),   priors = list(sensitivity_10[[2]]), inc_reference = FALSE)




#############################################################
#### Model averaging
#############################################################
# Get bayes factor for models with comparable likelihoods
bayes_f <- bayes_factor(SIR = list(sir_base[[1]],
                                   sensitivity_1[[1]], 
                                   sensitivity_2[[1]], 
                                   sensitivity_3[[1]],
                                   sensitivity_4[[1]],
                                   #sensitivity_5[[1]],
                                   #sensitivity_6[[1]],
                                   sensitivity_7[[1]],
                                   sensitivity_8[[1]],
                                   sensitivity_9[[1]],
                                   sensitivity_10[[1]]))


# Create a new model based on bayes factors
model_average <- weight_model(SIR = list(sir_base[[1]],
                                   sensitivity_1[[1]], 
                                   sensitivity_2[[1]], 
                                   sensitivity_3[[1]],
                                   sensitivity_4[[1]],
                                   #sensitivity_5[[1]],
                                   #sensitivity_6[[1]],
                                   sensitivity_7[[1]],
                                   sensitivity_8[[1]],
                                   sensitivity_9[[1]],
                                   sensitivity_10[[1]]), 
                        bayes_factor = bayes_f)

# For plotting make a vector of bayes factors, set NA for models that cant be compared (different likelihood)
bayes_vec <- round(c(bayes_f[1:5],  NA, NA, bayes_f[6:9], NA), 2)


# Compare Aposteriors of all
compare_posteriors(
  reference_sir = TRUE, 
  SIR = list(sir_base[[1]],
             sensitivity_1[[1]], 
             sensitivity_2[[1]], 
             sensitivity_3[[1]],
             sensitivity_4[[1]],
             sensitivity_5[[1]],
             sensitivity_6[[1]],
             sensitivity_7[[1]],
             sensitivity_8[[1]],
             sensitivity_9[[1]],
             sensitivity_10[[1]],
             model_average), 
  model_names = c( "B", paste0("S-", 1:10), "MA"), 
  bayes_factor = bayes_vec,
  file_name = paste0("Model runs/",file_names[11]),
  years = c(2021, 2030))

# Plot and get parameter values from Model Average
file_name <-paste0("Model runs/",file_names[11])
plot_trajectory(model_average, Reference = sir_base[[1]],  file_name = file_name)
plot_density(SIR = list(sir_base[[1]], model_average), priors = list(sir_base[[2]]),  file_name = file_name)
plot_ioa(model_average,  file_name = file_name, ioa_names = c("FG", "BG1"))
summary_table(model_average,  file_name = file_name)
save(model_average, file = paste0(file_name, ".Rdata"))
