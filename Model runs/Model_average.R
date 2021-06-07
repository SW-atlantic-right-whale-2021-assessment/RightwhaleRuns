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
  model_names = c( "R", paste0("D ", 1:7), paste0("C ", 1:7), paste0("G ", 1:2), paste0("M ", 1:2), "MA"), 
  bayes_factor = bayes_vec,
  file_name = "Cross scenario comparison/Figure_3_",
  years = c(2008, 2019))

# Plot and get parameter values from Model Average
file_name <- "Model runs/Model average/model_average"
plot_trajectory(new_mod, Reference = sir_reference[[1]],  file_name = file_name)
plot_density(SIR = list(sir_reference[[1]], new_mod), priors = list(sir_reference[[2]]),  file_name = file_name,  lower = c(NA, 20000, NA, NA, NA, 15000, NA, 21000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1))
plot_ioa(new_mod,  file_name = file_name, ioa_names = c("FG", "BG1"))
summary_table(new_mod,  file_name = file_name)
save(new_mod, file = paste0(file_name, ".Rdata"))