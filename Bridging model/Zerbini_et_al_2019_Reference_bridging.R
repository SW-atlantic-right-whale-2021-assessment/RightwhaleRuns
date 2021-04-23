# Code for Bridging HumpbackSIR to StateSpaceSIR

# Load data
source("Bridging model/InputData_HWAssessment_InitialAssessments_Oct2018.R")
library(StateSpaceSIR)

# Reference
Rel.Abundance.Pavanato$Index <- 2
rel_abund_ref <- rbind(Rel.Abundance.Branch, Rel.Abundance.Pavanato)
Rel.Abundance.Pavanato$Index <- 1

################################################################################
# Reference
################################################################################
file_name <- "Bridging model/Reference/Reference_bridging"
sir_reference <- list()
for(i in 1:2){
  sir_reference[[i]] <-  StateSpaceSIR(file_name = paste0(file_name, c("","prior")[i]),
                                      n_resamples = 1000,
                                      priors = make_prior_list(r_max = make_prior(runif, 0, 0.118),
                                                               N_obs = make_prior(runif, 500, 40000),
                                                               var_N = make_prior(0)),
                                      catch_multipliers = make_multiplier_list(
                                        make_prior(rnorm, 1.0185 , 0.0028)),
                                      premodern_catch_multipliers = make_multiplier_list(
                                        make_prior(rnorm, 1.71, 0.073)),
                                      target.Yr = 2008,
                                      num.haplotypes = 0,
                                      output.Yrs = c(2012, 2006, 2019, 2030),
                                      abs.abundance = rbind(Abs.Abundance.2008, Abs.Abundance.2012),
                                      rel.abundance = rel_abund_ref,
                                      rel.abundance.key = TRUE, # Indices of abundance
                                      count.data = Count.Data,
                                      count.data.key = FALSE,
                                      growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
                                      growth.rate.Yrs = c(1995, 1996, 1997, 1998),
                                      catch.data = Core.Catches2,
                                      premodern_catch_data = merge(PreModern.Catch.Min, PreModern.Catch.Max, by = "Year", all = T),
                                      control = sir_control(threshold = 0.01 * 1e-23, progress_bar = TRUE),
                                      realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_reference[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_reference[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_reference, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_reference[[1]],  file_name = file_name)
plot_trajectory(sir_reference[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_reference[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_reference[[2]]), inc_reference = FALSE)
plot_ioa(sir_reference[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_reference[[1]],  file_name = file_name)



################################################################################
# State-Space invgamma on process
################################################################################
file_name <- "Bridging model/Reference_state_space/Reference_state_space"
sir_state_space <- list()
for(i in 1:2){
  sir_state_space[[i]] <-  StateSpaceSIR(file_name = paste0(file_name, c("","prior")[i]),
                                       n_resamples = 1000,
                                       priors = make_prior_list(r_max = make_prior(runif, 0, 0.118),
                                                                N_obs = make_prior(runif, 500, 40000),
                                                                var_N = make_prior(rinvgamma, 4, 0.1)),
                                       catch_multipliers = make_multiplier_list(
                                         make_prior(rnorm, 1.0185 , 0.0028)),
                                       premodern_catch_multipliers = make_multiplier_list(
                                         make_prior(rnorm, 1.71, 0.073)),
                                       target.Yr = 2008,
                                       num.haplotypes = 0,
                                       output.Yrs = c(2012, 2006, 2019, 2030),
                                       abs.abundance = rbind(Abs.Abundance.2008, Abs.Abundance.2012),
                                       rel.abundance = rel_abund_ref,
                                       rel.abundance.key = TRUE, # Indices of abundance
                                       count.data = Count.Data,
                                       count.data.key = FALSE,
                                       growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
                                       growth.rate.Yrs = c(1995, 1996, 1997, 1998),
                                       catch.data = Core.Catches2,
                                       premodern_catch_data = merge(PreModern.Catch.Min, PreModern.Catch.Max, by = "Year", all = T),
                                       control = sir_control(threshold = 0.02 * 1e-23, progress_bar = TRUE),
                                       realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_state_space[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_state_space[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_state_space, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_state_space[[1]],  file_name = file_name)
plot_trajectory(sir_state_space[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_state_space[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_state_space[[2]]), inc_reference = FALSE)
plot_ioa(sir_state_space[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_state_space[[1]],  file_name = file_name)



