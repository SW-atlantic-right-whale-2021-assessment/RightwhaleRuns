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
  sir_reference[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.118),
                             N_obs = make_prior(runif, 500, 40000),
                             var_N = make_prior(0)),
    catch_multipliers = make_multiplier_list(
      make_prior(rnorm, 1.0185 , 0.0028), 
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
    catch.data = list(Core.Catches2, merge(PreModern.Catch.Min, PreModern.Catch.Max, by = "Year", all = T)),
    control = sir_control(threshold = 0.01 * 1e-22, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_reference[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_reference[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_reference, file = paste0(file_name, ".Rdata"))


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
  sir_state_space[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.118),
                             N_obs = make_prior(runif, 500, 40000),
                             var_N = make_prior(rinvgamma, 4, 0.1)),
    catch_multipliers = make_multiplier_list(
      make_prior(rnorm, 1.0185 , 0.0028), 
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
    catch.data = list(Core.Catches2, merge(PreModern.Catch.Min, PreModern.Catch.Max, by = "Year", all = T)),
    control = sir_control(threshold = 0.02 * 1e-22, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_state_space[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_state_space[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_state_space, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_state_space[[1]],  file_name = file_name)
plot_trajectory(sir_state_space[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_state_space[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_state_space[[2]]), inc_reference = FALSE)
plot_ioa(sir_state_space[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_state_space[[1]],  file_name = file_name)


################################################################################
# Read in and update catch
################################################################################
sw_right_data<-read.delim("Data/datosModeloBallenasmiles2020Miles.csv", sep=";",header=FALSE)   
names(sw_right_data)<- c("Year","CatchMin","CatchMax","Nt")

# Two periods of SLRs
# period 1: 1648-1770: struck and lost rate factor = 1
# period 2: 1771-1850: struck and lost rate factor ~ norm(1.60, 0.04^2 )
# period 3: 1851-1973: struck and lost rate factor ~ norm(1.09, 0.04^2 )
# period 4: 1974-2030: struck and lost rate factor = 1

catch_list <- list(sw_right_data[which(sw_right_data$Year < 1771),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1771 & sw_right_data$Year <= 1850),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1851 & sw_right_data$Year <= 1973),1:3],
                   sw_right_data[which(sw_right_data$Year > 1973),1:3])


file_name <- "Bridging model/Update_catch/Update_catch"
sir_update_catch <- list()
for(i in 1:2){
  sir_update_catch[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.118),
                             N_obs = make_prior(runif, 500, 40000),
                             var_N = make_prior(rinvgamma, 4, 0.1)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
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
    catch.data = catch_list,
    control = sir_control(threshold = 0.02 * 1e-22, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_update_catch[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_update_catch[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_update_catch, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_update_catch[[1]],  file_name = file_name)
plot_trajectory(sir_update_catch[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_update_catch[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_update_catch[[2]]), inc_reference = FALSE)
plot_ioa(sir_update_catch[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_update_catch[[1]],  file_name = file_name)



################################################################################
# Update priors
################################################################################
file_name <- "Bridging model/Update_priors/Update_priors"
sir_update_priors <- list()
for(i in 1:2){
  sir_update_priors[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.6, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 0,
    output.Yrs = c(2012, 2006, 2019, 2030),
    abs.abundance = rbind(Abs.Abundance.2008, Abs.Abundance.2012),
    rel.abundance = rel_abund_ref,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 0.06 * 1e-24, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_update_priors[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_update_priors[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_update_priors, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_update_priors[[1]],  file_name = file_name)
plot_trajectory(sir_update_priors[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_update_priors[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_update_priors[[2]]), inc_reference = FALSE)
plot_ioa(sir_update_priors[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_update_priors[[1]],  file_name = file_name)



################################################################################
# Update abundance
################################################################################
file_name <- "Bridging model/Update_abundance/Update_abundance"

Abs.Abundance.2009 <- data.frame(Year = 2009, N.obs = 4029, CV.obs = NA) # FIXME: not used as of 4/24/21
Abs.Abundance.2010 <- data.frame(Year = 2010, N.obs = 4245, CV.obs = 245/4245) # 2010: 4245 (SE: 245, 95% CI 3,765, 4,725). 


sir_update_abundance <- list()
for(i in 1:2){
  sir_update_abundance[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.6, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 0,
    output.Yrs = c(2021, 2030),
    abs.abundance = Abs.Abundance.2010,
    rel.abundance = rel_abund_ref,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 0.06 * 1e-19, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_update_abundance[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_update_abundance[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_update_abundance, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_update_abundance[[1]],  file_name = file_name)
plot_trajectory(sir_update_abundance[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_update_abundance[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_update_abundance[[2]]), inc_reference = FALSE)
plot_ioa(sir_update_abundance[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_update_abundance[[1]],  file_name = file_name)


################################################################################
# Update relative abundance
################################################################################
file_name <- "Bridging model/Update_relative_abundance/Update_relative_abundance"
sw_right_data.RelAbundance <- sw_right_data[which(sw_right_data$Nt>0),]
Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_data.RelAbundance)), Year = sw_right_data.RelAbundance$Year, IA.obs = sw_right_data.RelAbundance$Nt, CV.IA.obs = rep(0, nrow(sw_right_data.RelAbundance))) #Using 0.2 as a proxy

sir_update_relative_abundance <- list()
for(i in 1:2){
  sir_update_relative_abundance[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             add_VAR_IA = make_prior(rinvgamma, 2, 0.5),
                             Pmsy = make_prior(runif, 0.6, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 0,
    output.Yrs = c(2021, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 0.1 * 1e-64, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_update_relative_abundance[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_update_relative_abundance[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_update_relative_abundance, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_update_relative_abundance[[1]],  file_name = file_name)
plot_trajectory(sir_update_relative_abundance[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_update_relative_abundance[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_update_relative_abundance[[2]]), inc_reference = FALSE)
plot_ioa(sir_update_relative_abundance[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_update_relative_abundance[[1]],  file_name = file_name)


################################################################################
# Update old base run (no absolute abundance, estimate variance of IA, estimate q with prior, change catch)
################################################################################
file_name <- "Bridging model/Update_old_base/Update_old_base"
sw_right_data.RelAbundance <- sw_right_data[which(sw_right_data$Nt>0),]
Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_data.RelAbundance)), Year = sw_right_data.RelAbundance$Year, IA.obs = sw_right_data.RelAbundance$Nt, CV.IA.obs = rep(0, nrow(sw_right_data.RelAbundance))) #Using 0.2 as a proxy


catch_list_old <- list(sw_right_data[which(sw_right_data$Year < 1678),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1678 & sw_right_data$Year < 1901),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1901 & sw_right_data$Year <= 1973),1:3],
                   sw_right_data[which(sw_right_data$Year > 1973),1:3])

sir_old_base <- list()
for(i in 1:2){
  sir_old_base[[i]] <-  StateSpaceSIR(
    file_name = "NULL",
    n_resamples = 10000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             add_VAR_IA = make_prior(rinvgamma,2,0.5),
                             q_IA = make_prior(rnorm, 0.5, 0.4),
                             Pmsy = make_prior(0.6)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.5 , 0.03), 
      make_prior(rnorm, 1.0185, 0.0028),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 0,
    output.Yrs = c(2021, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = FALSE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list_old,
    control = sir_control(threshold = 0.1 * 1e-61, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_old_base[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_old_base[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_old_base, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_old_base[[1]],  file_name = file_name)
plot_trajectory(sir_old_base[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_old_base[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_old_base[[2]]), inc_reference = FALSE)
plot_ioa(sir_old_base[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_old_base[[1]],  file_name = file_name)



