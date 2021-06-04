library(StateSpaceSIR)



################################################################################
# Read in data
################################################################################
# -- Catch
sw_right_data<-read.delim("Data/datosModeloBallenasmiles2020Miles1648to2019.csv", sep=";",header=FALSE)   
names(sw_right_data)<- c("Year","CatchMin","CatchMax","Nt")

# Four periods of SLRs
# - Period 1: 1648-1770: SLR = 1
# - Period 2: 1771-1850: SLR ~ N(1.6, 0.04)
# - Period 3: 1851-1973: SLR ~ N(1.09, 0.04)
# - Period 4: 1974-Present: SLR = 1
catch_list <- list(sw_right_data[which(sw_right_data$Year < 1771),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1771 & sw_right_data$Year <= 1850),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1851 & sw_right_data$Year <= 1973),1:3],
                   sw_right_data[which(sw_right_data$Year > 1973),1:3])

# -- Absolute abundance
Abs.Abundance.2009 <- data.frame(Year = 2009, N.obs = 4029, CV.obs = NA) # FIXME: not used as of 4/24/21
Abs.Abundance.2010 <- data.frame(Year = 2010, N.obs = 4245, CV.obs = 245/4245) # 2010: 4245 (SE: 245, 95% CI 3,765, 4,725).

# -- Relative abundance
sw_right_rel_abundance<-read.csv("Data/Accumulated_n_whales_1999_to_2019.csv") 

Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_rel_abundance)), 
                                    Year = sw_right_rel_abundance$Year, 
                                    IA.obs = sw_right_rel_abundance$A_xy_mu_sim) #Using 0.2 as a proxy
Rel.Abundance.SWRight = cbind(Rel.Abundance.SWRight, sw_right_rel_abundance[,paste0("X",1:17)])

for(i in 1:8){
  dir.create(paste0("Model runs/Sensitivity_",i))
}

################################################################################
# Base model
################################################################################
file_name <- "Model runs/Base/Base"

sir_base <- list()
for(i in 1:2){
  sir_base[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_base[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_base[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_base, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base[[1]],  file_name = file_name)
plot_trajectory(sir_base[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_base[[1]]),  file_name = file_name,   priors = list(sir_base[[2]]), inc_reference = FALSE)
#plot_ioa(sir_base[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_base[[1]],  file_name = file_name)




################################################################################
# Sensitivity_1
################################################################################
# - Prior on rmax of lognormal(-2.65, 0.5)T(0, 0.11)
# - Using rlnormTrunc fron "EnvStats" package

file_name <- "Model runs/Sensitivity_1/Sensitivity_1"

sensitivity_1 <- list()
for(i in 1:2){
  sensitivity_1[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(rlnormTrunc, -2.65, 0.5, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_1[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_1[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_1, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_1[[1]],  file_name = file_name)
plot_trajectory(sensitivity_1[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_1[[1]]),  file_name = file_name,   priors = list(sensitivity_1[[2]]), inc_reference = FALSE)
#plot_ioa(sensitivity_1[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_1[[1]],  file_name = file_name)



################################################################################
# Sensitivity_2
################################################################################
file_name <- "Model runs/Sensitivity_2/Sensitivity_2"
# Invgamma prior on sigma (process error 2, 0.05)

sensitivity_2 <- list()
for(i in 1:2){
  sensitivity_2[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 2, 0.05),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_2[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_2[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_2, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_2[[1]],  file_name = file_name)
plot_trajectory(sensitivity_2[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_2[[1]]),  file_name = file_name,   priors = list(sensitivity_2[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_2[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_2[[1]],  file_name = file_name)



################################################################################
# Sensitivity_3
################################################################################
file_name <- "Model runs/Sensitivity_3/Sensitivity_3"
# Var_N (sigma) ~ invgamma(8, 0.2)

sensitivity_3 <- list()
for(i in 1:2){
  sensitivity_3[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 8, 0.2),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_3[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_3[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_3, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_3[[1]],  file_name = file_name)
plot_trajectory(sensitivity_3[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_3[[1]]),  file_name = file_name,   priors = list(sensitivity_3[[2]]), inc_reference = FALSE)
#plot_ioa(sensitivity_3[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_3[[1]],  file_name = file_name)



################################################################################
# Sensitivity_4
################################################################################
file_name <- "Model runs/Sensitivity_4/Sensitivity_4"
# No struck and loss rates

sensitivity_4 <- list()
for(i in 1:2){
  sensitivity_4[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(1), 
      make_prior(1),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_4[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_4[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_4, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_4[[1]],  file_name = file_name)
plot_trajectory(sensitivity_4[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_4[[1]]),  file_name = file_name,   priors = list(sensitivity_4[[2]]), inc_reference = FALSE)
#plot_ioa(sensitivity_4[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_4[[1]],  file_name = file_name)



################################################################################
# Sensitivity_5
################################################################################
file_name <- "Model runs/Sensitivity_5/Sensitivity_5"
# - Catch time series is only low

sensitivity_5 <- list()
for(i in 1:2){
  sensitivity_5[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             catch_sample = make_prior(0) # Set to 0 used low catch only
                             ),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_5[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_5[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_5, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_5[[1]],  file_name = file_name)
plot_trajectory(sensitivity_5[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_5[[1]]),  file_name = file_name,   priors = list(sensitivity_5[[2]]), inc_reference = FALSE)
#plot_ioa(sensitivity_5[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_5[[1]],  file_name = file_name)



################################################################################
# Sensitivity_6
################################################################################
file_name <- "Model runs/Sensitivity_6/Sensitivity_6"
# -- Catch time series is high

sensitivity_6 <- list()
for(i in 1:2){
  sensitivity_6[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             catch_sample = make_prior(1) # Set to 1 used high catch only
                             ),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_6[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_6[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_6, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_6[[1]],  file_name = file_name)
plot_trajectory(sensitivity_6[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_6[[1]]),  file_name = file_name,   priors = list(sensitivity_6[[2]]), inc_reference = FALSE)
#plot_ioa(sensitivity_6[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_6[[1]],  file_name = file_name)



################################################################################
# Sensitivity_7
################################################################################
file_name <- "Model runs/Sensitivity_7/Sensitivity_7"
# Include additional relative abundance information

sensitivity_7 <- list()
for(i in 1:2){
  sensitivity_7[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(rinvgamma, 4, 0.1),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold = 0.5 * 1e-2, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_7[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_7[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_7, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_7[[1]],  file_name = file_name)
plot_trajectory(sensitivity_7[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_7[[1]]),  file_name = file_name,   priors = list(sensitivity_7[[2]]), inc_reference = FALSE)
#plot_ioa(sensitivity_7[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_7[[1]],  file_name = file_name)



################################################################################
# Sensitivity_8 no process error
################################################################################
file_name <- "Model runs/Sensitivity_8/Sensitivity_8"

sensitivity_8 <- list()
for(i in 1:2){
  sensitivity_8[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(0),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
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
    control = sir_control(threshold =  1e-6, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sensitivity_8[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_8[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_8, file = paste0(file_name, ".Rdata"))



load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_8[[1]],  file_name = file_name)
plot_trajectory(sensitivity_8[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_8[[1]]),  file_name = file_name,   priors = list(sensitivity_8[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_8[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sensitivity_8[[1]],  file_name = file_name)

