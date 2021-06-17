library(StateSpaceSIR)
library(EnvStats)



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

for(i in 1:3){
  dir.create(paste0("Model runs/Test_",i))
}


################################################################################
# Test_1 model - parametric time-varying q
################################################################################
file_name <- "Model runs/Test_1/Test_1"

sir_test1 <- list()
for(i in 1:2){
  sir_test1[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 10000,
    priors = make_prior_list(r_max =  make_prior(rlnorm, log(0.06), 0.5),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             q_IA1 = make_prior(rlnorm, log(0.5), 0.4),
                             q_IA2 = make_prior(rnorm, 0, 0.5),
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
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_test1[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_test1[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_test1, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_test1[[1]],  file_name = file_name)
plot_trajectory(sir_test1[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_test1[[1]]),  file_name = file_name,   priors = list(sir_test1[[2]]), inc_reference = FALSE)
plot_ioa(sir_test1[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_test1[[1]],  file_name = file_name)





################################################################################
# Test_2 model - additional cv on IOA
################################################################################
file_name <- "Model runs/Test_2/Test_2"

sir_test2 <- list()
for(i in 1:2){
  sir_test2[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 10000,
    priors = make_prior_list(r_max =  make_prior(rlnorm, log(0.06), 0.5),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             add_VAR_IA = make_prior(rlnorm, log(0.2), 0.5),
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
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_test2[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_test2[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_test2, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_test2[[1]],  file_name = file_name)
plot_trajectory(sir_test2[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_test2[[1]]),  file_name = file_name,   priors = list(sir_test2[[2]]), inc_reference = FALSE)
plot_ioa(sir_test2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_test2[[1]],  file_name = file_name)




################################################################################
# Test_3 model - parametric time-varying q uninformative rmax
################################################################################
file_name <- "Model runs/Test_3/Test_3"

sir_test3 <- list()
for(i in 1:2){
  sir_test3[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 10000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             q_IA1 = make_prior(rlnorm, log(0.5), 0.4),
                             q_IA2 = make_prior(rnorm, 0, 0.5),
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
    control = sir_control(threshold = 1e-6, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_test3[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_test3[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_test3, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_test3[[1]],  file_name = file_name)
plot_trajectory(sir_test3[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_test3[[1]]),  file_name = file_name,   priors = list(sir_test3[[2]]), inc_reference = FALSE)
plot_ioa(sir_test3[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_test3[[1]],  file_name = file_name)
