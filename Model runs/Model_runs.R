library(StateSpaceSIR)



################################################################################
# Read in data
################################################################################
# -- Catch
sw_right_data<-read.delim("Data/datosModeloBallenasmiles2020Miles1648to2019.csv", sep=";",header=FALSE)   
names(sw_right_data)<- c("Year","CatchMin","CatchMax","Nt")

# Four periods of SLRs
catch_list <- list(sw_right_data[which(sw_right_data$Year < 1771),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1771 & sw_right_data$Year <= 1850),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1851 & sw_right_data$Year <= 1973),1:3],
                   sw_right_data[which(sw_right_data$Year > 1973),1:3])

# -- Absolute abundance
Abs.Abundance.2009 <- data.frame(Year = 2009, N.obs = 4029, CV.obs = NA) # FIXME: not used as of 4/24/21
Abs.Abundance.2010 <- data.frame(Year = 2010, N.obs = 4245, CV.obs = 245/4245) # 2010: 4245 (SE: 245, 95% CI 3,765, 4,725).

# -- Relative abundance
sw_right_rel_abundance<-read.csv("Data/Accumulated_n_whales_1999_to_2019.csv") 

sw_right_data.RelAbundance <- sw_right_data[which(sw_right_data$Nt>0),]
Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_data.RelAbundance)), Year = sw_right_data.RelAbundance$Year, IA.obs = sw_right_data.RelAbundance$Nt, CV.IA.obs = rep(0, nrow(sw_right_data.RelAbundance))) #Using 0.2 as a proxy



################################################################################
# Run models
################################################################################
file_name <- "Model runs/Base/Base"
sw_right_data.RelAbundance <- sw_right_data[which(sw_right_data$Nt>0),]
Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_data.RelAbundance)), Year = sw_right_data.RelAbundance$Year, IA.obs = sw_right_data.RelAbundance$Nt, CV.IA.obs = rep(0, nrow(sw_right_data.RelAbundance))) #Using 0.2 as a proxy

sir_base <- list()
for(i in 1:2){
  sir_base[[i]] <-  StateSpaceSIR(
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
    control = sir_control(threshold = 0.1 * 1e-63, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}
resample_summary_reference <- summary_sir(sir_base[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_base[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
#save(sir_base, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base[[1]],  file_name = file_name)
plot_trajectory(sir_base[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_base[[1]]),  file_name = file_name,  lower = c(NA, NA, NA, NA, NA, 15000, NA, 24000, NA, NA, NA, NA, 0.5, 0.85), upper = c(NA, NA, 2000, NA, 20500, NA, NA, NA,  0.06, NA, NA, NA, 1, 1), priors = list(sir_base[[2]]), inc_reference = FALSE)
plot_ioa(sir_base[[1]],  file_name = file_name, ioa_names = c("FG", "BG1") )
summary_table(sir_base[[1]],  file_name = file_name)
