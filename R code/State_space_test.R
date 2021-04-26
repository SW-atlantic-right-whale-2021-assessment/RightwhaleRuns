library(StateSpaceSIR)

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


Abs.Abundance.2009 <- data.frame(Year = 2009, N.obs = 4029, CV.obs = NA) # FIXME: not used as of 4/24/21
Abs.Abundance.2010 <- data.frame(Year = 2010, N.obs = 4245, CV.obs = 245/4245) # 2010: 4245 (SE: 245, 95% CI 3,765, 4,725). 


sw_right_data.RelAbundance <- sw_right_data[which(sw_right_data$Nt>0),]
Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_data.RelAbundance)), Year = sw_right_data.RelAbundance$Year, IA.obs = sw_right_data.RelAbundance$Nt, CV.IA.obs = rep(0, nrow(sw_right_data.RelAbundance))) #Using 0.2 as a proxy



file_name = "NULL"
n_resamples = 1000
priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                         N_obs = make_prior(runif, 100, 10000),
                         var_N = make_prior(rinvgamma, 4, 0.1),
                         z = make_prior(use = FALSE),
                         add_CV_IA = make_prior(rinvgamma, 2, 0.5),
                         Pmsy = make_prior(runif, 0.6, 0.8))
catch_multipliers = make_multiplier_list(
  make_prior(1),
  make_prior(rnorm, 1.60 , 0.04), 
  make_prior(rnorm, 1.09, 0.04),
  make_prior(1))
target.Yr = 2019
num.haplotypes = 10
output.Yrs = c(2021, 2030)
abs.abundance = Abs.Abundance.2010
rel.abundance = Rel.Abundance.SWRight
rel.abundance.key = TRUE # Indices of abundance
count.data = Count.Data # Not used
count.data.key = FALSE # Don't use count data
growth.rate.obs = c(0.074, 0.033, FALSE) # Do not include growth rate
growth.rate.Yrs = c(1995, 1996, 1997, 1998) # Not used
catch.data = catch_list
control = sir_control(threshold = 0.06 * 1e-67, progress_bar = TRUE)
realized_prior = ifelse(i == 1, "FALSE", "TRUE")
  begin.time <- Sys.time()
  
  ################################
  # Assigning variables
  ################################
  target.Yr <- target.Yr
  ## Use the first year of the projection is set as the first year in the
  ## catch series
  start_yr <- min(catch.data[[1]]$Year, na.rm = TRUE)
  for(i in 2:length(catch.data)){
    start_yr <- min(c(start_yr, catch.data[[i]]$Year), na.rm = TRUE)
  }
  
  ## The last year of the projection is set as the last year in the catch or
  ## abundance series, whichever is most recent
  end_yr <- max(max(abs.abundance$Year),
                max(rel.abundance$Year),
                output.Yrs)
  for(i in 1:length(catch.data)){
    end_yr <- max(c(end_yr, catch.data[[i]]$Year), na.rm = TRUE)
  }
  
  ## Setting the target year for the bisection method
  bisection.Yrs <- target.Yr-start_yr + 1
  ## Setting the years to project
  projection.Yrs <- end_yr-start_yr + 1
  Year <- seq(start_yr, end_yr, by = 1)
  
  # Expand the catch time series and fill missing values with 0
  catch_min <- list()
  catch_dif <- list()
  catch_max <- list()
  for(i in 1:length(catch.data)){
    catch.data[[i]] <- merge(data.frame(Year), catch.data[[i]], by="Year", all = TRUE)
    catch.data[[i]][is.na(catch.data[[i]])] <- 0
    
    # Get min and max of catch
    catch_original <- as.matrix(
      catch.data[[i]][,grep("catch", colnames(catch.data[[i]]), ignore.case = T)])
    catch_min[[i]] <- apply(catch_original, 1, FUN=min)
    catch_max[[i]] <- apply(catch_original, 1, FUN=max)
    catch_dif[[i]] <- catch_max[[i]] - catch_min[[i]]
  }
  
  ## Assigning the catch data
  n_catch_series <- length(catch.data)
  
  # Catch multiplier and parameters check
  if(length(catch_multipliers) != n_catch_series){
    stop("Number of catch multipliers (",
         length(catch_multipliers),
         ") does not equal number of catch periods (",
         n_catch_series,
         ")")
  }
  
  ## Determining the number of Indices of Abundance available
  num.IA <- max(rel.abundance$Index)
  
  ## Determining the number of Count Data sets available
  num.Count <- max(count.data$Index)
  
  ## Computing the value of sigma as in Zerbini et al. 2011
  rel.abundance$Sigma <- sqrt(log(1 + rel.abundance$CV.IA.obs^2))
  
  ## Computing the value of sigma for the count data as in Zerbini et al. (2011)
  count.data$Sigma <- sqrt(log(1 + count.data$CV.IA.obs^2))
  
  ## Computing the value of sigma as in Zerbini et al. 2011
  abs.abundance$Sigma <- sqrt(log(1 + abs.abundance$CV.obs^2))
  
  ## Computing the minimum viable population, if num.haplotypes=0, assumes no MVP
  MVP <- 3 * num.haplotypes
  
  
  ## Function to calculate Z if Pmsy is used
  NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z }
  if(priors$z$use & priors$Pmsy$use){
    warning("Priors were set on both Pmsy and Z, using the prior on Z")
  }
  
  ## Sample from prior for `z` or Pmsy (usually constant) if constant
  if (priors$z$use) {
    if(priors$z$class == "constant"){
      sample.z <- priors$z$rfn()
      sample.Pmsy <- uniroot(NmsyKz,z=sample.z,lower=0,upper=1)$root
    }
  } else {
    if(priors$Pmsy$class == "constant"){
      sample.Pmsy <- priors$Pmsy$rfn()
      sample.z <- uniroot(NmsyKz,NmsyK=sample.Pmsy,lower=1,upper=100)$root
    }
  }
  
  ## Start the loop
  i <- 0
  ## Keep track of number of draws
  draw <- 1
  Cumulative.Likelihood <- 0
  
  #Creating output vectors
  #-------------------------------------
  sir_names <- c("r_max", "K", "var_N", "z", "Pmsy", paste0("catch_multiplier_", 1:length(catch_multipliers)) , "catch_parameter",
                 "sample.N.obs", "add_CV", "add_CV_IA","Nmin", "YearMin",
                 "violate_MVP", paste0("N", target.Yr), paste0("N", output.Yrs),
                 paste0("ROI_IA", unique(rel.abundance$Index)),
                 paste0("q_IA", unique(rel.abundance$Index)),
                 paste0("ROI_Count", unique(count.data$Index)),
                 paste0("q_Count", unique(count.data$Index)),
                 "NLL.IAs", "NLL.Count", "NLL.N", "NLL.GR", "NLL", "Likelihood",
                 "Max_Dep",paste0("status", target.Yr), paste("status", output.Yrs, sep = ""), "draw", "save")
  
  proc_error_save <- matrix(NA, nrow = n_resamples, ncol = projection.Yrs-1)
  resamples_output <- matrix(NA, nrow = n_resamples, ncol = length(sir_names))
  resamples_trajectories <- matrix(NA, nrow = n_resamples, ncol = projection.Yrs)
  catch_trajectories <- matrix(NA, nrow = n_resamples, ncol = projection.Yrs)
  colnames(catch_trajectories) =  paste0("Catch_", Year)
  
  if (control$progress_bar) {
    pb <- txtProgressBar(min = 0, max = n_resamples, style = 3)
  }
  
  #Initiating the SIR loop
  while (i < n_resamples) {
    #Sampling from Priors
    #-------------------------------
    save <- FALSE #variable to indicate whether a specific draw is kept
    
    # Sampling for catch_multiplier and high-low sample
    sample.catch_multipliers <- sapply(catch_multipliers, function(x) x$rfn())
    sample.catch_parameter <- priors$catch_sample$rfn()
    
    catches <- rep(0, length(Year))
    
    for(p in 1:length(catch.data)){
      catches <- catches + (catch_min[[p]] + catch_dif[[p]] * sample.catch_parameter) * sample.catch_multipliers[p]
    }
    
    #Sampling for r_max
    sample.r_max <- priors$r_max$rfn()
    while (sample.r_max > 0.118) {
      sample.r_max <- priors$r_max$rfn()
    }
    
    ## Sampling from the N.obs prior
    sample.N.obs <- priors$N_obs$rfn()
    
    ## Prior on additional CV
    if (priors$add_CV$use) {
      sample.add_CV <- priors$add_CV$rfn()
    } else {
      sample.add_CV <- 0
    }
    
    if (priors$add_CV_IA$use) {
      sample.add_CV_IA <- priors$add_CV_IA$rfn()
    } else {
      sample.add_CV_IA <- 0
    }
    
    ## Sample from prior for variance of process error
    sample.var_N <- priors$var_N$rfn()
    sample.proc.error <- rlnorm(projection.Yrs-1, 0, sqrt(sample.var_N)) # Random process error
    
    ## Sample from prior for `z` or Pmsy (usually constant) if random
    if (priors$z$use) {
      if(priors$z$class == "function"){
        sample.z <- priors$z$rfn()
        sample.Pmsy <- uniroot(NmsyKz,z=sample.z,lower=0,upper=1)$root
      }
    } else {
      if(priors$Pmsy$class == "function"){
        sample.Pmsy <- priors$Pmsy$rfn()
        sample.z <- uniroot(NmsyKz,NmsyK=sample.Pmsy,lower=1,upper=100)$root
      }
    }
    
    ## Sampling from q priors if q.prior is TRUE; priors on q for indices of
    ## abundance
    if (priors$q_IA$use) {
      q.sample.IA <- replicate(num.IA, priors$q_IA$rfn())
    } else {
      ## FIXME: -9999 is probably not a good sentinel value here; NA?
      q.sample.IA <- rep(-9999, length(unique(rel.abundance$Index)))
    }
    
    ##priors on q for count data
    if (priors$q_count$use) {
      q.sample.Count <- replicate(num.Count, priors$q_count$rfn())
    } else {
      ## FIXME: Sentinel -9999 again
      q.sample.Count <- rep(-9999, length(unique(count.data$Index)))
    }
    
    ## Conduct logistic bisection
    sample.K <-  try(LOGISTIC.BISECTION.K(K.low = control$K_bisect_lim[1],
                                          K.high = 1e8,
                                          r_max = sample.r_max,
                                          z = sample.z,
                                          num_Yrs = bisection.Yrs,
                                          start_yr = start_yr,
                                          target.Pop = sample.N.obs,
                                          catches = catches,
                                          proc_error = sample.proc.error,
                                          MVP = MVP,
                                          tol = control$K_bisect_tol),
                     silent = TRUE)
    
    ## If population is too variable because of process error, give error and set likelihood to 0
    if(class(sample.K) == "try-error"){
      sample.K = 999
      K.error = TRUE
    } else{
      sample.K = sample.K
      K.error = FALSE
    }
    
    #Computing the predicted abundances with the samples from the priors
    #----------------------------------------
    Pred_N <- GENERALIZED_LOGISTIC(r_max = sample.r_max,
                                   K = sample.K,
                                   N1 = sample.K,
                                   z = sample.z,
                                   start_yr = start_yr,
                                   num_Yrs = projection.Yrs,
                                   catches = catches,
                                   proc_error = sample.proc.error,
                                   MVP = MVP)
    
    
    #Computing the predicted ROI for the IAs and Count data, if applicable
    #----------------------------------------
    #For IAs
    if (rel.abundance.key) {
      Pred.ROI.IA <- COMPUTING.ROI(data = rel.abundance,
                                   Pred_N = Pred_N$Pred_N,
                                   start_yr = start_yr)
    } else {
      Pred.ROI.IA <- rep(0, num.IA)
    }
    
    #For Count Data
    if (count.data.key) {
      Pred.ROI.Count <- COMPUTING.ROI(data = count.data,
                                      Pred_N = Pred_N$Pred_N,
                                      start_yr = start_yr)
    } else {
      Pred.ROI.Count <- rep(0, num.Count)
    }
    
    #Calculate Analytical Qs if rel.abundance.key is TRUE
    #---------------------------------------------------------
    if (rel.abundance.key) {
      if (!priors$q_IA$use) {
        q.sample.IA <- CALC.ANALYTIC.Q(rel.abundance,
                                       Pred_N$Pred_N,
                                       start_yr,
                                       sample.add_CV_IA,
                                       num.IA)
      } else {
        q.sample.IA <- q.sample.IA
      }
    }
    
    
    ## Calculate Analytical Qs if count.data.key is TRUE
    ## (NOT USED YET - AZerbini, Feb 2013)
    if (count.data.key) {
      if (!priors$q_count$use) {
        q.sample.Count <- CALC.ANALYTIC.Q(count.data,
                                          Pred_N$Pred_N,
                                          start_yr,
                                          sample.add_CV,
                                          num.Count)
      } else {
        q.sample.Count <- q.sample.Count
      }
    }
    
    if (control$verbose > 3) {
      message("r_max = ", sample.r_max,
              " N.obs = ", sample.N.obs,
              " K = ", sample.K,
              " Pred_N.target = ", Pred_N$Pred_N[bisection.Yrs],
              " q.IAs = ", q.sample.IA,
              " q.Count = ", q.sample.Count)
    }
    
    #Compute the likelihoods
    #--------------------------------
    # (1) relative indices (if rel.abundance.key is TRUE)
    if (rel.abundance.key) {
      lnlike.IAs <- LNLIKE.IAs(rel.abundance,
                               Pred_N$Pred_N,
                               start_yr,
                               q.sample.IA,
                               sample.add_CV_IA,
                               TRUE)
    } else {
      lnlike.IAs <- 0
    }
    
    # (2) count data (if count.data.key is TRUE)
    if (count.data.key) {
      lnlike.Count <- LNLIKE.IAs(count.data,
                                 Pred_N$Pred_N,
                                 start_yr,
                                 q.sample.Count,
                                 sample.add_CV,
                                 log=TRUE)
    } else {
      lnlike.Count <- 0
    }
    
    # (3) absolute abundance
    lnlike.Ns <- LNLIKE.Ns(abs.abundance,
                           Pred_N$Pred_N,
                           start_yr,
                           sample.add_CV,
                           log=TRUE)
    
    # (4) growth rate if applicable
    if (growth.rate.obs[3]) {
      Pred.GR <- PRED.GROWTH.RATE(growth.rate.Yrs=growth.rate.Yrs,
                                  Pred_N=Pred_N$Pred_N,
                                  start_yr=start_yr)
      lnlike.GR <- LNLIKE.GR(Obs.GR=growth.rate.obs[1],
                             Pred.GR=Pred.GR,
                             GR.SD.Obs=growth.rate.obs[2])
    } else {
      lnlike.GR <- 0
    }
    
    if (control$verbose > 2) {
      message("lnlike.IAs = ", lnlike.IAs,
              " lnlike.Count = ", lnlike.Count,
              " lnlike.Ns = ", lnlike.Ns,
              " lnlike.GR = ", lnlike.GR)
    }
    
    ## These use the likelihoods in Zerbini et al. (2011)
    NLL <- lnlike.IAs[[1]] + lnlike.Count[[1]] + lnlike.Ns[[1]] + lnlike.GR[[1]]
    Likelihood <- exp(-NLL)
    if (control$verbose > 1) {
      message("NLL = ", NLL,
              " Likelihood = ", Likelihood)
    }
    
    
    ## If population fell below minimum viable population size, set likelihood to 0
    if (Pred_N$Violate_Min_Viable_Pop) {
      Likelihood <- 0
      if (control$verbose > 0) {
        message("MVP violated on draw", draw)
      }
    }
    
    ## If population was too variable because of process error, set likelihood to 0
    if (K.error) {
      Likelihood <- 0
      if (control$verbose > 0) {
        message("Population dynamics too variable on draw", draw)
      }
    }
    
    
    Cumulative.Likelihood <- Cumulative.Likelihood + Likelihood
    
    # Trick to just extract realized prior
    if(realized_prior){
      Cumulative.Likelihood <- 2 * control$threshold
    }
    
    if (!Pred_N$Violate_Min_Viable_Pop) {
      
      while (Cumulative.Likelihood > control$threshold & i < n_resamples) {
        if (control$verbose > 0) {
          message("sample = ", i, " draw = ", draw)
        }
        if (control$verbose > 1) {
          message("draw = ", draw,
                  " Likelihood = ", Likelihood,
                  " Cumulative = ", Cumulative.Likelihood)
        }
        save <- TRUE
        Cumulative.Likelihood <- Cumulative.Likelihood-control$threshold
        resamples_trajectories[i+1,] <- Pred_N$Pred_N
        catch_trajectories[i+1,] <- catches
        proc_error_save[i+1,] <- sample.proc.error
        resamples_output[i+1,] <- c(sample.r_max,
                                    sample.K,
                                    sample.var_N,
                                    sample.z,
                                    sample.Pmsy,
                                    sample.catch_multipliers,
                                    sample.catch_parameter,
                                    sample.N.obs,
                                    sample.add_CV,
                                    sample.add_CV_IA,
                                    Pred_N$Min_Pop,
                                    ifelse(length(Pred_N$Min_Yr) == 1, Pred_N$Min_Yr, "Multiple"),
                                    Pred_N$Violate_Min_Viable_Pop,
                                    c(Pred_N$Pred_N[target.Yr - start_yr + 1]),
                                    c(Pred_N$Pred_N[output.Yrs - start_yr + 1]),
                                    Pred.ROI.IA,
                                    q.sample.IA,
                                    Pred.ROI.Count,
                                    q.sample.Count,
                                    lnlike.IAs[[1]],
                                    lnlike.Count[[1]],
                                    lnlike.Ns[[1]],
                                    lnlike.GR[[1]],
                                    NLL,
                                    Likelihood,
                                    Pred_N$Min_Pop / sample.K,
                                    c(Pred_N$Pred_N[target.Yr - start_yr + 1] /
                                        sample.K),
                                    c(Pred_N$Pred_N[output.Yrs - start_yr + 1] /
                                        sample.K),
                                    draw,
                                    save)
        i <- i+1
        if (control$progress_bar) {
          setTxtProgressBar(pb, i)
        }
      }
    }
    draw <- draw+1
  }
  
  # Save outputs
  resamples_output <- data.frame(resamples_output)
  resamples_output[] <- lapply(resamples_output, function(x) as.numeric(as.character(x)))
  names(resamples_output) <- sir_names
  if(!is.null(file_name)){
    write.csv(resamples_output,
              paste0(file_name, "_", "resamples_output.csv"))
  }
  
  resamples_trajectories <- data.frame(resamples_trajectories)
  resamples_trajectories[] <- lapply(resamples_trajectories, function(x) as.numeric(as.character(x)))
  names(resamples_trajectories) <- paste0("N_", Year)
  if(!is.null(file_name)){
    write.csv(resamples_trajectories,
              paste0(file_name, "_", "resamples_trajectories.csv"))
  }
  
  catch_trajectories <- data.frame(catch_trajectories)
  catch_trajectories[] <- lapply(catch_trajectories, function(x) as.numeric(as.character(x)))
  names(catch_trajectories) <- paste0("Catch_", Year)
  if(!is.null(file_name)){
    write.csv(catch_trajectories,
              paste0(file_name, "_", "catch_trajectories.csv"))
  }
  
  proc_error_save <- data.frame(proc_error_save)
  proc_error_save[] <- lapply(proc_error_save, function(x) as.numeric(as.character(x)))
  names(proc_error_save) <- paste0("Proc_error_", Year[1:(projection.Yrs-1)])
  if(!is.null(file_name)){
    write.csv(proc_error_save,
              paste0(file_name, "_", "proc_error.csv"))
  }
  
  resamples.per.samples <- draw / n_resamples
  if(resamples.per.samples < 3){
    warning("Number of resamples per sample is ",
            round(resamples.per.samples, 1),
            ", use higher threshold value.")
  } else if (resamples.per.samples > 20) {
    warning("Number of resamples per sample is ",
            round(resamples.per.samples, 1),
            ", use lower threshold value.")
  }
  
  end.time <- Sys.time()
  if (control$verbose > 0) {
    message("Time to Compute = ", (end.time-begin.time))
  }
  
  return_list <- list(call = call,
                      file_name = file_name,
                      Date.Time = Sys.time(),
                      Time.to.compute.in.minutes = paste((end.time-begin.time) / 60),
                      threshold = control$threshold,
                      Ratio.Resamples.per.Sample = paste("1 resample",
                                                         ":",
                                                         resamples.per.samples,
                                                         "samples"),
                      resamples_output = resamples_output,
                      resamples_trajectories = resamples_trajectories,
                      catch_trajectories = catch_trajectories,
                      inputs = list(draws = draw,
                                    n_resamples = n_resamples,
                                    prior_r_max = priors$r_max,
                                    catch_multipliers = catch_multipliers,
                                    priors_N_obs = priors$N_obs,
                                    target.Yr = target.Yr,
                                    start_yr = start_yr,
                                    MVP = paste("num.haplotypes = ",
                                                num.haplotypes,
                                                "MVP = ",
                                                3 * num.haplotypes),
                                    tolerance = control$K_bisect_tol,
                                    output.Years = output.Yrs,
                                    abs.abundance = abs.abundance,
                                    catch.data = catch.data,
                                    realized_prior = realized_prior))
  if(rel.abundance.key){ return_list$inputs$rel.abundance = rel.abundance}
  if(count.data.key){ return_list$inputs$count.data = count.data}
  
  class(return_list) <- "SIR" # Defines class for object
  