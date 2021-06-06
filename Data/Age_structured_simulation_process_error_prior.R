# Set up
# devtools::install_github("mcsiple/mmrefpoints")
library(mmrefpoints)

#' Generate one marine mammal population trajectory (Modified to adjust age-at-maturity)
#'
#' This function generates one trajectory for a marine mammal population, starting at a user-specified depletion level \code{InitDepl}.
#'
#' @details
#' The population model is a single-sex age-structured model in which the number of calves or pups born each year is density dependent, with the extent of density dependence a function of the number of mature adults \eqn{\tildeN}, the fecundity (pregnancy rate) at pre-exploitation equilibrium \eqn{f_0}, the maximum theoretical fecundity rate fmax, the degree of compensation \eqn{z}, and the abundance of individuals aged 1+ \eqn{N_{y+1}^{1+}} relative to carrying capacity \eqn{K^{1+}}. This function can be used alone but is intended to be used with \code{Projections()} to generate multiple simulations. NOTE: Either \code{ConstantCatch} or \code{ConstantF} can be specified, but not both.
#'
#' @param S0 Calf/pup survival, a numeric value between 0 and 1
#' @param S1plus Survival for animals age 1 year and older, a numeric value between 0 and 1
#' @param K1plus The pre-exploitation population size of individuals aged 1 and older.  If this value is unavailable, it can be approximated by using the initial depletion and the estimate of current abundance
#' @param AgeMat Age at maturity in years (assumed to be age at first parturition - 1).
#' @param InitDepl Starting depletion level
#' @param ConstantCatch Total bycatch each year, expressed as a vector of length \code{nyears}
#' @param ConstantF vector (length = \code{nyears}) rate of bycatch each year
#' @param z The degree of compensation.  The default value is \code{z = 2.39}.
#' @param nyears Number of years to project
#' @param nages "Maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity +2 years without losing accuracy. Must be greater than AgeMat.
#' @param lambdaMax Maximum steady rate of increase (population growth rate)
#'
#' @return A list containing a matrix \code{N} of numbers at age (dimensions \code{nyears} (rows) x \code{nages} (columns)) and one vector \code{TotalPop} (a vector of length \code{nyears}), containing the number of age 1+ individuals in the population.

# Note, nages = Plus Group Age, and Plus Group Age can = AgeMat+2 without losing accuracy (per AEP 11/30/18)
#'
#' @examples
#' # Generate a time series of abundance for a bowhead whale
#' dynamics(S0 = 0.944, S1plus = 0.99, K1plus = 9000, AgeMat = 17,
#'  InitDepl = 0.6, ConstantCatch = NA, ConstantF = rep(0.01, times = 100), 
#'  z = 2.39, nyears = 100, nages = 25, lambdaMax = 1.04)
#' @export
dynamics <- function(S0, S1plus, K1plus, AgeMat, InitDepl, 
                     ConstantCatch = NA, ConstantF = NA, z, 
                     nyears, nages, lambdaMax) {
  # Checks
  if (length(ConstantCatch) > 1 & length(ConstantF) > 1) {
    stop("Cannot have both constant F and constant catch- choose one and set the other to NA!")
  }
  
  if(AgeMat > nages){stop("Age at maturity cannot be larger than plus group age. Change AgeMat or nages.")}
  if(S0 < 0 | S0 >= 1){stop("Calf/pup survival must be between 0 and 1.")}
  if(S1plus < 0 | S1plus >= 1){stop("Adult survival must be between 0 and 1.")}
  if(K1plus < 0){stop("Carrying capacity K1plus must be greater than zero.")}
  
  if (InitDepl > 1) {
    InitDepl <- 1
  }
  
  nyrs <- nyears + 1
  AgePart <- AgeMat + 1 # Age at first parturition = age at maturity +1 (~gestation period)
  AgePartVec <- rep(0, nages + 1)
  AgePartVec[(ceiling(AgePart):(nages+1))] <- 1
  AgePartVec[floor(AgePart) - 1] <- AgePart - floor(AgePart)
  
  Neq <- Ninit <- vector(length = (nages + 1))
  N <- C <- matrix(0, nrow = nyrs, ncol = (nages + 1))
  
  Tot1P <- rep(0, length = nyrs)
  Nrep <- rep(0, length = nyrs) # number of reproductive individuals
  
  
  NPROut <- npr(S0 = S0, S1plus = S1plus, nages = nages, AgeMat = AgeMat, E = 0)
  N0 <- NPROut$npr # mature nums per recruit
  P0 <- NPROut$P1r # 1+ nums per recruit
  Neq <- NPROut$nvec # 1+ nums per recruit
  
  f0 = 1/N0
  #f0 <- (1 - S1plus) / (S0 * (S1plus)^(AgeMat - 1)) # analytical soln for f0
  fmax <- getfecmax(lambdaMax = lambdaMax, S1plus = S1plus, S0 = S0, AgeMat = AgeMat)
  
  # Equilibrium conditions (need outside of if() statement to get R0)
  Neq[1] <- 1 # Age 0
  Neq[2] <- S0 # Age 1
  for (a in 3:nages) {
    Neq[a] <- Neq[a - 1] * S1plus
  } # Age 2+
  Neq[nages + 1] <- (S0 * S1plus^(nages - 1)) / (1 - S1plus) # plus group
  R0 <- K1plus / sum(Neq[2:(nages + 1)]) # numerical soln
  
  
  # Initial conditions, equilibrium
  if (InitDepl == 1) { # pop starts at equilibrium
    N[1, ] <- Neq * R0
    Tot1P[1] <- K1plus
    Nrep[1] <- sum(N[1, ] * AgePartVec)
    
    # Initial conditions, non-equilibrium
  } else { # pop starts at InitDepl*K
    
    E <- get_f(
      f.start = 0.5,
      S0.w = S0,
      S1plus.w = S1plus,
      nages.w = nages,
      AgeMat.w = AgeMat,
      InitDepl.w = InitDepl,
      z.w = z,
      lambdaMax.w = lambdaMax,
      N0.w = N0,
      P0.w = P0
    )
    
    Ninit[1] <- 1 # N_exploited; Age 0
    Ninit[2] <- S0 # Age 1
    for (a in 3:nages) {
      Ninit[a] <- S0 * (S1plus * (1 - E))^(a - 2)
    }
    Ninit[nages + 1] <- (S0 * (S1plus * (1 - E))^(nages - 1)) / (1 - (S1plus * (1 - E)))
    
    #-----
    A <- (fmax - f0) / f0
    # Fec0 <- 1.0 / N0
    # A <- (fmax - Fec0) / Fec0
    
    RF <- get_rf(E_in = E, S0 = S0, S1plus = S1plus, nages = nages, AgeMat = AgeMat, z = z, A = A, P0 = P0, N0 = N0)
    InitNumsAtAge <- Ninit * RF # Initial nums at age
    PropsAtAge <- InitNumsAtAge / sum(InitNumsAtAge) # Proportions at age
    
    n0 <- PropsAtAge[1] / sum(PropsAtAge[-1]) # Proportion of the pop that is age 0
    N0.fished <- InitDepl * K1plus * n0 # Number of age 0 individuals @ the start
    
    N[1, ] <- N0.fished * Ninit
    
    Tot1P[1] <- sum(N[1, 2:(nages + 1)])
    Nrep[1] <- sum(N[1, ] * AgePartVec)
  } # end initial conditions
  
  if (length(ConstantCatch) > 1) {
    for (Yr in 1:nyears) {
      MortE <- min(ConstantCatch[Yr] / Tot1P[Yr], 0.99) # bycatch mortality rate
      
      N[Yr + 1, 2] <- N[Yr, 1] * S0
      N[Yr + 1, 3:(nages + 1)] <- N[Yr, 2:nages] * (1 - MortE) * S1plus
      N[Yr + 1, (nages + 1)] <- (N[Yr, nages] + N[Yr, nages + 1]) * (1 - MortE) * S1plus
      Tot1P[Yr + 1] <- sum(N[Yr + 1, 2:(nages + 1)])
      Nrep[Yr + 1] <- sum(N[Yr + 1,] * AgePartVec)
      N[Yr + 1, 1] <- Nrep[Yr + 1] * (f0 + (fmax - f0) * (1 - (Tot1P[Yr + 1] / K1plus)^z))
    }
  } else {
    #sel <- 1
    for (Yr in 1:nyears) {
      MortE <- ConstantF[Yr] #* sel
      
      N[Yr + 1, 2] <- N[Yr, 1] * S0
      N[Yr + 1, 3:(nages + 1)] <- N[Yr, 2:nages] * (1 - MortE) * S1plus
      N[Yr + 1, (nages + 1)] <- (N[Yr, nages] + N[Yr, nages + 1]) * (1 - MortE) * S1plus
      Tot1P[Yr + 1] <- sum(N[Yr + 1, 2:(nages + 1)])
      Nrep[Yr + 1] <- sum(N[Yr + 1, ] * AgePartVec)
      N[Yr + 1, 1] <- Nrep[Yr + 1] * (f0 + (fmax - f0) * (1 - (Tot1P[Yr + 1] / K1plus)^z)) # rec
    }
  }
  N <- N[-nyrs, ]
  Tot1P <- Tot1P[-nyrs]
  return(list(TotalPop = Tot1P, N = N))
}


NSim <- 10000
Nyears <- 200

# Sample juvenile survival
MeanS0 <- exp(-rnorm(NSim, 0.179, 0.027)) # Cooke 2013
AnnualS0dev <- rlnorm(NSim, 0, 0.097) # Cooke 2013 - Annual deviates
S0 <- MeanS0 * AnnualS0dev
hist(S0)

# Sample adult survival
logistic <- function(x){1/(1+exp(-x))}
MeanS1plus <- exp(-rnorm(NSim, 0.026, 0.003)) # Cooke 2013
AnnualS1plusdev <- rnorm(NSim, 0, 0.19423) # Price et al 2017 sigma of annual deviation in adult survival on logit
S1plus <- logistic(logit(MeanS1plus) + AnnualS1plusdev)
hist(S1plus)

# Sample lambda max
MeanlambdaMax <- rnorm(NSim, 1.065, 0.002) # 6.5% (S.E. 0.2%) (Cooke 2013)
AnnuallambdaMax <- exp(rnorm(NSim, 0, 0)) # No annual deviation
lambdaMax <- MeanlambdaMax * AnnuallambdaMax
hist(lambdaMax)

# Sample z
NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z } ## Function to calculate Z if Pmsy is used
sample.Pmsy <- runif(NSim, 0.5, 0.8)
sample.z <- sapply(sample.Pmsy, function(x) uniroot(NmsyKz,NmsyK=x,lower=1,upper=100)$root)

# Combine
RightWhale <- data.frame(S0 = S0,
                         S1plus = S1plus,
                         lambdaMax = lambdaMax,
                         AgePart = rnorm(NSim, 8.4, 0.4), # 9.1 yr (SE 0.4) (Cooke 2013)
                         z = sample.z)
RightWhale <- RightWhale[which(RightWhale$S0 < 1 & RightWhale$S1plus < 1 & RightWhale$lambdaMax > 1),] # Make sure survival is less than 1

# Apply function to each simulated draw
RightWhale$Nproj <- apply(RightWhale, 1, function(x) dynamics(S0 = x[1], S1plus = x[2], K1plus = 9000, AgeMat = round(x[4])-1,
                                                              InitDepl = 0.9, ConstantCatch = NA, ConstantF = rep(0, times = Nyears), 
                                                              z = x[5], nyears = Nyears, nages = 25, lambdaMax = x[3])$TotalPop[Nyears]) 

# Results
hist(RightWhale$Nproj, xlab = "Numbers") # - funky shape because of Z prior
var(log(RightWhale$Nproj)) # Estimate variance of log total numbers - use as prior


### DEMOGRAPHIC PARAMATERS
## 1. Adult survival
# -- Females
# - 0.992 (0.989-0.994) (Macarena) Females
# - Survival 0.990 (95% CL 0.985, 0.996) for adult female  (Brandão 2010) - South Africa
# - Mortality of 0.019 (SE.0.005) (Cooke 2001)
# - Mortality of 0.02 (SE.0.004) (Cooke 2003)
# - Survival 0.986 (95% CL 0.976, 0.999) (Best 2001) - South Africa
# - Survival 0.990 (95% CL 0.983, 0.997) (Best 2005) - South Africa

# -- Males
# - 0.957 (0.941-0.969) (Macarena) Males

# -- Both
# - 0.951 (0.937-0.963) (Macarena) Unknown sex
# - 0.99 (Tulloch et al 2018) Both sexes
# - Fixed 0.9 to 0.99  (Carroll 2011)
# - Estimated 0.75 and 0.91 (Carroll 2011)

## 2. First year survival
# - Survival 0.713 (95% CL 0.529, 0.896) for juveniles (Brandão 2010) - South Africa
# - 0.070, 0.081 (Tulloch et al 2018)
# - Survival 0.913 (95% CL 0.601, 0.994) (Best 2001) - South Africa
# - Survival 0.734 (95% CL 0.518, 0.95) (Best 2005) - South Africa
# - Mortality rate mean of 0.18 (SE 0.03), interannual SD of 0.01, yearly estimates ranging from 0.08 to 0.37 (Cooke 2013)
# 0 0.69 to 0.923

## 3. Sex Ratio
# - 1:1 sex ratio observed (Carrol et al 2011)
# - 0.470 Females juveniles (Tulloch et al 2018)

## 4.  Age at maturity
# 6 years - Maturity is considered here as adult females 1 year past the age at first parturition (Tullock 2018)
# Age at first calving 9.1 yr (SE 0.3) (Cooke 2001)
# Age at first calving 9.1 yr (SE 0.4) (Cooke 2003)
# Age at first parturition 7.88 (95% CI 7.17, 9.29) (Best 2001) - South Africa
# Age at first parturition 7.69 (95% CI 7.06, 8.32) (Best 2005) - South Africa
# Age at first parturition 7.74 (95% CI 7.15, 8.33) (Brandao 2010) - South Africa

## 5. Reproductive rates
# 0.321 and 0.275 (Tulloch et al 2018) - Annual number of offspring per female
# Calving interval of 3.16 (95% CL 3.13, 3.19) (Brandao 2010) - South Africa
# Calving interval of 3.12 (95% CL 3.07, 3.17) (Best 2001) - South Africa
# Maximum calving interval of 5 years (Brandao 2010) - South Africa
# Calving interval of 3.6 (95% CL 3.3, 4.1) (Payne 1990) 
# Calving interval of 3.5 (95% CL 3.11, 3.18) (Best 2005) - South Africa
# Calving interval of 3.35 (SE 0.05) (Cooke 2001)
# Calving interval of 3.42 (SE 0.11) (Cooke 2003)

## 6. Population growth rate
# - 0.071 per year (95% CI 0.059, 0.082) (Best 2001) - South Africa
# - 0.073 per year (95% CI 0.066, 0.079) (Best 2005) - South Africa
# - 0.070 per year (95% CI 0.065, 0.075) (Brandao 2010) - South Africa
# - Annual percentage rate of population increase 6.9% (SE = 0.7%) (Cooke 2001)
# - Annual percentage rate of population increase 6.8% (S.E. 0.5%) (Cooke 2003)
# - Annual percentage rate of population increase 7.6% (S.E. 1.7%) (Payne 1990)
# - %6.2 (SE 0.2) (Cooke 2015)

### CITATIONS
# Bannister, J., 2001. Status of southern right whales (Eubalaena australis) off Australia. J. Cetacean Res. Manag. 1991, 103–110. doi:10.47536/jcrm.vi.273
# Best, P.B., Brandão, A., Butterworth, D.S., 2001. Demographic parameters of southern right whales off South Africa. J. Cetacean Res. Manag. 161–169. doi:10.47536/jcrm.vi.296
# Best, P. B., Brandão, A., & Butterworth, D. S. (2005) Updated estimates of demographic parameters for southern right whales off South Africa. International Whaling Commission document: SC/57/BRG2, 1–17.
# Brandao A, Best P, Butterworth D (2010) Estimates of demographic parameters for southern right whales off South Africa from survey data 1979 to 2006. Unpublished report (SC/62/BRG30) presented to the Scientific Committee of the International Whaling Commission, Cambridge, UK
# Carroll, E.L., Patenaude, N.J., Childerhouse, S.J., Kraus, S.D., Fewster, R.M., Baker, C.S., 2011. Abundance of the New Zealand subantarctic southern right whale population estimated from photo-identification and genotype mark-recapture. Mar. Biol. 158, 2565–2575. doi:10.1007/s00227-011-1757-9
# Cooke J, Rowntree V, Payne R (2001) Estimates of demographic parameters for southern right whales (Eubalaena australis) observed off Peninsula Valdes, Argentina. J Cetacean Res Manage Special Issue 2:125–132
# Cooke J, Rowntree V, Payne R (2003) Analysis of inter-annual variation in reproductive success of South Atlantic right whales (Eubalaena australis) from photo-identification of calving females observed off Peninsula Valde ́s. Unpublished report (SC/55/O23) presented to the Scientific Committee of the International Whaling Commission, Cambridge, UK
# Payne, R., Rowntree, V., Perkins, J. S., Cooke, J. G., & Lankester, K. (1990). Population size, trends and reproductive parameters of right whales (Eubalaena australis) off Peninsula Valdes, Argentina. Rep. int. Whal. Commn, 271-8.
#Tulloch, V.J.D., Plagányi, É.E., Matear, R., Brown, C.J., Richardson, A.J., 2018. Ecosystem modelling to quantify the impact of historical whaling on Southern Hemisphere baleen whales. Fish Fish. 19, 117–137. doi:10.1111/faf.12241