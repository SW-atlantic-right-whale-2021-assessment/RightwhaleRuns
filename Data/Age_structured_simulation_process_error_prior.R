# Set up
# devtools::install_github("mcsiple/mmrefpoints")
library(mmrefpoints)

#' Calculate numbers per recruit  (Modified to adjust age-at-maturity)
#'
#' @description Calculate numbers-per-recruit as a function of the bycatch rate assuming that 1+ animals are subject to bycatch  \emph{E}.
#'
#' @param S0 Calf/pup survival, a numeric value between 0 and 1
#' @param S1plus Adult survival, a numeric value between 0 and 1
#' @param nages Plus group age in years
#' @param AgeMat Age at maturity in years (must be equal to or less than nages)
#' @param E Bycatch mortality rate, a numeric value between 0 and 1
#'
#' @return A list of numbers per recruit (\code{npr}), 1+ numbers per recruit (\code{P1r}), and numbers at age per recruit (\code{nvec})
#'
#' @examples
#' (unpr <- npr(S0 = 0.9, S1plus = 0.9, 
#' AgeMat = 11, nages = 13, E = 0)) # unfished nums per recruit
#' (nprf <- npr(S0 = 0.9, S1plus = 0.9, 
#' AgeMat = 11, nages = 13, E = 0.8)) # nums per recruit at bycatch rate E
#' @export
npr <- function(S0, S1plus, nages, AgeMat, E = 0) {
  if(AgeMat > nages){warning("Age at maturity cannot be larger than plus group age. Change AgeMat or nages.")}
  if(S0 < 0 | S0 >= 1){stop("Calf/pup survival must be between 0 and 1.")}
  if(S1plus < 0 | S1plus >= 1){stop("Adult survival must be between 0 and 1.")}
  
  AgePart <- AgeMat + 1 # Age at first parturition = age at maturity +1 (~gestation period)
  AgePartVec <- rep(0, nages + 1)
  AgePartVec[(ceiling(AgePart+1):(nages+1))] <- 1
  AgePartVec[floor(AgePart)+1] <- AgePart - floor(AgePart)
  
  N.vec <- vector(length = nages + 1) # Ages 0 thru nages --> vector 1:(nages+1)
  N.vec[1] <- 1 # Age 0
  N.vec[2] <- 1 * S0
  
  # Survival (due to natural and bycatch causes)
  OnePlusSurv <- (S1plus * (1 - E))
  
  # Calculate numbers-at-age
  for (a in 3:(nages)) {
    N.vec[a] <- S0 * (OnePlusSurv^(a - 2))
  }
  
  N.vec[nages + 1] <- (S0 * OnePlusSurv^(nages - 1)) / (1 - OnePlusSurv) # plus group age
  npr <- sum(N.vec * AgePartVec) # *Reproducing* animals (AgeMat+1 = age at first parturition)
  P1r <- sum(N.vec[2:(nages + 1)]) # 1+ whales/pinnipeds
  Outs <- list()
  Outs$npr <- npr
  Outs$P1r <- P1r
  Outs$nvec <- N.vec
  return(Outs)
}


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
  AgePartVec[(ceiling(AgePart+1):(nages+1))] <- 1
  AgePartVec[floor(AgePart)+1] <- AgePart - floor(AgePart)
  
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
  return(list(TotalPop = Tot1P, N = N, fmax = fmax, Neq = Neq, N0 = N0, P0 = P0))
}

# Set up
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

# Sample age-at-first pregnancy
AgeMat <- rnorm(NSim, 7.58, 0.18)
hist(AgeMat)

# Sample z
NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z } ## Function to calculate Z if Pmsy is used
sample.Pmsy <- runif(NSim, 0.5, 0.8)
sample.z <- sapply(sample.Pmsy, function(x) uniroot(NmsyKz,NmsyK=x,lower=1,upper=100)$root)

# Combine
RightWhale <- data.frame(S0 = S0,
                         S1plus = S1plus,
                         lambdaMax = lambdaMax,
                         AgeMat = AgeMat, # (Cooke 2013)
                         z = sample.z)
RightWhale <- RightWhale[which(RightWhale$S0 < 1 & RightWhale$S1plus < 1 & RightWhale$lambdaMax > 1),] # Make sure survival is less than 1

# Apply function to each simulated draw
RightWhale$Nproj <- apply(RightWhale, 1, function(x) 
  mmrefpoints::dynamics(S0 = x[1], S1plus = x[2], K1plus = 9000, AgeMat = round(x[4]),
           InitDepl = 0.9, ConstantCatch = NA, ConstantF = rep(0, times = Nyears), 
           z = x[5], nyears = Nyears, nages = 25, lambdaMax = x[3])$TotalPop[Nyears])

# Results
hist(RightWhale$Nproj, xlab = "Numbers") # - funky shape because of Z prior
var(log(RightWhale$Nproj)) # Estimate variance of log total numbers - use as prior

# Citations
# Cooke, Rowntree & Sironi (2015) Cooke J, Rowntree V, Sironi M. Southwest Atlantic right whales: interim updated population assessment from photo-id collected at Península Valdéz, Argentina. SC/66/IWC Southern Right Whale Assessment Workshop; 2015. p. 9.
# Pace, R.M., Corkeron, P.J., Kraus, S.D., 2017. State-space mark-recapture estimates reveal a recent decline in abundance of North Atlantic right whales. Ecol. Evol. 7, 8730–8741. doi:10.1002/ece3.3406