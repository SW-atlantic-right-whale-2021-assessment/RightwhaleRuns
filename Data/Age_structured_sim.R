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
# - Mortality rate mean of 0.18 (SE 0.03), interannual SD of 0.1, yearly estimates ranging from 0.08 to 0.37 (Cooke 2013)

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


library(mmrefpoints) # https://github.com/mcsiple/mmrefpoints

#' Code to produce estimates of numbers-at-age given prior parameter estimates and an age-structured model. Adapted from https://github.com/mcsiple/mmrefpoints
#'
#' @param S0 Age-0 to 1 survival
#' @param S1plus Age-1+ survival
#' @param lambdaMax Maximum rate of population increase
#' @param AgePart Age at first partruition
#' @param nages Number of ages (default = 20)
#' @param K1plus Carrying capacity of adults (defaults at 10,000)
#' @param InitDepl Initial depletion (defaults at 1)
#' @param z Shape parameter of logistic curve
#' @param nyrs Number of years to project
#'
#' @return Total numbers at age at equilibrium
#' @export
#' 
#' @examples 
#' example_results <- whale_fun( S0 = 0.71, S1plus = 0.99, lambdaMax = 1.07,  AgePart = 7.74, nages = 20, K1plus = 10000, InitDepl = 1, z = 2.39, nyrs = 1000)
whale_fun <- function( S0 = 0.71, 
                       S1plus = 0.99, 
                       lambdaMax = 1.07,  
                       AgePart = 7.74, 
                       nages = 20, 
                       K1plus = 10000, 
                       InitDepl = 1, 
                       z = 2.39, 
                       nyrs = 100){
  # Set up data
  AgeMat <- AgePart - 1 # Age at maturity = age at first partruition - 1 (~1 year gestation period)
  PartAtAge <- c(0, rep(0, nages))
  PartAtAge[(ceiling(AgePart):nages) + 1] <- 1
  PartAtAge[floor(AgePart) + 1] <- AgePart - floor(AgePart)
  Neq <- vector(length = (nages + 1)) # Numbers-at-age at equilibrium
  N <- C <- matrix(0, nrow = nyrs, ncol = (nages + 1)) # Numbers-at-age per year
  Tot1P <- rep(0, length = nyrs) # Total numbers per year
  Nrep <- rep(0, length = nyrs) # Numbers of reproductive individuals per year
  
  # Calculate numbers per-recruit
  NPROut <- npr(S0 = S0, S1plus = S1plus, nages = nages, AgeMat = AgeMat, E = 0)
  N0 <- NPROut$npr # mature nums per recruit
  P0 <- NPROut$P1r # 1+ nums per recruit
  Neq <- NPROut$nvec # 1+ nums per recruit
  
  # Calculate fecundity
  f0 = 1/N0 # Probability of reproducing this year
  fmax <- (lambdaMax^(AgeMat) - (S1plus * (lambdaMax^(AgeMat - 1)))) / (S0 * S1plus^(AgeMat - 1))
  
  # Equilibrium conditions (need outside of if() statement to get R0)
  Neq[1] <- 1 # Age 0
  Neq[2] <- S0 # Age 1
  for (a in 3:nages) {
    Neq[a] <- Neq[a - 1] * S1plus
  } # Age 2+
  Neq[nages + 1] <- (S0 * S1plus^(nages - 1)) / (1 - S1plus) # plus group
  R0 <- K1plus / sum(Neq[2:(nages + 1)]) # numerical soln
  
  # Project population
  # - Initial conditions
  PropsAtAge <- Neq / sum(Neq) # Proportions at age
  n0 <- PropsAtAge[1] / sum(PropsAtAge[-1]) # Proportion of the pop that is age 0
  N0.depl <- InitDepl * K1plus * n0 # Number of age 0 individuals @ the start
  N[1, ] <- N0.depl * Neq
  
  Tot1P[1] <- sum(N[1, 2:(nages + 1)])
  Nrep[1] <- sum(N[1, ] * PartAtAge)
  
  # - Project population through time
  for (Yr in 1:(nyrs-1)) {
    N[Yr + 1, 2] <- N[Yr, 1] * S0
    N[Yr + 1, 3:(nages + 1)] <- N[Yr, 2:nages] * S1plus
    N[Yr + 1, (nages + 1)] <- (N[Yr, nages] + N[Yr, nages + 1]) * S1plus
    Tot1P[Yr + 1] <- sum(N[Yr + 1, 2:(nages + 1)])
    Nrep[Yr + 1] <- sum(N[Yr + 1, ] * PartAtAge)
    N[Yr + 1, 1] <- Nrep[Yr + 1] * (f0 + (fmax - f0) * (1 - (Tot1P[Yr + 1] / K1plus)^z)) # rec
  }
  
  return(list(TotalN = Tot1P, N = N)) 
}


# Set up parameters
library(mmrefpoints) # https://github.com/mcsiple/mmrefpoints

ndraws <- 1000

NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z } ## Function to calculate Z if Pmsy is used
sample.Pmsy <- runif(ndraws, 0.5, 0.8)
sample.z <- c()
for(i in 1:ndraws){
  sample.z[i] <- uniroot(NmsyKz,NmsyK=sample.Pmsy[i],lower=1,upper=100)$root
}

whale_draws <- data.frame(S0 = runif(ndraws, 0.518, 0.95), 
                          S1plus = runif(ndraws, 0.937, 0.994),
                          lambdaMax = rnorm(ndraws, 1.062, 0.003), # 6.2% (S.E. 0.2%) (Cooke 2013)
                          AgePart = rnorm(ndraws, 8.4, 0.4), # 9.1 yr (SE 0.4) (Cooke 2013)
                          z = sample.z 
                          
)


whale_draws$Nproj <- NA
whale_draws$Nproj2 <- NA

nyears <- 100

# Run pop model for each draw draws
for(i in 1:ndraws){
  sol <- dynamics(S0 = whale_draws$S0[i], S1plus = whale_draws$S1plus[i], K1plus = 20000, AgeMat = round(whale_draws$AgePart[i] - 1),
                  InitDepl = 0.9, ConstantCatch = NA, ConstantF = rep(0, times = nyears), 
                  z = whale_draws$z[i], nyears = nyears, nages = 25, lambdaMax = whale_draws$lambdaMax[i])
  # K and initdepl should be arbitrary once the population reaches equilibriums

  whale_draws$Nproj[i] <- sol$TotalPop[nyears]
}

var(log(whale_draws$Nproj))
hist(whale_draws$Nproj)
