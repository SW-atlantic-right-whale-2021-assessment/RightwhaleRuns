# Set up
# devtools::install_github("mcsiple/mmrefpoints")
library(mmrefpoints)
NSim <- 1000
Nyears <- 200

# Sample juvenile survival
MeanS0 <- exp(-rnorm(NSim, 0.179, 0.027)) # Cooke 2013
AnnualS0dev <- rlnorm(NSim, 0, 0.097) # Cooke 2013 - Annual deviates
S0 <- MeanS0 * AnnualS0dev
hist(S0)

# Sample adult survival
MeanS1plus <- exp(-rnorm(NSim, 0.026, 0.003)) # Cooke 2013
AnnualS1plusdev <- exp(rnorm(NSim, 0, 0.01)) # This is a guess, probably less than juveniles  - Annual deviates
S1plus <- MeanS1plus * AnnualS1plusdev
hist(S1plus)

NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z } ## Function to calculate Z if Pmsy is used
sample.Pmsy <- runif(NSim, 0.5, 0.8)
sample.z <- sapply(sample.Pmsy, function(x) uniroot(NmsyKz,NmsyK=x,lower=1,upper=100)$root)

RightWhale <- data.frame(S0 = S0,
                         S1plus = S1plus,
                         lambdaMax = rnorm(NSim, 1.065, 0.002), # 6.5% (S.E. 0.2%) (Cooke 2013)
                         AgePart = rnorm(NSim, 8.4, 0.4), # 9.1 yr (SE 0.4) (Cooke 2013)
                         z = sample.z)

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