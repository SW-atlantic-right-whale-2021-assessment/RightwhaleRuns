# Set up
# devtools::install_github("mcsiple/mmrefpoints")
library(mmrefpoints)
NSim <- 1000
Nyears <- 200

# Sample priors
MeanS0 <- exp(-rnorm(NSim, 0.18, 0.03)) # Cooke 2013
AnnualS0dev <- exp(rnorm(NSim, 0, 0.01)) # Cooke 2013
S0 <- MeanS0 * AnnualS0dev

MeanS1plus <- exp(-rnorm(NSim, 0.019, 0.005)) # Cooke 2003
AnnualS1plusdev <- exp(rnorm(NSim, 0, 0.001)) # Cooke 2013
S1plus <- MeanS1plus * AnnualS1plusdev

NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z } ## Function to calculate Z if Pmsy is used
sample.Pmsy <- runif(NSim, 0.5, 0.8)
sample.z <- sapply(sample.Pmsy, function(x) uniroot(NmsyKz,NmsyK=x,lower=1,upper=100)$root)

RightWhale <- data.frame(S0 = S0,
                         S1plus = S1plus,
                         lambdaMax = rnorm(NSim, 1.062, 0.003), # 6.2% (S.E. 0.2%) (Cooke 2013)
                         AgePart = rnorm(NSim, 8.4, 0.4), # 9.1 yr (SE 0.4) (Cooke 2013)
                         z = sample.z)

# Apply function to each simulated draw
RightWhale$Nproj <- apply(RightWhale, 1, function(x) dynamics(S0 = x[1], S1plus = x[2], K1plus = 9000, AgeMat = round(x[4])-1,
                                                              InitDepl = 0.9, ConstantCatch = NA, ConstantF = rep(0, times = Nyears), 
                                                              z = x[5], nyears = Nyears, nages = 25, lambdaMax = x[3])$TotalPop[Nyears]) 

# Results
hist(RightWhale$Nproj, xlab = "Numbers") # - funky shape because of Z prior
var(log(RightWhale$Nproj)) # Estimate variance of log total numbers - use as prior