################################################################################
# Load libraries
################################################################################
rm(list=ls())
require(MASS)
library(glmmTMB)

################################################################################
# Read in and clean data
################################################################################
balle<- read.csv ( "Data/ballena3.csv" , sep=";" , dec=".")


balle<- balle [c(-11,-12),] #11 y 12 (2004) vuelos con el aerocomander// esta es la seleccion de datos que hice para la los paràmetros para el nùmero de ballenas que dan la vuelta anualmente por PV
balle <- subset (balle, Year<2020) # Subset years

####Variable Respuesta ############
balle$RTA <- (balle$T) # Total of observed whales


################################################################################
# First stage - regression model
################################################################################
#	Regresion model selected (up to  2019)
regresion.balle.nb.jul.cuad <- glm.nb(RTA ~ as.factor(Year) + Juliano + I(Juliano^2), data = balle, link = log)
regresion.balle.glmmTMB.jul.cuad <- glmmTMB(RTA ~ as.factor(Year) + Juliano + I(Juliano^2), data = balle, family = nbinom2)

summary (regresion.balle.nb.jul.cuad)
nb.Jul.cuad <- cbind(Estimate = coef(regresion.balle.nb.jul.cuad))
nb.Jul.cuad # Estimates
vcov.nb.Jul.cuad <- vcov(regresion.balle.nb.jul.cuad) # Variance-covariance matrix
vcov.nb.Jul.cuad


################################################################################
# Second Stage - accumulated number of whales
################################################################################
# -- Run for years of interest and using nb regression model
A_xy <- data.frame(Year = sort(unique(balle$Year)), B = c(0, nb.Jul.cuad[2:17])) # Years which we want to calculate the accumulated number of whales from years which we have data and the associated regression parameters. NOTE: 1990 is the intercept 

# Run function across years of interest
A_xy$A_xy <- NA
for(i in 1:nrow(A_xy)){
  A_xy$A_xy[i] <- accum_fun(a = nb.Jul.cuad[1], # Intercept
                            b = A_xy$B[i], # Year parameter
                            c = nb.Jul.cuad[18], # Julian day parameter
                            d = nb.Jul.cuad[19], # Julian day^2 parameter
                            mu = 60, # mu from manuscript
                            sigma = 8.66, # sigma from manuscript
                            x = 320) # Calculate until day 320
  
}


A_xy




