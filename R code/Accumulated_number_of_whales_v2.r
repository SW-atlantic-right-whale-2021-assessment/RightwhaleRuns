################################################################################
# Load libraries
################################################################################
rm(list=ls())
require(MASS)

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
regresion.balle.nb.jul.cuad<- glm.nb(RTA ~ Year + Juliano + I(Juliano^2), data = balle, link = log)
summary (regresion.balle.nb.jul.cuad)
nb.Jul.cuad <- cbind(Estimate = coef(regresion.balle.nb.jul.cuad))
nb.Jul.cuad # Estimates


################################################################################
# Second Stage - accumulated number of whales
################################################################################
#' Function to calculate the accumulated number of whales using previously estimated regression coefficients
#'
#' @param a Parameter from log-link GLM representing the intercept
#' @param b Parameter from log-link GLM representing the effect of year
#' @param c Parameter from log-link GLM representing the effect of day
#' @param d Parameter from log-link GLM representing the effect of day^2
#' @param mu Mean of normal distribution representing time whale remains in area
#' @param sigma Standard deviation of normal distribution representing time whale remains in area
#' @param year Year for calculating accumulated number of whales 
#' @param x Day of year to calculated accumulated number of whales: Assuming 320 because after 320 no new whales come in
#'
#' @return
#' @export
#'
#' @example 
#' accum_fun(year = 2014, x = 320)
accum_fun <- function(a = -77.41, 
                      b =  0.03234, 
                      c = 0.1605, 
                      d = -0.0003323, 
                      mu = 60, 
                      sigma = sqrt(60), 
                      year = 2014, 
                      x = 320){
  
  t <- 1:365 # Julian day
  W_t <- c(0, exp(a + b * year + c * t + d * t^2)) # Estimated number of whales on day t: W_0 = 0
  P_t <- dnorm(t, mu, sigma, FALSE) # Probability whale remains in area
  dW_t <- c(W_t[2:366] - W_t[1:365]) # delta w - 2:366 in the first term because W_t goes from day 0 to day 365
  A_x <- c(0) # Accumulated number of whales A_0 = 0
  
  # Loop through K
  # -- Should be as follows:
  # -- A_0 = 0
  # -- A_1 = dW_1 + p1 * A_0
  # -- A_2 = dW_2 + p1 * A_1 + p2 * A_0
  # -- A_3 = dW_3 + p1 * A_2 + p2 * A_1 + p3 * A_0
  for(k in 1:x){
    A_x[k+1] <- sum(dW_t[1:k]) + sum(P_t[1:k] * A_x[k:1])
  }
  
  return(A_x[x+1])
}
accum_fun(year = 2014, x = 113) # Example using parameter values from excel file matches 5.047502868 on dat 113 from 1 CONTEO BALLENAS.xls

# -- Run for years of interest and using nb regression model
Years <- 1998:2019 # Years which we want to calculate the accumulated number of whales
# Run function across years of interest
A_xy <- sapply(Years, function(x) accum_fun(a = nb.Jul.cuad[1], # Intercept
                                            b = nb.Jul.cuad[2], # Year parameter
                                            c = nb.Jul.cuad[3], # Julian day parameter
                                            d = nb.Jul.cuad[4], # Julian day^2 parameter
                                            year = x, # Year
                                            mu = 60, # mu from manuscript
                                            sigma = 8.66, # sigma from manuscript
                                            x = 320) # Calculate until day 320
               )
A_xy <- data.frame(Year = Years, A_x = A_xy) # Look good?
A_xy





