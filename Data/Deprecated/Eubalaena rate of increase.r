require (ordinal)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require (vcd)
library(AICcmodavg) 
library(pgirmess)
library (MuMIn)
library(xlsx)



rm(list=ls())
ls()
 notebook<-  c( "C:\\Users\\Cosca\\Dropbox\\Papers en preparación\\SRW polpulation trajectory\\ballena3.csv") 

 
balle<- read.csv ( notebook , sep=";" , dec=".")
balle<- balle [c(-11,-12),] #11 y 12 (2004) vuelos con el aerocomander// esta es la seleccion de datos que hice para la los paràmetros para el nùmero de ballenas que dan la vuelta anualmente por PV

#balle<- balle [c(-11,-12,-31,-32,-33,-44, -47),] #11 (octubre 2004);12 (otro vuelo con aerocomander); 33 (24 septiembre 2008)
 #balle<- balle [c(-14,-15,-35, -41),] # 26 (16 de septiembre de 2007, se contó menos que en los vuelos anterior y posterior), con el 31 (13 agosto 2008) pasó lo mismo.

 
 
 

 
attach(balle)
head(balle)
tail(balle)
list(Fecha)

balle<- subset (balle, Año<2020) #Change to get the yearly rate  of increase
attach(balle)

  ####Variable Respuesta ############
  
  RTA <- (T) #Total of observed whales
 
  Mes<- (balle$Mes)#Months as factor
  Juliano<- (balle$Juliano) #Jualian day
  Año<- (balle$Año) #Year
  
  length(Año); length (RTA); length (Juliano)
  
  

#	Regresion model selected (up to  2019)
  regresion.balle.nb.jul.cuad<- glm.nb(RTA ~ Año + Juliano + I(Juliano^2))
  summary (regresion.balle.nb.jul.cuad)
  nb.Jul.cuad <- cbind(Estimate = coef(regresion.balle.nb.jul.cuad), confint(regresion.balle.nb.jul.cuad))
  nb.Jul.cuad
  exp(nb.Jul.cuad)
  
	
	
	






#attributes(modelos)$models

###Lista de los coeficientes######
lista.coef<-list(
      nb.sinaño.Jul.cuad,     #2
      nb.Jul.cuad,            #1
      nb.mes.cuad,            #3 
      pois.jul.cuad,          #4
      pois.sinaño.jul.cuad,   #5
      pois.mes.cuad,          #6
      nb.mes,                 #7
      pois.mes,               #8
      nb.Jul,                 #9
      pois.jul,               #10
      nb.año,                 #11
      pois.año                #12
)

#Lista de los coeficientes par ala tabla de comparaciñnd e modelos 

lista.coef<-list(
  nb.sinaño.Jul.cuad,     #2
  nb.Jul.cuad,            #1
  nb.mes.cuad,            #3 
  pois.jul.cuad,          #4
  pois.sinaño.jul.cuad,   #5
  pois.mes.cuad,          #6
  nb.mes,                 #7
  pois.mes,               #8
  nb.Jul,                 #9
  pois.jul,               #10
  nb.año,                 #11
  pois.año                #12
)

 plot(RTA ~ Juliano , col="red", xlab="Dña", ylab="Nº of whales", type="p") 
	lines ( fitted.values( regresion.balle.nb.jul.cuad)~ Juliano, type="p", col="blue",data=balle)

	
	
	
	
	
	
	
	
