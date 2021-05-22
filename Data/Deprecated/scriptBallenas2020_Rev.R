library(openxlsx)
library("xlsx")
library("R2WinBUGS")
library("coda")
library("lattice")
library("mcmcplots")

datos1<-read.delim("datosModeloBallenasmiles2020.csv", sep=";",header=FALSE)   
names(datos1)<- c("Y","Htmin","Htmax","Nt")

############## model with density-dependence base-case##########

sink("ballenas.txt")  
cat("
    model {
    
    ######observation´s error######
    for (i in 1:N) { 
    Lmean[i] <- log(Z[i]*X)  
    L[i] ~ dlnorm(Lmean[i],itau3)
    
    ####replicated data
    L.rep[i] ~ dlnorm(Lmean[i],itau3) 
    }
    
    ######process´s error###### 
    for(j in 1:30) {        
    Zmean[j] <- 0   
    Z[j] ~ dlnorm(Zmean[j],isigma3)
    }
    
    for(i in 31:124) {          
    C[i-1]<- Cmin[i-1]+pi*(Cmax[i-1]- Cmin[i-1])
    Zmean[i] <- log(max(Z[i-1] + rl*Z[i-1]*(1-pow(Z[i-1],tita))-Kl*C[i-1],0.01)) 
    Z[i] ~ dlnorm(Zmean[i],isigma3)I(0,2)
    }
    
    for(i in 125:253) {          
    C[i-1]<- Cmin[i-1]+pi*(Cmax[i-1]- Cmin[i-1])
    Zmean[i] <- log(max(Z[i-1] + rl*Z[i-1]*(1-pow(Z[i-1],tita))-Kl*C[i-1]*SRL1,0.01)) 
    Z[i] ~ dlnorm(Zmean[i],isigma3)I(0,2)
    }
    
    for(i in 254:N) {
    C[i-1]<- Cmin[i-1]+pi*(Cmax[i-1]- Cmin[i-1])
    Zmean[i] <- log(max(Z[i-1] + rl*Z[i-1]*(1-pow(Z[i-1],tita))-Kl*C[i-1]*SRL2,0.01)) 
    Z[i] ~ dlnorm(Zmean[i],isigma3)I(0,2)
    }
    
    #### Priors on pi###########
    pi~ dunif(0,1)
    
    #### Priors on SRL1###########
    SRL1 ~ dnorm(1.5,33)
    
    #### Priors on SRL2###########
    SRL2 ~ dnorm(1.0185,357)
    
    #### Priors on tita (shape parameter)###########
    tita<- 2.39
    
    #######Priors on rl#####
    rl ~ dunif(0.01,0.11)
    
    #######Priors on X#####
    f~dnorm(0.5,0.25)I(0,1)
    X <- f*kl    
    
    #######Priors on kl#####
    kl ~ dunif(10,60)
    Kl<-1/kl
    
    #######Priors on isigma3 and itau3#####
    a1<-4;  b1<-0.1
    isigma3 ~ dgamma(a1,b1)
    Sigma3<-1/isigma3;
    
    c1<-2; d1<-0.446
    itau3 ~ dgamma(c1,d1)
    Tau3<-1/itau3
    
    
    #####Prediction: M-year extension to process equation#####
    for(i in (N+1):(N+5)) {
    Zmean[i] <- log(max(Z[i-1] + rl*Z[i-1]*(1-pow(Z[i-1],tita)),0.01))
    Z[i] ~ dlnorm(Zmean[i],isigma3)
    Lmean[i] <- log(Z[i]*X)  
    L.new[i] ~ dlnorm(Lmean[i],itau3)
    predictedL.new[i]<-exp(Lmean[i])
    predictedBl.new[i]<-exp(Lmean[i])/f
    }
    
    
    ############derived quantities##############
    for (i in 1:N) { 
    predictedL[i]<-exp(Lmean[i])
    predictedBl[i]<-exp(Lmean[i])/f
    
    residuall[i] <- log(L[i])-Lmean[i]
    sq[i] <- pow(residuall[i], 2)      # Squared residuals
    
    # replicated derived quantities
    residuall.rep[i] <- log(L.rep[i])-Lmean[i] 
    sq.new[i] <- pow(residuall.rep[i], 2)  # Squared residuals for new data
    }
    
    fit <- sum(sq[])              # Sum of squared residuals for actual data set
    fit.new <- sum(sq.new[])      # Sum of squared residuals for new data set
    test <- step(fit.new-fit) 		# Test whether new data set more extreme
    bpvalue <- mean(test) 		  	# Bayesian p-value
    
    }
    ",fill=TRUE)
sink()
# Bundle data

Cmin<-as.numeric(datos1$Htmin)

Cmax<-as.numeric(datos1$Htmax)

n<-length(Cmin)

N<-as.numeric(n)

L<- as.numeric(datos1$Nt)

win.data <- list("N","L","Cmin","Cmax")

# Inits function
inits <- function(){
  list( itau3=50, isigma3=50,f=0.4,rl=0.07,kl=30,SRL1=1.5,SRL2=1,pi=0.5)}  
  
#Parameters to monitor

params1 <- c("Sigma3","Tau3","f","rl","kl","SRL1","SRL2","pi","residuall","predictedL","predictedBl","L.rep","fit", "fit.new", "bpvalue","predictedL.new","predictedBl.new")

# MCMC settings

nc <- 3    # Number of chains
ni <- 1000000 # Number of draws from posterior for each chain
nb <- 10000  # Number of draws to discard as burn-in
nt <- 50    # Thinning rate

# Start Gibbs sampler
# 
out1 <- bugs(data = win.data, inits = inits, parameters = params1, model = "ballenas.txt", 
             n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE, working.directory = getwd())

print(out1, dig = 3)

save.image("C:/Users/romer/Documents/Resumenes-informes/2020/Paper modelo ballenas/SC - Revision/base_case.RData")

(pV<-var(out1$sims.list$deviance)/2)
(DICest<-out1$mean$deviance+pV)



###### Goodness-of-Fit assessment in Bayesian analyses####

datos1$res<-0
datos1$predictedL<-0

for (i in 1:length(datos1$Y)){
  if(is.na(datos1[i,"Nt"])) datos1[i,"res"]<-NA else datos1[i,"res"]<-out1$mean$residuall[i]
  if(is.na(datos1[i,"Nt"])) datos1[i,"predictedL"]<-NA else datos1[i,"predictedL"]<-out1$mean$predictedL[i]
}


plot(datos1$predictedL, datos1$res, las=1,xlab="Predicted values", ylab="Standardized Residuals ", pch=20, col="blue",ylim=c(-0.05,0.05))
abline(h=0)


par(mfrow=c(1, 2))

plot(datos1$predictedL, datos1$res, main="Residuals vs. predicted values", las=1, xlab="Predicted values", ylab="Residuals", pch=20, col="blue",ylim=c(-0.1,0.1))
abline(h=0)

plot(datos1$Nt~datos1$predictedL,ylab="Observed values",xlab="Predicted values")#,ylim=c(0,30),xlim=c(0,30))

##### Posterior Predictive Distributions and Bayesian p-Values ######

lim <- c(0, max(c(out1$sims.list$fit, out1$sims.list$fit.new)))
plot(out1$sims.list$fit, out1$sims.list$fit.new, main="Graphical posterior predictive check", las=1, xlab="SSQ for actual data set", ylab="SSQ for ideal (new) data sets", xlim=lim, ylim=lim, col="blue")
abline(0, 1)

mean(out1$sims.list$fit.new>out1$sims.list$fit)

################Convergence diagnostics with coda package#################

outCoda<-as.mcmc(out1)

outCoda<-as.mcmc.bugs(out1)

summary(outCoda)

summary(as.factor(outCoda))
str(outCoda)

xyplot(outCoda)#grafica los traceplot

qqmath(outCoda)

densplot(outCoda)

acfplot(outCoda)#grafica la autocorrelacion para todos los parametros

##paquete mcmcplot

mcmcplot(outCoda)
traplot(outCoda)#This function produces trace plots from an MCMC simulation on a single plot for all parameters
denplot(outCoda)#Creates a plot of densities for specified parameters from an MCMC simulation i
caterplot(outCoda)#Creates plots of credible intervals for parameters from an MCMC simulation

corplot(cor(as.matrix(outCoda)), cex.axis=0.75)#Creates an image plot of a correlation matrix where colors of different shades represent differing levels of correlation.

parcorplot(outCoda,col=heat.colors(13))#Creates an image plot of posterior correlations between model parameters from an MCMC
autplot1(outCoda)#Creates an autocorrelation or partial autocorrelation plot of MCMC output

rmeanplot(outCoda)

geweke.diag(outCoda, frac1=0.1, frac2=0.5)

geweke.plot(outCoda, frac1 = 0.1, frac2 = 0.5, nbins = 40,
            pvalue = 0.05, auto.layout = TRUE)

gelman.diag(outCoda, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

gelman.plot(outCoda, bin.width = 10, max.bins = 50,
            confidence = 0.95, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)

heidel.diag(outCoda, eps=0.1, pvalue=0.05)

autocorr.diag(outCoda)
crosscorr.plot(outCoda)

ParCor<-crosscorr(outCoda)
write.xlsx(ParCor, "ParCorBC.xlsx",col.names=TRUE,row.names=TRUE)


