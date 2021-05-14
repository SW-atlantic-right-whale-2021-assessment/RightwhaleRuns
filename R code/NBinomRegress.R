
library(TMB)
compile("NBinomRegress.cpp")
dyn.load(dynlib("NBinomRegress"))


X = model.matrix(regresion.balle.nb.jul.cuad)
data <- list(Y = balle$RTA, Xd=rep(1, length(balle$RTA)), X = X)
parameters <- list(betad=1, beta=rep(0, ncol(X)))
obj <- MakeADFun(data, parameters, DLL="NBinomRegress")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)