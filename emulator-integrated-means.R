## this is the one to use

library(rstan)
source("GPfunctionsOptim.R")
source("hetGPfunctions.R")
hetGP <- stan_model("hetGP.stan")
sml <- stan_model("het-SML.stan")
library("R.matlab")

## load data
#y data
y.test1 <- read.table(file = "emulator-data/y-val.txt",sep=",")[,1]
y.hetgp1<- read.table(file = "emulator-data/y-hetgp.txt",sep=",")[,1]
y.c1 <- read.table(file = "emulator-data/y-c.txt",sep=",")[,1]
y.e1 <- read.table(file = "emulator-data/y-e.txt",sep=",")[,1]
#x data
x.test <- as.matrix(read.table(file = "emulator-data/x-val.txt",sep=","))
x.hetgp <- as.matrix(as.matrix(read.table(file = "emulator-data/x-hetgp.txt",sep=",")))
x.c <- as.matrix(read.table(file = "emulator-data/x-c.txt",sep=","))
x.e <- as.matrix(read.table(file = "emulator-data/x-e.txt",sep=","))

probit <- function(p){
  qnorm(p)
}
inv.probit <- function(q){
  pnorm(q)
}

y.test <- probit(y.test1)
y.hetgp <- probit(y.hetgp1)
y.c <- probit(y.c1)
y.e <- probit(y.e1)

## plots
par(mfrow=c(3,3))
for(i in 1:9) plot(x.test[,i],y.test1,pch=20,cex=0.6)
for(i in 1:9) plot(x.hetgp[,i],y.hetgp1,pch=20,cex=0.6)
for(i in 1:9) plot(x.e[,i],y.e1,pch=20,cex=0.6)
for(i in 1:9) plot(x.c[,i],y.c1,pch=20,cex=0.6)

## define probit
## because availability in (0,1)



par(mfrow=c(3,3))
for(i in 1:9) plot(x.test[,i],y.test,pch=20,cex=0.6)
for(i in 1:9) plot(x.hetgp[,i],y.hetgp,pch=20,cex=0.6)
for(i in 1:9) plot(x.e[,i],y.e,pch=20,cex=0.6)
for(i in 1:9) plot(x.c[,i],y.c,pch=20,cex=0.6)

## define test metrics
MSE <- function(true, pred){
  mean((true-pred)^2)
}
Score <- function(y, m, v){
  # y = "true" value
  # m = emulator mean
  # v = emulator variance - including nugget term
  -((y-m)^2)/v - log(v)
}




## standardise inputs (x-mu)/sig

##training data
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.std <- scale(x.hetgp)

##scale validation data accordingly
x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))

#### fit emulators! ####

### first fit HetGP emulator ####
var.beta <- diag(1, 10)
mean.beta <- rep(0, 10)
data.hetGP <- list(
  m_p = 10, v_p = 10,
  N = length(y.hetgp), K = 9,
  x = x.std, 
  m_H = cbind(1, log(x.hetgp)), v_H = cbind(1, x.std),
  y = as.vector(y.hetgp),
  a = rep(1,length(y.hetgp)),
  ## prior
  m_beta_m = rep(0, 10), m_beta_s = var.beta,
  m_a_theta = rep(2,9), m_b_theta = rep(1,9),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget = 0,
  v_beta_m = rep(0, 10), v_beta_s = var.beta,
  v_a_theta = rep(2,9), v_b_theta = rep(1,9),
  v_a_sigma = 2, v_b_sigma = 2,
  v_nugget_a = 2, v_nugget_b = 2
)
temp <- list()

find.mode <- function(x){
  rstan::optimizing(hetGP, data = data.hetGP, verbose = F, as_vector = F)
}
st.hetgp <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 3)
en.hetgp <- Sys.time() - st.hetgp
en.hetgp
beepr::beep()
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))


het.fit <- temp[[best.emulator]]


pars <- het.fit$par


## predict response for hetgp
H <- cbind(1, log(x.hetgp))
H.new <- cbind(1, log(x.test))
covmat.hetgp <- H %*% var.beta %*% t(H) + cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(exp(pars$logLambda))
precmat.hetgp <- chol2inv(chol(covmat.hetgp))
x.v1 <- x.v1[1:100,1:9]
x.std <- x.std[1:100,1:9]
crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )
mean.hetgp <- crosscov %*% (precmat.hetgp %*% (y.hetgp - H %*% mean.beta))
var.hetgp <- H.new %*% var.beta %*% t(H.new) + cov.x1.x2(x.v1, x.v1, pars$m_sigma, pars$m_theta, 2)
var.hetgp <- var.hetgp - crosscov %*% precmat.hetgp %*% t(crosscov)
H.l <- cbind(1, x.std)
H.l.new <- cbind(1, x.v1)
covmat.hetgp.lambda <- H.l %*% var.beta %*% t(H.l) + cov.x1.x2(x.std, x.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.std)[1])
precmat.hetgp.lambda <- chol2inv(chol(covmat.hetgp.lambda))
crosscov.hetgp.lambda <- gp.crosscov(x.v1, x.std, pars$v_sigma, pars$v_theta,H.l.new, H.l, var.beta)
lambda.hetgp <- exp( crosscov.hetgp.lambda %*% precmat.hetgp.lambda %*%(pars$logLambda - as.vector(H.l%*%mean.beta) ))
var.hetgp.full <- var.hetgp + diag(lambda.hetgp)
plot(mean.hetgp, y.test)
mean((mean.hetgp - y.test)^2)
#### fit SML emulator ####

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std <- x.c.std[n.c:1,]
x.e.std <- x.e.std[n.ml:1,]

pairs(cbind(y.c2, x.c.std))
pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

m.h <- cbind(cbind(1, log(x.c[n.c:1,])), rbind(matrix(0, ncol = 10, nrow = n.c - n.e) ,cbind(1, log(x.e[n.e:1,]))))
tail(m.h)
var.beta <- diag(1, 20)
ml.data <- list(
  
  ## data ##
  
  m_p = 10, v_p = 10,
  N = n.e + n.c, K = 9,
  n_c = n.c, n_e = n.e,
  x_e = x.e.std, x_c = x.c.std,
  m_H = m.h, v_H = cbind(1, x.e.std),
  y_c = y.c2, y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, 10), m_beta_s = var.beta[1:10,1:10],
  m_a_theta = rep(2, 9), m_b_theta = rep(1, 9),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, 10), v_beta_s = var.beta[1:10,1:10],
  v_a_theta = rep(2, 9), v_b_theta = rep(1, 9),
  v_a_sigma = 2, v_b_sigma = 2,
  v_nugget_a = 2, v_nugget_b = 2,
  
  c_beta_m = rep(0, 10), c_beta_s = var.beta[1:10,1:10],
  c_a_theta = rep(2, 9), c_b_theta = rep(1, 9),
  c_a_sigma = 2, c_b_sigma = 2,
  c_nugget_a = 2, c_nugget_b = 2,
  m_rho = 1, s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 3)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par

designmat <- matrix(0, ncol=20, nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:10] <- cbind(1, log(x.c[n.c:1,]))
designmat[-(1:length(y.c)), (1:10)] <- pars$rho*cbind(1, log(x.e[n.e:1,]))
designmat[-(1:length(y.c)), -(1:10)] <- cbind(1, log(x.e[n.e:1,]))
dat.covmat <- dataCovmat(x.c.std, x.e.std, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat)); rm(dat.covmat)
H.c <- cbind(1, log(x.c[n.c:1,]))
H.e <- cbind(1, log(x.c[n.e:1,]))
H.ml.lambda <- cbind(1, x.e.std)
H <- cbind(1, log(x.test))
var.bc <- var.beta[1:10,1:10]; var.be <- var.beta[11:20,11:20]
crosscov.ml <- ml.crosscov(x.v2, x.c.std, x.e.std, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
covmat.ml.lambda <- H.ml.lambda %*% var.beta[1:10,1:10] %*% t(H.ml.lambda) + cov.x1.x2(x.e.std, x.e.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.e.std)[1])
precmat.ml.lambda <- chol2inv(chol(covmat.ml.lambda))
crosscov.ml.lambda <- gp.crosscov(x.v2, x.e.std, pars$v_sigma, pars$v_theta, cbind(1, x.v2), H.ml.lambda, var.beta[1:10,1:10])
lambda.ml <- exp( crosscov.ml.lambda %*%( precmat.ml.lambda %*% pars$logLambda) )
var.ml.full <- var.ml + diag(lambda.ml)

#### compare the fits ####

## score & MSE
#mse (probit)
MSE.hetgp <- MSE(y.test, mean.hetgp)
MSE.sml <- MSE(y.test, mean.ml)
MSE.hetgp; MSE.sml
## Rmse
sqrt(MSE.hetgp); sqrt(MSE.sml)
#mse (original scale)
MSE.hetgp0 <- MSE(y.test1, inv.probit(mean.hetgp))
MSE.sml0 <- MSE(y.test1, inv.probit(mean.ml))
sqrt(MSE.hetgp0); sqrt(MSE.sml0)
(MSE.hetgp0); (MSE.sml0)
#score
Score.hetgp <- Score(y.test, mean.hetgp, diag(var.hetgp.full))
Score.sml <- Score(y.test, mean.ml, diag(var.ml.full))
sum(Score.hetgp); sum(Score.sml)
## look at residuals
par(mfrow=c(1,2))

chol.het <- forwardsolve(t(chol(var.hetgp.full)), (y.test - mean.hetgp))
chol.ml <- forwardsolve(t(chol(var.ml.full)), (y.test - mean.ml))

## coverage plot


coverage <- function(resids){
  alpha <- seq(0, 1, length = length(resids))
  emp.cov <- rep(0, length(resids))
  for(i in 1:length(resids)){
    emp.cov[i] <- sum(abs(resids) < qnorm(1-alpha[i]/2))/length(resids)
  }
  list(y=emp.cov,x = 1 - alpha) 
}

par(mfrow=c(1,2), cex.axis=2, cex.lab=2,cex.main=2)
plot(coverage(chol.het), xlab = "Empirical Coverage", ylab="", main = "HetGP", pch = 20)
abline(0,1)
title(ylab = "Expected Covereage" , line = 2.4)
plot(coverage(chol.ml),  xlab = "Empirical Coverage", ylab="", main = "SML", pch = 20)
abline(0,1)
hist(chol.ml)
hist(chol.het)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(x.test[,i], xlab = i, chol.het); abline(h = c(-2,2))
}

par(mfrow=c(3,3))
for(i in 1:9){
  plot(x.test[,i], chol.ml,xlab = i); abline(h = c(-2,2))
}

par(mfrow=c(1,2),cex.axis=1, cex.lab=1,cex.main=1)
par(mar = c(5,2.5,3.5,0.05))
i=6
plot(x.test[,i], chol.het,pch=20,xlab = "Blade Lifetime", main  = "HetGP", ylim = c(-2.5,2.5)); abline(h = c(-1.96,1.96), col=2, lwd=1.3, lty=2)
plot(x.test[,i], chol.ml,pch=20,xlab = "Blade Lifetime",yaxt="n",ylab="",main="SML", ylim = c(-2.5,2.5)); abline(h = c(-1.96,1.96), col=2, lwd=1.3, lty=2)

## save some data ...

## save data ...
save.sml <- FALSE## change to save results

if(save.sml){
pars <- ml.fit$par
  x.v2 <- x.v2[1:100,1:9]
  
  designmat.new <- cbind( pars$rho*cbind(1, log(x.test)), cbind(1, log(x.test)) )
  
  crosscov <- ml.crosscov(x.v2,x.c.std, x.e.std, pars, cbind(1, log(x.c[n.c:1,])), cbind(1, log(x.e[n.e:1,])), cbind(1, log(x.test)), diag(1^2, 10), diag(1^2, 10))
  pred.mean <- function(y, y.hat, precmat, linear.fit, crosscov){
    ## y.hat is the "prior mean" of the observed data
    
    linear.fit + crosscov %*%( precmat %*% (y - y.hat))
    
  }
  pred.var <- function(prior.var, crosscov, precmat){
    prior.var - crosscov %*% precmat %*% t(crosscov)
  }
  var.beta <- diag(1, 20)
  designmat <- matrix(0, ncol=20, nrow=length(c(y.c, y.e)))
  designmat[1:length(y.c),1:10] <- cbind(1, log(x.c[n.c:1,]))
  designmat[-(1:length(y.c)), (1:10)] <- pars$rho*cbind(1, log(x.e[n.e:1,]))
  designmat[-(1:length(y.c)), -(1:10)] <- cbind(1, log(x.e[n.e:1,]))
  dat.covmat <- dataCovmat(x.c.std, x.e.std, pars, var.beta,  designmat)
  precmat <- chol2inv(chol(dat.covmat)) 

  
designmat.lambda <- cbind(1, x.e.std)
precmat.lambda <- cov.x1.x2(x.e.std, x.e.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, length(y.e))
precmat.lambda <- precmat.lambda + designmat.lambda%*%diag(1^2, 10)%*%(t(designmat.lambda))
precmat.lambda <- chol2inv(
  chol(precmat.lambda)
)
dat <- list()
dat$y.c <- y.c2
dat$y.e <- y.e2
dat$y <- c(dat$y.c,dat$y.e)
dat$pars <- pars
dat$b.prior <- rep(0, 20)
dat$B <- diag(1, 20)
dat$x.e <- x.e.std
dat$x.c <- x.c.std
dat$H.e <- m.h[-(1:(n.c-n.e)), 1:10]
dat$H.c <- m.h[,1:10]
dat$precmat <- precmat
dat$xCenter <- attr(scale(x.e), "scaled:center")
dat$xScale <- attr(scale(x.e), "scaled:scale")
dat$H <- designmat
dat$precmat.lambda <- precmat.lambda
dat$designmat.lambda <- cbind(1, x.e.std)
dat$B.lambda <- diag(1, 10)
dat$b.prior.lambda <- rep(0, 10)
designmat <- matrix(0, ncol=20, nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:10] <- cbind(1, log(x.c[n.c:1,]))
designmat[-(1:length(y.c)), (1:10)] <- pars$rho*cbind(1, log(x.e[n.e:1,]))
designmat[-(1:length(y.c)), -(1:10)] <- cbind(1, log(x.e[n.e:1,]))
dat$designmat <- designmat
saveRDS(dat,"sa-dat/ml-dat.RDS")
save.sml <- FALSE
}
save.hetgp <- FALSE
if(save.hetgp){
  dat.hgp <- list()
  pars <- het.fit$par
  dat.hgp$pars <- pars
  dat.hgp$y <- y.hetgp
  dat.hgp$x <- x.std
  dat.hgp$B <- diag(1, 10)
  dat.hgp$B.lambda <- diag(1, 10)
  dat.hgp$b.prior <- rep(0,10)
  dat.hgp$b.prior.lambda <- rep(0,10)
  dat.hgp$designmat <- cbind(1, log(x.hetgp))
  dat.hgp$designmat.lambda <- cbind(1, x.std)
  dat.hgp$xCenter<- colMeans(x.hetgp)
  dat.hgp$xScale <- apply(x.hetgp, 2, sd)
  dat.hgp$precmat <- precmat.hetgp
  dat.hgp$precmat.lambda <- precmat.hetgp.lambda
  saveRDS(dat.hgp,"sa-dat/hgp-dat.RDS")
  save.hetgp <- FALSE
}
