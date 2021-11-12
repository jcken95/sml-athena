rm(list=ls())
source("../GPfunctionsOptim.R")
source("../hetGPfunctions.R")
dat <- readRDS("../sa-dat/hgp-dat.RDS")
p0 <- dim(dat$x)[2]
Nsim <- 100
Nsobol <- 10^5
V.mat <- matrix(nrow=Nsim, ncol=p0)
tot.var.vec <- rep(0, Nsim)
S.T.eps.vec <- rep(0, Nsim)
pars <- dat$pars
gp.predict <- function(prior.mean,crosscov, E){
  as.vector(prior.mean) + as.vector(crosscov%*%E)
}
covmat <- cov.x1.x2(dat$x, dat$x, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(dat$x)[1])
precmat <- chol2inv(chol(covmat))

first.order <- function(N, dat, Beta, precmat, xrange = c(-1.7,1.7)){
  p <- dim(dat$x)[2]
  E.lambda <- precmat%*%(pars$logLambda - as.vector(dat$designmat.lambda%*%Beta))
  
  M <- lhs::randomLHS(N, p)%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  M.prime <- lhs::randomLHS(N, p)%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  crosscov.M <- cov.x1.x2(M, dat$x, pars$v_sigma, pars$v_theta, 2)
  mu.M <- cbind(1, M)%*%Beta
  y.hat <- as.vector( gp.predict(mu.M, crosscov.M, E.lambda))## predict from GP using design M
  
  E.y <- mean(y.hat)
  Var.y <- mean(y.hat^2)-E.y^2#(sum(y.hat^2)/(N-1)) - E.y^2##epistemic uncertainty in mean response
  total.var <- Var.y + pars$v_nugget
  ## now evaluate the sensitivy indices
  W <- rep(0,p)
  for(j in 1:p){
    M2 <- M.prime
    M2[,j] <- M[,j]
    crosscov.M2 <- cov.x1.x2(M2, dat$x, pars$v_sigma, pars$v_theta, 2)
    mu.M2 <- cbind(1, M2)%*%Beta
    y2 <- as.vector( gp.predict(mu.M2, crosscov.M2, E.lambda) )## predict from GP using design M2
    W[j] <- sum(y.hat * y2)/(N)
  }
  res <- list()
  res$V <- W-E.y^2
  res$total.var <- total.var
  res
}

b.lambda.post <- list()
b.lambda.post$mean <- dat$b.prior.lambda + dat$B.lambda%*%t(dat$designmat.lambda) %*% dat$precmat.lambda %*% (dat$pars$logLambda - as.vector(dat$designmat.lambda%*%dat$b.prior.lambda))
b.lambda.post$var <- dat$B.lambda - dat$B.lambda%*%t(dat$designmat.lambda) %*% dat$precmat.lambda %*% t(dat$B.lambda%*%t(dat$designmat.lambda))
k=0
while(k < Nsim){
  
  Beta <- MASS::mvrnorm(1, mu = b.lambda.post$mean, Sigma = b.lambda.post$var)
  ## test
  ##first.order does vectorised computation
  ## run mult times with same Beta to avoid memory overload
  tmprun <- first.order(Nsobol, dat, Beta,precmat, xrange = c(-1.7,1.7))
  k <- k+1
  V.mat[k,] <- tmprun$V
  tot.var.vec[k] <- tmprun$total.var
  print(k)
  
}
res <- list()
res$V <- V.mat
res$total.var <- tot.var.vec
res$S.T.eps <- S.T.eps.vec
saveRDS("hgp-var.RDS")
