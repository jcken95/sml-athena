## code up the SA using algorithm 8.4
## https://bookdown.org/rbg/surrogates/chap8.html#chap8ind
rm(list=ls())
source("../GPfunctionsOptim.R")
source("../hetGPfunctions.R")
dat <- readRDS("../sa-dat/hgp-dat.RDS")

#source("../ml-functions-new.R")
Nsim <- 100
Nsobol <- 10^5
first.order <- function(N, dat, Betas, Precmats, xrange = c(-1.7,1.7)){
  
  p <- dim(dat$x)[2] ## number input params
  pars <- dat$pars
  x.orig <- dat$x%*%diag(dat$xScale) + matrix(rep(dat$xCenter,dim(dat$x)[1]), ncol=p, byrow=T)
  mu.y <- as.vector(cbind(1, log(x.orig))%*%Betas$mean)
  mu.lambda <- as.vector(cbind(1, dat$x)%*%Betas$lambda)
  E <- Precmats$mean%*%(dat$y - mu.y)
  E.lambda <- Precmats$lambda%*%(pars$logLambda-mu.lambda)
  M <- diff(xrange)*lhs::randomLHS(N, p)  + min(xrange)
  M.prime <- diff(xrange)*lhs::randomLHS(N, p)  + min(xrange)
  M.orig <- M%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  M.prime.orig <- M.prime%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  mu.M <- cbind(1, log(M.orig)) %*% Betas$mean
  mu.M.prime <- cbind(1, log(M.prime.orig)) %*% Betas$mean  
  M.crosscov <- cov.x1.x2(M, dat$x, pars$m_sigma, pars$m_theta, 2)
  y.hat <- as.vector( gp.predict(mu.M, M.crosscov, E) )## predict from GP using design M
  E.y <- mean(y.hat)
  Var.y <- sum(y.hat^2)/(N-1) - E.y^2##epistemic uncertainty in mean response
  ## compute stochastic variation due to stochastic nature of simulations
  ## aleatory uncertainty
  lambda.crosscov <- cov.x1.x2(M, dat$x, pars$v_sigma, pars$v_theta,2)
  mu.lambda.new <-  as.vector(cbind(1, M)%*%Betas$lambda) 
  tmp <- lambda.crosscov%*% dat$chol.precmat.lambda
  log.var.diag <- diag(pars$v_sigma^2+pars$v_nugget,N) - diag( tmp%*%t(tmp) )
  E.lambdasq <- mean(exp(mu.lambda.new + 0.5 * log.var.diag))
  total.var <- E.lambdasq + Var.y
  print(E.lambdasq)
  
  ## now evaluate the sensitivity indices
  W <- rep(0,p)
  if(T){
    for(j in 1:p){
      M2 <- M.prime
      M2[,j] <- M[,j]
      M2.orig <- M2%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
      M2.crosscov <-cov.x1.x2(M2, dat$x, pars$m_sigma, pars$m_theta,2)
      mu2 <- cbind(1, log(M2.orig))%*%(Betas$mean)
      y2 <- as.vector( gp.predict(mu2, M2.crosscov,  E) )## predict from GP using design M2
      W[j] <- sum(y.hat * y2)/(N-1)
      }
  }
  
  
  V <- W-E.y^2
  
  res <- list()
  
  res$V <- V
  res$total.var <- total.var
  res$E.lambdasq<- E.lambdasq
  
  res
}
p0 <- dim(dat$x)[2]
V.mat <- matrix(nrow=Nsim, ncol=p0)
tot.var.vec <- rep(0, Nsim)
S.T.eps.vec <- rep(0, Nsim)
pars <- dat$pars
dat$chol.precmat.lambda <- chol2inv(
  chol(cov.x1.x2(dat$x, dat$x, pars$v_sigma, pars$v_theta, 2)) + diag(pars$v_nugget, length(dat$y))
) 
gp.predict <- function(prior.mean,crosscov, E){
  as.vector(prior.mean) + as.vector(crosscov%*%E)
}
Precmats <- list()
Betas <- list()

b.lambda.post <- list()
b.lambda.post$mean <- dat$b.prior.lambda + dat$B.lambda%*%t(dat$designmat.lambda) %*% dat$precmat.lambda %*% (dat$pars$logLambda - as.vector(dat$designmat.lambda%*%dat$b.prior.lambda))
b.lambda.post$var <- dat$B.lambda - dat$B.lambda%*%t(dat$designmat.lambda) %*% dat$precmat.lambda %*% t(dat$B.lambda%*%t(dat$designmat.lambda))
k = 0
while(k < Nsim){
  Betas <- list()
  Betas$lambda <- MASS::mvrnorm(1, mu = b.lambda.post$mean, Sigma = b.lambda.post$var)
  log.lam.post <- list()
  C.lambda <- cov.x1.x2(dat$x,dat$x,pars$v_sigma, pars$v_theta)
  log.lam.precmat <- chol2inv(
    chol(C.lambda + diag(pars$v_nugget, dim(dat$x)[1]))
  )
  E.lambda <- log.lam.precmat%*%(pars$logLambda - as.vector(cbind(1, dat$x.e)%*%Betas$lambda))
  log.lam.mean <- cbind(1, dat$x)%*%Betas$lambda + C.lambda%*%E.lambda
  log.lam.var <- C.lambda - C.lambda%*% log.lam.precmat %*%t(C.lambda)
  log.lam <- MASS::mvrnorm(1, mu = log.lam.mean, Sigma = log.lam.var)
  
  ## draw values of b.lambda via an MVN generator
  
  ## now draw values of beta for the mean fn

  covmat <- cov.x1.x2(dat$x, dat$x, pars$m_sigma, pars$m_theta,2) + diag(exp(log.lam))
  covmat <- covmat + dat$designmat%*%dat$B%*%t(dat$designmat)
  
  Precmats.tmp <- chol2inv(chol(covmat))
  
  b.post <- list()
  
  b.post$mean <- dat$b.prior + dat$B %*% t(dat$designmat) %*% (Precmats.tmp ) %*% (dat$y - dat$designmat%*%dat$b.prior)
  b.post$var <- dat$B - dat$B %*% t(dat$designmat) %*% Precmats.tmp %*% t(dat$B %*% t(dat$designmat))
  
  Betas$mean <- MASS::mvrnorm(n=1, mu = b.post$mean, Sigma = b.post$var)

  Precmats$lambda <- chol2inv(
    chol(
      cov.x1.x2(dat$x, dat$x, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(dat$x)[1])
    )
  )
  Precmats$mean <- chol2inv(chol(
    covmat
  ))
  
  ## test
  tmprun <- first.order(Nsobol, dat, Betas, Precmats, xrange = c(-1.7,1.7))
  
  k <- k+1
  V.mat[k,] <- tmprun$V
  tot.var.vec[k] <- tmprun$total.var
  S.T.eps.vec[k] <- tmprun$E.lambdasq
  print(k)
  
  
}
k=0
## draw betas from the posterior
res <- list()
res$V <- V.mat
res$total.var <- tot.var.vec
res$S.T.eps <- S.T.eps.vec
saveRDS("hgp-mean.RDS")
