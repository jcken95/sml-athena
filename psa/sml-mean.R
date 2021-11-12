## code up the SA using algorithm 8.4
## https://bookdown.org/rbg/surrogates/chap8.html#chap8ind
rm(list=ls())
source("../GPfunctionsOptim.R")
source("../hetGPfunctions.R")
dat <- readRDS("../sa-dat/ml-dat.RDS")
#source("../ml-functions-new.R")
Nsim <- 100
Nsobol <- 10^5
first.order <- function(N, dat, Betas, Precmats, xrange = c(-1.7,1.7)){
  
  p <- dim(dat$x.e)[2] ## number input params
  pars <- dat$pars
  
  x.c.orig <- dat$x.c%*%diag(dat$xScale) + matrix(rep(dat$xCenter,dim(dat$x.c)[1]), ncol=p, byrow=T)
  x.e.orig <- dat$x.e%*%diag(dat$xScale) + matrix(rep(dat$xCenter,dim(dat$x.e)[1]), ncol=p, byrow=T)
  mu.y <- as.vector(c(
    cbind(1, log(x.c.orig))%*%Betas$cheap, cbind(1, log(x.e.orig))%*%(pars$rho*Betas$cheap + Betas$exp)
  ))
  print( length((dat$y - mu.y)) )
  print(dim(Precmats$mean))
  mu.lambda <- as.vector(cbind(1, dat$x.e)%*%Betas$lambda)
  E <- Precmats$mean%*%(dat$y - mu.y)
  E.lambda <- Precmats$lambda%*%(pars$logLambda-mu.lambda)
  M <- diff(xrange)*lhs::randomLHS(N, p)  + min(xrange)
  print(range(M))
  M.prime <- diff(xrange)*lhs::randomLHS(N, p)  + min(xrange)
  M.orig <- M%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  M.prime.orig <- M.prime%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  mu.M <- cbind(1, log(M.orig)) %*% (pars$rho*Betas$cheap + Betas$exp)
  mu.M.prime <- cbind(1, log(M.prime.orig)) %*% (pars$rho*Betas$cheap + Betas$exp)  
  M.crosscov <- ml.crosscov.mle(M, dat$x.e, dat$x.c, pars)
  M.prime.crosscov <- ml.crosscov.mle(M, dat$x.e, dat$x.c, pars)
  y.hat <- as.vector( gp.predict(mu.M, M.crosscov, E) )## predict from GP using design M
  E.y <- mean(y.hat)
  Var.y <- sum(y.hat^2)/(N-1) - E.y^2##epistemic uncertainty in mean response
  print(Var.y)
  ## compute stochastic variation due to stochastic nature of simulations
  ## aleatory uncertainty
  lambda.crosscov <- cov.x1.x2(M, dat$x.e, pars$v_sigma, pars$v_theta,2)
  mu.lambda.new <-  as.vector(cbind(1, M)%*%Betas$lambda) 
  tmp <- lambda.crosscov%*%dat$chol.precmat.lambda
  log.var.mean <- as.vector( gp.predict(mu.lambda.new, lambda.crosscov, E.lambda) )## predict (log) GP using design M
  log.var.var <- diag(pars$v_sigma^2+pars$v_nugget, N) - tmp%*%t(tmp)
  E.lambdasq <- mean(exp(log.var.mean + 0.5*log.var.var))
  total.var <- E.lambdasq + Var.y
  print(E.lambdasq)
  
  ## now evaluate the sensitivity indices
  W <- rep(0,p)
  if(T){
    for(j in 1:p){
      
      M2 <- M.prime
      M2[,j] <- M[,j]
      M2.orig <- M2%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
      M2.crosscov <- ml.crosscov.mle(M2, dat$x.e, dat$x.c, pars)
      mu2 <- cbind(1, log(M2.orig))%*%(pars$rho*Betas$cheap + Betas$exp)
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
p0 <- dim(dat$x.e)[2]
V.mat <- matrix(nrow=Nsim, ncol=p0)
tot.var.vec <- rep(0, Nsim)
S.T.eps.vec <- rep(0, Nsim)
ml.crosscov.mle <- function(x.new, x.exp, x.cheap, pars){
  cov.cheap <- cov.x1.x2(x.new, x.cheap, pars$c_sigma*pars$rho, pars$c_theta, 2)
  cov.exp <- cov.x1.x2(x.new, x.exp, pars$m_sigma, pars$m_theta, 2)
  
  cbind(cov.cheap, cov.exp)
}
pars <- dat$pars
gp.predict <- function(prior.mean,crosscov, E){
  as.vector(prior.mean) + as.vector(crosscov%*%E)
}
Precmats <- list()
Betas <- list()
dat$chol.precmat.lambda <- chol2inv(
  chol(cov.x1.x2(dat$x.e, dat$x.e, pars$v_sigma, pars$v_theta, 2)) + diag(pars$v_nugget, length(dat$y.e))
) 
k = 0
while(k < Nsim){
  b.lambda.post <- list()
  b.lambda.post$mean <- dat$b.prior.lambda + dat$B.lambda%*%t(dat$designmat.lambda) %*% dat$precmat.lambda %*% (dat$pars$logLambda - as.vector(dat$designmat.lambda%*%dat$b.prior.lambda))
  b.lambda.post$var <- dat$B.lambda - dat$B.lambda%*%t(dat$designmat.lambda) %*% dat$precmat.lambda %*% t(dat$B.lambda%*%t(dat$designmat.lambda))
  Betas <- list()
  Betas$lambda <- MASS::mvrnorm(1, mu = b.lambda.post$mean, Sigma = b.lambda.post$var)
  
  
  
  
  log.lam.post <- list()
  C.lambda <- cov.x1.x2(dat$x.e,dat$x.e,pars$v_sigma, pars$v_theta)
  log.lam.precmat <- chol2inv(
    chol(C.lambda + diag(pars$v_nugget, dim(dat$x.e)[1]))
  )
  E.lambda <- log.lam.precmat%*%(pars$logLambda - as.vector(cbind(1, dat$x.e)%*%Betas$lambda))
  log.lam.mean <- cbind(1, dat$x.e)%*%Betas$lambda + C.lambda%*%E.lambda
  log.lam.var <- C.lambda - C.lambda%*% log.lam.precmat %*%t(C.lambda)
  log.lam <- MASS::mvrnorm(1, mu = log.lam.mean, Sigma = log.lam.var)
  
  ## draw values of b.lambda via an MVN generator
  
  ## now draw values of beta for the mean fn
  var.yc <- cov.x1.x2(dat$x.c, dat$x.c, pars$c_sigma, pars$c_theta,2) + diag(pars$c_nugget, dim(dat$x.c)[1])
  var.ye <- cov.x1.x2(dat$x.e, dat$x.e, pars$c_sigma * pars$rho, pars$c_theta,2)
  var.ye <- var.ye + cov.x1.x2(dat$x.e, dat$x.e, pars$m_sigma, pars$m_theta,2) + diag(exp(log.lam))
  cov.ye.yc <- cov.x1.x2(dat$x.e, dat$x.c, pars$c_sigma, pars$c_theta, 2)*pars$rho
  covmat <- rbind(cbind(var.yc, t(cov.ye.yc)), cbind(cov.ye.yc, var.ye)) + dat$designmat%*%dat$B%*%t(dat$designmat)
  
  Precmats.tmp <- chol2inv(chol(covmat))
  
  b.post <- list()
  
  b.post$mean <- dat$b.prior + dat$B %*% t(dat$designmat) %*% (Precmats.tmp ) %*% (dat$y - dat$designmat%*%dat$b.prior)
  b.post$var <- dat$B - dat$B %*% t(dat$designmat) %*% Precmats.tmp %*% t(dat$B %*% t(dat$designmat))
  
  beta.draw <- MASS::mvrnorm(n=1, mu = b.post$mean, Sigma = b.post$var)
  
  p <- dim(dat$x.c)[2]+1
  Betas$cheap <- beta.draw[1:p]
  Betas$exp <- beta.draw[-(1:p)]
  
  Precmats$lambda <- chol2inv(
    chol(
      cov.x1.x2(dat$x.e, dat$x.e, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(dat$x.e)[1])
    )
  )
  Precmats$mean <- chol2inv(chol(
    rbind(cbind(var.yc, t(cov.ye.yc)), cbind(cov.ye.yc, var.ye)) 
  ))
  
  ## test
  tmprun <- first.order(Nsobol, dat, Betas, Precmats, xrange = c(-1.7,1.7))

  k <- k+1
  V.mat[k,] <- tmprun$V
  tot.var.vec[k] <- tmprun$total.var
  S.T.eps.vec[k] <- tmprun$E.lambdasq
  print(k)

  
}

res <- list()
res$V <- V.mat
res$total.var <- tot.var.vec
res$S.T.eps <- S.T.eps.vec
saveRDS("sml-mean.RDS")
