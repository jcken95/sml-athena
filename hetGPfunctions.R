## het GP functions


cond.mean2 <- function(y.obs,x.obs,x.new,sigma,omega,mu.Y,mu.Yobs,nugget=0,nu=2)
{
  
  ## computes the conditional mean at points x.new
  ## based on observations y.obs at points x.obs  
  
  ## assuming independent measurement error with variance = nugget^2 
  
  Sigma.Yobs.Yobs <- cov.x1.x2(x.obs,x.obs,sigma,omega,nu)
  Sigma.Y.Yobs <- cov.x1.x2(x.new,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs+diag(nugget)
  
  ## find solve(Sigma.Yobs.Yobs) carefully...
  
  R <- chol(Sigma.Yobs.Yobs) ##sqrt s.yobs.yobs
  
  R.inv <- backsolve(R, diag(1, dim(R)[1]))
  
  mu.Y+Sigma.Y.Yobs %*%  ((R.inv) %*% (t(R.inv) %*%(y.obs-mu.Yobs)))
}

cond.var.het <- function(y.obs.v, x.obs, x.new, sigma, omega, nugget, mu.Y.v, mu.Yobs.v, nu = 2){
  ## all I need to do it to "predict" a new variance ...

  ## (sigma, omega, nugget) are parameters of latent log-variance GP
  
  exp(cond.mean(as.vector(y.obs.v), x.obs, x.new, sigma, omega, mu.Y.v, as.vector(mu.Yobs.v), nugget^0.5 , nu = 2))
  
}

cond.covmat.het <- function(var.obs, var.new, x.obs, x.new, sigma, omega, nu=2){
  ## computes the conditional covariance matrix at points x.new
  ## based on observations y.obs at points x.obs  
  
  ## allows for non-constant variance in response
  
  ## predict the "new noise" with a seperate GP ...
  ## e.g. cond.var.het
  ## N.B. cond.var.het only predicts the aleatory uncertainty, does not account for additional epistemic uncertainty
  
  
  Sigma.Y.Y <- cov.x1.x2(x.new,x.new,sigma,omega,nu) + diag(as.vector(var.new))
  Sigma.Y.Yobs <- cov.x1.x2(x.new,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- cov.x1.x2(x.obs,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs+diag(as.vector(var.obs))
  
  L <- t(chol(Sigma.Yobs.Yobs)) ##sqrt s.yobs.yobs
  
  
  L.inv.S.O.Y <- forwardsolve(L, t(Sigma.Y.Yobs)) ## inv sqrt s.yobs.yobs
  
  
  
  Sigma.Y.Y - t(L.inv.S.O.Y)  %*% (L.inv.S.O.Y) ## syy - s.y.yobs * inv s.yobs.yobs * s.yobs.y
  
}


cond.mean.het <- function(y.obs, x.obs, x.new, sigma, omega, var.Y, nugget, mu.Y, mu.Yobs, nu = 2){
  
  
  Sigma.Y.Y <- diag(var.Y)+ diag(nugget, length(var.Y))
  
  Sigma.Y.Y <- Sigma.Y.Y + cov.x1.x2(x.obs, x.obs, sigma, omega, 2)
  
  L.Syy <- t(chol(Sigma.Y.Y)) # lower cholesky facotr
  
  
  Sigma.Yobs.Y <- cov.x1.x2(x.obs, x.new, sigma, omega, nu)
  
  ## below is the prediction
  
  mu.Y + t(forwardsolve(L.Syy, Sigma.Yobs.Y))%*%forwardsolve(L.Syy, y.obs - mu.Yobs)
  
}

#Sigma.Yobs.Yobs <- exp(maps.test$par[11:40])+ diag(1e-1, 30)

#Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs + cov.x1.x2(x/1.345622, x/1.345622, 1, 1, 2)

#L.Syy <- t(chol(Sigma.Yobs.Yobs)) # lower cholesky facotr


#Sigma.Yobs.Y <- cov.x1.x2(x/1.345622, xx/1.345622, 1, 1, 2)

## below is the prediction

#rep(0, 100) + t(forwardsolve(L.Syy, (Sigma.Yobs.Y)))%*%forwardsolve(L.Syy, rep(1, 30))

## adapt these for SML!


cond.var.hetSML <- function(y.obs.v, x.obs, x.new, sigma, omega, nugget, mu.Y.v, mu.Yobs.v, nu = 2){
  ## all I need to do it to "predict" a new variance ...
  ## just the same as het mean really
  exp(cond.mean(y.obs.v, x.obs, x.new, sigma, omega, mu.Y.v, mu.Yobs.v, nugget^0.5 , nu = 2))
  
}

cond.covmat.hetSML <- function(var.obs, var.new, x.obs, x.new, sigma.e, omega.e, sigma.c, omega.c, lambda.c, rho, N.c, N, nu=2){
  
  ## N = total code runs
  ## N.c = cheap code runs
  ## NB I assume that lambda.c is already squared!
  ## likewise for var.obs/var.new
  
  ## x.obs = (x.cheap, x.expensive)
  
  x.obs <- as.matrix(x.obs)
  x.new <- as.matrix(x.new)
  
  Sigma.Y.Y <- cov.x1.x2(x.new, x.new, sigma.e, omega.e, nu) + diag(as.vector(var.new)) + rho*rho*cov.x1.x2(x.new, x.new, sigma.c, omega.c, 2)
  
  Sigma.Y.Yobs <- cbind(rho*cov.x1.x2(x.new,x.obs[1:N.c,],sigma.c,omega.c,nu), 
                        rho*rho*cov.x1.x2(x.new, x.obs[(1+N.c):N,], sigma.c, omega.c, nu) + cov.x1.x2(x.new,x.obs[(1+N.c):N,],sigma.e,omega.e,nu))
  
  V.EE <- cov.x1.x2(x.obs[(1+N.c):N,], x.obs[(1+N.c):N,], sigma.e, omega.e, nu) + diag(as.vector(var.obs)) + 
    rho*rho*( cov.x1.x2(x.obs[(1+N.c):N,], x.obs[(1+N.c):N,], sigma.c, omega.c, 2) )
  V.CC <- cov.x1.x2(x.obs[1:N.c,], x.obs[1:N.c,], sigma.c, omega.c, nu) + diag(lambda.c, N.c)
  V.CE <- rho*cov.x1.x2(x.obs[1:N.c,], x.obs[(1+N.c):N,], sigma.c, omega.c, nu)
  
  t.x <- cbind( rho*cov.x1.x2(x.new, x.obs[1:N.c,], sigma.c, omega.c, nu), 
           rho*rho*cov.x1.x2(x.new, x.obs[(1+N.c):N,], sigma.c, omega.c, nu) +  cov.x1.x2(x.new, x.obs[(1+N.c):N,], sigma.e, omega.e, nu))
  
  V <- rbind( cbind(V.CC, V.CE)  , cbind(t(V.CE), V.EE))
  
  L <- t(chol(V))
  
  tmp <- forwardsolve(L, t(t.x))       
         
  Sigma.Y.Y - t(tmp)%*%tmp
  
}
cond.covmat.hetSML(rnorm(10)^2, rnorm(5)^2, rnorm(50) ,rnorm(5), 1, 1, 1, 1, 1, 1, 40, 50, 2)


cond.mean.hetSML <- function(y.obs, x.obs, x.new, xi.C, xi.E, var.Y, rho,  mu.Y, mu.Yobs, N.c, nu = 2, cov = F){
  
  ## xi.C = (sigma_c, theta_c, lambda_c)
  ## xi.E = (simga_e, theta_e, lambda_e = fixed const)
  p <- length(xi.C)
  N <- length(y.obs)
  x.obs <- as.matrix(x.obs)
  x.new <- as.matrix(x.new)
  
  
  V11 <- cov.x1.x2(x.obs[1:N.c,], x.obs[1:N.c,], xi.C[1], xi.C[2:(p-1)], nu) + diag(xi.C[p], N.c)
  
  V12 <- rho * cov.x1.x2(x.obs[1:N.c, ], x.obs[(1+N.c):N, ], xi.C[1], xi.C[2:(p-1)], nu)
  
  V22 <- diag(var.Y) + diag(xi.E[p]^2, length(var.Y))
  
  V22 <- V22 + cov.x1.x2(x.obs[(1+N.c):N, ], x.obs[(1+N.c):N, ], xi.E[1], xi.E[2:(p-1)], 2) + rho*rho*cov.x1.x2(x.obs[(1+N.c):N, ], x.obs[(1+N.c):N, ], xi.C[1], xi.C[2:(p-1)], 2)
  
  Sigma.Y.Y <- rbind(cbind(V11, V12), cbind(t(V12), V22))
  
  eigenvalues <- eigen(Sigma.Y.Y)$values
  
  #if(min(eigenvalues)<0){
   # L.Syy <- t( chol( Sigma.Y.Y - diag( min(eigenvalues) - 1e-4, length(eigenvalues) ) ) ) # lower cholesky facotr
  #}else{
    L.Syy <- t(chol(Sigma.Y.Y))
  #}
  t.x1 <- rho * cov.x1.x2(x.new, x.obs[1:N.c,], xi.C[1], xi.C[2:(p-1)], nu )
  t.x2 <- rho*rho*cov.x1.x2(x.new, x.obs[(1+N.c):N,], xi.C[1], xi.C[2:(p-1)], nu ) + cov.x1.x2(x.new, x.obs[(1+N.c):N,], xi.E[1], xi.E[2:(p-1)], nu )
  
  t.x <- cbind(t.x1, t.x2) ## sigma_Y_Yobs^T
  ## below is the prediction
  
  res <- mu.Y + t(forwardsolve(L.Syy, t(t.x)))%*%forwardsolve(L.Syy, y.obs - mu.Yobs)
  if(cov){
    res <- list(pred = res, cov = Sigma.Y.Y)
  }
  res
}
## new fns

dataCovmat <- function(x.c,x.e, gp.hypers, var.beta, H){
  V.c <- cov.x1.x2(x.c, x.c, gp.hypers$c_sigma, gp.hypers$c_theta, 2) + diag(pars$c_nugget, dim(x.c)[1])
  V.e <- cov.x1.x2(x.e, x.e, gp.hypers$m_sigma, gp.hypers$m_theta, 2)
  V.e <- V.e + (gp.hypers$rho^2)*cov.x1.x2(x.e, x.e, gp.hypers$c_sigma, gp.hypers$c_theta, 2) + diag(exp(gp.hypers$logLambda))
  C.c.e <- gp.hypers$rho*cov.x1.x2(x.c, x.e, gp.hypers$c_sigma, gp.hypers$c_theta, 2)
  V <- rbind(cbind(V.c, C.c.e), cbind(t(C.c.e), V.e))
  V + H%*%var.beta%*%t(H)
}
betaMoments <- function(y, H, prior.mean, prior.var, data.precmat, cov.beta.data){
  post.mean <- prior.mean + cov.beta.data%*%data.precmat%*%(y-H%*%prior.mean)
  post.var <- prior.var - cov.beta.data%*%data.precmat%*%t(cov.beta.data)
  res <- list(
    mean=post.mean, var=post.var
  )
  res 
}
gp.crosscov <- function(x, x.obs, sigma, theta, H, H.obs, var.beta){
  H %*% var.beta %*% t(H.obs) + cov.x1.x2(x, x.obs, sigma, theta, 2)
}
ml.crosscov <- function(x, x.c, x.e, pars, H.c, H.e, H, var.bc, var.be){
  ## the variances for the betas are the prior variances
  cov.y.yc <- pars$rho*H %*% var.bc %*% t(H.c) + pars$rho*cov.x1.x2(x, x.c, pars$c_sigma, pars$c_theta, 2)
  cov.y.ye <- H %*% (pars$rho^2*var.bc + var.be) %*% t(H.e) + (pars$rho^2)*cov.x1.x2(x, x.e, pars$c_sigma, pars$c_theta, 2)+cov.x1.x2(x, x.e, pars$m_sigma, pars$m_theta, 2)
  cbind(cov.y.yc, cov.y.ye)
}

post.cov.x1.x2 <- function(x1, x2, x1.scaled, x2.scaled, x.e, x.c, H.e, H.c, pars, precmat){
  ## for multilevel emulator
  ## the scaled variables are on the original scale
  ## unscaled are standardised for the covariance function
  m.x1.x2 <- as.vector(x1 - x2)/pars$m_theta
  m.x1.x2 <- -sum(m.x1.x2^2)
  c.x1.x2 <- as.vector(x1 - x2)/pars$c_theta
  c.x1.x2 <- -sum(c.x1.x2^2)
  h1 <- c(1, log(x1.scaled))
  h2 <- c(1, log(x2.scaled))
  h.both <- rbind(h1, h2)
  
  ## first the GP variance
  prior.cov <- pars$m_sigma^2 * exp(m.x1.x2) + (pars$rho * pars$c_sigma)^2*exp(c.x1.x2)
  
  ## now regression varaince
  prior.cov <- prior.cov + sum(h1*h2*diag(dat$B)[1:10]) + (pars$rho^2)*sum(h1*h2*diag(dat$B)[11:20])
  x.both <- rbind(x1, x2)
  k.x.exp <- cov.x1.x2(x.both, x.e, pars$m_sigma, pars$m_theta, 2) + cov.x1.x2(x.both, x.e, pars$c_sigma * pars$rho, pars$c_theta, 2) + (h.both %*% dat$B[11:20,11:20] %*% t(H.e)) + (pars$rho^2) * h.both %*% dat$B[1:10,1:10] %*% t(H.e)
  k.x.cheap <-  pars$rho*cov.x1.x2(x.both, x.c, pars$c_sigma, pars$c_theta, 2) + pars$rho * h.both %*% dat$B[1:10,1:10] %*% t(H.c)
  k.x <- cbind(k.x.cheap, k.x.exp)
  #print(mean(k.x))
  #post.cov <- prior.cov - sum(t((k.x[1,]) * precmat) * (k.x[2,]))
  post.cov <- prior.cov - (t(k.x[1,])) %*% precmat %*% t(t(k.x[2,]))
  #print(prior.cov - k.x %*% precmat %*% t(k.x))
  post.cov
  
}
