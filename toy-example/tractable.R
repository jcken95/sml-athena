library("rstan")
source("../GPfunctionsOptim.R")
source("../hetGPfunctions.R")
het.uni <- stan_model("hetGP-univariate.stan")
SML <- stan_model("SML-hetGP-univariate.stan")
## lets define some functions ...


f1 <- function(x){
  4*sin(7*pi*x)
}
f2 <- function(x){
  f1(x) + 5*(2*x+1) + 3*log(x+0.01)
}
plot(f1(x.test), f2(x.test))
f1.rand <- function(x){
  
  f1(x) + 4*rnorm(length(x))
  
}

f2.rand <- function(x){
  
  f2(x) + (5*x+2)*rnorm(length(x))
  
}
NN <- 100 ## number reps for sim study

x <- seq(0 , 1, length = NN)
plot(x, f1(x), type = "l")
plot(x, f2.rand(x), type = "l")
plot(f1(x), f2(x)); cor(f1(x), f2(x))



############################
############################
############################
######## start here ########
############################
############################
############################
NN <- 100
score.emp1 <- mse.emp1 <- matrix(0, ncol = 2, nrow = NN)
set.seed(1234)
st1 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  
  ## prior on params
  
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  #het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  ## construct prior ...
  # log variance
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat
  dmat.logv <- cbind(1, x.std)
  dmat <- cbind(1, x.std)
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x.std, x.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.std)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x.std, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # mean
  dmat.new <- cbind(1, x.new)
  prior.mean.new <- dmat.new%*%m_beta_m
  prior.var.new <-  dmat.new%*%(m_beta_s^2)%*%t(dmat.new) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2) + diag(exp(log.lambdasq.hat))
  # data
  prior.mean <- dmat%*%m_beta_m
  var.dat <-  dmat%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(exp(pars$logLambda))
  cov.mean <- dmat.new%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.new, x.std, pars$m_sigma, pars$m_theta, 2)
  dat.precmat <- chol2inv(chol(var.dat))
  # posterior on observable quantity
  post.var.hgp <- diag(prior.var.new - cov.mean %*% dat.precmat %*% t(cov.mean))
  post.mean.hgp <- prior.mean.new + cov.mean %*% (dat.precmat %*% (y1 - prior.mean))
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  n.c <- 100
  n.e <- length(x1) - 1
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  c_beta_m <- rep(0,2); c_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    v_a_theta = 2, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = c_beta_m, c_beta_s = diag(c_beta_s),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 0, s_rho = 1/3
    
  )

  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  #sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  # log variance --- this is essentially the same as HetGP
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat (of log lambdasq)
  dmat.logv <- cbind(1, x2)

  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x2, x2, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x2)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x2, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # construct data covmat
  dmat <- rbind(cbind(cbind(1, x.c), matrix(0, ncol=2,nrow=length(x.c))), cbind( pars$rho*cbind(1, x2), cbind(1,x2)  ))
  dmat.new <- cbind(1, x.new) ## only predicting on the "expensive scale"
  ### now construct prior (and then post.) for observables
  #means
  prior.mean <- dmat%*%c(c_beta_m, m_beta_m)
  prior.mean.new <- dmat.new %*% c(pars$rho*c_beta_m + m_beta_m)
  #vars
  dat.covmat <- dmat %*% (diag(c(diag(c_beta_s),diag(m_beta_s)))^2) %*% t(dmat)
  cov.exp.cheap <- pars$rho*cov.x1.x2(x.c, x2, pars$c_sigma, pars$c_theta,2)
  gp.cov.bit <- rbind(
    cbind( cov.x1.x2(x.c, x.c, pars$c_sigma, pars$c_theta,2)+diag(pars$c_nugget,length(x.c)), cov.exp.cheap ),
    cbind(t(cov.exp.cheap), cov.x1.x2(x2, x2, pars$rho*pars$c_sigma, pars$c_theta,2)+cov.x1.x2(x2, x2, pars$m_sigma, pars$m_theta,2)+diag(exp(pars$logLambda)))
    )
  dat.covmat <- dat.covmat + gp.cov.bit ## dat cov mat
  prior.var.new <- dmat.new %*% ((pars$rho*pars$rho)*c_beta_s^2 + m_beta_s^2) %*% t(dmat.new)
  prior.var.new <- prior.var.new + cov.x1.x2(x.new, x.new, pars$c_sigma*pars$rho, pars$c_theta, 2) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2)
  prior.var.new <- prior.var.new + diag(exp(as.vector(log.lambdasq.hat)))
  cov.mean <- rbind(pars$rho*cov.x1.x2(x.c, x.new, pars$c_sigma, pars$c_theta, 2),
                    cov.x1.x2(x2, x.new, pars$rho*pars$c_sigma, pars$c_theta, 2) + cov.x1.x2(x2, x.new, pars$m_sigma, pars$m_theta, 2))
  cov.mean <- t(cov.mean + dmat %*% diag( c(diag(c_beta_s), diag(m_beta_s))^2 ) %*% t(cbind(pars$rho*dmat.new, dmat.new)))
  chol.var <- t(chol(dat.covmat))
  L.c <- forwardsolve(chol.var, t(cov.mean))
  # posterior on observable quantity
  post.var.sml <- diag(prior.var.new) - diag(t(L.c) %*% L.c)
  post.mean.sml <- prior.mean.new + (t(L.c) %*%  (forwardsolve(chol.var,c(y.c, y2) - dmat%*%c(c_beta_m, m_beta_m)  )))

  test.data <- f2.rand(x.test)
  
  mse.emp1[K, ] <-   c(  mean((post.mean.hgp - test.data)^2), mean((post.mean.sml - test.data)^2) )
  score.emp1[K, ] <- c(sum( (-(post.mean.hgp - test.data)^2)/post.var.hgp - log(post.var.hgp) ), sum( (-(post.mean.sml - test.data)^2)/post.var.sml - log(post.var.sml)))
  
  
  print(K)
  
}
en1 <- Sys.time()
tot1 <- en1 - st1

tot1
boxplot(mse.emp1)
boxplot(score.emp1)
#### simulation 2 ####
par(mfrow=c(1,1))
plot(x.new, post.mean.hgp)
points(x.new, post.mean.sml)
score.emp2 <- mse.emp2 <- matrix(0, ncol = 2, nrow = NN)

st2 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  
  ## prior on params
  
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  #het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  ## construct prior ...
  # log variance
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat
  dmat.logv <- cbind(1, x.std)
  dmat <- cbind(1, x.std)
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x.std, x.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.std)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x.std, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # mean
  dmat.new <- cbind(1, x.new)
  prior.mean.new <- dmat.new%*%m_beta_m
  prior.var.new <-  dmat.new%*%(m_beta_s^2)%*%t(dmat.new) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2) + diag(exp(log.lambdasq.hat))
  # data
  prior.mean <- dmat%*%m_beta_m
  var.dat <-  dmat%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(exp(pars$logLambda))
  cov.mean <- dmat.new%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.new, x.std, pars$m_sigma, pars$m_theta, 2)
  dat.precmat <- chol2inv(chol(var.dat))
  # posterior on observable quantity
  post.var.hgp <- diag(prior.var.new - cov.mean %*% dat.precmat %*% t(cov.mean))
  post.mean.hgp <- prior.mean.new + cov.mean %*% (dat.precmat %*% (y1 - prior.mean))
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  b <- 2
  n.c <- 100*b
  n.e <- length(x1) - b
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  c_beta_m <- rep(0,2); c_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    v_a_theta = 2, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = c_beta_m, c_beta_s = diag(c_beta_s),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 0, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  
  #sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  # log variance --- this is essentially the same as HetGP
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat (of log lambdasq)
  dmat.logv <- cbind(1, x2)
  
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x2, x2, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x2)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x2, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # construct data covmat
  dmat <- rbind(cbind(cbind(1, x.c), matrix(0, ncol=2,nrow=length(x.c))), cbind( pars$rho*cbind(1, x2), cbind(1,x2)  ))
  dmat.new <- cbind(1, x.new) ## only predicting on the "expensive scale"
  ### now construct prior (and then post.) for observables
  #means
  prior.mean <- dmat%*%c(c_beta_m, m_beta_m)
  prior.mean.new <- dmat.new %*% c(pars$rho*c_beta_m + m_beta_m)
  #vars
  dat.covmat <- dmat %*% (diag(c(diag(c_beta_s),diag(m_beta_s)))^2) %*% t(dmat)
  cov.exp.cheap <- pars$rho*cov.x1.x2(x.c, x2, pars$c_sigma, pars$c_theta,2)
  gp.cov.bit <- rbind(
    cbind( cov.x1.x2(x.c, x.c, pars$c_sigma, pars$c_theta,2)+diag(pars$c_nugget,length(x.c)), cov.exp.cheap ),
    cbind(t(cov.exp.cheap), cov.x1.x2(x2, x2, pars$rho*pars$c_sigma, pars$c_theta,2)+cov.x1.x2(x2, x2, pars$m_sigma, pars$m_theta,2)+diag(exp(pars$logLambda)))
  )
  dat.covmat <- dat.covmat + gp.cov.bit ## dat cov mat
  prior.var.new <- dmat.new %*% ((pars$rho*pars$rho)*c_beta_s^2 + m_beta_s^2) %*% t(dmat.new)
  prior.var.new <- prior.var.new + cov.x1.x2(x.new, x.new, pars$c_sigma*pars$rho, pars$c_theta, 2) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2)
  prior.var.new <- prior.var.new + diag(exp(as.vector(log.lambdasq.hat)))
  cov.mean <- rbind(pars$rho*cov.x1.x2(x.c, x.new, pars$c_sigma, pars$c_theta, 2),
                    cov.x1.x2(x2, x.new, pars$rho*pars$c_sigma, pars$c_theta, 2) + cov.x1.x2(x2, x.new, pars$m_sigma, pars$m_theta, 2))
  cov.mean <- t(cov.mean + dmat %*% diag( c(diag(c_beta_s), diag(m_beta_s))^2 ) %*% t(cbind(pars$rho*dmat.new, dmat.new)))
  chol.var <- t(chol(dat.covmat))
  L.c <- forwardsolve(chol.var, t(cov.mean))
  # posterior on observable quantity
  post.var.sml <- diag(prior.var.new) - diag(t(L.c) %*% L.c)
  post.mean.sml <- prior.mean.new + (t(L.c) %*%  (forwardsolve(chol.var,c(y.c, y2) - dmat%*%c(c_beta_m, m_beta_m)  )))

  test.data <- f2.rand(x.test)
  
  mse.emp2[K, ] <-   c(  mean((post.mean.hgp - test.data)^2), mean((post.mean.sml - test.data)^2) )
  score.emp2[K, ] <- c(sum( (-(post.mean.hgp - test.data)^2)/post.var.hgp - log(post.var.hgp) ), sum( (-(post.mean.sml - test.data)^2)/post.var.sml - log(post.var.sml)))
  
  
  print(K)
  
  
}
en2 <- Sys.time()
tot2 <- en2 - st2

tot2
boxplot(mse.emp2, ylim=c(9,15))
boxplot(score.emp2, ylim= c(-4000,-3000))

boxplot(mse.emp1, ylim=c(9,15))
boxplot(score.emp1, ylim= c(-4000,-3000))

score.emp3 <- mse.emp3 <- matrix(0, ncol = 2, nrow = NN)

st3 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  
  
  ## first generate "purely exp" emulator
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  
  ## prior on params
  
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  #het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  ## construct prior ...
  # log variance
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat
  dmat.logv <- cbind(1, x.std)
  dmat <- cbind(1, x.std)
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x.std, x.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.std)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x.std, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # mean
  dmat.new <- cbind(1, x.new)
  prior.mean.new <- dmat.new%*%m_beta_m
  prior.var.new <-  dmat.new%*%(m_beta_s^2)%*%t(dmat.new) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2) + diag(exp(log.lambdasq.hat))
  # data
  prior.mean <- dmat%*%m_beta_m
  var.dat <-  dmat%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(exp(pars$logLambda))
  cov.mean <- dmat.new%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.new, x.std, pars$m_sigma, pars$m_theta, 2)
  dat.precmat <- chol2inv(chol(var.dat))
  # posterior on observable quantity
  post.var.hgp <- diag(prior.var.new - cov.mean %*% dat.precmat %*% t(cov.mean))
  post.mean.hgp <- prior.mean.new + cov.mean %*% (dat.precmat %*% (y1 - prior.mean))
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  b <- 3
  n.c <- 100*b
  n.e <- length(x1) - b
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  c_beta_m <- rep(0,2); c_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    v_a_theta = 2, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = c_beta_m, c_beta_s = diag(c_beta_s),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 0, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  
  #sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  # log variance --- this is essentially the same as HetGP
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat (of log lambdasq)
  dmat.logv <- cbind(1, x2)
  
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x2, x2, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x2)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x2, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # construct data covmat
  dmat <- rbind(cbind(cbind(1, x.c), matrix(0, ncol=2,nrow=length(x.c))), cbind( pars$rho*cbind(1, x2), cbind(1,x2)  ))
  dmat.new <- cbind(1, x.new) ## only predicting on the "expensive scale"
  ### now construct prior (and then post.) for observables
  #means
  prior.mean <- dmat%*%c(c_beta_m, m_beta_m)
  prior.mean.new <- dmat.new %*% c(pars$rho*c_beta_m + m_beta_m)
  #vars
  dat.covmat <- dmat %*% (diag(c(diag(c_beta_s),diag(m_beta_s)))^2) %*% t(dmat)
  cov.exp.cheap <- pars$rho*cov.x1.x2(x.c, x2, pars$c_sigma, pars$c_theta,2)
  gp.cov.bit <- rbind(
    cbind( cov.x1.x2(x.c, x.c, pars$c_sigma, pars$c_theta,2)+diag(pars$c_nugget,length(x.c)), cov.exp.cheap ),
    cbind(t(cov.exp.cheap), cov.x1.x2(x2, x2, pars$rho*pars$c_sigma, pars$c_theta,2)+cov.x1.x2(x2, x2, pars$m_sigma, pars$m_theta,2)+diag(exp(pars$logLambda)))
  )
  dat.covmat <- dat.covmat + gp.cov.bit ## dat cov mat
  prior.var.new <- dmat.new %*% ((pars$rho*pars$rho)*c_beta_s^2 + m_beta_s^2) %*% t(dmat.new)
  prior.var.new <- prior.var.new + cov.x1.x2(x.new, x.new, pars$c_sigma*pars$rho, pars$c_theta, 2) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2)
  prior.var.new <- prior.var.new + diag(exp(as.vector(log.lambdasq.hat)))
  cov.mean <- rbind(pars$rho*cov.x1.x2(x.c, x.new, pars$c_sigma, pars$c_theta, 2),
                    cov.x1.x2(x2, x.new, pars$rho*pars$c_sigma, pars$c_theta, 2) + cov.x1.x2(x2, x.new, pars$m_sigma, pars$m_theta, 2))
  cov.mean <- t(cov.mean + dmat %*% diag( c(diag(c_beta_s), diag(m_beta_s))^2 ) %*% t(cbind(pars$rho*dmat.new, dmat.new)))
  chol.var <- t(chol(dat.covmat))
  L.c <- forwardsolve(chol.var, t(cov.mean))
  # posterior on observable quantity
  post.var.sml <- diag(prior.var.new) - diag(t(L.c) %*% L.c)
  post.mean.sml <- prior.mean.new + (t(L.c) %*%  (forwardsolve(chol.var,c(y.c, y2) - dmat%*%c(c_beta_m, m_beta_m)  )))
  
  test.data <- f2.rand(x.test)
  
  mse.emp3[K, ] <-   c(  mean((post.mean.hgp - test.data)^2), mean((post.mean.sml - test.data)^2) )
  score.emp3[K, ] <- c(sum( (-(post.mean.hgp - test.data)^2)/post.var.hgp - log(post.var.hgp) ), sum( (-(post.mean.sml - test.data)^2)/post.var.sml - log(post.var.sml)))
  
  
  print(K)
  
  
}
en3 <- Sys.time()
tot3 <- en3 - st3



score.emp4 <- mse.emp4 <- matrix(0, ncol = 2, nrow = NN)

st4 <- Sys.time()
for(K in 1:NN){
  
  
  ## first generate "purely exp" emulator
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(40, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  
  ## prior on params
  
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  #het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  ## construct prior ...
  # log variance
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat
  dmat.logv <- cbind(1, x.std)
  dmat <- cbind(1, x.std)
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x.std, x.std, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.std)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x.std, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # mean
  dmat.new <- cbind(1, x.new)
  prior.mean.new <- dmat.new%*%m_beta_m
  prior.var.new <-  dmat.new%*%(m_beta_s^2)%*%t(dmat.new) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2) + diag(exp(log.lambdasq.hat))
  # data
  prior.mean <- dmat%*%m_beta_m
  var.dat <-  dmat%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(exp(pars$logLambda))
  cov.mean <- dmat.new%*%(m_beta_s^2)%*%t(dmat) + cov.x1.x2(x.new, x.std, pars$m_sigma, pars$m_theta, 2)
  dat.precmat <- chol2inv(chol(var.dat))
  # posterior on observable quantity
  post.var.hgp <- diag(prior.var.new - cov.mean %*% dat.precmat %*% t(cov.mean))
  post.mean.hgp <- prior.mean.new + cov.mean %*% (dat.precmat %*% (y1 - prior.mean))
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  b <- 4
  n.c <- 100*b
  n.e <- length(x1) - b
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  m_beta_m <- rep(0,2); m_beta_s <- diag(1, 2)
  c_beta_m <- rep(0,2); c_beta_s <- diag(1, 2)
  v_beta_m <- rep(0,2); v_beta_s <- diag(1, 2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = m_beta_m, m_beta_s = diag(m_beta_s),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = v_beta_m, v_beta_s = diag(v_beta_s),
    v_a_theta = 2, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = c_beta_m, c_beta_s = diag(c_beta_s),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 4,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 1, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  
  #sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  # log variance --- this is essentially the same as HetGP
  dmat.logv.new <- cbind(1, x.new)
  prior.mean.logv.new <- dmat.logv.new%*%v_beta_m
  prior.var.logv.new <-  dmat.logv.new%*%(v_beta_s^2)%*%t(dmat.logv.new) + cov.x1.x2(x.new, x.new, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.new)[1])
  # data variance mat (of log lambdasq)
  dmat.logv <- cbind(1, x2)
  
  prior.var.logv <-  dmat.logv%*%(v_beta_s^2)%*%t(dmat.logv) + cov.x1.x2(x2, x2, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x2)[1])
  prior.mean.logv <- dmat.logv%*%v_beta_m
  precmat.logv <- chol2inv(chol(prior.var.logv))
  # covar btwn data and observed
  cov.logv <- dmat.logv.new  %*% v_beta_s^2 %*% t(dmat.logv) + cov.x1.x2(x.new, x2, pars$v_sigma, pars$v_theta, 2)
  # update variance for x.new
  ## nb. this is redundant for a median prediction of lambda^2(x) --- comment included for completeness
  #post.var.logv.new <- prior.var.logv.new - cov.logv %*% precmat.logv %*% t(cov.logv)
  # predict log lambda2(x.new)
  log.lambdasq.hat <- prior.mean.logv.new + cov.logv %*% (precmat.logv%*%(pars$logLambda - as.vector(prior.mean.logv)))
  # construct data covmat
  dmat <- rbind(cbind(cbind(1, x.c), matrix(0, ncol=2,nrow=length(x.c))), cbind( pars$rho*cbind(1, x2), cbind(1,x2)  ))
  dmat.new <- cbind(1, x.new) ## only predicting on the "expensive scale"
  ### now construct prior (and then post.) for observables
  #means
  prior.mean <- dmat%*%c(c_beta_m, m_beta_m)
  prior.mean.new <- dmat.new %*% c(pars$rho*c_beta_m + m_beta_m)
  #vars
  dat.covmat <- dmat %*% (diag(c(diag(c_beta_s),diag(m_beta_s)))^2) %*% t(dmat)
  cov.exp.cheap <- pars$rho*cov.x1.x2(x.c, x2, pars$c_sigma, pars$c_theta,2)
  gp.cov.bit <- rbind(
    cbind( cov.x1.x2(x.c, x.c, pars$c_sigma, pars$c_theta,2)+diag(pars$c_nugget,length(x.c)), cov.exp.cheap ),
    cbind(t(cov.exp.cheap), cov.x1.x2(x2, x2, pars$rho*pars$c_sigma, pars$c_theta,2)+cov.x1.x2(x2, x2, pars$m_sigma, pars$m_theta,2)+diag(exp(pars$logLambda)))
  )
  dat.covmat <- dat.covmat + gp.cov.bit ## dat cov mat
  prior.var.new <- dmat.new %*% ((pars$rho*pars$rho)*c_beta_s^2 + m_beta_s^2) %*% t(dmat.new)
  prior.var.new <- prior.var.new + cov.x1.x2(x.new, x.new, pars$c_sigma*pars$rho, pars$c_theta, 2) + cov.x1.x2(x.new, x.new, pars$m_sigma, pars$m_theta, 2)
  prior.var.new <- prior.var.new + diag(exp(as.vector(log.lambdasq.hat)))
  cov.mean <- rbind(pars$rho*cov.x1.x2(x.c, x.new, pars$c_sigma, pars$c_theta, 2),
                    cov.x1.x2(x2, x.new, pars$rho*pars$c_sigma, pars$c_theta, 2) + cov.x1.x2(x2, x.new, pars$m_sigma, pars$m_theta, 2))
  cov.mean <- t(cov.mean + dmat %*% diag( c(diag(c_beta_s), diag(m_beta_s))^2 ) %*% t(cbind(pars$rho*dmat.new, dmat.new)))
  chol.var <- t(chol(dat.covmat))
  L.c <- forwardsolve(chol.var, t(cov.mean))
  # posterior on observable quantity
  post.var.sml <- diag(prior.var.new) - diag(t(L.c) %*% L.c)
  post.mean.sml <- prior.mean.new + (t(L.c) %*%  (forwardsolve(chol.var,c(y.c, y2) - dmat%*%c(c_beta_m, m_beta_m)  )))
  
  test.data <- f2.rand(x.test)
  
  mse.emp4[K, ] <-   c(  mean((post.mean.hgp - test.data)^2), mean((post.mean.sml - test.data)^2) )
  score.emp4[K, ] <- c(sum( (-(post.mean.hgp - test.data)^2)/post.var.hgp - log(post.var.hgp) ), sum( (-(post.mean.sml - test.data)^2)/post.var.sml - log(post.var.sml)))
  
  
  print(K)
  
  

## plot HGP emulator ...
x0 <- sort(x.new,index.return=T)
xx <- x0$x; ind <- x0$ix
X <- x.test[ind]
  plot(X, f2(X), type="l", ylim = c(-10,30))
## varaince bands
lwr <- f2(X) - 1.96*abs(5*X+2)
upr <- f2(X) + 1.96*abs(5*X+2)
polygon(c(rev(X), (X)), c(rev(lwr) , upr), col=rgb(0,0.9,0.9,alpha=0.5),border=NA)
#x.new
lines(X, f2(X), lwd=3)

lines(X, post.mean.hgp[ind])
lines(X, (post.mean.hgp - 1.96*sqrt(post.var.hgp))[ind], lty=2)
lines(X, (post.mean.hgp + 1.96*sqrt(post.var.hgp))[ind], lty=2)
points(x1, y1)
points(X, post.mean.sml[ind])
en4 <- Sys.time()
tot4 <- en4 - st4
boxplot(mse.emp4)
boxplot(cbind(score.emp1, score.emp2, score.emp3, score.emp4), ylim = c(-5000, -3800))
boxplot(cbind(mse.emp1, mse.emp2, mse.emp3, mse.emp4))
Scores <- vector("list", 4)
Scores[[1]] <- score.emp1; Scores[[2]] <- score.emp2
Scores[[3]] <- score.emp3; Scores[[4]] <- score.emp4

Mses <- vector("list", 4)
Mses[[1]] <- mse.emp1; Mses[[2]] <- mse.emp2
Mses[[3]] <- mse.emp3; Mses[[4]] <- mse.emp4
saveRDS(Scores, "res/score.RDS")
saveRDS(Mses, "res/mse.RDS")
