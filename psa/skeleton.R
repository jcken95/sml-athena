
first.order <- function(){
  M <- lhs::randomLHS(N, p)%*%diag(dat$xScale) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  M.prime <- lhs::randomLHS(N, p) + matrix(rep(dat$xCenter,N), ncol=p, byrow=T)
  
  y.hat <- as.vector( ... )## predict from GP using design M
  
  E.y <- mean(y.hat)
  Var.y <- mean(y.hat^2)-E.y^2##epistemic uncertainty in mean response
  
  ## compute stochastic variation due to stochastic nature of simulations
  ## aleatory uncertainty
  log.var <- as.vector( ... )## predict (log) GP using design M
  E.lambdasq <- mean(exp(log.var))
  ## now evaluate the sensitivy indices
  W <- rep(0,p)
  for(j in 1:p){
    M2 <- M.prime
    M2[,j] <- M[,j]
    y2 <- as.vector( ... )## predict from GP using design M2
    W[j] <- sum(y.hat * y2)/(N-1)
  }
}