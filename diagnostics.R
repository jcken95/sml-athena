##within sample plot


## diagnostics
Xlabs <- c(
  "gearbox", "generator","freq. conv.", "transformer","MSB","blades","tower","foundations","catch all"
)
par(mfrow=c(3,3))
Smoother=FALSE#set to true to add a loess line to resids
for(i in 1:9){
  plot(x.test[,i],  chol.het, ylim=c(-3,3),xlab=Xlabs[i], pch=19,ylab=""); abline(h = c(-1.96,1.96), col="red",lwd=1.1)
  xx <- sort(x.test[,i],index.return=T)
  if(Smoother){
  mod.tmp <- loess(chol.het[xx$ix]~xx$x, span=2/3)
  lines(x=seq(0.01,5,length=100),y=predict(mod.tmp, newdata = as.data.frame(seq(0.01,5,length=100))), col="blue")
  }
}
#save as 6x9 pdf
par(mfrow=c(3,3))
for(i in 1:9){
  plot(x.test[,i], chol.ml, ylim=c(-3,3),xlab=Xlabs[i], pch=19,ylab=""); abline(h = c(-1.96,1.96), col="red",lwd=1.1)
  xx <- sort(x.test[,i],index.return=T)
  if(Smoother){
  mod.tmp <- loess(chol.ml[xx$ix]~xx$x, span=2/3)
  lines(x=seq(0.01,5,length=100),y=predict(mod.tmp, newdata = as.data.frame(seq(0.01,5,length=100))), col="blue")
  }
}

## qq plots

par(mfrow=c(1,1))
nn <- length(chol.ml)
Z <- matrix(rnorm(nn^2),ncol=nn)
qqnorm(chol.het,main="QQ plot HetGP",ylim=range(Z))
Z.sort <- Z
for(i in 1:nn){
  Z.sort[,i] <- sort(Z[,i])
}
for(i in 1:nn){
  points(qnorm(ppoints(nn)), Z.sort[,i], pch=20,col=rgb(0,0,1,alpha=0.1))
}
points(qnorm(ppoints(nn)), sort(chol.het),pch=16)
lines(qnorm(ppoints(nn)), apply(Z.sort, 1, min))
lines(qnorm(ppoints(nn)), apply(Z.sort, 1, max))
abline(0,1)

par(mfrow=c(1,1))
qqnorm(chol.ml,main="QQ plot SML",ylim=range(Z))
for(i in 1:100){
  points(qnorm(ppoints(nn)), Z.sort[,i], pch=20,col=rgb(0,0,1,alpha=0.1))
}
points(qnorm(ppoints(100)), sort(chol.ml),pch=16)
# lines(qnorm(ppoints(nn)), apply(Z.sort, 1, min))
# lines(qnorm(ppoints(nn)), apply(Z.sort, 1, max))
abline(0,1)

## coverage
par(mfrow=c(1,2))
plot(coverage(chol.het), pch=19,main="HetGP",xlab="Expected Coverage", ylab = "Sample Coverage")
for(i in 1:nn) points(coverage(Z[,i]),col=rgb(1,0,0,alpha=0.1),pch=20)
points(coverage(chol.het), pch=19)
abline(0,1,lwd=1.5)


plot(coverage(chol.ml), pch=19,main="SML",xlab="Expected Coverage",ylab="")

for(i in 1:nn) points(coverage(Z[,i]),col=rgb(1,0,0,alpha=0.1),pch=20)
points(coverage(chol.ml), pch=19)
abline(0,1,lwd=1.5)

## other plots ...
## predict the training data
pars <- ml.fit$par
crosscov.within <- cbind(
  pars$rho*cov.x1.x2(x.e.std, x.c.std, pars$c_sigma, pars$c_theta, 2)+pars$rho*H.e%*%dat$B[1:10,1:10]%*%t(H.c),
  cov.x1.x2(x.e.std, x.e.std, pars$c_sigma*pars$rho, pars$c_theta, 2)+cov.x1.x2(x.e.std, x.e.std, pars$m_sigma, pars$m_theta, 2)+H.e%*%(dat$B[1:10,1:10] + pars$rho^2*dat$B[-(1:10),-(1:10)])%*%t(H.e)
)#cov(eta(x), cheap) and cov(eta(x), expensive)
## used a prior mean of zero
M <- crosscov.within%*%precmat.ml%*%(c(y.c2, y.e2))
## just find pointwise variance
V <- diag(H.e%*%(dat$B[1:10,1:10] + pars$rho^2*dat$B[-(1:10),-(1:10)])%*%t(H.e))
V <- V + (pars$rho*pars$c_sigma)^2 + pars$m_sigma^2 - diag(
  crosscov.within%*%precmat.ml%*%t(crosscov.within)
)
par(mfrow=c(1,1))
plot(y.e2,M, ylim=c(0.5,2.3),xlab="Observed Probit Availability",ylab="Predictive Density",main="Predictive Density for Observed Availability")
abline(0,1)
for(i in 1:length(y.e2)){
  points(rep(y.e2[i],100),rnorm(100,M[i], sqrt(V[i])), col=rgb(1, 0, 0, alpha=0.2), cex=0.7)
}
points(y.e2,M,pch=19)

