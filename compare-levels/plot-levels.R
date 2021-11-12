## plot levels
library(R.matlab)
cheap.mat <- exp.mat <- matrix(ncol=10, nrow=100)

for(i in 1:100){
  for(j in 1:10){
    filname.cheap <- paste("cheap/out",i,"-",j,".mat",sep="")
    filename.exp <- paste("expensive/out",i,"-",j,".mat",sep="")
    cheap.mat[i,j] <- mean(unlist(R.matlab::readMat(filname.cheap)$y[[j]]))
    exp.mat[i,j] <- mean(unlist(R.matlab::readMat(filename.exp)$y[[j]]))
  }
}

probit <- function(p) qnorm(p)
y.c <- rowMeans(probit(cheap.mat))[1:95]
y.e <- rowMeans(probit(exp.mat))[1:95]
plot(y.c, y.e,xlab = "Probit Availability (Cheap)", ylab = "Probit Availability (Expensive)",pch=16,main = "Cheap and Expensive Runs")
cor(y.e, y.c)^2
####
# extra analysis

abline(lm(y.e~y.c))
linmod <- lm(y.e~y.c + I(y.c^2) + I(y.c^3)+ I(y.c^4))
logmod <- lm(y.e ~ log(y.c))
summary(linmod)
x0 <- seq(from=1.4, to=1.8, length=50)
lines(x0,cbind(1, x0, x0^2, x0^3, x0^4)%*%linmod$coefficients,col=2)
lines(x0, cbind(1, log(x0))%*%logmod$coefficients,col=3)
colVars <- function(X){
  apply(X, 1, var)
}
plot(colVars(probit(cheap.mat)), colVars(probit(exp.mat)))
