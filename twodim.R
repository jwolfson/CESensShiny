S0 <- rnorm(100)
S1 <- rnorm(100,mean=S0)
S0.pts <- seq(min(S0),max(S0),length.out=100)
S1.pts <- seq(min(S1),max(S1),length.out=100)
x <- S0.pts
y <- S1.pts
cor.sens <- 0.3
sigma <- toeplitz(c(var(S0),cor.sens))
sigma[2,2] <- var(S1)
z <- matrix(dmvnorm(expand.grid(S0.pts,S1.pts),mean=c(mean(S0),mean(S1)),sigma=sigma),ncol=length(S0.pts),byrow=TRUE)
filled.contour(x,y,z,color.palette=heat.colors)
