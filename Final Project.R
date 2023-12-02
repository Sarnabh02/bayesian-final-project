library(readxl)
hdata <- read_xlsx("healthdata.xlsx", sheet =2)
hdefs <- read_xlsx("healthdata.xlsx", sheet =1)



library(tidyverse)
hdata1 <- hdata %>% select(SEF28_medHHInc, HB10_RxAbuse, HO18_LifeExpect)
cor(hdata1)

fit <- lm(HO18_LifeExpect ~ SEF28_medHHInc + HB10_RxAbuse, data = hdata1)
Xmat<-model.matrix(fit)
XtX<-t(Xmat)%*%Xmat
sig2hat<-summary(fit)$sigma^2
p<- ncol(Xmat)-1

library(rjags)

cat('
data {
  dim.X <- dim(X)
  n <- dim.X[1]
}

model {
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], tausq)
    mu[i] <- alpha + X[i,]%*%beta[]
  }
  tausq ~ dgamma(nu0/2,nu0*sig20/2)
  sigma2 <- 1/tausq

  alpha ~ dnorm(0,.0001)
  beta ~ dmnorm(B0,ginv*tausq*Sig0)
  ginv ~ dgamma(.5,n/2)
}
',file = {final.proj.model1 = tempfile()})

d <- list(y = hdata1$HO18_LifeExpect,
          X = Xmat[,-1],
          B0 = rep(0,p),
          Sig0 = XtX[-1,-1],
          nu0 =1,
          sig20=sig2hat
)


inits <- list(list(tausq=1, alpha = 1, beta=rnorm(p) ),
              list(tausq=1, alpha=-1, beta=rnorm(p,sd=3) ),
              list(tausq=1, alpha=0, beta=rnorm(p,sd=.5) )
)

m <- jags.model(final.proj.model1, d, inits, n.chains=3)

### Make a preliminary run of 1000 iterations, with monitoring

x <- coda.samples(m, c("alpha","beta","sigma2",'ginv'), n.iter=1000)


### Assess convergence

gelman.diag(x, autoburnin=FALSE, multivariate=FALSE)


### Run 10000 more iterations

x <- coda.samples(m, c("alpha","beta","sigma2",'ginv'), n.iter=10000)

### Assess convergence

gelman.diag(x, autoburnin=FALSE, multivariate=FALSE)

### Check stats after burn-in

summary(x)

