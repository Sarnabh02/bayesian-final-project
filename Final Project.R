library(readxl)
hdata <- read_xlsx("data/healthdata.xlsx", sheet =2)
hdefs <- read_xlsx("data/healthdata.xlsx", sheet =1)



library(tidyverse)
hdata1 <- hdata %>% select(SEF28_medHHInc, HB10_RxAbuse, HO18_LifeExpect)
cor(hdata1)

fit <- lm(HO18_LifeExpect ~ SEF28_medHHInc + HB10_RxAbuse, data = hdata1)
summary(fit)
Xmat<-model.matrix(fit)
XtX<-t(Xmat)%*%Xmat
sig2hat<-summary(fit)$sigma^2
p<- ncol(Xmat)-1

library(rjags)


# Model 1 : Mixture of g Priors -----------------------------------------------


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

m <- jags.model(final.proj.model1, d, inits, n.chains=3, n.adapt =0)

### Make a preliminary run of 1000 iterations, with monitoring

x <- coda.samples(m, c("alpha","beta","sigma2",'ginv'), n.iter=1000)


### Assess convergence

gelman.diag(x, autoburnin=FALSE, multivariate=FALSE)
gelman.plot(x, autoburnin=FALSE)
plot(x)

### Run 10000 more iterations

x <- coda.samples(m, c("alpha","beta","sigma2",'ginv'), n.iter=10000)

### Assess convergence

gelman.diag(x, autoburnin=FALSE, multivariate=FALSE)
gelman.plot(x, autoburnin=FALSE)
plot(x)

### Check stats after burn-in

summary(x)

### Replication and DIC check

update(m, 1000)
dic.samples(m, 10000, type="pD")

# Model 2: Unit Information Prior (No g)  ----------------------------------------------------------

XtXdn<-t(Xmat)%*%Xmat/nrow(Xmat)
bhat <- fit$coefficients

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
  beta ~ dmnorm(B0,tausq*Sig0)
}
',file = {final.proj.model1 = tempfile()})

d <- list(y = hdata1$HO18_LifeExpect,
          X = Xmat[,-1],
          B0 = bhat[-1],
          Sig0 = XtXdn[-1,-1],
          nu0 =1,
          sig20=sig2hat
)


inits <- list(list(tausq=1, alpha = 1, beta=c(0,0)),
              list(tausq=1, alpha=-1, beta=c(1,1)),
              list(tausq=1, alpha=0, beta=c(-1,-1))
)

m <- jags.model(final.proj.model1, d, inits, n.chains=3, n.adapt=0)

### Make a preliminary run of 1000 iterations, with monitoring

x <- coda.samples(m, c("alpha","beta","sigma2"), n.iter=1000)

### Assess convergence

gelman.diag(x, autoburnin=FALSE, multivariate=FALSE)
gelman.plot(x, autoburnin=FALSE)
plot(x)

### Check stats after burn-in

summary(x)

### Replication and DIC check

update(m, 1000)
dic.samples(m, 10000, type="pD")








