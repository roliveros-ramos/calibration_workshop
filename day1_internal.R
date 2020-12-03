logistic = function(r, K, C, B0=10, T=100) {

  time = seq_len(T)
  B = rep(NA, T)
  B[1] = B0
  for(t in time[-T]) {
    B[t+1] = max(0, B[t] + r*B[t]*(1-B[t]/K) - C[t])
  }
  return(B)
}


# simulator of the logistic model
logistic2 = function(r, K, B0, h, T=100) {
  time = seq_len(T)
  B = rep(NA, T)
  C = rep(NA, T)
  B[1] = B0
  h = runif(T, min=h[1], max=h[2])
  for(t in time[-T]) {
    B[t+1] = B[t] + r*B[t]*(1-B[t]/K)
    C[t] = min(h[t]*B[t], 0.99*B[t+1])
    B[t+1] = B[t+1] - C[t]
  }
  return(list(B=B, C=C, h=h))
}

C = rep(0, 100)
C = sort(runif(n=100, min=5, max=20))
obs = logistic2(r=0.18, K=1200, B0=500, h=c(0.0, 0.2))
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(1,1,1,1))
plot(obs$B, type="b")
plot(obs$C, type="h")

dat = data.frame(obs)
dat$h = NULL
dat = round(dat, 3)
names(dat) = c("biomass", "catch")

ind = sample(1:100, size=70)

dat$biomass[ind] = NA
write.csv(dat, file="day1_data.csv")

x = logistic(r=0.18, K=1200, B0=500, C=dat$catch)

# error function
RSS = function(par, obs, ...) {
  # translating the 'par' vector to be our model
  r = par[1]
  K = par[2]
  # run the model
  sim = logistic(r=r, K=K, ...)
  # estimation method
  errors = obs - sim
  rss = sum(errors^2, na.rm=TRUE)
  return(rss)
}

RSS(par=c(0.28, 1200), B0=500, obs=dat$biomass, C=dat$catch)

optim(par=c(0.31, 1000), fn=RSS, B0=500, obs=dat$biomass, C=dat$catch)


optx = calibrate(par=c(0.1, 500), fn=RSS, B0=500, obs=dat$biomass, C=dat$catch,
                 method="BFGS")

opt0 = optim(par=c(0.1, 500), fn=RSS, B0=500, obs=dat$biomass, C=dat$catch)
opt1 = optim(par=opt0$par, fn=RSS, B0=500, obs=dat$biomass, C=dat$catch)
opt2 = optim(par=opt1$par, fn=RSS, B0=500, obs=dat$biomass, C=dat$catch)

opt0$par
opt1$par
opt2$par

optb0 = optim(par=c(0.8, 500), fn=RSS, obs=obs, C=C,
             method="BFGS")

optc0 = calibrate(par=c(0.8, 500), fn=RSS, obs=obs, C=C,
              method="BFGS")

rs = seq(from=0.01, to=0.5, by=0.01)
Ks = seq(from=50, to=500, by=10)

pars = expand.grid(r=rs, K=Ks)

# looping over all the parameters combinations ?apply
out = apply(pars, 1, FUN=RSS, obs=obs, C=C)
dim(out) = c(length(rs), length(Ks))

library(fields)
image.plot(rs, Ks, log(out+1e-9))
points(0.3, 100, col="red",pch=1)



