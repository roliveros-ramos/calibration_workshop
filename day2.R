
# reading the data
dat = read.csv("day1_data.csv", row.names = 1)

# simulator of the logistic model
logistic = function(r, K, C, B0=10, T=100) {

  time = seq_len(T)
  B = rep(NA, T)
  B[1] = B0
  for(t in time[-T]) {
    B[t+1] = B[t] + r*B[t]*(1-B[t]/K) - C[t]
    B[t+1] = max(0, B[t+1])
  }
  return(B)
}


C = dat$catch
obs = dat$biomass
# obs = obs*rlnorm(length(obs), mean=0, sd=0.1)
plot(obs, type="b")

# error function
RSS = function(par, obs, ...) {
  # translating the 'par' vector to be our model
  r = par[1]
  K = par[2]
  # run the model
  sim = logistic(r=r, K=K, ...)
  # estimation method
  errors = obs - sim
  rss = sum(errors^2, na.rm = TRUE)
  return(rss)
}

# error function
RSS2 = function(par, obs, ...) {
  # translating the 'par' vector to be our model
  r = par[1]
  K = par[2]
  B0 = par[3]
  # run the model
  sim = logistic(r=r, K=K, B0=B0, ...)
  # estimation method
  errors = obs - sim
  rss = sum(errors^2, na.rm = TRUE)
  return(rss)
}

# error function: maximum likelihood
RSS3 = function(par, obs, ...) {
  # translating the 'par' vector to be our model
  r = par[1]
  K = par[2]
  B0 = par[3]
  sigma = par[4]
  # run the model
  sim = logistic(r=r, K=K, B0=B0, ...)
  # estimation method
  errors = obs - sim
  ll = dnorm(x=errors, mean=0, sd=sqrt(sigma),
             log=TRUE)
  nll = -sum(ll, na.rm=TRUE)
  return(nll)
}

# reduced maximum likelihood
RSS4 = function(par, obs, ...) {
  # translating the 'par' vector to be our model
  r = par[1]
  K = par[2]
  B0 = par[3]
  # run the model
  sim = logistic(r=r, K=K, B0=B0, ...)
  # estimation method
  # errors = obs - sim # normal distributed errors
  # errors = log(obs) - log(sim) # lognormal distributed errors
  # errors = log((obs + 1e-8)/(sim + 1e-8) ) # lognormal distributed errors
  # sigma = var(errors, na.rm=TRUE)
  # ll = dlnorm(x=errors, mean=0, sd=sqrt(sigma),
  #            log=TRUE)
  # nll = -sum(ll, na.rm=TRUE)
  nll = calibrar:::lnorm2(obs, sim)
  return(nll)
}

RSS4(par=c(0.18, 1201, 500), obs=obs, C=C)

# using optim
calibrate(par=c(0.31, 101, 200, 1), fn=RSS3, obs=obs, C=C,
          method="BFGS")
calibrate(par=c(0.31, 800, 200, 1), fn=RSS3, obs=obs, C=C,
          method="BFGS")

calibrate(par=c(0.31, 101, 200), fn=RSS4, obs=obs, C=C,
          method="BFGS")
calibrate(par=c(0.31, 800, 200), fn=RSS4, obs=obs, C=C,
          method="BFGS")
calibrate(par=c(0.20, 800, 200), fn=RSS4, obs=obs, C=C,
          method="BFGS")

calibrate(par=c(0.18, 1000, 500), fn=RSS4, obs=obs, C=C,
          method="BFGS")

calibrate(par=c(0.18, 1000, 500), fn=RSS4, obs=obs, C=C)




