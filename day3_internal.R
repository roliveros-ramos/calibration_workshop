library(calibrar)

# simulator of the logistic model
runModel = function(par, T=100) {
  # translating the 'par' vector to be our model
  r  = par[1] # par["r"]
  K  = par[2] # par["K"]
  h  = par[3] # par["h"]
  B0 = par[4] # par["B0"]
  time = seq_len(T)
  B = rep(NA, T)
  C = rep(NA, T)
  B[1] = B0
  for(t in time[-T]) {
    B[t+1] = B[t] + r*B[t]*(1-B[t]/K)
    C[t]   = min(h*B[t], 0.99*B[t+1])
    B[t+1] = B[t+1] - C[t]
  }
  return(list(biomass=B, catch=C))
}

errorFunction = function(par, obs, ...) {
  # run the model
  sim = runModel(par, ...)
  nll_B = calibrar:::lnorm2(obs$biomass, sim$biomass)
  nll_C = calibrar:::lnorm2(obs$catch, sim$catch)
  nll = nll_B + nll_C
  return(nll)
}

real_par = c(r=0.2, K=1024, h=0.1, B0=901)

obs = runModel(real_par, T=110)
obs$biomass = tail(obs$biomass, 100)
obs$catch = tail(obs$catch, 100)
sd = 0.1
obs2 = obs
obs2$biomass = obs2$biomass*rlnorm(100, sd=10*sd)
obs2$catch = obs2$catch*rlnorm(100, sd=0.1*sd)

dat = as.data.frame(obs2)
write.csv(dat, file="day3_data.csv")

par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(obs2$biomass)
lines(obs$biomass, col="blue")

plot(obs2$catch)
lines(obs$catch, type="h", col="blue")

guess = c(r=0.1, K=800, h=0.0, B0=500)
guess = c(r=0.1, K=1000, h=0.0, B0=800)

errorFunction(guess, obs)

calibrate(par=guess, fn=errorFunction, obs=obs2,
          lower=c(0, 0, 0, 0), upper=c(4, 1e4, 20, 1000),
          replicates=3)

calibrate(par=guess, fn=errorFunction, obs=obs,
          lower=c(0, 0, 0, 0), upper=c(4, 1e4, 20, 1000),
          phases=c(1,2,1,2))


















