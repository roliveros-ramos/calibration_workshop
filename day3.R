
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

real_par = c(r=0.2, K=1024, h=0.1, B0=301)

obs = runModel(real_par)

plot(obs$biomass)
plot(obs$catch)

guess = c(r=0.1, K=800, h=0.2, B0=500)

calibrate(par=guess, fn=errorFunction, obs=obs)











