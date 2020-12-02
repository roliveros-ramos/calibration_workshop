
# simulator of the logistic model
logisticH = function(par, T=100) {
  # translating the 'par' vector to be our model
  r  = par[1]
  K  = par[2]
  h  = par[3]
  B0 = par[4]
  time = seq_len(T)
  B = rep(NA, T)
  C = rep(NA, T)
  B[1] = B0
  for(t in time[-T]) {
    B[t+1] = B[t] + r*B[t]*(1-B[t]/K)
    C[t]   = min(exp(h)*B[t], 0.99*B[t+1])
    B[t+1] = B[t+1] - C[t]
  }
  return(list(biomass=B, catch=C))
}

errorFunction = function(par, obs, ...) {
  # run the model
  sim = logisticH(par, ...)
  nll_B = calibrar:::lnorm2(obs$biomass, sim$biomass)
  nll_C = calibrar:::lnorm2(obs$catch, sim$catch)
  nll = nll_B + nll_C
  return(nll)
}













