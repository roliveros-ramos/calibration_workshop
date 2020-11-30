
# simulator of the logistic model
logistic = function(r, K, C, B0=10, T=100) {

  time = seq_len(T)
  B = rep(NA, T)
  B[1] = B0
  for(t in time[-T]) {
    B[t+1] = B[t] + r*B[t]*(1-B[t]/K) - C[t]
  }
  return(B)
}

C = rep(0, 100)
obs = logistic(r=0.3, K=100, C=C)
plot(sim)

# error function
RSS = function(par, obs, ...) {
  # translating the 'par' vector to be our model
  r = par[1]
  K = par[2]
  # run the model
  sim = logistic(r=r, K=K, ...)
  # estimation method
  errors = obs - sim
  rss = sum(errors^2)
  return(rss)
}

RSS(par=c(0.4, 120), obs=obs, C=C)





