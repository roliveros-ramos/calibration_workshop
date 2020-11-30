

logistic = function(r, K, C, B0=10, T=100) {

  time = seq_len(T)
  B = rep(NA, T)
  B[1] = B0

  for(t in time[-T]) {
    B[t+1] = B[t] + r*B[t]*(1-B[t]/K) - C[t]
  }

  return(B)

}

sim = logistic(r=0.3, K=100, C=rep(0, 100))
plot(sim)







