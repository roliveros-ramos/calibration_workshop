
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

RSS(par=c(0.31, 101), obs=obs, C=C)

optim(par=c(0.31, 101), fn=RSS, obs=obs, C=C)

optim(par=c(0.8, 500), fn=RSS, obs=obs, C=C)

opt0 = optim(par=c(0.8, 500), fn=RSS, obs=obs, C=C,
      control=list(maxit=500))

opt1 = optim(par=opt0$par, fn=RSS, obs=obs, C=C,
      control=list(maxit=500))

opt2 = optim(par=opt1$par, fn=RSS, obs=obs, C=C,
      control=list(maxit=500))


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



