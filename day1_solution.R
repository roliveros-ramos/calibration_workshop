
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

RSS(par=c(0.31, 800), obs=obs, C=C, B0=500)
RSS2(par=c(0.31, 800, 500), obs=obs, C=C)

# using optim
optim(par=c(0.31, 101), fn=RSS, obs=obs, C=C, B0=500)
optim(par=c(0.31, 800), fn=RSS, obs=obs, C=C, B0=500)

B0s = seq(from=300, to=1000, by=10)

errors = rep(NA, length(B0s))
for(i in seq_along(B0s)) {
  B0 = B0s[i]
  errors[i] = optim(par=c(0.31, 800), fn=RSS, obs=obs, C=C, B0=B0)$value
}

plot(B0s, errors, type="l")


# using calibrar
calibrate(par=c(0.31, 101), fn=RSS, obs=obs, C=C, B0=500)
calibrate(par=c(0.31, 800), fn=RSS, obs=obs, C=C, B0=500)

B0s = seq(from=300, to=1000, by=10)

errors = rep(NA, length(B0s))
for(i in seq_along(B0s)) {
  B0 = B0s[i]
  errors[i] = calibrate(par=c(0.31, 800), fn=RSS, obs=obs, C=C, B0=B0)$value
}

plot(B0s, errors, type="l")


# estimating B0


# using optim
optim(par=c(0.31, 101, 200), fn=RSS2, obs=obs, C=C)
optim(par=c(0.31, 800, 200), fn=RSS2, obs=obs, C=C)

calibrate(par=c(0.31, 800, 200), fn=RSS2, obs=obs, C=C)
calibrate(par=c(0.31, 800, 200), fn=RSS2, obs=obs, C=C,
          method="BFGS")
calibrate(par=c(0.31, 800, 200), fn=RSS2, obs=obs, C=C,
          method="CG")
calibrate(par=c(0.31, 800, 200), fn=RSS2, obs=obs, C=C,
          method="SANN")
calibrate(par=c(0.31, 800, 200), fn=RSS2, obs=obs, C=C,
          method="LINB")

