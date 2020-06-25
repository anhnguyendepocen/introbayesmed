# Exercise 1 -----

### 1) ----
sample_50 <- rnorm(50, mean = 2, sd = 3)
mean(sample_50)

sample_20000 <- rnorm(20000, mean = 2, sd = 3)
mean(sample_20000)

### 2a) ----

roulette_coord <- function(ngrid = 35){
  x <- sample(x = 0:ngrid, size = 1)
  y <- sample(x = 0:ngrid, size = 1)
  return(c(x, y))
}
roulette_coord()

### 2b) ----

inside_disk_fun <- function(p, ngrid = 35){
  d <- sqrt((p[1] - ngrid/2)^2 + (p[2] - ngrid/2)^2)
  return(d <= ngrid/2)
}

inside_disk_fun(p = c(16,17))
inside_disk_fun(p = c(0,0))

### 2c) ----

piMC <- function(in_disk){
  return(mean(4 * in_disk))
}

piMC(c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE))

### 2d) ----

# Grid size (resolution)
ngrid <- 1000

# Monte Carlo sample size
npoints <- 1000

# Points generation
pp <- matrix(NA, ncol = 2, nrow = npoints)
for(i in 1:nrow(pp)){
  pp[i, ] <- roulette_coord(ngrid)
}

# Estimate pi
in_disk <- apply(X = pp, MARGIN = 1, FUN = inside_disk_fun, ngrid = ngrid)
piMC(in_disk)

# Plot
## first we initialise an empty plot with the right size 
## using argument
plot(x = pp[, 1], y = pp[, 2], 
     xlim = c(0, ngrid), ylim = c(0, ngrid), 
     axes = 0, xlab = "x", ylab = "y", 
     type="n")
## we tick the x and then y axes from 1 to ngrid
axis(1, at=c(0:ngrid))
axis(2, at=c(0:ngrid))
## we add a square around the plot
box() 
## we plot the grid (using dotted lines thanks to the argument `lty = 3`)
## onto which the points are sample
for(i in 0:ngrid){
  abline(h=i, lty = 3)
  abline(v=i, lty = 3)
}
## we add the sampled points
lines(x = pp[, 1], y = pp[, 2], 
      xlim = c(0, ngrid), ylim = c(0, ngrid), 
      xlab = "x", ylab = "y",
      type= "p", pch=16)
## we add the circle display
x.cercle <- seq(0, ngrid, by = 0.1)
y.cercle <- sqrt((ngrid/2)^2 - (x.cercle-ngrid/2)^2)
lines(x.cercle, y = y.cercle + ngrid/2, col = "red")
lines(x.cercle, y = - y.cercle + ngrid/2, col = "red")
## finally we color in red the points sampled inside the disk
lines(x = pp[in_disk, 1], y = pp[in_disk, 2], 
      xlim = c(0, ngrid), ylim = c(0, ngrid), 
      xlab = "x", ylab = "y",
      type= "p", pch=16, col="red", cex=0.7)





























# Exercise 2 -----

my_sampler_expodist <- function(n, lambda){
  u <- runif(n)
  e <- -1/lambda * log(1-u)
  return(e)
}

nsamp <- 10000
my_samp <- my_sampler_expodist(n = nsamp, lambda = 10)
r_samp <- rexp(n = nsamp, rate = 10)
plot(density(my_samp), main = "", ylab="Probability density")
lines(density(r_samp), col = "red", lty = 2)
legend("topright", 
       legend = c("Inverse transform", "R rexp()"), 
       lty=c(1,2), 
       col=c("black", "red"))


# Exercise 3 ----

### 1) ----
post_num_hist <- function(theta, log = FALSE) {
  
  n <- 493472 # the data
  S <- 241945 # the data
  
  if (log) {
    num <- S * log(theta) + (n - S) * log(1 - theta)
  } else {
    num <- theta^S * (1 - theta)^(n - S) # the numerator of the posterior
  }
  return(num) # the output of the funtion
}

post_num_hist(0.2, log=TRUE)
post_num_hist(0.6, log=TRUE)
post_num_hist(0.49, log=TRUE)


### 2) ----
myMH <- function(niter, post_num) {
  x_save <- numeric(length = niter)  #create a vector of 0s of length niter to store the sampled values
  alpha <- numeric(length = niter)  #create a vector of 0s of length niter to store the acceptance probabilities
  # initialise x0
  x <- runif(n = 1, min = 0, max = 1)
  # accpetance-rejection loop
  for (t in 1:niter) {
    # sample y from the proposal (here uniform prior)
    y <- runif(n = 1, min = 0, max = 1)
    # compute the acceptance-rejection probability
    alpha[t] <- min(1, exp(post_num(y, log = TRUE) - post_num(x, log = TRUE)))
    # accept or reject
    u <- runif(1)
    if (u <= alpha[t]) {
      x_save[t] <- y
    } else {
      x_save[t] <- x
    }
    # update the current value
    x <- x_save[t]
  }
  return(list(theta = x_save, alpha = alpha))
}

### 3) ----

sampleMH <- myMH(2000, post_num = post_num_hist)

par(mfrow = c(2, 2))
plot(density(sampleMH$theta[-c(1:500)]), col = "red", xlim = c(0, 1), ylab = "Posterior probability density", 
     xlab = expression(theta), main = "")
curve(dbeta(x, 241945 + 1, 251527 + 1), from = 0, to = 1, n = 10000, add = TRUE)
legend("topright", c("M-H", "theory"), col = c("red", "black"), lty = 1)

plot(density(sampleMH$theta[-c(1:500)]), col = "red", ylab = "Posterior probability density", 
     xlab = expression(theta), main = "Zoom")
curve(dbeta(x, 241945 + 1, 251527 + 1), from = 0, to = 1, n = 10000, add = TRUE)
legend("topright", c("M-H", "theory"), col = c("red", "black"), lty = 1)

plot(sampleMH$alpha, type = "h", xlab = "Iteration", ylab = "Probabilité d'acceptation", 
     ylim = c(0, 1), col = "springgreen")
plot(sampleMH$theta, type = "l", xlab = "Iteration", ylab = expression(paste("Sampled value for ", 
                                                                             theta)), ylim = c(0, 1))

### 4) ----
post_num_beta <- function(theta, a = 3, b = 3, log = TRUE) {
  
  n <- 100  #number of trials (births)
  S <- 49  #number of success (feminine births)
  
  if (log) {
    num <- (a + S - 1) * (log(theta)) + (b + n - S - 1) * log(1 - theta)
  } else {
    num <- theta^(a + S - 1) * (1 - theta)^(b + n - S - 1)
  }
  return(num)
}

myMH_betaprior <- function(niter, post_num, a = 3, b = 3) {
  
  theta_save <- numeric(length = niter)  # vector of theta values ready to be saved
  alpha <- numeric(length = niter)  # vector of alpha values ready to be saved
  
  # initialise theta
  theta <- runif(n = 1, min = 0, max = 1)
  for (t in 1:niter) {
    
    # sample a value from the proposal (beta prior)
    theta_prop <- rbeta(n = 1, a, b)
    
    # compute acceptance-rejection probability
    alpha[t] <- min(1, exp(post_num(theta_prop, a = a, b = b, log = TRUE) - post_num(theta, 
                                                                                     a = a, b = b, log = TRUE)))
    # acceptance-rejection step
    u <- runif(1)
    if (u <= alpha[t]) {
      theta <- theta_prop  # acceptance of theta_prop as new current value
    }
    # saving the current value of theta
    theta_save[t] <- theta
  }
  return(list(theta = theta_save, alpha = alpha))
}

sampleMH <- myMH_betaprior(10000, post_num = post_num_beta)


par(mfrow = c(2, 2))
plot(density(sampleMH$theta[-c(1:500)]), col = "red", xlim = c(0, 1), ylab = "Posterior probability density", 
     xlab = expression(theta), main = "")
curve(dbeta(x, 49 + 1, 51 + 1), from = 0, to = 1, add = TRUE)
legend("topright", c("M-H", "theory"), col = c("red", "black"), lty = 1)
plot.new()
plot(sampleMH$alpha, type = "h", xlab = "Iteration", ylab = "Acceptance probability", 
     ylim = c(0, 1), col = "springgreen")
plot(sampleMH$theta, type = "l", xlab = "Iteration", ylab = expression(paste("Sampled value for ", 
                                                                             theta)), ylim = c(0, 1))
### 5) ----

post_num_beta_hist <- function(theta, a = 3, b = 3, log = TRUE) {
  
  n <- 493472  #number of trials (births)
  S <- 241945  #number of success (feminine births)
  
  if (log) {
    num <- (a + S - 1) * log(theta) + (b + n - S - 1) * log(1 - theta)
  } else {
    num <- theta^(a + S - 1) * (1 - theta)^(b + n - S - 1)
  }
  return(num)
}

myMH_betaprior_randomwalk <- function(niter, post_num, a = 3, b = 3) {
  
  theta_save <- numeric(length = niter)
  alpha <- numeric(length = niter)
  
  # initialise theta
  theta <- runif(n = 1, min = 0, max = 1)
  
  for (t in 1:niter) {
    # sample a value from the proposal (random walk)
    theta_prop <- theta + runif(1, -0.01, 0.01)  #rnorm(1, mean = 0, sd = 0.02)
    
    # compute acceptance-rejection probability
    alpha[t] <- min(1, exp(post_num(theta_prop, a = a, b = b, log = TRUE) - post_num(theta, 
                                                                                     a = a, b = b, log = TRUE)))
    
    # acceptance-rejection step
    u <- runif(1)
    if (u <= alpha[t]) {
      theta <- theta_prop  # accept theta_prop and update current value
    }
    
    # save current value
    theta_save[t] <- theta
  }
  
  return(list(theta = theta_save, alpha = alpha))
}

sampleMH <- myMH_betaprior_randomwalk(20000, post_num = post_num_beta_hist)


par(mfrow = c(2, 2))
plot(density(sampleMH$theta[-c(1:1000)]), col = "red", ylab = "Posterior probability density", 
     xlab = expression(theta), main = "")
curve(dbeta(x, 241945 + 1, 251527 + 1), from = 0, to = 1, n = 10000, add = TRUE)
legend("topright", c("M-H", "theory"), col = c("red", "black"), lty = 1)
plot(sampleMH$alpha, type = "h", xlab = "Iteration", ylab = "Acceptance probability", 
     ylim = c(0, 1), col = "springgreen")
plot(sampleMH$theta, type = "l", xlab = "Iteration", ylab = expression(paste("Sampled value for ", 
                                                                             theta)), ylim = c(0, 1))
plot(sampleMH$theta, type = "l", xlab = "Iteration", main = "Zoom", ylab = expression(paste("Sampled value for ", 
                                                                                            theta)), ylim = c(0.45, 0.55))























# Exercise 4 ----

### 0) ----
library(rjags)

### 1) ----
N <- 50 # the number of observations
obs <- rnorm(n = N, mean = 2, sd = 3) # the (fake) observed data

### 3) ----
myfirstjags <- jags.model("normalBUGSmodel.txt", 
                          data = list('obs' = obs, 'N' = length(obs)))

### 4) ----
res <- coda.samples(model = myfirstjags, 
                    variable.names = c("mu", "sigma"), 
                    n.iter = 2000)

### 5) ----
plot(res)
resum <- summary(res)
resum
resum$statistics["mu", "Mean"]
resum$statistics["sigma", "Mean"]
resum$quantiles["mu", "50%"]
resum$quantiles["sigma", "50%"]
resum$quantiles["mu", c(1,5)]
resum$quantiles["sigma", c(1,5)]

### 6) ----
library(coda)

### 7) ----
myjags2 <- jags.model("normalBUGSmodel.txt", data = list('obs' = obs, 'N' = N),
                      n.chains = 3, inits = list(list("mu" = 0, "tau" = 1), 
                                                 list("mu" = -10, "tau" = 1/100), 
                                                 list("mu" = 100, "tau" = 1/10)))


res2 <- coda.samples(model = myjags2, 
                     variable.names = c('mu', 'sigma'), 
                     n.iter = 5000)
plot(res2)
library(lattice)
densityplot(res2)
summary(res2)

#### 8) ----
gelman.plot(res2)

### 9) ----
acfplot(res2)
#par(mfrow=c(3,2))
#autocorr.plot(res2, ask = FALSE, auto.layout = FALSE) # one graph per parameter and per chain

### 10) ----
cumuplot(res2)

### 11) ----
crosscorr.plot(res2)

### 12) ----
library(HDInterval)
hdCI <- HDInterval::hdi(res2)
hdCI
symCI <- summary(res2)$quantiles[ ,c(1,5)]
symCI

symCI[, 2] - symCI[, 1]
hdCI[2, ] - hdCI[1, ]



# Exercise 5 ----

### 3) ----

#Individual data
ycontrol <- c(rep(0, 125-57), rep(1, 57))
yecmo <- c(rep(0, 124-44), rep(1, 44))

#sampling
library(rjags)
goligher_jags_indiv <- jags.model(file = "goligherBUGSmodel_indiv.txt", 
                                  data = list("ycontrol" =  ycontrol, 
                                              "ncontrol" = length(ycontrol),
                                              "yecmo" =  yecmo, 
                                              "necmo" = length(yecmo)
                                  ), 
                                  n.chains = 3)
res_goligher_indiv <- coda.samples(model = goligher_jags_indiv, 
                                   variable.names = c('pc', 'RR'), 
                                   n.iter = 20000)

#postprocessing
res_goligher_burnt_indiv <- window(res_goligher_indiv, start=5001) # remove burn-in for Markov chain convergence

#Population data
zcontrol <- 57
zecmo <- 44
#sampling
goligher_jags_pop <- jags.model(file = "goligherBUGSmodel_pop.txt", 
                                data = list("zcontrol" =  zcontrol, 
                                            "ncontrol" = 125,
                                            "zecmo" =  zecmo, 
                                            "necmo" = 124
                                ), 
                                n.chains = 3)
res_goligher_pop <- coda.samples(model = goligher_jags_pop, 
                                 variable.names = c('pc', 'RR'), 
                                 n.iter = 20000)

#post-processing
res_goligher_burnt_pop <- window(res_goligher_pop, start=5001) # remove burn-in for Markov chain convergence

### 4) ----
effectiveSize(res_goligher_burnt_pop)
plot(res_goligher_burnt_pop)
gelman.plot(res_goligher_burnt_pop)
acfplot(res_goligher_burnt_pop)
par(mfrow=c(3, 2))
cumuplot(res_goligher_burnt_pop, ask=FALSE, auto.layout = FALSE)
par(mfrow=c(1, 1))

summary(res_goligher_burnt_pop)
summary(res_goligher_burnt_indiv)


# shortest 95% Credibility interval:
HDInterval::hdi(res_goligher_burnt_pop) 

# posterior porbability of RR <1:
mean(c(sapply(res_goligher_burnt_pop, "[", , 1))<1)


# Exercise 6 ----

### 0) ----
library(bayesmeta)
data(CrinsEtAl2014)

### 2) ----
library("metafor")
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
crins.es[,c("publication", "yi", "vi")]

### 3) ----
res_crins_bayesmeta <- bayesmeta(y=crins.es$yi, sigma=crins.es$vi, 
                                 labels=crins.es$publication,
                                 tau.prior = function(t){dunif(t, max=4)},
                                 mu.prior = c(0, 4),
                                 interval.type = "central")
summary(res_crins_bayesmeta)
plot(res_crins_bayesmeta)

### 5) ----
# Sampling
library(rjags)
crins_jags_res <- jags.model(file = "crinsBUGSmodel.txt", 
                             data = list("logOR" =  crins.es$yi,
                                         "sigma" = sqrt(crins.es$vi),
                                         "N" = length(crins.es$yi)
                             ), 
                             n.chains = 3)
res_crins_jags_res <- coda.samples(model = crins_jags_res, 
                                   variable.names = c("mu", "tau"), 
                                   n.iter = 20000)

# Postprocessing
res_crins_jags_res <- window(res_crins_jags_res, start=5001) # remove burn-in for Markov chain convergence
summary(res_crins_jags_res)
HDInterval::hdi(res_crins_jags_res)
plot(res_crins_jags_res)


