source("~/Desktop/Research/scalemixture_new/scalemix_utils.R")
source("~/Desktop/Research/scalemixture_new/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture_new/scalemix_priors.R")
source("~/Desktop/Research/scalemixture_new/generic_samplers.R")
source("~/Desktop/Research/scalemixture_new/scalemix_sampler_02.R")

library(fields)   # For rdist

# ------------ 1. Simulation settings -------------
n.s <- 100        # Number of sites
n.t <- 40         # Number of time points
tau <- 4          # Nugget SD
delta <- 0.7      # For R
range <- 1        # Matern range
nu <-  3/2        # Matern smoothness



# -------------- 2. Generate fake data -----------------
S     <- cbind(runif(n.s, 0, 1), runif(n.s, 0, 1))
# Cor   <- corr.fn(rdist(S), lambda = lambda, gamma = gamma)
Cor   <- corr.fn(rdist(S), theta=c(range,nu))
C.Cor <- chol(Cor)

# set.seed(3333)
u<-rep(NA,n.t)
R<-rep(NA,n.t)
for(t in 1:n.t){
  u[t] <-runif(1,0,1)
  R[t]<-pow(1/(1-u[t]),delta/(1-delta))
}


X <- matrix(NA, n.s, n.t)
X.s <- matrix(NA, n.s, n.t)
for(t in 1:n.t) {
  Z.t<-crossprod(C.Cor, rnorm(n.s))
  Z.to.W.s<-1/(1-pnorm(Z.t))
  X.s[ ,t] <- R[t]*Z.to.W.s
  X[ ,t] <- X.s[ ,t] + sqrt(tau)*rnorm(n.s)
  
}


prob.below <- 0.95
# Threshold for fitting
thresh <- 11
theta.gpd <- c(thresh, 1, 0)

thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau, delta = delta) 
sum(X < thresh.X) / length(X)
cen <- X < thresh.X  



## ------------ 3. Marginal transformation -----------------
library(evd)

Y <- X
Y[cen] <- NA
Y[!cen] <- scalemix.me.2.gpd(x = X[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below)




## --------------- 4. Saving initial values -------------------
initial.values <- list(delta = delta-0.05, tau=tau+1, theta.gpd=c(thresh,theta.gpd[2],theta.gpd[3]+0.2), theta.c=c(range,nu), prob.below=prob.below, X.s=X.s, R=R, S=S)
n.updates <- 50000
thin <- 10
echo.interval <- 50
true.params <- list(delta = delta, tau=tau, theta.gpd=theta.gpd, theta.c=c(range,nu), prob.below=prob.below, X.s=X.s, R=R, S=S)
save(Y, X, thresh.X, cen, true.params, file="Initial.RData")




## --------------- 5. Running Metropolis -------------------
scalemix.sampler.02(Y=Y, S=S, cen=cen, thresh=thresh,
                      initial.values=initial.values,
                      n.updates=n.updates, thin=thin,
                      experiment.name="Huser-wadsworth-bigdelta",
                      echo.interval=echo.interval,
                      sigma.m=NULL, prop.Sigma=NULL,
                      true.params=true.params, sd.ratio=NULL, lower.prob.lim=0.5)

