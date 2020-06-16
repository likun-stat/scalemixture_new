Rcpp::sourceCpp('~/Desktop/Research/scalemixture/integration.cpp')
Rcpp::sourceCpp('~/Desktop/Research/scalemixture/likelihood.cpp')
source("~/Desktop/Research/scalemixture/scalemix_utils.R")
source("~/Desktop/Research/scalemixture/scalemix_likelihoods.R")

## ------------- 1. Marginal CDF ---------------
## Plot the CDF function
delta <- 0.8
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 20, length = 10000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
plot(x_vals, cdf_vals, type="l")
grid()
legend("topleft",pch=20,legend=expression(delta==0.8))

## Compare with asymptotic estimate
x_vals <- seq(1.01, 50, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.277
equiv<-function(x)   (delta/(2*delta-1))*x^{(delta-1)/delta} - (1-delta)*x^{-1}/(2*delta-1)
cdf_est<-equiv(x_vals)
plot(x_vals, 1-cdf_vals, type="l", col="grey",lwd=3, ylab="1-F(x)", 
     xlab=expression(x), main=expression(delta==0.8))
lines(x_vals, cdf_est,col="red")
grid(col="grey50")
legend('topright', lwd=c(3,1), col=c("grey","red"),legend = c("1-F(x)","Asymptotic estimate"))





## ------------- 2. Marginal PDF ---------------
## Compute the PDF function
delta <- 0.8
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 50, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
system.time(pdf_vals<-dmixture_me(x_vals, tau, delta))  #0.183

## Compare with the numerical results from CDF values
library(pspline)
deriv<-predict(sm.spline(x_vals, cdf_vals), x_vals, 1)


## Compare with the asymptotic results 
equiv.pdf<-function(x) {(1-delta)/(2*delta-1)*(x^{-1/delta}-x^{-2})}
pdf_est<-equiv.pdf(x_vals)


## Plot all the results together
layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
par(mai=c(0.1,0.62,0.82,0.32))
plot(x_vals, pdf_vals, type="l",col="grey30",lwd=3,ylim=c(min(pdf_vals),max(pdf_est)),
     xlab="",ylab="Marginal density", main=expression(delta==0.8), xaxt="n")
grid(col="grey40")
lines(x_vals,deriv,col="red")
lines(x_vals,pdf_est, col="blue")
legend('topright', lwd=c(3,1,1), col=c("grey30","red","blue"),legend = c("f(x)","Numerical derivative","Asymptotic estimate"))

par(mai=c(0.6,0.62,0.1,0.32))
plot(x_vals,pdf_est-pdf_vals,type='n',ylab="Diffrence",xlab=expression(x))
grid(ny=0,col="grey40")
abline(h=0)
lines(x_vals,pdf_est-pdf_vals,lwd=2,col="darkred")

layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
par(mai=c(1.02,0.82,0.82,0.42))



## ------------- 3. Find xrange ---------------
delta <- 0.3
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 40, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
plot(x_vals, cdf_vals, type="l")
grid()
legend("bottomright",pch=20,legend=expression(delta==0.3))

left<-which(cdf_vals<0.61&cdf_vals>0.59)[1]
right<-which(cdf_vals<0.82&cdf_vals>0.81)[1]
lines(c(x_vals[left],x_vals[left]),c(0,cdf_vals[left]),lty=2)
lines(c(x_vals[right],x_vals[right]),c(0,cdf_vals[right]),lty=2)

fd<-find_xrange_pmixture_me(cdf_vals[left], cdf_vals[right], c(20,30), tau, delta)
lines(c(fd[1],fd[1]),c(0,pmixture_me(fd[1], tau, delta)),col="red")
lines(c(fd[2],fd[2]),c(0,pmixture_me(fd[2], tau, delta)), col="red")




## ------------- 4. Quantile function for marginal x ---------------
delta <- 0.3
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 40, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
plot(x_vals, cdf_vals, type="l")
legend("bottomright",pch=20,legend=expression(delta==0.3))

p.vals<-seq(0.51,0.91,by=0.1)
q.vals<-qmixture.me.interp(p=p.vals, tau_sqd = tau, delta=delta)
abline(h=p.vals,lty=2,col="grey")
abline(v=q.vals,lty=3,col="blue")




## ------------- 5. Marginal transformation ---------------
Y<-scalemix.me.2.gpd(c(20,30),1,0.3,c(11,1,0))
gpd.2.scalemix.me(Y, tau_sqd=1, delta=0.3, theta.gpd=c(11,1,0))




## ------------- 6. Fix issues for large quantile --------------
pmixture_me_uni(4, tau, delta)
pmixture_me_old(4, tau, delta)
pmixture_me(4, tau, delta)



## ------------- 7. Use interpolation for PDF --------------
delta <- 0.3
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 50, length = 20000)
system.time(pdf_vals<-dmixture_me(x_vals, tau, delta))  #0.089
system.time(pdf_interp<-dmixture.me(x_vals, tau, delta))  #0.002


## Plot all the results together
plot(x_vals, pdf_vals, type="l",col="grey30",lwd=3,
     xlab="",ylab="Marginal density", main=expression(delta==0.3), xaxt="n")
grid(col="grey40")
lines(x_vals,pdf_interp,col="red")
legend('topright', lwd=c(3,1), col=c("grey30","red"),legend = c("f(x)","Numerical derivative"))



## ------------- 8. Use RcppArmadillo for likelihood / multivariate update --------------
eig2logdet_c(1:4)
eig2logdet(1:4)

library(fields)
S     <- cbind(seq(0, 1, length=200), rep(1, 200))
Cor   <- corr.fn(rdist(S), lambda = 0.5, gamma = 1)
eig.Sigma <- eigen(Cor, symmetric=TRUE)
V <- eig.Sigma$vectors
d <- eig.Sigma$values

eig2inv.quadform.vector(V, 1/d, 1:200)
eig2inv_quadform_vector(V, 1/d, 1:200)

eig2inv.times.vector(V, 1/d, 1:200)
eig2inv_times_vector(V, 1/d, 1:200)

M<-matrix(rnorm(600),ncol=3)
dmvn.eig(M, V, 1/d)
dmvn_eig(M, V, 1/d)




n.s<-200; n.t<-2; tau<-2; delta<-0.6
X <- matrix(NA, n.s, n.t)
X.s <- matrix(NA, n.s, n.t)
C.Cor <- chol(Cor)
u<-rep(NA,n.t)
R<-rep(NA,n.t)
for(t in 1:n.t){
  u[t] <-runif(1,0,1)
  R[t]<-pow(1/(1-u[t]),delta/(1-delta))
}
for(t in 1:n.t) {
  Z.t<-crossprod(C.Cor, rnorm(n.s))
  Z.to.W.s<-1/(1-pnorm(Z.t))
  X.s[ ,t] <- R[t]*Z.to.W.s
  X[ ,t] <- X.s[ ,t] + sqrt(tau)*rnorm(n.s)
}

X_s_likelihood_conditional(X.s[,1], R[1], V, d)
X.s.likelihood.conditional(X.s[,1], R[1], V, d)

X_s_likelihood_conditional_on_X(X.s[,1], X[,1], R[1], V, d, 0.04)
X.s.likelihood.conditional.on.X(X.s[,1], X[,1], R[1], V, d, 0.04)

system.time(X_s<-var_at_a_time_update_X_s(X[,1], R[1], V, d, 0.04,n_chain = 100)) #3.352  6.562  9.819 
system.time(X_s<-var.at.a.time.update.X.s(X[,1], R[1], V, d, 0.04,n.chain = 100)) #4.431  8.627  12.811


#########################################
## Update using the false normal proposal
#########################################

t<-1
X_s<-var_at_a_time_update_X_s(X[,t], R[t], V, d, 0.04,n_chain = 100)
plot(X[,t],xlab="s",ylab="X", main=sprintf("R=%.4f", R[t]))
grid()
lines(1:200, X_s$X.s, col="red", lwd=2)

regularized.d.inv <- 1/d * 1/R[t]^2 + 1/tau
prop.mean <- eig2inv.times.vector(V, 1/regularized.d.inv, X[ ,t]) / tau
prop.X.s <- prop.mean + eig2inv.times.vector(V, 1/sqrt(regularized.d.inv), rnorm(n.s))
lines(1:200, prop.X.s, col="blue", lwd=2)

thresh.X <- qmixture.me.interp(0.8, tau_sqd = tau, delta = delta) 
abline(h=thresh.X, lty=2, col="grey50")
legend("topright",lwd=2, col=c("red","blue"),legend=c("metropolis","normal"))
abline(h=R[t])



## ------------- 9. Test if OpenMP actually works --------------
system.time(s4 <- long_computation_omp(5000, 4))  # 0.043
system.time(tmp <- tryfor(1000000))  # 0.043
system.time(tmp <- tryfor_omp(1000000, 2))  # 0.023 
system.time(tmp <- tryfor_omp(1000000, 4))  # 0.013 
system.time(tmp <- tryfor_omp(1000000, 6))  # 0.008 
system.time(tmp <- tryfor_omp(1000000, 7))  # 0.0065
system.time(tmp <- tryfor_omp(1000000, 8))  # 0.008 

prob.below <- 0.8
theta.gpd <- c(11, 1, 0)
thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau, delta = delta) 
sum(X < thresh.X) / length(X)
cen <- X < thresh.X  

library(evd)
Y <- X
Y[cen] <- NA
Y[!cen] <- scalemix.me.2.gpd(x = X[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below)

X.imp <- X
X.star <- X.s

# OpenMP failure
#system.time(res <- var_at_a_time_update_X_s_cols(R, X.imp, tau, V, d, v_q=0.5, n_chain=100, threads=2))  # 8.998




## ------------- 10. Use foreach for parallelizaztion instead --------------
ptm<-proc.time()
for(t in 1:2){
  # Proposal
  # Treat the X process as the response, ignoring the marginal transformation
  # i.e. prop.X.s ~ X.s | X, R, other.params
  
  metropolis <- var.at.a.time.update.X.s(X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau, v.q = 0.5, n.chain = 100)
  metropolis$X.s
}
proc.time()-ptm


ptm<-proc.time()
library(doParallel)
library(foreach)
registerDoParallel(cores=1)
res<-foreach(t = 1:2, .combine="cbind") %dopar% {
  # Proposal
  # Treat the X process as the response, ignoring the marginal transformation
  # i.e. prop.X.s ~ X.s | X, R, other.params
  
  metropolis <- var.at.a.time.update.X.s(X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau, v.q = 0.5, n.chain = 100)
  metropolis$X.s
}
proc.time()-ptm

system.time(Xes <- X.s.update.mixture.me(R, Y, X.imp, X.star, cen, 
                                            prob.below, theta.gpd, delta,
                                            tau, V, d, v.q=0.5, n.chain = 100, thresh.X = thresh.X)$X.s) #8.706   215.607 

system.time(Xes <- X.s.update.mixture.me.par(R, Y, X.imp, X.star, cen, 
                                            prob.below, theta.gpd, delta,
                                            tau, V, d, v.q=0.5, n.chain = 100, thresh.X = thresh.X)$X.s) #4.371   110.953 88.217


## ------------- 11. Apply generic sampler to update each parameter --------------
# Run the imputation.animation.R before the MCMC run to generate data
# 1. delta
Res <- static.metr(z = R, starting.theta = delta, 
    likelihood.fn = delta.update.mixture.me.likelihood, prior.fn = interval.unif,
    hyper.params = 1, n.updates = 90000, prop.Sigma = NULL, sigma.m=NULL, verbose=FALSE, 
    Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, theta.gpd = theta.gpd, tau_sqd = tau)

Res <- adaptive.metr(z = R, starting.theta = delta,
    likelihood.fn = delta.update.mixture.me.likelihood, prior.fn = interval.unif,
    hyper.params = 1, n.updates = 90000, prop.Sigma = NULL, adapt.cov = FALSE, return.prop.Sigma.trace = TRUE,
    r.opt = .234, c.0 = 10, c.1 = .8, K = 10,
    Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, theta.gpd = theta.gpd, tau_sqd = tau)
plot(1:90000, Res$trace[,1], type="l", xlab = "Index", ylab = expression(delta))
abline(h=0.7, lty=2, col="red")


# 2. lambda, gamma (failed)
Res <- static.metr(z = R, starting.theta = c(lambda, gamma), 
    likelihood.fn = lam.gam.update.mixture.me.likelihood, prior.fn = lam.gam.prior,
    hyper.params = 2, n.updates = 90000, prop.Sigma = NULL, sigma.m=NULL, verbose=FALSE, 
    X.s = X.s, R = R, S = S)

Res <- adaptive.metr(z = R, starting.theta = c(1, 1.5),
    likelihood.fn = lam.gam.update.mixture.me.likelihood, prior.fn = lam.gam.prior,
    hyper.params = 2, n.updates = 10000, prop.Sigma = NULL, adapt.cov = TRUE, return.prop.Sigma.trace = FALSE,
    r.opt = .234, c.0 = 1, c.1 = .8, K = 100,
    X.s = X.s, R = R, S = S)


plot(1:10000, Res$trace[,1], type="l", xlab = "Index", ylab = expression(lambda))
plot(1:10000, Res$trace[,2], type="l", xlab = "Index", ylab = expression(gamma))


# 2a. rho (NA likelihood)
Res <- adaptive.metr(z = R, starting.theta = rho,
    likelihood.fn = rho.update.mixture.me.likelihood, prior.fn = rho.prior,
    hyper.params = 2, n.updates = 30000, prop.Sigma = NULL, adapt.cov = FALSE, return.prop.Sigma.trace = FALSE,
    r.opt = .234, c.0 = 10, c.1 = .8, K = 100,
    X.s = X.s, R = R, S = S)

plot(1:30000, Res$trace[,1], type="l", xlab = "Index", ylab = expression(rho))
abline(h=log(rho), lty=2, col="red")


# 3. tau
Res <- static.metr(z = R, starting.theta = tau, 
    likelihood.fn = tau.update.mixture.me.likelihood, prior.fn = tau.sqd.prior,
    hyper.params = c(0.1, 0.1), n.updates = 90000, prop.Sigma = NULL, sigma.m=NULL, verbose=FALSE, 
    Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, theta.gpd = theta.gpd)

Res <- adaptive.metr(z = R, starting.theta = tau,
    likelihood.fn = tau.update.mixture.me.likelihood, prior.fn = tau.sqd.prior,
    hyper.params = c(0.1, 0.1), n.updates = 30000, prop.Sigma = NULL, adapt.cov = FALSE, return.prop.Sigma.trace = FALSE,
    r.opt = .234, c.0 = 10, c.1 = .8, K = 10,
    Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, theta.gpd = theta.gpd)

plot(1:30000, Res$trace[,1], type="l", xlab = "Index", ylab = expression(tau^2))
abline(h=tau,lty=2,col="red")


# 4. theta.gpd: scale, shape
Res <- adaptive.metr(z = R, starting.theta = theta.gpd[2:3],
                     likelihood.fn = theta.gpd.update.mixture.me.likelihood, prior.fn = gpd.prior,
                     hyper.params = 0, n.updates = 9000, prop.Sigma = NULL, adapt.cov = TRUE, return.prop.Sigma.trace = FALSE,
                     r.opt = .234, c.0 = 20, c.1 = .8, K = 10,
                     Y =Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, tau_sqd = tau, loc = thresh, thresh.X=thresh.X)
tmp <- cov(Res$trace[4000:9000,])

Res <- adaptive.metr(z = R, starting.theta = theta.gpd[2:3],
    likelihood.fn = theta.gpd.update.mixture.me.likelihood, prior.fn = gpd.prior,
    hyper.params = 0, n.updates = 9000, prop.Sigma = tmp, adapt.cov = TRUE, return.prop.Sigma.trace = FALSE,
    r.opt = .234, c.0 = 20, c.1 = .8, K = 10,
    Y =Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, tau_sqd = tau, loc = thresh, thresh.X=thresh.X)

Res <- adaptive.metr(z = R, starting.theta = 1,
                     likelihood.fn = scale.gpd.update.mixture.me.likelihood, prior.fn = half.cauchy,
                     hyper.params = 1, n.updates = 9000, prop.Sigma = NULL, adapt.cov = FALSE, return.prop.Sigma.trace = FALSE,
                     r.opt = .234, c.0 = 20, c.1 = .8, K = 10,
                     shape=0, Y =Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, tau_sqd = tau, loc = thresh, thresh.X=thresh.X)
Res <- adaptive.metr(z = R, starting.theta = 1,
                     likelihood.fn = shape.gpd.update.mixture.me.likelihood, prior.fn = log.rho.prior,
                     hyper.params = 1, n.updates = 9000, prop.Sigma = NULL, adapt.cov = FALSE, return.prop.Sigma.trace = FALSE,
                     r.opt = .234, c.0 = 20, c.1 = .8, K = 10,
                     scale=1, Y =Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, tau_sqd = tau, loc = thresh, thresh.X=thresh.X)


plot(Res$trace[,1], type="l", xlab = "Index", ylab = "scale")
plot(Res$trace[-(1:100),2], type="l", xlab = "Index", ylab = "shape")


# 5. Rt
t <-2
Res <- adaptive.metr(z = R, starting.theta = R[t],
    likelihood.fn = Rt.update.mixture.me.likelihood, prior.fn = huser.wadsworth.prior,
    hyper.params = delta, n.updates = 30000, prop.Sigma = NULL, adapt.cov = FALSE, return.prop.Sigma.trace = FALSE,
    r.opt = .234, c.0 = 10, c.1 = .8, K = 10,
    X.s = X.s[,t], delta = delta, V = V, d = d)

plot(1:30000, Res$trace[1:30000,1], type="l", xlab = "Index", ylab = sprintf("R%i", t))
abline(h=R[t], lty=2, col="red")



## ------------- 12. Test giant scalemix_sampler that updates all parameter --------------
# Save the input for function arguments
initial.values <- list(delta = delta, rho=rho, tau=tau, theta.gpd=theta.gpd, prob.below=prob.below, X.s=X.s, R=R)
n.updates <- 2
thin <- 10
echo.interval <- 50
true.params <- list(delta = delta, rho=rho, tau=tau, theta.gpd=theta.gpd, prob.below=prob.below, X.s=X.s, R=R)
sigma.m=NULL; prop.Sigma=NULL
experiment.name="Huser-wadsworth"

save(Y, S, cen, thresh, initial.values, n.updates, thin, echo.interval, sigma.m, prop.Sigma, true.params, 
     experiment.name, file="input.RData")
load("input.RData")



## ------------- 13. Update X.s in Rcpp --------------
dgpd(c(12,13,14), loc=11, scale=1, shape=0, log=TRUE)
d_gpd(c(12,13,14), loc=11, scale=1, shape=0, islog=TRUE)

dgpd(c(123,135,146), loc=11, scale=1, shape=2, log=TRUE)
d_gpd(c(123,135,146), loc=11, scale=1, shape=2, islog=TRUE)

ptm=proc.time()
marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
                                                      theta.gpd, delta,
                                                      tau, thresh.X=thresh.X)
proc.time()-ptm

ptm=proc.time()
marg_transform_data_mixture_me_likelihood(as.vector(Y), as.vector(X), as.vector(X.s), as.vector(cen), prob.below,
                                              theta.gpd, delta,
                                              tau, thresh_X=thresh.X)
proc.time()-ptm

t<-2; tau_sqd<-9
res<- update_X_s_onetime(Y[,t], X[,t], X.s[,t], cen[,t], prob.below, theta.gpd, delta, tau_sqd, thresh.X, v.q[,t], R[t], V, d)
res

plot(X.s[,t],type="b",col="red")
points(res$X.s, type='b', col="blue")

ptm=proc.time()
X.s.res<-X.s.update.mixture.me.update.par.once.without.X.par(R, Y, X, X.s, cen,
                                                         prob.below, theta.gpd, delta,
                                                         tau, V, d, v.q=sigma.m$X.s, thresh.X=thresh.X)
proc.time()-ptm   #0.362
colSums(X.s.res$accepted)

ptm=proc.time()
X.s.res<-X.s.update.mixture.me.update.par.once.without.X(R, Y, X, X.s, cen,
                                                            prob.below, theta.gpd, delta,
                                                            tau, V, d, v.q=sigma.m$X.s, thresh.X=thresh.X)
proc.time()-ptm   #16.651
colSums(X.s.res$accepted)
