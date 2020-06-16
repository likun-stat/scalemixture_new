################################################################################
## The log likelihood of the data, where the data comes from a scale mixture
## of Gaussians, transformed to GPD (matrix/vector input)
##
## NOT ACTUALLY depending on X, just calculated Fx^-1(Fy(Y)) ahead of time
##
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X ................................. a matrix of Y, transformed to HOT 
##                                     scale mixture of Gaussian
## X.s ............................... X, but without the measurement error
## cen 
## prob.below
## theta.gpd
## delta
## tau_sqd
##

marg.transform.data.mixture.me.likelihood <- function(Y,  X, X.s, cen, prob.below, theta.gpd, delta, tau_sqd, 
                                                        thresh.X=NULL) {
  if (!is.matrix(Y)) Y <- matrix(Y, ncol=1)
  if (!is.matrix(X)) X <- matrix(X, ncol=1)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  
  ll <- matrix(0, nrow(Y), ncol(Y))
  
  loc <- theta.gpd[1]
  scale <- theta.gpd[2]
  shape <- theta.gpd[3]
  
  if (is.null(thresh.X)) thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau_sqd, delta = delta)
  
  if (sum(cen) > 0)  {
    ll[cen] <- pnorm(thresh.X, mean=X.s[cen], sd=sqrt(tau_sqd), log.p=TRUE)
  }
  
  if (sum(!cen) > 0) {
    ll[!cen] <- dnorm(X[!cen], mean=X.s[!cen], sd=sqrt(tau_sqd), log=TRUE) +
      dgpd(Y[!cen], loc=loc, scale=scale, shape=shape, log=TRUE)  -
      dmixture.me(X[!cen], tau_sqd = tau_sqd, delta = delta, log=TRUE)
  }
  
  which <- is.na(ll)
  if(any(which)) ll[which] <- -Inf  # Normal density larger order than marginal density of scalemix
  return(sum(ll)) 
}

#                                                                              #
################################################################################



################################################################################
## Updates the X, just by drawing from the log likelihood of the
## data, censored at the threshold (on the Y scale), where the data comes
## from a Huser-wadsworth scale mixture
## of Gaussians, transformed to GPD
##
## X.s ............................... X, but without the measurement error
## cen 
## prob.below
## theta.mix
## theta.gaussian
##


# Only for Y's uncensored
gpd.2.scalemix.me <- function(y, tau_sqd, delta, theta.gpd, prob.below=0) {
  require(evd)
  
  thresh <- theta.gpd[1]
  scale <- theta.gpd[2]
  shape <- theta.gpd[3]
  
  # unifs <- pgpd(y, loc=thresh, scale=scale, shape=shape)
  unifs <- (1-prob.below) * pgpd(y, loc=thresh, scale=scale, shape=shape) + prob.below
  scalemixes <- qmixture.me.interp(unifs, tau_sqd = tau_sqd, delta = delta, n.x=500)

  return(scalemixes)
}


X.update<-function(Y, cen, X.s, delta, tau_sqd, theta.gpd, prob.below){
  X <- matrix(0, nrow(X.s), ncol(X.s))
  sd <- sqrt(tau_sqd)
  
  ## Get the threshold
  thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau_sqd, delta = delta)
  
  # Update X[cen]
  if (sum(cen) > 0)  {
    X[cen] <- truncnorm::rtruncnorm(1, mean=X.s[cen], sd=sd, a=-Inf, b=thresh.X)
  }
  
  if (sum(!cen) > 0) {
    X[!cen] <- gpd.2.scalemix.me(y=Y[!cen], tau_sqd=tau_sqd, delta=delta, 
                                 theta.gpd=theta.gpd, prob.below=prob.below)
}
  
  return(X)
}

#                                                                              #
################################################################################







################################################################################
## The log likelihood of X.s, when it is conditioned on scaling factor R and Σ(λ,γ)
##
## X.s ............................... X, but without the measurement error
## R 
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

X.s.likelihood.conditional<-function(X.s, R, V, d){
  if(any(X.s<R)) return(-Inf) else{
    X.s.to.Z <- qnorm(1-R/X.s)
    loglik <- -0.5*as.vector(eig2inv.quadform.vector(V, 1/d, X.s.to.Z))+0.5*sum(X.s.to.Z^2)-2*sum(log(X.s))-0.5*sum(log(d))+length(X.s)*log(R)
    return(loglik)
  }
}


# X.s.likelihood.conditional.on.X<-function(X.s, X, R, V, d, tau_sqd){
#   if(any(X.s<R)) return(-Inf) else{
#     X.s.to.Z <- qnorm(1-R/X.s)
#     loglik <- -0.5*as.vector(eig2inv.quadform.vector(V, 1/d, X.s.to.Z))+0.5*sum(X.s.to.Z^2)-2*sum(log(X.s))-0.5*sum(log(d))+length(X.s)*log(R)
#     loglik <- loglik - 0.5*sum((X.s-X)^2)/tau_sqd
#     return(loglik)
#   }
# }
# 
# var.at.a.time.update.X.s <- function(X, R, V, d, tau_sqd, v.q=0.5, n.chain=100){
#   X <- as.vector(X)
#   n.s <- length(X)
#   X.s <- X
#   X.s[which(X.s<R)] <- R + abs(rnorm(length(which(X.s<R)),sd=0.5)) # X.s has to be larger than R
#   accept <- rep(0, n.s)
#   
#   for(i in 1:n.chain){
#     for(iter in 1:n.s){
#       X.s.update<-X.s
#       X.s.update[iter]<-X.s[iter]+rnorm(1,0,sd=v.q)
#       log.num <- X.s.likelihood.conditional.on.X(X.s.update, X=X, R=R, V=V, d=d, tau_sqd=tau_sqd)
#       log.denom <- X.s.likelihood.conditional.on.X(X.s, X=X, R=R, V=V, d=d, tau_sqd=tau_sqd)
#       
#       r <- exp(log.num - log.denom)
#       
#       if(runif(1) < r){    
#         X.s <- X.s.update
#         accept[iter]<-accept[iter]+1
#       }
#     }
#   }
#   
#   return(list(X.s=X.s, accept=accept))
# }

#                                                                              #
################################################################################





################################################################################
##  Updates the X.s, for the generic Metropolis sampler
## Samples from the scaled Gaussian process (update the smooth process).
## The mixing distribution comes from from the Huser-wadsworth scale mixing distribution.
## The PROPOSAL will be the conjugate update from the model that treats the
## X process (i.e. X.s, but with measurement error) as the response, ignoring
## the marginal transformation.  Then the Metropolis part either rejects or
## accepts the entire draw.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters (theta in dhuser.thibaud):
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## theta.mix
## tau_sqd
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

# X.s.update.mixture.me <- function(R, Y, X, X.s, cen, 
#                                   prob.below, theta.gpd, delta,
#                                   tau_sqd, V, d, v.q=0.5, n.chain=100,
#                                   thresh.X=NULL) {
#   
#   n.s <- nrow(Y)
#   n.t <- ncol(Y)
#   
#   accepted <- rep(FALSE, n.t)  
#   
#   
#   if (is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
#   
#   for(t in 1:n.t){
#     # Proposal
#     # Treat the X process as the response, ignoring the marginal transformation
#     # i.e. prop.X.s ~ X.s | X, R, other.params
#     
#     metropolis <- var.at.a.time.update.X.s(X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd, v.q = v.q, n.chain = n.chain)
#     prop.X.s <- metropolis$X.s
#     
#     
#     # M-H ratio    
#     log.rat <- 
#       X.s.likelihood.conditional.on.X(X.s[,t], X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd)+ # Prop density of current value
#       marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], prop.X.s,               # Likelihood of proposal
#                                                 cen[ ,t], prob.below,
#                                                 theta.gpd, delta,
#                                                 tau_sqd, thresh.X=thresh.X) +
#       X.s.likelihood.conditional(prop.X.s, R[t], V, d) -                                # Likelihood of proposal
#       X.s.likelihood.conditional.on.X(prop.X.s, X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd)-# Prop density of proposal
#       marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], X.s[ ,t],               # Likelihood of current value
#                                                 cen[ ,t], prob.below,
#                                                 theta.gpd, delta,
#                                                 tau_sqd, thresh.X=thresh.X) -
#       X.s.likelihood.conditional(X.s[ ,t], R[t], V, d)                                  # Likelihood of current value
#     
#     
#     
#     
#     if (runif(1) < exp(log.rat)) {
#       X.s[ ,t] <- prop.X.s
#       accepted[t] <- TRUE
#     }
#   }
#   
#   return(list(X.s=X.s, accepted=accepted))
# }
# 
# 
# 
# X.s.update.mixture.me.par <- function(R, Y, X, X.s, cen, 
#                                   prob.below, theta.gpd, delta,
#                                   tau_sqd, V, d, v.q=0.5, n.chain=100,
#                                   thresh.X=NULL) {
#   
#   library(doParallel)
#   library(foreach)
#   n.s <- nrow(Y)
#   n.t <- ncol(Y)
#   registerDoParallel(cores=n.t)
#   
#   accepted <- rep(FALSE, n.t)  
#   
#   
#   if (is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
#   
#   Res <- foreach (t = 1:n.t, .combine = "cbind") %dopar% {
#     # Proposal
#     # Treat the X process as the response, ignoring the marginal transformation
#     # i.e. prop.X.s ~ X.s | X, R, other.params
#     
#     metropolis <- var_at_a_time_update_X_s(X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd, v_q = v.q, n_chain =  n.chain)
#     prop.X.s <- metropolis$X.s
#     res <- c(X.s[,t],0)
#     
#     # M-H ratio    
#     log.rat <- 
#       X.s.likelihood.conditional.on.X(X.s[,t], X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd)+ # Prop density of current value
#       marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], prop.X.s,               # Likelihood of proposal
#                                                 cen[ ,t], prob.below,
#                                                 theta.gpd, delta,
#                                                 tau_sqd, thresh.X=thresh.X) +
#       X.s.likelihood.conditional(prop.X.s, R[t], V, d) -                                # Likelihood of proposal
#       X.s.likelihood.conditional.on.X(prop.X.s, X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd)-# Prop density of proposal
#       marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], X.s[ ,t],               # Likelihood of current value
#                                                 cen[ ,t], prob.below,
#                                                 theta.gpd, delta,
#                                                 tau_sqd, thresh.X=thresh.X) -
#       X.s.likelihood.conditional(X.s[ ,t], R[t], V, d)                                  # Likelihood of current value
#     
#     
#     
#     
#     if (runif(1) < exp(log.rat)) {
#       res<-c(prop.X.s,1)
#     }
#     res
#   }
#   
#   return(list(X.s=Res[1:n.s, ], accepted=Res[n.s+1, ]))
# }
# 
# 
# X.s.update.mixture.me.update.par.once.without.X <- function(R, Y, X, X.s, cen,
#                                                             prob.below, theta.gpd, delta,
#                                                             tau_sqd, V, d, v.q=NULL, 
#                                                             thresh.X=NULL){
#   library(doParallel)
#   library(foreach)
#   n.s <- nrow(Y)
#   n.t <- ncol(Y)
#   registerDoParallel(cores=n.t)
#   
#   accepted <- matrix(0,nrow=n.s, ncol=n.t) 
#   if(is.null(v.q)) v.q<-matrix(2.4^2,nrow=n.s, ncol=n.t)
#   
#   if(is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
#   
#   for(iter in 1:n.s){
#     Res<-foreach(t = 1:n.t, .combine = "cbind")%dopar%{
#       prop.X.s<-X.s[ ,t]
#       prop.X.s[iter]<-X.s[iter,t]+rnorm(1,0,sd=v.q[iter,t])
#       res<-c(X.s[,t],0)
#       
#       # M-H ratio    
#       log.rat <- 
#         marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], prop.X.s,               # Likelihood of proposal
#                                                   cen[ ,t], prob.below,
#                                                   theta.gpd, delta,
#                                                   tau_sqd, thresh.X=thresh.X) +
#         X_s_likelihood_conditional(prop.X.s, R[t], V, d) -                                # Likelihood of proposal
#         marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], X.s[ ,t],               # Likelihood of current value
#                                                   cen[ ,t], prob.below,
#                                                   theta.gpd, delta,
#                                                   tau_sqd, thresh.X=thresh.X) -
#         X_s_likelihood_conditional(X.s[ ,t], R[t], V, d)                                  # Likelihood of current value
#       
#       if(is.na(log.rat)) log.rat<- -Inf
#       if(runif(1) < exp(log.rat)){
#         res<-c(prop.X.s,1)
#       }
#       res
#     }
#     
#     X.s<-Res[1:n.s, ]
#     accepted[iter,]<-Res[n.s+1,]
#   }
#   
#   return(list(X.s=X.s, accepted=accepted))
# }
# 

X.s.update.mixture.me.update.par.once.without.X.par <- function(R, Y, X, X.s, cen,
                                                            prob.below, theta.gpd, delta,
                                                            tau_sqd, V, d,
                                                            thresh.X=NULL){
  library(doParallel)
  library(foreach)
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  registerDoParallel(cores=n.t)
  
  accepted <- matrix(0,nrow=n.s, ncol=n.t) 

  if(is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
  
  Res<-foreach(t = 1:n.t, .combine = "cbind")%dopar%{
    res<- update_X_s_onetime(Y = Y[,t], X = X[,t], X_s = X.s[,t], cen = cen[,t], prob_below = prob.below, 
                      theta_gpd = theta.gpd, delta = delta, tau_sqd = tau_sqd, thresh_X = thresh.X, R = R[t], V = V, d = d)
    c(res$X.s, res$accept, t)
  }
  
  X.s <- Res[1:n.s, ]
  accepted <-Res[(n.s+1):(2*n.s),]
  O<-order(Res[2*n.s+1,])
  X.s <- X.s[,O]
  accepted <- accepted[,O]
  
  return(list(X.s=X.s, accepted=accepted))
}


#                                                                              #
################################################################################








################################################################################
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. delta (from Huser and Wadsworth 2017)
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## tau_sqd
##

delta.update.mixture.me.likelihood <- function(data, params, Y, X.s, cen, 
                                                   prob.below, theta.gpd,
                                                   tau_sqd) {
  R <- data
  delta <- params
  if(delta < 0 || delta > 1) return(-Inf)
  
  X <- NA * Y
  X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd=tau_sqd, delta=delta, 
                               theta.gpd=theta.gpd, prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
               theta.gpd, delta, tau_sqd) + dhuser.wadsworth(R, delta, log=TRUE)

  return(ll)
}
# delta.update.mixture.me.likelihood(R, delta, Y, X.s, cen, prob.below, theta.gpd, tau)

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the measurement error variance (on the X scale), for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
## Just a wrapper for marg.transform.data.mixture.me.likelihood
##
##   *********** If we do end up updating the prob.below parameter, this
##   *********** is a good place to do it.
##
## data............................... a n.t vector of scaling factors
## params............................. tau_sqd
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## delta
##
tau.update.mixture.me.likelihood <- function(data, params, Y, X.s, cen, 
                                                 prob.below, delta, 
                                                 theta.gpd) {
  
  R <- data
  tau_sqd <- params

  X <- NA * Y
  X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd=tau_sqd, delta=delta, 
                               theta.gpd=theta.gpd, prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
                                                  theta.gpd, delta, tau_sqd)
  
  return(ll)
}
#tau.update.mixture.me.likelihood(R, tau, Y, X.s, cen, prob.below, delta, theta.gpd)

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the parameters of the GPD response distribution, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
## Just a wrapper for marg.transform.data.mixture.me.likelihood
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters:
##                                     theta[1] = GPD scale
##                                     theta[2] = GPD shape
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## delta
## tau_sqd
##
theta.gpd.update.mixture.me.likelihood <- function(data, params, Y, X.s, cen, 
                                                   prob.below, delta,
                                                   tau_sqd, loc, thresh.X=NULL) {
  
  R <- data
  theta.gpd <- c(loc, params)
  
  scale <- params[1]
  shape <- params[2]
  if (shape >= 0) max.support <- Inf  else max.support <- loc - scale/shape
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if (max(Y, na.rm=TRUE) > max.support) return(-Inf)
  
  X <- NA * Y
  X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd=tau_sqd, delta=delta, 
                               theta.gpd=theta.gpd, prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
                                  theta.gpd, delta, tau_sqd, thresh.X = thresh.X)
  
  return(ll)
}


# scale.gpd.update.mixture.me.likelihood <- function(data, params, shape, Y, X.s, cen, 
#                                                    prob.below, delta,
#                                                    tau_sqd, loc, thresh.X=NULL) {
#   
#   R <- data
#   theta.gpd <- c(loc, params, shape)
#   
#   scale <- params[1]
#   if (shape >= 0) max.support <- Inf  else max.support <- loc - scale/shape
#   
#   # If the parameters imply support that is not consistent with the data,
#   # then reject the parameters.
#   if (max(Y, na.rm=TRUE) > max.support) return(-Inf)
#   
#   X <- NA * Y
#   X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd=tau_sqd, delta=delta, 
#                                theta.gpd=theta.gpd, prob.below = prob.below)
#   
#   ll <- matrix(0, nrow(Y), ncol(Y))
#   ll[!cen] <- dnorm(X[!cen], mean=X.s[!cen], sd=sqrt(tau_sqd), log=TRUE) +
#     dgpd(Y[!cen], loc=loc, scale=scale, shape=shape, log=TRUE)  -
#     dmixture.me(X[!cen], tau_sqd = tau_sqd, delta = delta, log=TRUE)+log(1-prob.below)
#   
#   return(sum(ll))
# }
# 
# 
# shape.gpd.update.mixture.me.likelihood <- function(data, params, scale, Y, X.s, cen, 
#                                                    prob.below, delta,
#                                                    tau_sqd, loc, thresh.X=NULL) {
#   
#   R <- data
#   theta.gpd <- c(loc, scale, params)
#   
#   shape <- params[1]
#   if (shape >= 0) max.support <- Inf  else max.support <- loc - scale/shape
#   
#   # If the parameters imply support that is not consistent with the data,
#   # then reject the parameters.
#   if (max(Y, na.rm=TRUE) > max.support) return(-Inf)
#   
#   X <- NA * Y
#   X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd=tau_sqd, delta=delta, 
#                                theta.gpd=theta.gpd, prob.below = prob.below)
#   
#   ll <- matrix(0, nrow(Y), ncol(Y))
#   ll[!cen] <- dnorm(X[!cen], mean=X.s[!cen], sd=sqrt(tau_sqd), log=TRUE) +
#     dgpd(Y[!cen], loc=loc, scale=scale, shape=shape, log=TRUE)  -
#     dmixture.me(X[!cen], tau_sqd = tau_sqd, delta = delta, log=TRUE)+log(1-prob.below)
#   
#   return(sum(ll))
# }
# theta.gpd.update.mixture.me.likelihood(R, c(1,0), Y, X.s, cen, prob.below, delta,
#                                               tau, loc=thresh, thresh.X=thresh.X)

#                                                                              #
################################################################################



################################################################################
## Update covariance parameters. For the generic Metropolis sampler
## Samples from the parameters of the underlying Gaussian process, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters:
##                                     lambda = params[1]
##                                     gamma  = params[2]
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## R
## S
## V, d
##
theta.c.update.mixture.me.likelihood <- function(data, params, X.s, S, 
                                                 V=NULL, d=NULL) {
  
  library(doParallel)
  library(foreach)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  n.t <- ncol(X.s)
  registerDoParallel(cores=n.t)
  
  R <- data
  range <- params[1]
  nu <- params[2]
  # if(lambda<0 || gamma<0 || gamma>2)  return(-Inf)
  
  if(is.null(V)){
    Cor   <- corr.fn(rdist(S), theta=c(range,nu))
    eig.Sigma <- eigen(Cor, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values
  }
  
  ll<-foreach(i = 1:n.t, .combine = "c") %dopar% {
    X.s.likelihood.conditional(X.s[,i], R[i], V, d)
    # dmvn.eig(qnorm(1-R[i]/X.s[, i]), V = V, d.inv = 1/d)
  }
  
  return(sum(ll))
}


#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the scaling factors, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. R[t]
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error (vector)
## V, d
## delta
##
Rt.update.mixture.me.likelihood <- function(data, params, X.s, delta, 
                                                 V=NULL, d=NULL) {
  R <- data
  Rt <- params
  if(Rt < 1)  return(-Inf)
  
  ll <- X.s.likelihood.conditional(X.s, Rt, V, d)
  return(as.vector(ll))
}

# Rt.update.mixture.me.likelihood(R, R[1], X.s[,1], delta, V, d)

#                                                                              #
################################################################################


