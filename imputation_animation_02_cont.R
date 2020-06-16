source("~/Desktop/Research/scalemixture/scalemix_utils.R")
source("~/Desktop/Research/scalemixture/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture/scalemix_priors.R")
source("~/Desktop/Research/scalemixture/generic_samplers.R")
source("~/Desktop/Research/scalemixture/scalemix_sampler_02.R")

library(fields)   # For rdist

# ---------------- 1. Load the current states --------------------
tmp1<-file.mtime("Huser-wadsworth-bigdelta_progress_F.RData")
tmp2<-file.mtime("Huser-wadsworth-bigdelta_progress_T.RData")
Tdiff<-as.numeric(difftime(tmp1,tmp2), units="secs")
if(Tdiff>0) {load("Huser-wadsworth-bigdelta_progress_F.RData")} else {load("Huser-wadsworth-bigdelta_progress_T.RData")}

Y <- state$Y
cen <- state$cen
S <- state$S
thresh <- state$thresh
initial.values <- list(delta=state$delta, rho=state$rho, tau=state$tau, theta.gpd=state$theta.gpd, 
                       prob.below=state$prob.below, X.s=state$X.s, R=state$R)
i_prev <- state$i
n.updates<-10000

sigma.m <- state$sigma.m
prop.Sigma <- state$prop.Sigma
sd.ratio <- state$sd.ratio.trace[i_prev]

true.params <- list(delta = out.obj$delta.trace[1], rho=out.obj$rho.trace[1], tau=out.obj$tau.trace[1],
                    theta.gpd=c(thresh, out.obj$theta.gpd.trace[1,]), prob.below=state$prob.below, X.s=out.obj$X.s.trace[1,,],
                    R=out.obj$R.trace[1,])

rm(state)


## --------------- 2. Running Metropolis -------------------
scalemix.sampler.02.cont(Y=Y, S=S, cen=cen, thresh=thresh,
                          initial.values=initial.values , i_prev=i_prev,
                          n.updates=n.updates, thin=10,
                          experiment.name="Huser-wadsworth-bigdelta",
                          echo.interval=50,
                          sigma.m=sigma.m, prop.Sigma=prop.Sigma, 
                          true.params=true.params, sd.ratio=sd.ratio, lower.prob.lim=0.5, out.obj=out.obj)


