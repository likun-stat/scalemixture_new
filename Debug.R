## ---------- Make sure full conditionals attain maximum around the true values ----------
tmp <- X.update(Y=Y, cen=cen, X.s=true.params$X.s, delta=true.params$delta, tau_sqd=true.params$tau, theta.gpd=true.params$theta.gpd, prob.below=prob.below)

which<-4; plot(true.params$X.s[,which]);points(tmp[,which],pch=20,col='red');abline(h=thresh.X)


Cor   <- corr.fn(rdist(true.params$S), c(range,nu))
eig.Sigma <- eigen(Cor, symmetric=TRUE)
V <- eig.Sigma$vectors
d <- eig.Sigma$values

X.s.up <- X.s.update.mixture.me.update.par.once.without.X.par(R=true.params$R, Y=Y, X=tmp, X.s=true.params$X.s, 
                                    cen=cen, prob.below=prob.below, theta.gpd=true.params$theta.gpd, delta=true.params$delta, 
                                    tau_sqd=true.params$tau, V=V, d=d, thresh.X=NULL)
which<-6; plot(tmp[,which]);points(X.s.up$X.s[,which],pch=20,col='red');abline(h=thresh.X)
X.s.up$accepted[,which]

fun<- function(delta) delta.update.mixture.me.likelihood(data=true.params$R, params=delta, Y=Y, X.s=X.s.up$X.s, cen=cen, 
                                               prob.below=prob.below, theta.gpd=true.params$theta.gpd,
                                               tau_sqd=true.params$tau)
fun<-Vectorize(fun)
Delta <- seq(0.52,0.55,length.out = 50)
Lik <- fun(Delta)
plot(Delta,Lik,type='l',ylab = 'Full Conditional for Delta parameter')
abline(v=0.7,lty=2,col='red')
abline(v=Delta[which.max(Lik)],lty=2,col='blue')
grid()


fun<- function(tau) tau.update.mixture.me.likelihood(data=true.params$R, params=tau, Y=Y, X.s=X.s.up$X.s, cen=cen, 
                                                                 prob.below=prob.below, delta=true.params$delta, 
                                                                 theta.gpd=true.params$theta.gpd)
fun<-Vectorize(fun)
Tau <- seq(3,5,length.out = 50)
Lik <- fun(Tau)
plot(Tau,Lik,type='l',ylab = 'Full Conditional for Tau parameter')
abline(v=4,lty=2,col='red')
abline(v=Tau[which.max(Lik)],lty=2,col='blue')
grid()


fun<- function(shape) theta.gpd.update.mixture.me.likelihood(data = true.params$R, 
                  params = c(true.params$theta.gpd[2],shape),Y = Y, X.s = X.s.up$X.s, 
                  cen = cen, prob.below = true.params$prob.below,delta = true.params$delta,
                  tau_sqd = true.params$tau, loc = true.params$theta.gpd[1], thresh.X = thresh.X)
fun<-Vectorize(fun)
Shape <- seq(-0.005,0.005,length.out = 50)
Lik <- fun(Shape)
plot(Shape,Lik,type='l',ylab = 'Full Conditional for Shape parameter')
abline(v=0,lty=2,col='red')
grid()


fun<- function(scale) theta.gpd.update.mixture.me.likelihood(data = true.params$R, 
                          params = c(scale,true.params$theta.gpd[3]),Y = Y, X.s = true.params$X.s, 
                          cen = cen, prob.below = true.params$prob.below,delta = true.params$delta, 
                          tau_sqd = true.params$tau, loc = true.params$theta.gpd[1], thresh.X = thresh.X)
fun<-Vectorize(fun)
Scale <- seq(0.995,1.005,length.out = 50)
Lik <- fun(Scale)
plot(Scale,Lik,type='l',ylab = 'Full Conditional for Scale parameter')
abline(v=1,lty=2,col='red')
grid()


fun<- function(range) theta.c.update.mixture.me.likelihood(data=R, params=c(range, true.params$theta.c[2]), X.s=true.params$X.s, S=true.params$S, 
                                                                       V=NULL, d=NULL)
fun<-Vectorize(fun)
Range <- seq(0.9,1.1,length.out = 50)
Lik <- fun(Range)
plot(Range,Lik,type='l',ylab = 'Full Conditional for Range parameter')
abline(v=Range[which.max(Lik)],lty=2,col='red')
grid()


fun<- function(nu) theta.c.update.mixture.me.likelihood(data=R, params=c(true.params$theta.c[1], nu), X.s=true.params$X.s, S=true.params$S, 
                                                           V=NULL, d=NULL)
fun<-Vectorize(fun)
Nu <- seq(1.3,1.6,length.out = 50)
Lik <- fun(Nu)
plot(Nu,Lik,type='l',ylab = 'Full Conditional for Nu parameter')
abline(v=Nu[which.max(Lik)],lty=2,col='red')
grid()


## ---------- Check the trace plots & the current full conditionals ----------
plot(out.obj$delta.trace[1:state$i],type='l',ylab=expression(delta))
state$sigma.m$delta  # current proposal variance
fun<- function(delta) delta.update.mixture.me.likelihood(data=state$R, params=delta, Y=state$Y, X.s=state$X.s, cen=state$cen, 
                                                         prob.below=state$prob.below, theta.gpd=state$theta.gpd,
                                                         tau_sqd=state$tau)
fun<-Vectorize(fun)
Delta <- seq(0.69,0.71,length.out = 50)
Lik <- fun(Delta)
plot(Delta,Lik,type='l',ylab = 'Full Conditional for Delta parameter')
abline(v=0.65,lty=2,col='red')
abline(v=Delta[which.max(Lik)],lty=2,col='blue')
grid()

plot(out.obj$tau.trace[1:state$i],type='l',ylab=expression(tau^2))
state$sigma.m$tau # current proposal variance
fun<- function(tau) tau.update.mixture.me.likelihood(data=state$R, params=tau, Y=state$Y, X.s=state$X.s, cen=state$cen, 
                                                     prob.below=state$prob.below, delta=state$delta, 
                                                     theta.gpd=state$theta.gpd)
fun<-Vectorize(fun)
Tau <- seq(3,6,length.out = 50)
Lik <- fun(Tau)
plot(Tau,Lik,type='l',ylab = 'Full Conditional for Tau parameter')
abline(v=4,lty=2,col='red')
abline(v=Tau[which.max(Lik)],lty=2,col='blue')
grid()



plot(out.obj$theta.c.trace[1:state$i,1],type='l',ylab='range')
plot(out.obj$theta.c.trace[1:state$i,2],type='l',ylab=expression(nu))
state$sigma.m$theta.c
state$prop.Sigma$theta.c

fun<- function(range) theta.c.update.mixture.me.likelihood(data=R, params=c(range, true.params$theta.c[2]), X.s=true.params$X.s, S=true.params$S, 
                                                           V=NULL, d=NULL)
fun<-Vectorize(fun)
Range <- seq(0.9,1.1,length.out = 50)
Lik <- fun(Range)
plot(Range,Lik,type='l',ylab = 'Full Conditional for Range parameter')
abline(v=Range[which.max(Lik)],lty=2,col='red')
grid()


fun<- function(nu) theta.c.update.mixture.me.likelihood(data=R, params=c(true.params$theta.c[1], nu), X.s=true.params$X.s, S=true.params$S, 
                                                        V=NULL, d=NULL)
fun<-Vectorize(fun)
Nu <- seq(1.3,1.6,length.out = 50)
Lik <- fun(Nu)
plot(Nu,Lik,type='l',ylab = 'Full Conditional for Nu parameter')
abline(v=Nu[which.max(Lik)],lty=2,col='red')
grid()




delta=0.62
X<-state$X
X[!cen] <- gpd.2.scalemix.me(state$Y[!cen], tau_sqd=state$tau, delta=delta, 
                             theta.gpd=state$theta.gpd, prob.below =state$prob.below)

ll <- matrix(0, nrow(state$Y), ncol(state$Y))

loc <- state$theta.gpd[1]
scale <- state$theta.gpd[2]
shape <- state$theta.gpd[3]
thresh.X <- qmixture.me.interp(state$prob.below, tau_sqd = state$tau, delta = delta)

if (sum(cen) > 0)  {
  ll[state$cen] <- pnorm(thresh.X, mean=state$X.s[state$cen], sd=sqrt(state$tau), log.p=TRUE)
}
cen <-state$cen
if (sum(!cen) > 0) {
  ll[!cen] <- dnorm(X[!cen], mean=state$X.s[!cen], sd=sqrt(state$tau), log=TRUE) +
    dgpd(state$Y[!cen], loc=loc, scale=scale, shape=shape, log=TRUE)  -
    dmixture.me(X[!cen], tau_sqd = state$tau, delta = delta, log=TRUE)
}

