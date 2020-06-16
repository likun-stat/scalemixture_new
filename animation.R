library(scales)
setwd("~/Desktop")
j=9
out.obj$X.s.trace[1,,j]
thresh<-qmixture.me.interp(p = 0.8, tau_sqd = 9, delta = 0.7)

png(file=sprintf("./plots/ani_%4.4i-00.png", 1), height=5, width=6, units="in", res=150)
plot(c(state$S[1 ,1],state$S[ ,1]), c(thresh,out.obj$X.trace[1, ,j]), type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
matplot(state$S[ ,1], out.obj$X.s.trace[1,,j], col=alpha("gray30", 0.5), type="l", add=TRUE)
points(state$S[ ,1], out.obj$X.trace[1, ,j], pch=c(20, 4)[state$cen[ ,j]+1], col=2*(state$cen[ ,j]+1), cex=0.5)
abline(h=thresh, lty=2, col="gray80", lwd=3)
dev.off()

cen.X<-out.obj$X.trace[1, ,j]
cen.X[state$cen[,j]]<-thresh
P<-890
count=0
for(t in 1:100){
  if (t > 1) {
    points(state$S[ ,1], cen.X, pch=c(20, 4)[state$cen[ ,j]+1], col=2*(state$cen[ ,j]+1), cex=0.5)
    dev.off()
  }
  
  #true smooth process in grey + new updated Y[!cen]
  count=count+1
  png(file=sprintf("./plots/trace%i.png", count), height=5, width=6, units="in", res=150)
  plot(c(state$S[1 ,1],state$S[ ,1]), c(thresh,out.obj$X.trace[1, ,j]), type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(state$S[ ,1], out.obj$X.s.trace[1,,j], col=alpha("purple", 0.3), type="l", add=TRUE)
  abline(h=thresh, lty=2, col="gray80", lwd=3)
  points(state$S[ ,1], cen.X, pch=c(20, 4)[state$cen[ ,j]+1], col=2*(state$cen[ ,j]+1), cex=0.5)
  dev.off() 
  
  #true smooth process in grey + new updated Y[!cen] + update X.s
  count=count+1
  png(file=sprintf("./plots/trace%i.png", count), height=5, width=6, units="in", res=150)
  plot(c(state$S[1 ,1],state$S[ ,1]), c(thresh,out.obj$X.trace[1, ,j]), type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(state$S[ ,1], out.obj$X.s.trace[1,,j], col=alpha("purple", 0.3), type="l", add=TRUE)
  abline(h=thresh, lty=2, col="gray80", lwd=3)
  lines(state$S[ ,1], out.obj$X.s.trace[t+P+1,,j], col="gray50", lwd=2)
  points(state$S[ ,1], cen.X, pch=c(20, 4)[state$cen[ ,j]+1], col=2*(state$cen[ ,j]+1), cex=0.5)
  dev.off() 
  
  #true smooth process in grey + new updated X.s
  count=count+1
  png(file=sprintf("./plots/trace%i.png", count), height=5, width=6, units="in", res=150)
  plot(c(state$S[1 ,1],state$S[ ,1]), c(thresh,out.obj$X.trace[1, ,j]), type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(state$S[ ,1], out.obj$X.s.trace[1,,j], col=alpha("purple", 0.3), type="l", add=TRUE)
  abline(h=thresh, lty=2, col="gray80", lwd=3)
  lines(state$S[ ,1], out.obj$X.s.trace[t+P+1,,j], col="gray50", lwd=2)
  points(state$S[state$cen[,j] ,1], cen.X[state$cen[,j]], pch=4, col=4, cex=0.5)
  dev.off() 
  
  #true smooth process in grey + new updated X.s + update Y[!cen]
  count=count+1
  png(file=sprintf("./plots/trace%i.png", count), height=5, width=6, units="in", res=150)
  cen.X[!state$cen[,j]]<-out.obj$X.trace[t+P+2,!state$cen[,j],j]
  plot(c(state$S[1 ,1],state$S[ ,1]), c(thresh,out.obj$X.trace[1, ,j]), type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(state$S[ ,1], out.obj$X.s.trace[1,,j], col=alpha("purple", 0.3), type="l", add=TRUE)
  abline(h=thresh, lty=2, col="gray80", lwd=3)
  lines(state$S[ ,1], out.obj$X.s.trace[t+P+1,,j], col="gray50", lwd=2)
}
points(state$S[ ,1], cen.X, pch=c(20, 4)[state$cen[ ,j]+1], col=2*(state$cen[ ,j]+1), cex=0.5)
dev.off()

library(magick)
img<-image_read(sprintf("~/Desktop/plots/ani_%4.4i-01.png", t))
for(t in 1:50){
  tmp<-image_read(sprintf("~/Desktop/plots/ani_%4.4i-01.png", t))
  img<-c(img,tmp)
  tmp<-image_read(sprintf("~/Desktop/plots/ani_%4.4i-02.png", t))
  img<-c(img,tmp)
  tmp<-image_read(sprintf("~/Desktop/plots/ani_%4.4i-03.png", t))
  img<-c(img,tmp)
  tmp<-image_read(sprintf("~/Desktop/plots/ani_%4.4i-04.png", t))
  img<-c(img,tmp)
}


ani<-image_animate(img, fps = 4, dispose = "previous",loop=1)
image_write(ani, "~/Desktop/plots/sim.gif")

i=5
plot(out.obj$X.trace[1:6000,which(!state$cen[,j])[i],j],type='l',xlab="index",ylab="Uncensored X")
abline(h=out.obj$X.trace[1,which(!state$cen[,j])[i],j],lty=2,col="red")
