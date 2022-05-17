library(data.table)
library(latex2exp)
dt.org <- data.frame(read.table('gauge.txt',header = T))
dt.org

#Q1: -----
#raw plot
tiff("Fig 1: raw regression.tiff", width = 9, height = 5, units = "in", res = 300)
par(mfrow=c(1,2))
plot(dt.org,
     pch = 20,
     cex = 0.6,
     cex.axis = 0.6)
abline(reg.raw,
       col = "blue",
       lwd = 1.5)
legend(x = "topright",
       legend = c("linear-linear regression"),
       col =c("blue"),
       lty = 1,
       cex = 0.6,
       lwd = 1.5)
plot(residuals(reg.raw)~dt.org$density,
     pch = 20,
     cex = 0.6,
     cex.axis = 0.6,
     ylab = "gain: residual",
     xlab = "density")
abline(h=0,
       lty = "dashed",
       col = "blue",
       lwd = 1.5)
legend(x = "topright",
       legend = c("0 residual reference"),
       col =c("blue"),
       lty = 2,
       cex = 0.6,
       lwd = 1.5)
dev.off()

qqnorm(residuals(reg.raw),pch = 20,cex = 0.7)
abline(0,1)

#linear-linear regression
reg.raw <- lm(data = dt.org, gain~density)
reg.raw
plot(dt.org,
     pch = 20,
     cex = 0.6)
abline(reg.raw )

tiff("raw reg: residual.tiff", width = 6, height = 5, units = "in", res = 300)

dev.off()

#Q2: ----
dt.trans <- data.table(dt.org)
dt.trans[,gain.log:= log(gain)]
reg.log_linear <- lm(data = dt.trans,
                     gain.log~density)

tiff("log-linear.tiff", width = 6, height = 5, units = "in", res = 300)
plot(data = dt.trans,
     gain.log ~ density,
     pch = 20,
     cex = 0.6,
     ylab = "log(gain)")
abline(reg.log_linear)
dev.off()

tiff("log-linear reg: residual.tiff", width = 6, height = 5, units = "in", res = 300)
plot(residuals(reg.log_linear)~dt.trans$density,
     pch = 20,
     cex = 0.6,
     xlab = "density",
     ylab = "residual: log(gain)")
abline(h=0,
       lty = "dashed",
       col = "dark blue")
dev.off()

residual.sd <- sd(residuals(reg.log_linear))

#Q3: ----
dt.trans[,id:=1:90]
tuning.n = 0.01  #adjust normal simulation SD
dt.trans[density!=0.001,
         den.Nsim1:= rnorm(1, mean = density, sd = tuning.n),by = id]
reg.log_linear.N1 <- lm(data = dt.trans[density!=0.001,],
                        gain.log ~ den.Nsim1)
plot(data = dt.trans,
     gain.log ~ den.Nsim1,
     pch = 20,
     cex = 0.6,
     ylab = "log(gain)")
abline(reg.log_linear.N1)
summary(reg.log_linear.N1)
summary(reg.log_linear)
dt.trans
tuning.u = 0.01  #adjust uniform simulation width
dt.trans[density!=0.001,
         den.Usim1:= runif(1,min = density - tuning.u,max = density + tuning.u),by = id]
reg.log_linear.U1 <- lm(data = dt.trans[density!=0.001,],
                        gain.log ~ den.Usim1)

plot(data = dt.trans,
     gain.log ~ den.Usim1,
     pch = 20,
     cex = 0.6,
     ylab = "log(gain)")
abline(reg.log_linear.U1)

summary(reg.log_linear.U1)
summary(reg.log_linear) #slope = -4.60594 intercept = 5.99727

#noise level <- N(density,0.005)
dt.trans$density
densityindex <- 1:90
means <- 1:length(seq(from = 0.001, to = 0.01, by= 0.001))

density_data
rnorm(n=10,mean = 0.001,sd = 0.01)
dist(density_data)
distance <- c(0.082,0.096,0.096,0.094,0.095,0.075,0.068,0.079)
mean(distance)/2
noise.set <- seq(from = 0.001, to = 0.01, by = 0.001)


holder <- rep(NA,1000)
slope.set <- data.frame(s1=holder,s2=holder,s3=holder,s4=holder,s5=holder,s6=holder,s7=holder,s8=holder,s9=holder,s10=holder)
intercept.set <- data.frame(s1=holder,s2=holder,s3=holder,s4=holder,s5=holder,s6=holder,s7=holder,s8=holder,s9=holder,s10=holder)
Rsq.set <- data.frame(s1=holder,s2=holder,s3=holder,s4=holder,s5=holder,s6=holder,s7=holder,s8=holder,s9=holder,s10=holder)
df.set <- data.frame(s1=holder,s2=holder,s3=holder,s4=holder,s5=holder,s6=holder,s7=holder,s8=holder,s9=holder,s10=holder)

for (trial in 1:10) {
  noise <- noise.set[trial]
  slopes <- rep(NA, 1000)
  intercepts <- rep(NA, 1000)
  sizes <- rep(NA,1000)
  Rsq <- rep(NA,1000)
  for (i in 1:1000) {
    dens.sim <- rep(NA, 90)
    for (idx in 1:90) {
      dens.sim[idx] <- rnorm(n = 1 ,
                             mean = dt.trans$density[idx],
                             sd = noise)
    }
    data.sim <- data.table(data.frame(g=dt.trans$gain.log,den =dens.sim))
    reg.curr <- lm(data = data.sim[dens.sim>0,],
       g~den)
    slopes[i] <- reg.curr$coefficients[2]
    intercepts[i] <- reg.curr$coefficients[1]
    Rsq[i] <- summary(reg.curr)$r.squared
    sizes[i] <- data.sim[dens.sim>0,.N]
  }
  slope.set[,trial] <- slopes
  intercept.set[,trial] <- intercepts
  df.set[,trial] <- sizes
  Rsq.set[,trial] <- Rsq
}

intercept.set
Rsq.set

for(i in 1:10){
  print(mean(Rsq.set[,i]))
}

summary(reg.log_linear)
b1.ctr <- -4.60594
b1.se <- 0.03182
b0.ctr <- 5.997
b0.se <- 0.01274
b1.lwr <- (-4.60594 - 0.03182)
b1.upr <- (-4.60594 + 0.03182)

cover.1sd <- rep(NA,10)
sign95 <- rep(NA,10)
sign99 <- rep(NA,10)
sign90 <- rep(NA,10)
sign95.b0 <- rep(NA,10)
sign99.b0 <- rep(NA,10)
sign90.b0 <- rep(NA,10)

for(col in 1:10){
  t.stats <- rep(NA,1000)
  t.stats.b0 <- rep(NA,1000)
  t.crit95 <- rep(NA,1000)
  t.crit99 <- rep(NA,1000)
  t.crit90 <- rep(NA,1000)
  
  for(i in 1:1000){
    t.stats[i] <- (slope.set[i,col]-b1.ctr)/b1.se
    t.stats.b0[i] <- (intercept.set[i,col]-b0.ctr)/b0.se
    t.crit95[i] <- qt(p = 0.95, df = df.set[i,col]-2)
    t.crit99[i] <- qt(p = 0.99, df = df.set[i,col]-2)
    t.crit90[i] <- qt(p = 0.90, df = df.set[i,col]-2)
  }
  
  sign95[col] <- mean(abs(t.stats[1:1000]) <= t.crit95[1:1000])
  sign99[col] <- mean(abs(t.stats[1:1000]) <= t.crit99[1:1000])
  sign90[col] <- mean(abs(t.stats[1:1000]) <= t.crit90[1:1000])
  sign95.b0[col] <- mean(abs(t.stats.b0[1:1000]) <= t.crit95[1:1000])
  sign99.b0[col] <- mean(abs(t.stats.b0[1:1000]) <= t.crit99[1:1000])
  sign90.b0[col] <- mean(abs(t.stats.b0[1:1000]) <= t.crit90[1:1000])
  cover.1sd[col] <- mean(slope.set[,col]>b1.lwr & slope.set[,col]<b1.upr)
}
sign90
sign95
sign90.b0
sign95.b0
Rsq
reg.log_linear

mean(slope.set[,1]>4)

slope.set
means
reg.log_linear

Rsq.means <- rep(NA,10)

for(i in 1:10){
  Rsq.means[i] <- mean(Rsq.set[,i])
}
Rsq.means
round(Rsq.means,4)
plot(1:10,sign90*100,
     type = "b",
     pch = 20,
     ylim = c(70,100),
     cex = 0.8,
     xaxt = 'n',
     ylab = 'Performance indicators (in %)',
     xlab = expression(paste("noise level: ", italic("s"))))
lines(1:10,sign90.b0*100,type = "b",pch = 20,col = "blue",cex = 0.8)
lines(1:10,sign95*100,type = "b",pch = 4,cex = 0.8)
lines(1:10,sign95.b0*100,type = "b",pch = 4,col = "blue",cex = 0.8)
lines(1:10,Rsq.means*100,type ="b",cex = 0.8)
axis(side = 1, labels = seq(from = 0.001, to = 0.01, length.out = 10),at = 1:10)
grid()

plot(1:10,xlab = TeX('$\\e^\\Y_i$'))
#Q4 ----
den.org <- dt.org$density
den.factor <- unique(den.org)
pred <- predict(reg.log_linear, newdata = data.frame(density = den.factor), interval = 'prediction')
plot(data = dt.trans,
     gain.log ~ density,
     pch = 20,
     cex = 0.6,
     ylab = "log(gain)")
abline(reg.log_linear)
lines(pred[,2] ~ den.factor,lty = 2, col = "blue")
lines(pred[,3] ~ den.factor,lty = 2, col = "blue")

dt.pred <- data.table(pred)
dt.pred[,width:= (upr - lwr)/2]
dt.pred[width==min(width),] # most accurate at log(gain) of 4.5325
dt.pred[width==max(width),] # least accurate at log(gain) of 2.837

#specific prediction at 0.508 and 0.001
den.point <- c(0.508,0.001)
pred.points <- predict(reg.log_linear,newdata = data.frame(density = den.point), interval = 'prediction')

fullden <- seq(from = min(0.001) ,to = max(den.factor) ,length.out = 500)
fullprediction <- predict(reg.log_linear,newdata = data.frame(density = fullden), interval = 'prediction')
fullprediction <- data.table(fullprediction)
fullprediction[,width:= (upr-lwr)/2]
fullprediction/sd(residuals(reg.log_linear))

pred.points <- data.table(pred.points)
pred.points[,width:=(upr - lwr)/2]
pred.points
point0.001 <- c(3.65745, 3.521268, 3.793631)
ponit0.508 <- c(5.99266, 5.855343, 6.129976)
exp(point0.001)
exp(point0.508)
dist(exp(ponit0.508))


dt.pred

dt.trans[density==0.508|density == 0.001,]
dt.org
#helper funct
ion
pred.reverse <- function(point){
  return(-0.217*point+1.302)
}

pred.reverse(3.6574)

#Q5
reg.log_linear
beta0 = reg.log_linear$coefficients[1]
beta1 = reg.log_linear$coefficients[2]

#Reverse
beta0_Rev=as.numeric(-beta0/beta1)
beta0_Rev
beta1_Rev=as.numeric(1/beta1)
beta1_Rev

pred.points
pred.points.rev<- beta0_Rev+pred.points*beta1_Rev
pred.points.rev

points <- c(3.521268,3.65745, 3.793631)
pred.reverse <- function(point){
  return(-0.217*point+1.302)
}
pred.reverse(3.6574)
pred.reverse(points)
dist(pred.reverse(points))

beta0_Rev+log(33.827)*beta1_Rev

dt.pred
X_Predict<- beta0_Rev+dt.pred*beta1_Rev
X_Predict<- as.matrix(X_Predict)
X_Predict
pred

plot(data = dt.trans,
     density~gain.log,
     pch = 20,
     cex = 0.6,
     ylab = "density")
abline(a=beta0_Rev,b=beta1_Rev)

lines(X_Predict[,2] ~ pred[,1],lty = 2, col = "blue")
lines(X_Predict[,3] ~ pred[,1],lty = 2, col = "blue")

gain.point <- c(38.6, 426.7)
gain.point.log<- log(gain.point)
gain.point.log
X <- beta0_Rev+gain.point.log*beta1_Rev
X


#Q6
dt.sub1 <- dt.trans[density!= 0.508,c("density","gain","gain.log","id")]
reg.sub1 <- lm(data = dt.sub1,gain.log~density)
summary(reg.sub1)

pred.sub<-predict(reg.sub1,data.frame(density = den.factor),interval = 'prediction')

plot(y=dt.sub1[,gain.log],x=dt.sub1[,density])
#abline(reg.log_linear,col = "red")
abline(reg.sub1,col = "blue")
lines(pred.sub[,2]~den.factor,col = "blue",lty = 2)
#lines(pred[,2]~den.factor,col = "red",lty = 2)
lines(pred.sub[,3]~den.factor,col = "blue",lty = 2)

sub1.lower <- lm(pred.sub[,2]~den.factor)
summary(sub1.lower) #5.853688
sub1.lower$coefficients[1]
rev.sub1.lower <- function(loggain){
  return((loggain-sub1.lower$coefficients[1])/sub1.lower$coefficients[2])
}
rev.sub1.lower(log(38.6))
sub1.upper <- lm(pred.sub[,3]~den.factor)
rev.sub1.upper <- function(loggain){
  return((loggain-sub1.upper$coefficients[1])/sub1.upper$coefficients[2])
}
rev.sub1.upper(log(38.6))

rev.sub1 <- function(loggain){
  return((loggain-reg.sub1$coefficients[1])/reg.sub1$coefficients[2])
}

(log(38.6)-reg.sub1$coefficients[1])/reg.sub1$coefficients[2]

dt.sub2 <- dt.trans[density!= 0.001,c("density","gain","gain.log","id")]
dt.sub2
reg.sub2 <- lm(data = dt.sub2,gain.log~density)
reg.sub2
pred.sub2<-predict(reg.sub2,newdata = data.frame(density = den.factor),interval = 'prediction')
pred.sub2
exp(2.851833)
sub2.lower <- lm(pred.sub2[,2]~den.factor)
summary(sub2.lower) 
sub2.lower$coefficients[1]
rev.sub2.lower <- function(loggain){
  return((loggain-sub2.lower$coefficients[1])/sub2.lower$coefficients[2])
}
rev.sub2.lower(log(426.7))
sub2.upper <- lm(pred.sub2[,3]~den.factor)
rev.sub2.upper <- function(loggain){
  return((loggain-sub2.upper$coefficients[1])/sub2.upper$coefficients[2])
}
rev.sub2.upper(log(426.7))
rev.sub2 <- function(loggain){
  return((loggain-reg.sub2$coefficients[1])/reg.sub2$coefficients[2])
}


(log(426.7)-reg.sub2$coefficients[1])/reg.sub2$coefficients[2]
reg.sub2


dt.trans
plot(x=dt.sub2[,density],y=dt.sub2[,gain.log])
abline(reg.sub2)

gain426<- dt.trans[density == 0.001,gain]
gain426
gain38 <- dt.trans[density == 0.508,gain]
gain38

rev.sub1(log(gain38))
rev.sub1.upper(log(gain38))
rev.sub1.lower(log(gain38))

pred38<- data.frame(gain38 = gain38,
                    d.center38=rev.sub1(log(gain38)),
                    d.lower38=rev.sub1.lower(log(gain38)),
                    d.upper38=rev.sub1.upper(log(gain38)))

rev.sub2(log(gain426))
rev.sub2.upper(log(gain426))
rev.sub2.lower(log(gain426))

pred426 <- data.frame(gain426 = gain426,
                      d.center426 = rev.sub2(log(gain426)),
                      d.lower426 =rev.sub2.lower(log(gain426)),
                      d.upper426 =rev.sub2.upper(log(gain426)) )

pred38  <- data.table(pred38)
pred426 <-data.table(pred426)

pred38[,distToUpr:= d.upper38-0.508]
pred38[,distToLwr:= -d.lower38+0.508]

pred426[,distToUpr:= d.upper426-0.001]
pred426[,distToLwr:= -d.lower426+0.001]
pred38
pred426

dt.sub1
dt.sub2
-sd(residuals(reg.log_linear)) + 5.997

dt.trans[,.N,by=density]

#adv
dt.adv <- dt.trans[density!= 0.508 & density!= 0.001,c("density","gain","gain.log")]
dt.adv

K = 10
n = dt.adv[,.N]
n.vld <- n/K
fold <- sample(rep(1:K, length = n))
fold
fold

loggain.v <- dt.adv[,gain.log]
density.v <- dt.adv[,density]

plot(density.v,loggain.v)
reg.adv <- lm(loggain.v~density.v)

ROOS1 <- rep(NA,10)
ROOS2 <- rep(NA,10)

for (run in 1:10) {
  id.vld <- which(fold == run)
  g <- loggain.v[-id.vld]
  d <- density.v[-id.vld]
  reg.train <- lm(g~d)
  g.hat <- predict(reg.train,data.frame(d=density.v[id.vld]))
  #print(mean((loggain.v[id.vld]-g.hat)^2))
  benchmark <- mean(g)
  g.hat2 <- predict(reg.adv,data.frame(density.v=density.v[id.vld]))
  MSE_valid <- mean((loggain.v[id.vld]-g.hat)^2)
  MSE_benchmark <- mean((loggain.v[id.vld]-benchmark)^2)
  MSE_benchmark2 <- mean((loggain.v[id.vld]-g.hat2)^2)
  ROOS1[run]<-(1-MSE_valid/MSE_benchmark)
  print("---")
  ROOS2[run]<-(1-MSE_valid/MSE_benchmark2)
  #print(MSE_valid)
  #print(density.v[id.vld])
  #print(g.hat)
}
ROOS1
ROOS2
mean(ROOS1)
mean(ROOS2)


dt.pred[,4]/sd(residuals(reg.log_linear))
sd(residuals(reg.log_linear))



dt.pred

