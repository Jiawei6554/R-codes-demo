library(data.table)

setwd("/Users/jiawei/Documents/R files/MATH 189/HW3")
dt.org <- data.frame(read.table('hcmv.txt',header = T))

N <-length(dt.org$location)
pos.max <- max(dt.org$location)
gene <- seq(1,pos.max)

id <- 1:N
DT <- data.table(id,dt.org)

#medians ----
sim <- data.table(data.frame(id))

for (run in 1:500) {
  lable <- paste("sim", run, sep = "")
  set.seed(run^2+1)
  sim.curr <- sort(sample.int(pos.max,N),
              decreasing = FALSE)
  #assign(lable,sim.curr)
  sim<-cbind(sim,sim.curr)
}

sim<-sim[,-1]
sim.df<-data.frame(t(sim))
medians.sim <- rep(-1,296)
firstQ.sim <- rep(-1,296)
thirdQ.sim <- rep(-1,296)
for (idx in 1:296) {
  firstQ.sim[idx] <- quantile(sim.df[,idx],0.25,names = FALSE)
  medians.sim[idx] <- median(sim.df[,idx])
  thirdQ.sim[idx] <- quantile(sim.df[,idx],0.75,names = FALSE)
}
firstQ.sim
medians.sim
thirdQ.sim

medians.dt <- data.frame(id,medians.sim)

tiff("median-middle50.tiff", units="in", width=10, height=7, res=300)
plot(medians.sim,
     pch = 20,
     cex = 0.2,
     col = "blue",
     ylab = "location",,
     cex.axis = 0.6)
#abline(v=c(115,126),col = "#80808080",lty ="dashed")
abline(h=93000,col = "#80808080",lty = "dashed")
#abline(v=c(91,96),col = "#80808080",lty ="dashed")
abline(h=76000,col = "#80808080",lty = "dashed")
abline(h=195000,col = "#80808080",lty = "dashed")
for (idx2 in 1:296) {
  IQR = thirdQ.sim[idx2] - firstQ.sim[idx2]
  segments(x0 = idx2,y0=firstQ.sim[idx2],
           x1 = idx2,y1=thirdQ.sim[idx2],
           col = "#4682B450")
}
points(DT,cex = 0.2)
legend(x = "bottomright",
       legend = c("CMV data","median group"),
       col = c("black","blue"),
       pch = c(20,20),
       lty = c(0,1),
       cex = 0.8)
dev.off()

tiff("median-iqr.tiff", units="in", width=10, height=7, res=300)
plot(medians.sim,
     pch = 20,
     cex = 0.2,
     col = "blue",
     ylab = "location",
     cex.axis = 0.6)

for (idx2 in 1:296) {
  IQR = thirdQ.sim[idx2] - firstQ.sim[idx2]
  segments(x0 = idx2,y0=firstQ.sim[idx2]-1.5*IQR,
           x1 = idx2,y1=thirdQ.sim[idx2]+1.5*IQR,
           col = "#3cb37150")
}
points(DT,cex = 0.2)
legend(x = "bottomright",
       legend = c("CMV data","median group"),
       col = c("black","blue"),
       pch = c(20,20),
       lty = c(0,1),
       cex = 0.8)
dev.off()


#
#Q1----
set.seed(1999)
sim.group1 <- sort(sample.int(pos.max,N))
set.seed(2021)
sim.group2 <- sort(sample.int(pos.max,N))
set.seed(1234)
sim.group3 <- sort(sample.int(pos.max,N))
set.seed(9999)
sim.group4 <- sort(sample.int(pos.max,N))


tiff("ecdf.tiff", units="in", width=10, height=4, res=300)
plot(ecdf(DT[,location]), do.points = FALSE,
     xlab = "location",
     cex = 0.4,
     main ="")
plot(ecdf(sim.group1), do.points = FALSE, add = TRUE, col = "blue",
     cex = 0.2)
plot(ecdf(sim.group2), do.points = FALSE, add = TRUE, col = "red",
     cex = 0.2)
plot(ecdf(sim.group3), do.points = FALSE, add = TRUE, col = "pink",
     cex = 0.2)
plot(ecdf(sim.group4), do.points = FALSE, add = TRUE, col = "green",
     cex = 0.2)
plot(ecdf(DT[,location]), do.points = FALSE,add = TRUE,
     cex = 0.4)
abline(v = 80000, lty = "dashed")
abline(v = 93000, lty = "dashed")
legend(x = "bottomright",
       legend = c("group 1","group 2","group 3","group 4","CMV data"),
       col = c("blue","red","pink","green","black"),
       lty = c(1,1,1,1,1),
       cex = 0.7)
dev.off()

#region split----
regionsplit <- function(n.region, gene, site){
  count.int <- table(cut(site, breaks = seq(1, length(gene), length.out=n.region+1), include.lowest=TRUE))
  count.vector <- as.vector(count.int)
  count.tab <- table(factor(count.vector,levels=0:max(count.vector)))
  return (count.tab)
}

#Q3 ----
split100 <- regionsplit(100,gene,DT[,location])
split100
split100.counts <- c(rep(0,9),rep(1,11),rep(2,26),rep(3,21),rep(4,13),rep(5,13),rep(6,3),rep(7,3),
                     rep(8,0),rep(9,0),rep(10,0),rep(11,0),rep(12,0),rep(13,0),rep(14,1))
group1.split100 <- regionsplit(100,gene,sim.group1)
group1.split100
group1.split100.count <- c(rep(0,6),rep(1,13),rep(2,23),rep(3,23),rep(4,15),rep(5,14),rep(6,4),rep(7,2))

group2.split100 <- regionsplit(100,gene,sim.group2)
group2.split100
group2.split100.count <- c(rep(0,3),rep(1,17),rep(2,18),rep(3,21),rep(4,28),rep(5,10),rep(6,3))

group3.split100 <- regionsplit(100,gene,sim.group3)
group3.split100
group3.split100.count <- c(rep(0,7),rep(1,16),rep(2,20),rep(3,22),rep(4,14),rep(5,14),rep(6,3),rep(7,2),rep(8,2))

group4.split100 <- regionsplit(100,gene,sim.group4)
group4.split100
group4.split100.count <-c(rep(0,8),rep(1,10),rep(2,21),rep(3,28),rep(4,22),rep(5,4),rep(6,3),
                    rep(7,0),rep(8,3),rep(9,0),rep(10,1))

plot(density(split100.counts),
     ylim = c(0,0.3),
     main = "",
     xlab = "N = 100, Palindromes in a unit interval",
     xaxt = "n")
axis(side = 1, at = seq(0,15,by=1),cex.axis = 0.8)
lines(density(group1.split100.count),lty = 2,col = "blue")
lines(density(group2.split100.count),lty = 2,col = "red")
lines(density(group3.split100.count),lty = 2,col = "pink")
lines(density(group4.split100.count),lty = 2,col = "dark green")
legend(x = "topright",
       legend = c("CMV data","group1","group2","group3","group4"),
       lty = c(1,2,2,2,2),
       col = c("black","blue","red","pink","dark green"),
       cex = 0.8)


split60 <- regionsplit(60,gene,DT[,location])
split60
split60.counts <- c(rep(0,0),rep(1,4),rep(2,5),rep(3,8),rep(4,14),rep(5,9),rep(6,3),rep(7,9),
                     rep(8,5),rep(9,0),rep(10,1),rep(11,1),rep(12,0),rep(13,0),rep(14,0),rep(15,1))
group1.split60 <- regionsplit(60,gene,sim.group1)
group1.split60
group1.split60.count <- c(rep(0,2),rep(1,1),rep(2,3),rep(3,6),rep(4,13),rep(5,13),rep(6,10),rep(7,6),
                           rep(8,3),rep(9,2),rep(10,1))

group2.split60 <- regionsplit(60,gene,sim.group2)
group2.split60
group2.split60.count <- c(rep(0,0),rep(1,3),rep(2,5),rep(3,5),rep(4,10),rep(5,14),rep(6,10),rep(7,8),
                          rep(8,4),rep(9,0),rep(10,1))

group3.split60 <- regionsplit(60,gene,sim.group3)
group3.split60
group3.split60.count <- c(rep(0,0),rep(1,5),rep(2,3),rep(3,10),rep(4,11),rep(5,10),rep(6,8),rep(7,2),
                          rep(8,4),rep(9,5),rep(10,0),rep(11,2))

group4.split60 <- regionsplit(60,gene,sim.group4)
group4.split60
group4.split60.count <- c(rep(0,1),rep(1,0),rep(2,5),rep(3,11),rep(4,10),rep(5,11),rep(6,10),rep(7,4),
                          rep(8,5),rep(9,2),rep(10,0),rep(11,0),rep(12,1))

tiff("counts.tiff", units="in", width=12, height=5, res=300)
par(mfrow=c(1,2))
plot(density(split100.counts),
     ylim = c(0,0.3),
     main = "",
     xlab = "N = 100, Palindromes in a unit interval",
     xaxt = "n",
     cex.axis = 0.8,
     lwd = 1.5)
axis(side = 1, at = seq(0,15,by=1),cex.axis = 0.5)
abline(v = (0:15),col = "#80808080",lty ="dotted")
abline(h = seq(0,0.3, by =0.05),col = "#80808080",lty ="dotted")
lines(density(group1.split100.count),lty = 2,col = "blue",
      lwd = 1.5)
lines(density(group2.split100.count),lty = 2,col = "red",
      lwd = 1.5)
lines(density(group3.split100.count),lty = 2,col = "pink",
      lwd = 1.5)
lines(density(group4.split100.count),lty = 2,col = "dark green",
      lwd = 1.5)
lines(dpois(x=0:15,lambda=2.96),lwd = 2, col = "gray")
legend(x = "topright",
       legend = c("CMV data","group1","group2","group3","group4","Pois(2.96)"),
       lty = c(1,2,2,2,2),
       col = c("black","blue","red","pink","dark green","grey"),
       cex = 0.7)

plot(density(split60.counts),
     ylim = c(0,0.3),
     main = "",
     xlab = "N = 60, Palindromes in a unit interval",
     xaxt = "n",
     cex.axis = 0.8,
     lwd = 1.5)
axis(side = 1, at = seq(0,16,by=1),cex.axis = 0.5)
abline(v = (0:16),col = "#80808080",lty ="dotted")
abline(h = seq(0,0.3, by =0.05),col = "#80808080",lty ="dotted")
lines(density(group1.split60.count),lty = 2,col = "blue",
      lwd = 1.5)
lines(density(group2.split60.count),lty = 2,col = "red",
      lwd = 1.5)
lines(density(group3.split60.count),lty = 2,col = "pink",
      lwd = 1.5)
lines(density(group4.split60.count),lty = 2,col = "dark green",
      lwd = 1.5)
legend(x = "topright",
       legend = c("CMV data","group1","group2","group3","group4","Pois(4.93)"),
       lty = c(1,2,2,2,2,1),
       col = c("black","blue","red","pink","dark green","grey"),
       cex = 0.7)
lines(dpois(x=0:15,lambda=300/60),lwd = 2, col = "gray")
dev.off()




split296 <- regionsplit(296,gene,DT[,location])
split296
split296.counts <- c(rep(0,119),rep(1,108),rep(2,36),rep(3,22),rep(4,8),rep(5,1),rep(6,1),rep(7,1))
hist(split296.counts,breaks = 7,probability = 1,ylim = c(0,1))
lines(dpois(x=0:7,lambda=296/296),type="b")
lines(density(split296.counts))

group1.split296 <- regionsplit(296,gene,sim.group1)
group1.split296 
group1.split296.count <-c(rep(0,113),rep(1,106),rep(2,47),rep(3,25),rep(4,5),rep(5,1))
hist(group1.split296.count,breaks = 7,probability = 1,add = TRUE)



split50 <- regionsplit(50,gene,DT[,location]) 
split50
split50.counts <- c(rep(0,1),rep(1,1),rep(2,1),rep(3,5),rep(4,8),rep(5,6),rep(7,12),rep(8,4),rep(9,3),
                    rep(10,1),rep(11,1),rep(16,1))
hist(split50.counts)


hist(DT[,location],breaks = 100)

?hist()
a <- c(rep(0,44),rep(1,80),rep(2,36),rep(3,26),rep(4,10),rep(5,0),rep(6,2),rep(7,2))
hist(a)
hist(split1,probability = FALSE)
plot(density(split1))
split1.c <- regionsplit(200,gene,sim.group2)
#medians.sim

split1.c 
a1 <- dpois(seq(0,7,by =1),lambda = (N/200))
plot(a1, type = "h")

plot(density(a1))
hist(split1.c,breaks = 100)

plot(ecdf(medians.sim),cex = 0.2)
plot(ecdf(sim.group1),cex = 0.2,add = TRUE)
plot(ecdf(sim.group2),cex = 0.2,add = TRUE)
plot(ecdf(sim.group3),cex = 0.2,add = TRUE)     
plot(ecdf(sim.group4),cex = 0.2,add = TRUE)

split100.counts
split60.counts
group1.split60.count
length(split60.counts)
chisq.test(split60.counts,rep(296/60,60))



#------
library(car)
library(tidyverse)
library(lmtest)
library(sandwich)
library(zoo)
library(lmtest)
library(data.table)
setwd("~/Desktop")
data = read.table('hcmv.txt', header=T)
DT <- data.table(data)

Countofpalindrome<- as.numeric(count(DT))
Countofpalindrome
Palindrome <- DT$location
Palindrome
sumofpairs <- c(0)
sumoftriplets <- c(0)
sumofrandompairs <- c(0)
sumofrandomtriplets <- c(0)
L = max(DT$location)
L

#Q1 Random Generate
set.seed(1999)
sim.group1 <- runif(Countofpalindrome, min = 1, max = L)
sim.group1
set.seed(2021)
sim.group2 <- runif(Countofpalindrome, min = 1, max = L)
sim.group2
set.seed(1234)
sim.group3 <- runif(Countofpalindrome, min = 1, max = L)
sim.group3
set.seed(9999)
sim.group4 <- runif(Countofpalindrome, min = 1, max = L)
sim.group4

#Q1 Medians
#medians ----
dt.org <- data.frame(read.table('hcmv.txt',header = T))
N <-length(dt.org$location)
pos.max <- max(dt.org$location)
#medians ----
id <- 1:296
DT <- data.table(id,dt.org)
sim <- data.table(data.frame(id))

for (run in 1:500) {
  lable <- paste("sim", run, sep = "")
  set.seed(run^2+1)
  sim.curr <- sort(sample.int(pos.max,N),
                   decreasing = FALSE)
  #assign(lable,sim.curr)
  sim<-cbind(sim,sim.curr)
}

sim<-sim[,-1]
sim.df<-data.frame(t(sim))
medians.sim <- rep(-1,296)
for (idx in 1:296) {
  medians.sim[idx] <- median(sim.df[,idx])
}
medians.sim

#Q2 
#Spaces  CMV
SpaceConsecutivePal<-c(0)

for (i in 1:Countofpalindrome-1){
  SpaceConsecutivePal[i] <-sum(Palindrome[i+1] - Palindrome[i])
}
scatterplot(seq(from=1,to=295,by=1),SpaceConsecutivePal, main="Spacings between Consecutive Palindromes (CMV)")
regressionofSpace <- lm(SpaceConsecutivePal~(seq(from=1,to=295,by=1)))
regressionofSpace
abline(regressionofSpace, lwd=2)
cov((seq(from=1,to=295,by=1)),SpaceConsecutivePal)
summary(SpaceConsecutivePal)
sd(SpaceConsecutivePal)

#Spacings CMV Sum of pairs & triplets
for (i in 1:Countofpalindrome-1){
  sumofpairs[i] <-sum((Palindrome[i] + Palindrome[i+1]))
}

for (i in 1:Countofpalindrome-2){
  sumoftriplets[i] <- sum((Palindrome[i] + Palindrome[i+1] + Palindrome[i+2]))
}

summary(sumofpairs)
summary(sumoftriplets)

SpacePairsCMV <- c(0)
SpaceTripletsCMV <- c(0)

for (i in 1:Countofpalindrome-2){
  SpacePairsCMV[i] <-sum(sumofpairs[i+1] - sumofpairs[i])
}
for (i in 1:Countofpalindrome-3){
  SpaceTripletsCMV[i] <-sum(sumoftriplets[i+1] - sumoftriplets[i])
}

summary(SpacePairsCMV)
sd(SpacePairsCMV)
summary(SpaceTripletsCMV)
sd(SpaceTripletsCMV)

#Spacings Random Medium Group
SpaceConsecutiveRand <- c(0)
medians.sim <- sort(medians.sim, decreasing=FALSE)
for (i in 1:Countofpalindrome-1){
  SpaceConsecutiveRand[i] <-sum(medians.sim[i+1] - medians.sim[i])
}
SpaceConsecutiveRand
scatterplot(seq(from=1,to=295,by=1),SpaceConsecutiveRand, main="Spacings between Consecutive Palindromes (Group1Random)")
regressionofSpaceRandom <- lm(SpaceConsecutiveRand~(seq(from=1,to=295,by=1)))
regressionofSpaceRandom
abline(regressionofSpaceRandom, lwd=2)
cov((seq(from=1,to=295,by=1)),SpaceConsecutiveRand)
summary(SpaceConsecutiveRand)
sd(SpaceConsecutiveRand)


#Random pairs/triplets (Medium)
for (i in 1:Countofpalindrome-1){
  sumofrandompairs[i] <-sum(medians.sim[i] + medians.sim[i+1])
}

for (i in 1:Countofpalindrome-2){
  sumofrandomtriplets[i] <- sum((medians.sim[i] + medians.sim[i+1] + medians.sim[i+2]))
}

sumofrandompairs <- sort(sumofrandompairs, decreasing = F)
sumofrandomtriplets <- sort(sumofrandomtriplets, decreasing = F)
summary(sumofrandompairs)
summary(sumofrandomtriplets)

SpacePairsRan <- c(0)
SpaceTripletsRan <- c(0)

for (i in 1:Countofpalindrome-2){
  SpacePairsRan[i] <-sum(sumofrandompairs[i+1] - sumofrandompairs[i])
}
for (i in 1:Countofpalindrome-3){
  SpaceTripletsRan[i] <-sum(sumofrandomtriplets[i+1] - sumofrandomtriplets[i])
}

summary(SpacePairsRan)
sd(SpacePairsRan)
summary(SpaceTripletsRan)
sd(SpaceTripletsRan)

#Group 1 Spacings
SpaceConsecutiveRand1 <- c(0)
sim.group1 <- sort(sim.group1, decreasing=FALSE)
for (i in 1:Countofpalindrome-1){
  SpaceConsecutiveRand1[i] <-sum(sim.group1[i+1] - sim.group1[i])
}
SpaceConsecutiveRand1

#sum of pairs triplets&Spacings
sumofrandompairs1 <- c(0)
sumofrandomtriplets1 <- c(0)
for (i in 1:Countofpalindrome-1){
  sumofrandompairs1[i] <-sum(sim.group1[i] + sim.group1[i+1])
}
for (i in 1:Countofpalindrome-2){
  sumofrandomtriplets1[i] <- sum((sim.group1[i] + sim.group1[i+1] + sim.group1[i+2]))
}
sumofrandompairs1 <- sort(sumofrandompairs1, decreasing = F)
sumofrandomtriplets1 <- sort(sumofrandomtriplets1, decreasing = F)
summary(sumofrandompairs1)
summary(sumofrandomtriplets1)

SpacePairsRan1 <- c(0)
SpaceTripletsRan1 <- c(0)

for (i in 1:Countofpalindrome-2){
  SpacePairsRan1[i] <-sum(sumofrandompairs1[i+1] - sumofrandompairs1[i])
}

for (i in 1:Countofpalindrome-3){
  SpaceTripletsRan1[i] <-sum(sumofrandomtriplets1[i+1] - sumofrandomtriplets1[i])
}

summary(SpacePairsRan1)
sd(SpacePairsRan1)
summary(SpaceTripletsRan1)
sd(SpaceTripletsRan1)

#Group 2 Spacings
SpaceConsecutiveRand2 <- c(0)
sim.group2 <- sort(sim.group2, decreasing=FALSE)
for (i in 1:Countofpalindrome-1){
  SpaceConsecutiveRand2[i] <-sum(sim.group2[i+1] - sim.group2[i])
}
SpaceConsecutiveRand2

#Group 2 sum of pairs triplets
sumofrandompairs2 <- c(0)
sumofrandomtriplets2 <- c(0)
for (i in 1:Countofpalindrome-1){
  sumofrandompairs2[i] <-sum(sim.group2[i] + sim.group2[i+1])
}
for (i in 1:Countofpalindrome-2){
  sumofrandomtriplets2[i] <- sum((sim.group2[i] + sim.group2[i+1] + sim.group2[i+2]))
}
sumofrandompairs2 <- sort(sumofrandompairs2, decreasing = F)
sumofrandomtriplets2 <- sort(sumofrandomtriplets2, decreasing = F)
summary(sumofrandompairs2)
summary(sumofrandomtriplets2)

SpacePairsRan2 <- c(0)
SpaceTripletsRan2 <- c(0)

for (i in 1:Countofpalindrome-2){
  SpacePairsRan2[i] <-sum(sumofrandompairs2[i+1] - sumofrandompairs2[i])
}
for (i in 1:Countofpalindrome-3){
  SpaceTripletsRan2[i] <-sum(sumofrandomtriplets2[i+1] - sumofrandomtriplets2[i])
}

summary(SpacePairsRan2)
sd(SpacePairsRan2)
summary(SpaceTripletsRan2)
sd(SpaceTripletsRan2)

#Group 3 Spacings
SpaceConsecutiveRand3 <- c(0)
sim.group3 <- sort(sim.group3, decreasing=FALSE)
for (i in 1:Countofpalindrome-1){
  SpaceConsecutiveRand3[i] <-sum(sim.group3[i+1] - sim.group3[i])
}
SpaceConsecutiveRand3

#Group 3 sum of pairs triplets
sumofrandompairs3 <- c(0)
sumofrandomtriplets3 <- c(0)
for (i in 1:Countofpalindrome-1){
  sumofrandompairs3[i] <-sum(sim.group3[i] + sim.group3[i+1])
}
for (i in 1:Countofpalindrome-2){
  sumofrandomtriplets3[i] <- sum((sim.group3[i] + sim.group3[i+1] + sim.group3[i+2]))
}
sumofrandompairs3 <- sort(sumofrandompairs3, decreasing = F)
sumofrandomtriplets3 <- sort(sumofrandomtriplets3, decreasing = F)
summary(sumofrandompairs3)
summary(sumofrandomtriplets3)

SpacePairsRan3 <- c(0)
SpaceTripletsRan3 <- c(0)

for (i in 1:Countofpalindrome-2){
  SpacePairsRan3[i] <-sum(sumofrandompairs3[i+1] - sumofrandompairs3[i])
}
for (i in 1:Countofpalindrome-3){
  SpaceTripletsRan3[i] <-sum(sumofrandomtriplets3[i+1] - sumofrandomtriplets3[i])
}

summary(SpacePairsRan3)
sd(SpacePairsRan3)
summary(SpaceTripletsRan3)
sd(SpaceTripletsRan3)

#Group 4 Spacings
SpaceConsecutiveRand4 <- c(0)
sim.group4 <- sort(sim.group4, decreasing=FALSE)
for (i in 1:Countofpalindrome-1){
  SpaceConsecutiveRand4[i] <-sum(sim.group4[i+1] - sim.group4[i])
}
SpaceConsecutiveRand4

#Group 4 sum of pairs triplets
sumofrandompairs4 <- c(0)
sumofrandomtriplets4 <- c(0)
for (i in 1:Countofpalindrome-1){
  sumofrandompairs4[i] <-sum(sim.group4[i] + sim.group4[i+1])
}
for (i in 1:Countofpalindrome-2){
  sumofrandomtriplets4[i] <- sum((sim.group4[i] + sim.group4[i+1] + sim.group4[i+2]))
}
sumofrandompairs4 <- sort(sumofrandompairs4, decreasing = F)
sumofrandomtriplets4 <- sort(sumofrandomtriplets4, decreasing = F)
summary(sumofrandompairs4)
summary(sumofrandomtriplets4)

SpacePairsRan4 <- c(0)
SpaceTripletsRan4 <- c(0)

for (i in 1:Countofpalindrome-2){
  SpacePairsRan4[i] <-sum(sumofrandompairs4[i+1] - sumofrandompairs4[i])
}
for (i in 1:Countofpalindrome-3){
  SpaceTripletsRan4[i] <-sum(sumofrandomtriplets4[i+1] - sumofrandomtriplets4[i])
}

summary(SpacePairsRan4)
sd(SpacePairsRan4)
summary(SpaceTripletsRan4)
sd(SpaceTripletsRan4)







#ecdf
plot(ecdf(SpaceConsecutivePal),col="blue",main="Empirical Cumulative Distribution Function: Spacings")
plot(ecdf(SpaceConsecutiveRand),col="red",add=TRUE)
plot(ecdf(SpaceConsecutiveRand1),col="green",add=TRUE)
plot(ecdf(SpaceConsecutiveRand2),col="purple",add=TRUE)
plot(ecdf(SpaceConsecutiveRand3),col="yellow",add=TRUE)
plot(ecdf(SpaceConsecutiveRand4),col="cyan",add=TRUE)
legend("bottomright", legend=c("CMV Data", "Random(Medium)","Random(Group1)","Random(Group2)","Random(Group3)","Random(Group4)"),horiz=FALSE, cex=0.8,
       fill=c("blue","red","green","purple","yellow","cyan"))


#ecdf
plot(ecdf(SpacePairsCMV),col="blue",main="Ecdf: Spacings between sum of pairs")
plot(ecdf(SpacePairsRan),col="red",add=TRUE)
plot(ecdf(SpacePairsRan1),col="green",add=TRUE)
plot(ecdf(SpacePairsRan2),col="purple",add=TRUE)
plot(ecdf(SpacePairsRan3),col="yellow",add=TRUE)
plot(ecdf(SpacePairsRan4),col="cyan",add=TRUE)
legend("bottomright", legend=c("CMV Data", "Random(Medium)","Random(Group1)","Random(Group2)","Random(Group3)","Random(Group4)"),horiz=FALSE, cex=0.8,
       fill=c("blue","red","green","purple","yellow","cyan"))

plot(ecdf(SpaceTripletsCMV),col="blue",main="Ecdf: Spacings between sum of triplets")
plot(ecdf(SpaceTripletsRan),col="red",add=TRUE)
plot(ecdf(SpaceTripletsRan1),col="green",add=TRUE)
plot(ecdf(SpaceTripletsRan2),col="purple",add=TRUE)
plot(ecdf(SpaceTripletsRan3),col="yellow",add=TRUE)
plot(ecdf(SpaceTripletsRan4),col="cyan",add=TRUE)
legend("bottomright", legend=c("CMV Data", "Random(Medium)","Random(Group1)","Random(Group2)","Random(Group3)","Random(Group4)"),horiz=FALSE, cex=0.8,
       fill=c("blue","red","green","purple","yellow","cyan"))

#Q4


# Generate data
X = sample(1:L, size = 296, replace = F)

# Histogram 1 when interval length is 10
M1 = L/10   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M1),  xlab="Location",main="Frequency Distribution with 10 Intervals")
m1 = L/M1    # Number of categories
m1
expected = 296/m1
abline(h = expected, col='blue', lwd = 2)

# Histogram 2 when interval length is 50
M2 = L/50   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M2), xlab="Location",main="Frequency Distribution with 50 Intervals")
m2 = L/M2   # Number of categories
m2
expected = 296/m2
abline(h = expected, col='blue', lwd = 2)

# Histogram 3 when interval length is 100
M3 = L/100   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M3),  xlab="Location",main="Frequency Distribution with 100 Intervals")
m3 = L/M3    # Number of categories
m3
expected = 296/m3
abline(h = expected, col='blue', lwd = 2)

M4 = L/500   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M4),xlab="Location", main="Frequency Distribution with 500 Intervals")
m4 = L/M4    # Number of categories
m4
expected = 296/m4
abline(h = expected, col='blue', lwd = 2)
M4

M5 = L/1000   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M5),  xlab="Location",main="Frequency Distribution with 1000 Intervals")
m5 = L/M5    # Number of categories
m5
expected = 296/m5
abline(h = expected, col='blue', lwd = 2)
M5

hist(DT$location,xlim=c(90000,95000),breaks = 100,xlab="Location",main="Histogram with Interval between 90,000th and 95,000th")

hist(DT$location,xlim=c(195000,200000),breaks = 100,xlab="Location",main="Histogram with Interval between 195,000th and 200,000th")



#Advanced K.S. Test
ks.test(medians.sim,DT$location)
ks.test(sim.group1,DT$location)
ks.test(sim.group2,DT$location)
ks.test(sim.group3,DT$location)
ks.test(sim.group4,DT$location)


