library(car)
library(tidyverse)
library(lmtest)
library(sandwich)
library(zoo)
library(lmtest)
library(data.table)

data <- data.frame(read.table('hcmv.txt',header = T))
DT <- data.table(data)

Countofpalindrome<- as.numeric(count(DT))
Countofpalindrome
Palindrome <- DT$location
Palindrome
sumofpairs <- c(0)
sumoftriplets <- c(0)

#Q2 TBD
for (i in 1:Countofpalindrome-1){
  append(sumofpairs,(Palindrome[i] + Palindrome[i+1]))
}

for (i in 1:Countofpalindrome-2){
  append(sumoftriplets,(Palindrome[i] + Palindrome[i+1] + Palindrome[i+2]))
}


sumofpairs

sumoftriplets


#Q4

L = max(DT$location)
L

# Generate data
X = sample(1:L, size = 296, replace = F)

# Histogram 1 when interval length is 10
M1 = L/10   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M1))
m1 = L/M1    # Number of categories
m1
expected = 296/m1
abline(h = expected, col='blue', lwd = 2)

#Chi-square for interval=10
chi.square.stat = sum((DT$location - expected)^2/expected)
chi.square.stat
d.f = m1 - 1
?pchisq
1 - pchisq(chi.square.stat, df = d.f)


# Histogram 2 when interval length is 20
M2 = L/20  # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M2))
m2 = L/M2    # Number of categories
m2
expected = 296/m2
abline(h = expected, col='blue', lwd = 2)

# Histogram 3 when interval length is 50
M3 = L/50   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M3))
m3 = L/M3    # Number of categories
m3
expected = 296/m3
abline(h = expected, col='blue', lwd = 2)

# Histogram 4 when interval length is 100
M4 = L/100   # interval length
X.hist = hist(DT$location, breaks = seq(0, L, by = M4))
m4 = L/M4    # Number of categories
m4
expected = 296/m4
abline(h = expected, col='blue', lwd = 2)






