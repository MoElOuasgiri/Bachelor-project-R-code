##############################PARAMETRIC##############################
#Initialization parameters
step = 2            #change size(mean, variance, normal)
step1 = 1.1           #change size(mean, student-t)
step2 = 100            #change size(degrees of freedom)
size = 10000          #Monte-carlo sample size
point = 50            #change-point location
length = 100           #length of data
alpha = 0.05          #significance level
P = 5                 #boundary cost-function(depends on data/assumptions)

countmcvn = 0         #used for computation power
countmcvu = 0
countvcl = 0
countsicvar = 0

#####################################################################################################
#Phase 1(empirical critical values)
statcritmcvn <- matrix(NA, nrow = length - 3, ncol = 1)
stat1critmcvn <- matrix(NA, nrow = size, ncol = 1)
statcritmcvu <- matrix(NA, nrow = length - 3, ncol = 1)
stat1critmcvu <- matrix(NA, nrow = size, ncol = 1)
statcritvcl <- matrix(NA, nrow = length - 3, ncol = 1)
stat1critvcl <- matrix(NA, nrow = size, ncol = 1)
for (i in 1:size){
  chg1 <- c(rnorm(length, 0, 1))          #normal data
  #chg1 <- c(rt(length, 3))               #student-t data
  for (j in 2:(length-2)){
    if(abs(chg1[j]-median(chg1)) < P){
      point = j
      #statcritmcvn[j - 1] <- sum((chg1 - mean(chg1))^2) - (sum((chg1[1:point] - mean(chg1[1:point]))^2) + sum((chg1[(point+1):length] - mean(chg1[(point+1):length]))^2))      #mean change, var known
      statcritmcvn[j - 1] <- sum((chg1 - median(chg1))^2) - (sum((chg1[1:point] - median(chg1[1:point]))^2) + sum((chg1[(point+1):length] - median(chg1[(point+1):length]))^2))
      #statcritmcvu[j - 1] <- (point*(mean(chg1[1:point]) - mean(chg1))^{2} + (length - point)*(mean(chg1[(point+1):length]) - mean(chg1))^{2})/(sum((chg1 - mean(chg1))^{2}) - (point*(mean(chg1[1:point]) - mean(chg1))^{2} + (length - point)*(mean(chg1[(point+1):length]) - mean(chg1))^{2}))    #mean change, var unknown
      statcritmcvu[j - 1] <- (point*(median(chg1[1:point]) - median(chg1))^{2} + (length - point)*(median(chg1[(point+1):length]) - median(chg1))^{2})/(sum((chg1 - median(chg1))^{2}) - (point*(median(chg1[1:point]) - median(chg1))^{2} + (length - point)*(median(chg1[(point+1):length]) - median(chg1))^{2}))
      statcritvcl[j - 1] <- length*log((sum((chg1 - 0)^2))/(length)) - point*log((sum((chg1[1:point] - 0)^2))/(point)) - (length - point)*log((sum((chg1[(point+1):length] - 0)^2))/(length - point))       #variance change
      } else {
      point = j
      statcritmcvn[j - 1] <- 0
      statcritmcvu[j - 1] <- 0
      statcritvcl[j - 1] <- 0
    }
  }
  stat1critmcvn[i] <- max(abs(statcritmcvn))
  stat1critmcvu[i] <- max(abs(statcritmcvu))
  stat1critvcl[i] <- max(abs(statcritvcl))
}

critmcvn <- quantile(stat1critmcvn, 1 - alpha/(length-3), type=1)     #critical values
critmcvu <- quantile(stat1critmcvu, 1 - alpha/(length-3), type=1)
critvcl <- quantile(stat1critvcl, 1 - alpha/(length-3), type=1)
alog <- (2*log(log(length)))^{0.5}
blog <- 2*log(log(length)) + 0.5*log(log(log(length))) - log(gamma(0.5))
critsicvar <- (-(log(log((1 - (alpha/(length-3)) + exp(-2*exp(blog)))^{-0.5}))/alog) + blog/alog)^{2} - log(length)
#####################################################################################
point = 50              #reinitialisize parameter
#####################################################################################
#Phase 2(determine power and cp-location estimate)
statmcvn <- matrix(NA, nrow = length - 3, ncol = 1)
stat1mcvn <- matrix(NA, nrow = size, ncol = 1)
statmcvu <- matrix(NA, nrow = length - 3, ncol = 1)
stat1mcvu <- matrix(NA, nrow = size, ncol = 1)
statvcl <- matrix(NA, nrow = length - 3, ncol = 1)
stat1vcl <- matrix(NA, nrow = size, ncol = 1)
statsicvark <- matrix(NA, nrow = length - 3, ncol = 1)
stat1sicvark <- matrix(NA, nrow = size, ncol = 1)
cpemcvn <- matrix(NA, nrow = size, ncol = 1)
cpemcvu <- matrix(NA, nrow = size, ncol = 1)
cpevcl <- matrix(NA, nrow = size, ncol = 1)
cpesicvar <- matrix(NA, nrow = size, ncol = 1)
for (i in 1:size){
  chg1 <- c(rnorm(point, 0, 1), rnorm(length - point, step, 1))             #normal mean change
  chg2 <- c(rnorm(point, 0, 1), rnorm(length - point, 0, (step*1)^{0.5}))   #normal var change
  #chg1 <- c(rt(point, 3), (rt(length - point, 3) + step1))                 #student t mean change
  #chg2 <- c(rt(point, 3), rt(length - point, 3 + step2))                   #student t var change 
  sicvarT <- length*log(2*pi) + length*log((sum((chg2 - 0)^2))/(length)) + length + log(length)   #MSB method
  for (j in 2:(length-2)){
    if(abs(chg1[j]-median(chg1)) < P){
      #statmcvn[j - 1] <- sum((chg1 - mean(chg1))^2) - (sum((chg1[1:j] - mean(chg1[1:j]))^2) + sum((chg1[(j+1):length] - mean(chg1[(j+1):length]))^2))
      statmcvn[j - 1] <- sum((chg1 - median(chg1))^2) - (sum((chg1[1:j] - median(chg1[1:j]))^2) + sum((chg1[(j+1):length] - median(chg1[(j+1):length]))^2))
      #statmcvu[j - 1] <- (j*(mean(chg1[1:j]) - mean(chg1))^{2} + (length - j)*(mean(chg1[(j+1):length]) - mean(chg1))^{2})/(sum((chg1 - mean(chg1))^{2}) - (j*(mean(chg1[1:j]) - mean(chg1))^{2} + (length - j)*(mean(chg1[(j+1):length]) - mean(chg1))^{2}))
      statmcvu[j - 1] <- (j*(median(chg1[1:j]) - median(chg1))^{2} + (length - j)*(median(chg1[(j+1):length]) - median(chg1))^{2})/(sum((chg1 - median(chg1))^{2}) - (j*(median(chg1[1:j]) - median(chg1))^{2} + (length - j)*(median(chg1[(j+1):length]) - median(chg1))^{2}))
      } else {
      statmcvn[j - 1] <- 0
      statmcvu[j - 1] <- 0
    }
    if(abs(chg2[j]-median(chg2)) < P){
      statvcl[j - 1] <- (length*log((sum((chg2 - 0)^2))/(length)) - j*log((sum((chg2[1:j] - 0)^2))/(j)) - (length - j)*log((sum((chg2[(j+1):length] - 0)^2))/(length - j)))^{0.5}
    } else {
      statvcl[j - 1] <- 0
    }
  }
  for (j in 2:(length-2)){
    statsicvark[j - 1] <- length*log(2*pi) + j*log((sum((chg2[1:j] - 0)^2))/(j)) + (length - j)*log((sum((chg2[(j+1):length] - 0)^2))/(length - j)) + length + 2*log(length)
  }
  stat1mcvn[i] <- max(abs(statmcvn))
  stat1mcvu[i] <- max(abs(statmcvu))
  stat1vcl[i] <- max(abs(statvcl))
  stat1sicvark[i] <- min(statsicvark)
  if(stat1mcvn[i] > critmcvn){          #check significance and determine cp-location
    countmcvn = countmcvn + 1
    cpemcvn[i] <- which.max(statmcvn) + 1
  } else {
    countmcvn = countmcvn
    cpemcvn[i] <- 0
  }
  if(stat1mcvu[i] > critmcvu){
    countmcvu = countmcvu + 1
    cpemcvu[i] <- which.max(statmcvu) + 1
  } else {
    countmcvu = countmcvu
    cpemcvu[i] <- 0
  }
  if(stat1vcl[i] > critvcl){
    countvcl = countvcl + 1
    cpevcl[i] <- which.max(statvcl) + 1
  } else {
    countvcl = countvcl
    cpevcl[i] <- 0
  }
  if(sicvarT > stat1sicvark[i] + critsicvar){
    countsicvar = countsicvar + 1
    cpesicvar[i] <- which.min(statsicvark) + 1
  } else {
    countsicvar = countsicvar
    cpesicvar[i] <- 0
  }
} 
cpemcvn <- cpemcvn[cpemcvn != 0]
cpemcvu <- cpemcvu[cpemcvu != 0]
cpevcl <- cpevcl[cpevcl != 0]
cpesicvar <- cpesicvar[cpesicvar != 0]
powermcvn = countmcvn / size 
powermcvu = countmcvu / size
powervcl = countvcl / size
powersicvar = countsicvar / size

powermcvn;powermcvu;powervcl;powersicvar
mean(cpemcvn);mean(cpemcvu);mean(cpevcl);mean(cpesicvar)
sd(cpemcvn);sd(cpemcvu);sd(cpevcl);sd(cpesicvar)