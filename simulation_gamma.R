###########################PARAMETRIC(GAMMA)##############################
#Initialization parameters
step1 = 0.5           #change size scale
step2 = 0.1           #change size shape
size = 10000          #Monte-carlo sample size
point = 15            #change-point location
length = 30           #length of data
alpha = 0.05          #significance level
Psc = 12              #boundary cost-function(depends on data/assumptions)
Psh = 5
T1 = 0
T2 = 0
KS = 0
val = 0
Q = 0
beta = 0

countsc = 0           #used for computation power
countsh = 0
countkssc = 0
countkssh = 0
countcvmsc = 0
countcvmsh = 0

#####################################################################################################
#Phase 1(empirical critical values)
statcritsc <- matrix(NA, nrow = length - 3, ncol = 1)       #Param. scale
stat1critsc <- matrix(NA, nrow = size, ncol = 1)
statcritsh <- matrix(NA, nrow = length - 3, ncol = 1)       #Param. shape
stat1critsh <- matrix(NA, nrow = size, ncol = 1)
statcritks <- matrix(NA, nrow = length - 3, ncol = 1)       #Kolmogorov Smirnoff
stat1critks <- matrix(NA, nrow = size, ncol = 1)
statcritcvm <- matrix(NA, nrow = length - 3, ncol = 1)      #Cramer von Mises
stat1critcvm <- matrix(NA, nrow = size, ncol = 1)
for (i in 1:size){
  chg1 <- c(rgamma(length, shape = 1, scale = 1))
  vect <- 1:(length - 1)
  vect2 <- 1:length
  vect3 <- 2*vect2 - length - 1
  chg1stat <- vect*chg1[2:length]
  chg1stat1 <- vect3*log(chg1)
  for (j in 2:(length-2)){
    point = j
    chg1statright <- vect[1:(length - point)]*chg1[(point+1):(length)]
    statcritsc[j - 1] <- (length*1/sum(chg1))*sum(chg1stat)-((point*1/sum(chg1[1:point]))*sum(chg1stat[1:point]) + ((length - point)*1/sum(chg1[(point+1):length]))*sum(chg1statright))
    vect3left <- 2*vect2[1:point]-point-1
    chg1stat1left <- vect3left*log(chg1[1:point])
    vect3right <- 2*vect2[1:(length-point)]-(length-point)-1
    chg1stat1right <- vect3right*log(chg1[(point+1):length])
    statcritsh[j - 1] <- abs(sum(chg1stat1)) - (abs(sum(chg1stat1left))+abs(sum(chg1stat1right)))
    V <- chg1[1:point] 
    W <- chg1[(point+1):length]
    T1 <- length(V)
    T2 <- length(W)
    P = ecdf(V)
    O = ecdf(W)
    z = seq(min(chg1), max(chg1), by=0.1)
    lengthz = length(z)
    dif <- matrix(NA, nrow = lengthz, ncol = 1)
    for (k in 1:lengthz){
      dif[k] <- abs(P(k) - O(k))  
    }
    KS <- max(dif)
    if(T1 > 2*T2){
      beta = 1/(2*sqrt(T1))
    } else {
      if((T2 <= T1 & T1 <= 2*T2) | (T1/T2 == round(T1/T2)) == TRUE){
        beta = 2/(3*sqrt(T1))
      } else {
        beta = 2/(5*sqrt(T1))
      }
    }
    val = KS*sqrt(T1*T2/(T1 + T2)) + beta
    Q = 2*(exp(-2*val^{2}) - exp(-8*val^{2}))
    statcritks[j - 1] <- 1 - Q
    dif1 <- sum((P(chg1) - O(chg1))^{2})
    statcritcvm[j - 1] <- ((T1*T2)/(T1 + T2))*dif1
  }
  stat1critsc[i] <- max(abs(statcritsc))
  stat1critsh[i] <- max(abs(statcritsh))
  stat1critks[i] <- max(abs(statcritks))
  stat1critcvm[i] <- max(abs(statcritcvm))
} 
critsc <- quantile(stat1critsc, 1 - alpha/(length-3), type=1)
critsh <- quantile(stat1critsh, 1 - alpha/(length-3), type=1)
critks <- quantile(stat1critks, 1 - alpha/(length-3), type=1)
critcvm <- quantile(stat1critcvm, 1 - alpha/(length-3), type=1)
########################################################################
point = 15         #reinitialisize parameter
########################################################################
#Phase 2(determine power and cp-location estimate)
statsc <- matrix(NA, nrow = length - 3, ncol = 1)         #Param. scale
stat1sc <- matrix(NA, nrow = size, ncol = 1)
statsh <- matrix(NA, nrow = length - 3, ncol = 1)         #Param. shape
stat1sh <- matrix(NA, nrow = size, ncol = 1)
statkssc <- matrix(NA, nrow = length - 3, ncol = 1)       #KS, scale
stat1kssc <- matrix(NA, nrow = size, ncol = 1)
statkssh <- matrix(NA, nrow = length - 3, ncol = 1)       #KS, shape
stat1kssh <- matrix(NA, nrow = size, ncol = 1)
statcvmsc <- matrix(NA, nrow = length - 3, ncol = 1)      #CVM, scale
stat1cvmsc <- matrix(NA, nrow = size, ncol = 1)
statcvmsh <- matrix(NA, nrow = length - 3, ncol = 1)      #CVM, shape
stat1cvmsh <- matrix(NA, nrow = size, ncol = 1)
cpesc <- matrix(NA, nrow = size, ncol = 1)
cpesh <- matrix(NA, nrow = size, ncol = 1)
cpekssc <- matrix(NA, nrow = size, ncol = 1)
cpekssh <- matrix(NA, nrow = size, ncol = 1)
cpecvmsc <- matrix(NA, nrow = size, ncol = 1)
cpecvmsh <- matrix(NA, nrow = size, ncol = 1)
for (i in 1:size){
  chg1 <- c(rgamma(point, shape = 1, scale = 1), rgamma(length - point, shape = 1, scale = 1+step1))      #scale change
  vect <- 1:(length - 1)
  chg1stat <- vect*chg1[2:length]
  chg2 <- c(rgamma(point, shape = 1, scale = 1), rgamma(length - point, shape = 1 + step2, scale = 1))    #shape change
  vect2 <- 1:length
  vect3 <- 2*vect2 - length - 1
  chg2stat1 <- vect3*log(chg2)
  for (j in 2:(length-2)){
    if(abs(chg1[j]) < Psc){
      chg1statright <- vect[1:(length - j)]*chg1[(j+1):(length)]
      statsc[j - 1] <- (length*1/sum(chg1))*sum(chg1stat)-((j*1/sum(chg1[1:j]))*sum(chg1stat[1:j]) + ((length - j)*1/sum(chg1[(j+1):length]))*sum(chg1statright))
      Vsc = chg1[1:j]
      Wsc = chg1[(j+1):length]
      T1 <- length(Vsc)
      T2 <- length(Wsc)
      Psc = ecdf(Vsc)
      Osc = ecdf(Wsc)
      zsc = seq(min(chg1), max(chg1), by=0.1)
      lengthzsc = length(zsc)
      difsc <- matrix(NA, nrow = lengthzsc, ncol = 1)
      for (k in 1:lengthzsc){
        difsc[k] <- abs(Psc(k) - Osc(k))  
      }
      KSsc <- max(difsc)
      if(T1 > 2*T2){
        beta = 1/(2*sqrt(T1))
      } else if((T2 <= T1 & T1 <= 2*T2) | ((T1/T2 == round(T1/T2)) == TRUE )){
        beta = 2/(3*sqrt(T1))
      } else {
        beta = 2/(5*sqrt(T1))
      }
      valsc = KSsc*sqrt(T1*T2/(T1 + T2)) + beta
      Qsc = 2*(exp(-2*valsc^{2}) - exp(-8*valsc^{2}))
      statkssc[j - 1] <- 1 - Qsc
      dif1sc <- sum(abs(Psc(chg1) - Osc(chg1))^{2})
      statcvmsc[j - 1] <- ((T1*T2)/(T1 + T2))*dif1sc
    } else {
      statsc[j - 1] <- 0
      statkssc[j - 1] <- 0
      statcvmsc[j - 1] <- 0
    }
    if(abs(chg2[j]) < Psh){
      vect3left <- 2*vect2[1:j]-j-1
      chg2stat1left <- vect3left*log(chg2[1:j])
      vect3right <- 2*vect2[1:(length-j)]-(length-j)-1
      chg2stat1right <- vect3right*log(chg2[(j+1):length])
      statsh[j - 1] <- abs(sum(chg2stat1))-(abs(sum(chg2stat1left))+abs(sum(chg2stat1right)))
      Vsh = chg2[1:j]
      Wsh = chg2[(j+1):length]
      T1 = length(Vsh)
      T2 = length(Wsh)
      Psh = ecdf(Vsh)
      Osh = ecdf(Wsh)
      zsh = seq(min(chg2), max(chg2), by=0.1)
      lengthzsh = length(zsh)
      difsh <- matrix(NA, nrow = lengthzsh, ncol = 1)
      for (k in 1:lengthzsh){
        difsh[k] <- abs(Psh(k) - Osh(k))  
      }
      KSsh <- max(difsh)
      if(T1 > 2*T2){
        beta = 1/(2*sqrt(T1))
      } else if((T2 <= T1 & T1 <= 2*T2) | ((T1/T2 == round(T1/T2)) == TRUE )){
        beta = 2/(3*sqrt(T1))
      } else {
        beta = 2/(5*sqrt(T1))
      }
      valsh = KSsh*sqrt(T1*T2/(T1 + T2)) + beta
      Qsh = 2*(exp(-2*valsh^{2}) - exp(-8*valsh^{2}))
      statkssh[j - 1] <- 1 - Qsh
      dif1sh <- sum(abs(Psh(chg2) - Osh(chg2))^{2})
      statcvmsh[j - 1] <- ((T1*T2)/(T1 + T2))*dif1sh
    } else {
      statsh[j - 1] <- 0
      statkssh[j - 1] <- 0
      statcvmsh[j - 1] <- 0
    }
  }
  stat1sc[i] <- max(abs(statsc))
  stat1sh[i] <- max(abs(statsh))
  stat1kssc[i] <- max(abs(statkssc))
  stat1kssh[i] <- max(abs(statkssh))
  stat1cvmsc[i] <- max(abs(statcvmsc))
  stat1cvmsh[i] <- max(abs(statcvmsh))
  if(stat1sc[i] > critsc){                     #check significance and determine cp-location
    countsc = countsc + 1
    cpesc[i] <- which.max(abs(statsc))
  } else {
    countsc = countsc
    cpesc[i] <- 0
  }
  if(stat1sh[i] > critsh){
    countsh = countsh + 1
    cpesh[i] <- which.max(abs(statsh))
  } else {
    countsh = countsh
    cpesh[i] <- 0
  }
  if(stat1kssc[i] > critks){
    countkssc = countkssc + 1
    cpekssc[i] <- which.max(abs(statkssc)) + 1
  } else {
    countkssc = countkssc
    cpekssc[i] <- 0
  }
  if(stat1kssh[i] > critks){
    countkssh = countkssh + 1
    cpekssh[i] <- which.max(abs(statkssh)) + 1
  } else {
    countkssh = countkssh
    cpekssh[i] <- 0
  }
  if(stat1cvmsc[i] > critcvm){
    countcvmsc = countcvmsc + 1
    cpecvmsc[i] <- which.max(abs(statcvmsc)) + 1
  } else {
    countcvmsc = countcvmsc
    cpecvmsc[i] <- 0
  }
  if(stat1cvmsh[i] > critcvm){
    countcvmsh = countcvmsh + 1
    cpecvmsh[i] <- which.max(abs(statcvmsh)) + 1
  } else {
    countcvmsh = countcvmsh
    cpecvmsh[i] <- 0
  }
} 
cpesc <- cpesc[cpesc != 0]
cpesh <- cpesh[cpesh != 0]
cpekssc <- cpekssc[cpekssc != 0]
cpekssh <- cpekssh[cpekssh != 0]
cpecvmsc <- cpecvmsc[cpecvmsc != 0]
cpecvmsh <- cpecvmsh[cpecvmsh != 0]
powersc = countsc / size 
powersh = countsh / size
powerkssc = countkssc / size
powerkssh = countkssh / size
powercvmsc = countcvmsc / size
powercvmsh = countcvmsh / size

powersc;powersh;powerkssc;powerkssh;powercvmsc;powercvmsh
mean(cpesc);mean(cpesh);mean(cpekssc);mean(cpekssh);mean(cpecvmsc);mean(cpecvmsh)
sd(cpesc);sd(cpesh);sd(cpekssc);sd(cpekssh);sd(cpecvmsc);sd(cpecvmsh)

