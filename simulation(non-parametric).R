##############################NON-PARAMETRIC##############################
#Initialization parameters
step1 = 1.1        #change size(mean, student-t)
step2 = -2         #change size(degrees of freedom)
size = 1000        #Monte-carlo sample size
point = 50         #change-point location
length = 100       #length of data
alpha = 0.05       #significance level
T1 = 0
T2 = 0
KS = 0
val = 0
Q = 0
beta = 0
ranks <- matrix(NA, nrow = length, ncol = 1)

countmw = 0                  #used for computation power
countmood = 0
countdavid = 0
countab = 0
countksmean = 0
countkssigma = 0
countcvmmean = 0
countcvmsigma = 0

#####################################################################################################
#Phase 1(empirical critical values)
statcritmw <- matrix(NA, nrow = length - 3, ncol = 1)           #Mann-Whitney
stat1critmw <- matrix(NA, nrow = size, ncol = 1)
statcritmood <- matrix(NA, nrow = length - 3, ncol = 1)         #Mood's test
stat1critmood <- matrix(NA, nrow = size, ncol = 1)   
statcritdavid <- matrix(NA, nrow = length - 3, ncol = 1)        #David's test
stat1critdavid <- matrix(NA, nrow = size, ncol = 1)
statcritab <- matrix(NA, nrow = length - 3, ncol = 1)           #Ansari-Bradley
stat1critab <- matrix(NA, nrow = size, ncol = 1)
statcritks <- matrix(NA, nrow = length - 3, ncol = 1)           #Kolmogorov-Smirnoff
stat1critks <- matrix(NA, nrow = size, ncol = 1)
statcritcvm <- matrix(NA, nrow = length - 3, ncol = 1)          #Cramer von Mises
stat1critcvm <- matrix(NA, nrow = size, ncol = 1)
for (i in 1:size){
  chg1 <- c(rt(length, 3))
  for (j in 2:(length-2)){
      point = j
      ranks <- rank(chg1)
      V <- chg1[1:point]
      W <- chg1[(point+1):length]
      T1 <- length(V)
      T2 <- length(W)
      rankV <- ranks[1:point]
      rankW <- ranks[(point+1):length]
      rankmean <- sum(rank(V))/point
      statcritmw[j - 1] <- (sum(ranks[1:point] - (length + 1)/2) - (T1*T2/2))/(sqrt(T1*T2*(T1 + T2 + 1)/12))
      statcritdavid[j - 1] <- (sum((ranks[1:point] - rankmean)^{2})/(point - 1) - (length*(length + 1))/12)/(sqrt(length*(length - point)*(length  + 1)*(3*(length + 1)*(point + 1) - length*point)/(360*point*(point - 1))))
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
  stat1critmw[i] <- max(abs(statcritmw))
  stat1critdavid[i] <- max(abs(statcritdavid))
  stat1critks[i] <- max(abs(statcritks))
  stat1critcvm[i] <- max(abs(statcritcvm))
} 
critmw <- quantile(stat1critmw, 1 - alpha/(length-3), type=1)         #critical values
critdavid <- quantile(stat1critdavid, 1 - alpha/(length-3), type=1)
critks <- quantile(stat1critks, 1 - alpha/(length-3), type=1)
critcvm <- quantile(stat1critcvm, 1 - alpha/(length-3), type=1)
########################################################################
point = 50         #reinitialisize parameter
########################################################################
#Phase 2(determine power and cp-location estimate)
statmw <- matrix(NA, nrow = length - 3, ncol = 1)
stat1mw <- matrix(NA, nrow = size, ncol = 1)
statmood <- matrix(NA, nrow = length - 3, ncol = 1)
stat1mood <- matrix(NA, nrow = size, ncol = 1)
statdavid <- matrix(NA, nrow = length - 3, ncol = 1)
stat1david <- matrix(NA, nrow = size, ncol = 1)
statab <- matrix(NA, nrow = length - 3, ncol = 1)
stat1ab <- matrix(NA, nrow = size, ncol = 1)
statksmean <- matrix(NA, nrow = length - 3, ncol = 1)
stat1ksmean <- matrix(NA, nrow = size, ncol = 1)
statkssigma <- matrix(NA, nrow = length - 3, ncol = 1)
stat1kssigma <- matrix(NA, nrow = size, ncol = 1)
statcvmmean <- matrix(NA, nrow = length - 3, ncol = 1)
stat1cvmmean <- matrix(NA, nrow = size, ncol = 1)
statcvmsigma <- matrix(NA, nrow = length - 3, ncol = 1)
stat1cvmsigma <- matrix(NA, nrow = size, ncol = 1)
cpemw <- matrix(NA, nrow = size, ncol = 1)
cpemood <- matrix(NA, nrow = size, ncol = 1)
cpedavid <- matrix(NA, nrow = size, ncol = 1)
cpeab <- matrix(NA, nrow = size, ncol = 1)
cpeksmean <- matrix(NA, nrow = size, ncol = 1)
cpekssigma <- matrix(NA, nrow = size, ncol = 1)
cpecvmmean <- matrix(NA, nrow = size, ncol = 1)
cpecvmsigma <- matrix(NA, nrow = size, ncol = 1)
for (i in 1:size){
  chg1 <- c(rt(point, 3), (rt(length - point, 3) + step1))        #mean change
  chg2 <- c(rt(point, 3), (rt(length - point, 3 + step2)))        #df change
  for (j in 2:(length-2)){
    ranks1 <- rank(chg1)
    ranks2 <- rank(chg2)
    Vmean = chg1[1:j]
    Wmean = chg1[(j+1):length]
    Vsigma = chg2[1:j]
    Wsigma = chg2[(j+1):length]
    T1 <- length(Vmean)
    T2 <- length(Wmean)
    rankV <- ranks2[1:j]
    rankW <- ranks2[(j + 1):length]
    rankmean <- sum(rank(Vsigma))/j
    Pmean = ecdf(Vmean)
    Psigma = ecdf(Vsigma)
    Omean = ecdf(Wmean)
    Osigma = ecdf(Wsigma)
    zmean = seq(min(chg1), max(chg1), by=0.1)
    lengthzmean = length(zmean)
    difmean <- matrix(NA, nrow = lengthzmean, ncol = 1)
    for (k in 1:lengthzmean){
      difmean[k] <- abs(Pmean(k) - Omean(k))
    }
    KSmean <- max(difmean)
    zsigma = seq(min(chg2), max(chg2), by=0.1)
    lengthzsigma = length(zsigma)
    difsigma <- matrix(NA, nrow = lengthzsigma, ncol = 1)
    for (k in 1:lengthzsigma){
      difsigma[k] <- abs(Psigma(k) - Osigma(k))
    }
    KSsigma <- max(difsigma)
    if(T1 > 2*T2){
      beta = 1/(2*sqrt(T1))
    } else if((T2 <= T1 & T1 <= 2*T2) | ((T1/T2 == round(T1/T2)) == TRUE )){
      beta = 2/(3*sqrt(T1))
    } else {
      beta = 2/(5*sqrt(T1))
    }
    valmean = KSmean*sqrt(T1*T2/(T1 + T2)) + beta
    Qmean = 2*(exp(-2*valmean^{2}) - exp(-8*valmean^{2}))
    dif1mean <- sum(abs(Pmean(chg1) - Omean(chg1))^{2})
    valsigma = KSsigma*sqrt(T1*T2/(T1 + T2)) + beta
    Qsigma = 2*(exp(-2*valsigma^{2}) - exp(-8*valsigma^{2}))
    dif1sigma <- sum(abs(Psigma(chg2) - Osigma(chg2))^{2})
    statmw[j - 1] <- (sum(ranks1[1:j] - (length + 1)/2) - (T1*T2/2))/(sqrt(T1*T2*(T1 + T2 + 1)/12))
    statksmean[j - 1] <- 1 - Qmean
    statcvmmean[j - 1] <- ((T1*T2)/(T1 + T2))*dif1mean
    statmood[j - 1] <- mood.test(Vsigma, Wsigma, alternative = c("two.sided"))$p.value
    statdavid[j - 1] <- (sum((ranks2[1:j] - rankmean)^{2})/(j - 1) - (length*(length + 1))/12)/(sqrt(length*(length - j)*(length  + 1)*(3*(length + 1)*(j + 1) - length*j)/(360*j*(j - 1))))
    statab[j - 1] <- ansari.test(Vsigma, Wsigma, alternative = c("two.sided"), exact = NULL, conf.level = 0.95)$p.value
    statkssigma[j - 1] <- 1 - Qsigma
    statcvmsigma[j - 1] <- ((T1*T2)/(T1 + T2))*dif1sigma
  }
  stat1mw[i] <- max(abs(statmw))
  stat1mood[i] <- min(statmood)
  stat1david[i] <- max(abs(statdavid))
  stat1ab[i] <- min(statab)
  stat1ksmean[i] <- max(abs(statksmean))
  stat1kssigma[i] <- max(abs(statkssigma))
  stat1cvmmean[i] <- max(abs(statcvmmean))
  stat1cvmsigma[i] <- max(abs(statcvmsigma))
  if(stat1mw[i] > critmw){             #check significance and determine cp-location
    countmw = countmw + 1
    cpemw[i] <- which.max(abs(statmw))
  } else {
    countmw = countmw
    cpemw[i] <- 0
  }
  if(stat1mood[i] < alpha){
    countmood = countmood + 1
    cpemood[i] <- which.min(abs(statmood))
  } else {
    countmood = countmood
    cpemood[i] <- 0
  }
  if(stat1david[i] > critdavid){
    countdavid = countdavid + 1
    cpedavid[i] <- which.max(abs(statdavid))
  } else {
    countdavid = countdavid
    cpedavid[i] <- 0
  }
  if(stat1ab[i] <= alpha){    
    countab = countab + 1
    cpeab[i] <- which.min(statab) + 1
  } else {
    countab = countab
    cpeab[i] <- 0
  }
  if(stat1ksmean[i] > critks){
    countksmean = countksmean + 1
    cpeksmean[i] <- which.max(abs(statksmean)) + 1
  } else {
    countksmean = countksmean
    cpeksmean[i] <- 0
  }
  if(stat1kssigma[i] > critks){
    countkssigma = countkssigma + 1
    cpekssigma[i] <- which.max(abs(statkssigma)) + 1
  } else {
    countkssigma = countkssigma
    cpekssigma[i] <- 0
  }
  if(stat1cvmmean[i] > critcvm){
    countcvmmean = countcvmmean + 1
    cpecvmmean[i] <- which.max(abs(statcvmmean)) + 1
  } else {
    countcvmmean = countcvmmean
    cpecvmmean[i] <- 0
  }
  if(stat1cvmsigma[i] > critcvm){
    countcvmsigma = countcvmsigma + 1
    cpecvmsigma[i] <- which.max(abs(statcvmsigma)) + 1
  } else {
    countcvmsigma = countcvmsigma
    cpecvmsigma[i] <- 0
  }
} 
cpemw <- cpemw[cpemw != 0]
cpemood <- cpemood[cpemood != 0]
cpedavid <- cpedavid[cpedavid != 0]
cpeab <- cpeab[cpeab != 0]
cpeksmean <- cpeksmean[cpeksmean != 0]
cpekssigma <- cpekssigma[cpekssigma != 0]
cpecvmmean <- cpecvmmean[cpecvmmean != 0]
cpecvmsigma <- cpecvmsigma[cpecvmsigma != 0]
powermw = countmw / size 
powermood = countmood / size
powerdavid = countdavid / size
powerab = countab / size
powerksmean = countksmean / size
powerkssigma = countkssigma / size
powercvmmean = countcvmmean / size
powercvmsigma = countcvmsigma / size

powermw;powermood;powerdavid;powerab;powerksmean;powerkssigma;powercvmmean;powercvmsigma
mean(cpemw);mean(cpemood);mean(cpedavid);mean(cpeab);mean(cpeksmean);mean(cpekssigma);mean(cpecvmmean);mean(cpecvmsigma)
sd(cpemw);sd(cpemood);sd(cpedavid);sd(cpeab);sd(cpeksmean);sd(cpekssigma);sd(cpecvmmean);sd(cpecvmsigma)

