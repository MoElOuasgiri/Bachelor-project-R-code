wbs=function(data,K=5,M=1000,pen=0){
                                              #cost functions:
  normalmean=function(xsq,x,len){             #normal model, mean change
    return(-0.5*(xsq-(x^2)/len))
  }
  
  # normalvar=function(x,len){                #normal model, variance change
  #   cond = x<=0
  #   x[cond==TRUE]=0.00000000001
  #   return(-0.5*len*(log(2*pi)+log(x/len)+1))
  # }
  
  L = length(data)
  ysq = c(0,cumsum(data^2))
  ysum = c(0,cumsum(data))
  xsq = c(cumsum(data^2))
  xsum = c(cumsum(data))
  tau = c(0,L)
  cpt = matrix(0,nrow=2,ncol=K)
  
  for(q in 1:K){
    discr <- rep(0,L-1)
    i=1
    s=tau[1]+1;e=tau[2]
    null=normalmean(ysq[e+1]-ysq[s],ysum[e+1]-ysum[s],e-s+1)
    #null=normalvar(ysq[e+1]-ysq[s],e-s+1)
    s1 <- sample(s:e, M, replace=T)
    e1 <- sample(s:e, M, replace=T)
    discr2 <- matrix(NA, nrow=2, ncol=M)
    discr3 <- matrix(NA, nrow=2, ncol=length(tau) - 1)
    for(g in 1:M){
      if((e1[g] - s1[g]) < 1){ 
        discr2[1,g]=0; discr2[2,g]=0
        } else {
          null1 = normalmean(xsq[e1[g]]-xsq[s1[g]],xsum[e1[g]]-xsum[s1[g]],e1[g]-s1[g]+1)
          #null1 = normalvar(xsq[e1[g]]-xsq[s1[g]],e1[g]-s1[g]+1)
          discr1 = rep(0, e1[g]-s1[g])
          for(t in 1:(e1[g]-s1[g])){
            discr1[t] <- normalmean(xsq[t+s1[g]]-xsq[s1[g]],xsum[t+s1[g]]-xsum[s1[g]],t+1)+normalmean(xsq[e1[g]]-xsq[t+s1[g]],xsum[e1[g]]-xsum[t+s1[g]],e1[g]+1-s1[g]-t)-null1
            #discr1[t] <- normalvar(xsq[t+s1[g]]-xsq[s1[g]], t+1)+normalvar(xsq[e1[g]]-xsq[t+s1[g]],e1[g]+1-s1[g]-t)-null1
          }
        }
      discr2[1,g]=max(discr1);discr2[2,g]=which.max(discr1) + s1[g]
    }
    discr3[1,i]=max(discr2[1,]);loc <- which.max(discr2[1,]);discr3[2,i]=discr2[2,loc]
    for(j in 1:(L-1)){
      if(j==e){
        s=e+1;i=i+1;e=tau[i+1]
        null=normalmean(ysq[e+1]-ysq[s],ysum[e+1]-ysum[s],e-s+1)
        #null=normalvar(ysq[e+1]-ysq[s],e-s+1)
        s1 <- sample(s:e, M, replace=T)
        e1 <- sample(s:e, M, replace=T)
        for(g in 1:M){
          if((e1[g] - s1[g]) < 1){ 
            discr2[1,g]=0; discr2[2,g]=0
          } else {
            null1 = normalmean(xsq[e1[g]]-xsq[s1[g]],xsum[e1[g]]-xsum[s1[g]],e1[g]-s1[g]+1)
            #null1 = normalvar(xsq[e1[g]]-xsq[s1[g]], e1[g]-s1[g]+1)
            discr1 = rep(0, e1[g]-s1[g])
            for(t in 1:(e1[g]-s1[g])){
              discr1[t] <- normalmean(xsq[t+s1[g]]-xsq[s1[g]],xsum[t+s1[g]]-xsum[s1[g]],t+1)+normalmean(xsq[e1[g]]-xsq[t+s1[g]],xsum[e1[g]]-xsum[t+s1[g]],e1[g]+1-s1[g]-t)-null1
              #discr1[t] <- normalvar(xsq[t+s1[g]]-xsq[s1[g]],t+1)+normalmean(xsq[e1[g]]-xsq[t+s1[g]],e1[g]+1-s1[g]-t)-null1
            }
          }
          discr2[1,g]=max(discr1);discr2[2,g]=which.max(discr1) + s1[g]
        }
        discr3[1,i]=max(discr2[1,]);loc <- which.max(discr2[1,]);discr3[2,i]=discr2[2,loc]
      }else{
        discr[j]=normalmean(ysq[j+1]-ysq[s],ysum[j+1]-ysum[s],j-s+1)+normalmean(ysq[e+1]-ysq[j+1],ysum[e+1]-ysum[j+1],e-j)-null
        #discr[j]=normalvar(ysq[j+1]-ysq[s],j-s+1)+normalvar(ysq[e+1]-ysq[j+1],e-j)-null
      }
    }
    discrwild <- max(discr3[1,]);k1 = which.max(discr3[1,])
    discrmax <- max(discr)
    if(max(discrwild,discrmax) == discrwild){
      k=discr3[2,k1]
      cpt[2,q]=discrwild  
      } else {
      k=which.max(discr)[1]
      cpt[2,q]=max(discr)
    }
    cpt[1,q]=k
    tau=sort(c(tau,k))
  }
  signcps=NULL                                             #significant change-points
  for(i in 1:length(pen)){
    thres=(2*cpt[2,])>=pen[i]                              #threshold for significance using penalty value
    if(sum(thres)==0){
      signcps=0
    }
    else{
      signcps=c(signcps,max(which((thres)==TRUE)))
    }
  }
  return(list(cpt=cpt,signcps=signcps,pen=pen))
}

