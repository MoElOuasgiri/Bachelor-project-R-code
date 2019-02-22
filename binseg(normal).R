bs=function(data,K=5,pen=0){
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
  tau = c(0,L)
  cpt = matrix(0,nrow=2,ncol=K)

  for(i in 1:K){
    discr <- rep(0,L-1)
    l=1
    s=tau[1]+1;e=tau[2]
    null=normalmean(ysq[e+1]-ysq[s],ysum[e+1]-ysum[s],e-s+1)
    #null=normalvar(ysq[e+1]-ysq[s],e-s+1)
    for(j in 1:(L-1)){
      if(j==e){
        s=e+1;l=l+1;e=tau[l+1]
        null=normalmean(ysq[e+1]-ysq[s],ysum[e+1]-ysum[s],e-s+1)
        #null=normalvar(ysq[e+1]-ysq[s],e-s+1)
      }else{
        discr[j]=normalmean(ysq[j+1]-ysq[s],ysum[j+1]-ysum[s],j-s+1)+normalmean(ysq[e+1]-ysq[j+1],ysum[e+1]-ysum[j+1],e-j)-null
        #discr[j] <- normalvar(ysq[j+1]-ysq[s],j-s+1)+normalvar(ysq[e+1]-ysq[j+1],e-j)-null
      }
    }
    k=which.max(discr)[1]
    cpt[1,i]=k;cpt[2,i]=max(discr)
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

