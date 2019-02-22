segn=function(data,K=5,pen=0){
  L = length(data)
  segmconf=matrix(0,ncol=L,nrow=L)
  for(i in 1:L){                          #mean change
    ssq=0
    sumx=0
    for(j in i:L){
      len=j-i+1
      sumx=sumx+data[j]
      ssq=ssq+data[j]^2
      segmconf[i,j]=-0.5*(ssq-(sumx^2)/len)
    }
  }
  # for(i in 1:L){                          #variance change
  #   ssq=0
  #   for(j in i:L){
  #     len=j-i+1
  #     ssq=ssq+(data[j]-mean(data))^2
  #     if(ssq <= 0){
  #       varest=0.00000000001/len
  #     } else{
  #       varest=ssq/len
  #     }
  #     segmconf[i,j]=-(len/2)*(log(2*pi)+log(varest)+1)
  #   }
  # }
  likemat=matrix(0,ncol=L,nrow=K)
  likemat[1,]=segmconf[1,]
  cp=matrix(0,ncol=L,nrow=K)
  for(k in 2:K){
    for(j in k:L){
      discr=NULL                                            ######################
      intval=(k-1):(j-1)                                    #mean change
      discr=likemat[k-1,intval]+segmconf[intval+1,j]        ######################
      # if((j-2-k)<0){                                      #
      #   discr=-Inf                                        #  
      # }                                                   #variance change
      # else{                                               #
      #   intval=(k):(j-2)                                  #
      #   discr=likemat[k-1,intval]+segmconf[intval+1,j]    #
      # }                                                   ######################
      likemat[k,j]= max(discr,na.rm=TRUE)
      cp[k,j]=which(discr==max(discr,na.rm=TRUE))[1]+(k-2)
    }
  }
  cps=matrix(0,ncol=K,nrow=K)
  for(k in 2:K){
    cps[k,1]=cp[k,L]
    for(i in 1:(k-1)){
      cps[k,(i+1)]=cp[(k-i),cps[k,i]]
    }
  }
  
  signcps=NULL
  for(i in 1:length(pen)){
    thres=(-2*likemat[,L]) >= pen[i]
    if(sum(thres)==0){
      signcps=0
    }
    else{
      signcps=c(signcps,max(which((thres)==TRUE)))
    }
  }
  return(list(cps=cps,signcps=signcps,pen=pen,discr=thres[signcps+1]))
}
