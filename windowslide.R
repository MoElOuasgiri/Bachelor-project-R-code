windowslide=function(data, width = 5){
                                                    #cost function:
  normalmean=function(xsq,x,len){                   #normal model, mean change         
    return(-0.5*(xsq-(x^2)/len))
  }
  
  # normalvar=function(x,len){                      #normal model, variance change
  #   cond = x<=0
  #   x[cond==TRUE]=0.00000000001
  #   return(-0.5*len*(log(2*pi)+log(x/len)+1))
  # }
  
  L = length(data)
  discr <- matrix(NA, nrow = 1, ncol = L - 2*width + 2)
  for (t in ((width+1) : (L - width + 2))){                #Compute discrepancy curve
    left <- data[(t - width):(t - 1)]
    right <- data[(t - 1):(t + width - 2)]
    window <- data[(t - width):(t + width - 2)]
    ysq=c(cumsum(window^2))
    ysum=c(cumsum(window))
    discr[t - width] <- normalmean(ysq[width]-ysq[1],ysum[width]-ysum[1],width)+normalmean(ysq[2*width - 1]-ysq[width],ysum[2*width - 1]-ysum[width],width)-normalmean(ysq[2*width - 1]-ysq[1],ysum[2*width - 1]-ysum[1],2*width-1)
    #discr[t - width] <- normalvar(ysq[width]-ysq[1],width)+normalvar(ysq[2*width - 1]-ysq[width],width)-normalvar(ysq[2*width - 1]-ysq[1],2*width-1)
  }
  return(discr)
}


