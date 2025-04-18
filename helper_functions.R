library(MASS)
library(nlme)
library(lme4)
library(dplyr)

# Design matrix for time effects
bmf = function(lt,m){

  ns=lt-1 # lt: number of time periods, 5 or 11
  bm=kronecker(rep(1,ns*m),diag(1,lt)) # Computes the generalized Kronecker product of two arrays
  colnames(bm)=paste0("b",1:lt)	

  return(bm)
}

# Indicator matrix for intervention
xmf = function(lt,m){

  ns=lt-1
  xm=c(lower.tri(diag(1,lt))[,-lt]*1)
  xm=kronecker(diag(m),xm)
  colnames(xm)=paste0("x",1:m)

  return(xm)
}

# Design matrix for exposure time
zmf = function(lt,m){

  ns=lt-1
  zm=rbind(0,diag(1,ns))

  for(k in 2:ns){
    m01=matrix(0,k,ns)
    m02=matrix(0,lt-k,k-1)
    zmk=rbind(m01,cbind(diag(1,lt-k),m02))
    zm=rbind(zm,zmk)
  }
  
  zm=kronecker(diag(m),zm)
  colnames(zm)=paste0("z",1:(m*ns))
  return(zm)
}

fmatrix = function(time,exind){

  m=2
  ns=time-1
  nts=time*ns

  cond=((2<=exind)&(exind<=ns))
  if(!cond) print("error")

  xm1=c(lower.tri(diag(1,time))[,-time]*1)

  zm=rbind(0,diag(1,ns))

  for(k in 2:ns){
    m01=matrix(0,k,ns)
    m02=matrix(0,time-k,k-1)
    zmk=rbind(m01,cbind(diag(1,time-k),m02))
    zm=rbind(zm,zmk)
  }
  
  zm1=zm
  eseq=exind:ns
  le=length(eseq)

  xm2=rep(0,nts)
  zm2=matrix(0,nts,ns)

  for(k in 1:le){
    zind=which(zm1[,eseq[k]]==1)
    zm2[zind,k]=1
    xm2[zind]=1
  }

  zm=rbind(cbind(zm1,zm2),cbind(zm2,zm1))
  colnames(zm)=paste0("z",1:(m*ns))

  xm=rbind(cbind(xm1,xm2),cbind(xm2,xm1))
  colnames(xm)=paste0("x",1:m)

  return(list(xm=xm,zm=zm))
}


supp_matrix = function(time,exind){
  
  m=2
  ns=time-1
  nts=time*ns
  
  cond=((2<=exind)&(exind<=ns))
  if(!cond) print("error")
  
  xm1=c(lower.tri(diag(1,time))[,-time]*1)
  
  zm=rbind(0,diag(1,ns))
  
  for(k in 2:ns){
    m01=matrix(0,k,ns)
    m02=matrix(0,time-k,k-1)
    zmk=rbind(m01,cbind(diag(1,time-k),m02))
    zm=rbind(zm,zmk)
  }
  
  zm1=zm
  eseq=exind:ns
  le=length(eseq)
  xm2=rep(0,nts)
  zm2=matrix(0,nts,ns)
  
  for(k in 1:le){
    zind=which(zm1[,eseq[k]]==1)
    zm2[zind,k]=1
    xm2[zind]=1
  }
  
  zm=rbind(cbind(zm1,zm2),cbind(zm1,zm2))
  colnames(zm)=paste0("z",1:(m*ns))
  
  xm=rbind(cbind(xm1,xm2),cbind(xm1,xm2))
  colnames(xm)=paste0("x",1:m)
  
  return(list(xm=xm,zm=zm))
}


model = function(time,m,bd,dd){

  ns=time-1 

  ldeltav=seq(bd,dd,length=m*ns)
  ldeltam=matrix(ldeltav,nrow=ns,ncol=m) 

  cons=round(colMeans(ldeltam),2) 
  condeltav=rep(cons,each=ns) 

  lagdeltav1=c(rep(0,ns/2),rep(2*cons[1],ns/2))
  lagdeltav2=c(rep(0,ns/2),rep(2*cons[2],ns/2))
  lagdeltav=c(lagdeltav1,lagdeltav2)

  lagdeltav1_my=c(0,rep(ns*cons[1]/(ns-1),ns-1))
  lagdeltav2_my=c(0,rep(ns*cons[2]/(ns-1),ns-1))
  lagdeltav_my=c(lagdeltav1_my,lagdeltav2_my)

  con1 = cons[1]
  con2 = cons[2]
  fv1=con1*log(seq(ns/2,2*ns,length=ns))
  nldeltav1=con1+fv1-mean(fv1)
  fv2=con2*exp(seq(-ns/2,0.1,length=ns))
  nldeltav2=con2+fv2-mean(fv2)
  nldeltav=c(nldeltav1,nldeltav2)

  return(list(condeltav=condeltav, # constant effect
            ldeltav=ldeltav, # linearly increasing effect
            lagdeltav=lagdeltav, # half-period lagged effect
            lagdeltav_my =lagdeltav_my, # one-period lagged effect
            nldeltav = nldeltav, # curved effect
            cons=cons))
}



# Model fitting 

## Constant treatment effect model
fit_A=function(sdas,yt,lt,k){
  
  # column names
  bf=paste0("b",1:lt)
  xf=paste0("x",1:k)
  
  fa=as.formula(paste0(yt,"~",paste0(paste0(c(bf,xf),collapse="+"),"-1")))
  co=lme(fa,random=~1|cid,data=sdas,method="ML")
  
  eind=-(1:lt)
  xf=co$coef$fixed[eind]
  
  return(xf)
}

## Time-varying fixed treatment effect model
fit_B=function(sdas,yt,lt,k){
  
  bf=paste0("b",1:lt)
  zf=paste0("z",1:(k*(lt-1)))
  
  fo=lme(as.formula(paste0(yt,"~",paste0(paste0(c(bf,zf),collapse="+"),"-1"))),
         random=~1|cid,data=sdas,method="ML")
  
  ns=lt-1
  eind=-(1:lt)
  zf=fo$coef$fixed[eind] 
  x1_est_tv = mean(zf[1:ns])
  x2_est_tv = mean(zf[(ns+1): (2*ns)])
  
  return(c(x1_est_tv, x2_est_tv))
}

## Random treatment effect model
rfit3=function(sdas,yt,lt,k){
  if(lt == 5){
   sdas = sdas %>%
    mutate(expt1 = ifelse(z1 == 1, 1, ifelse(z2 ==1, 2, ifelse(z3 == 1, 3, ifelse(z4 == 1, 4, 0)))),
         expt2 = ifelse(z5 == 1, 1, ifelse(z6 ==1, 2, ifelse(z7 == 1, 3, ifelse(z8 == 1, 4, 0)))))
  }else if(lt == 11){
    sdas = sdas %>%
      mutate(expt1 = ifelse(z1 == 1, 1, 
                            ifelse(z2 ==1, 2, 
                                   ifelse(z3 == 1, 3, 
                                          ifelse(z4 == 1, 4, 
                                                 ifelse(z5 == 1, 5, 
                                                        ifelse(z6 == 1, 6, 
                                                               ifelse(z7 == 1, 7, 
                                                                      ifelse(z8 == 1, 8, 
                                                                             ifelse(z9 == 1, 9, 
                                                                                    ifelse(z10 == 1, 10, 0)))))))))),
             expt2 = ifelse(z11 == 1, 1, 
                            ifelse(z12 ==1, 2, 
                                   ifelse(z13 == 1, 3, 
                                          ifelse(z14 == 1, 4, 
                                                 ifelse(z15 == 1, 5, 
                                                        ifelse(z16 == 1, 6, 
                                                               ifelse(z17 == 1, 7, 
                                                                      ifelse(z18 == 1, 8, 
                                                                             ifelse(z19 == 1, 9, 
                                                                                    ifelse(z20 == 1, 10, 0)))))))))))
  }
   
  sdas = sdas %>%
    mutate(x1_str= paste0(expt1, cid),
           x2_str= paste0(expt2, cid))
  
   ff = as.formula(paste0(yt, "~ as.factor(time) + x1 + x2 + (0 + x1 | x1_str) + (0 + x2 | x2_str) + (1|cid)"))
   ro = lmer(ff,data = sdas)
   
   eind=-(1:lt)
   rf=summary(ro)$coef[eind, "Estimate"]
   return(rf)
}


