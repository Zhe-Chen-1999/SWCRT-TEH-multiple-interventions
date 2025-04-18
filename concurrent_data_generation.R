source("helper_functions.R")

icc=0.05 
tt=3 
sigma2_a=tt*icc 
sigma2_e=tt*(1-icc) 
nsim=500
nbp=500
m=2
times=c(5,11)
lts=length(times)
bbeta=0.1
dbeta=0.5
trsatv=rep(0.2,2)

sobj=list(n=n, icc=icc,sigma2_a=sigma2_a,sigma2_e=sigma2_e,nsim=nsim,times=times, 
          lts=lts,bd=bd,dd=dd,bbeta=bbeta,dbeta=dbeta,m=m,nbp=nbp,trsatv=trsatv)


for (bd in c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)){
  
  dd = bd + 0.5
  
  for (n in c(30, 100, 500)) {
  
    onev=rep(1,n*nsim)

    for(ts in 1:lts){

      lt=times[ts] # T
      ns=lt-1 

      bm=bmf(lt,m)
      zm=zmf(lt,m)
      nz=ncol(zm)
      xm=xmf(lt,m)
      nx=ncol(xm)

      fcovy=sigma2_a*matrix(1,lt,lt)+diag(sigma2_e,lt) 

      # Fixed time effects
      betav=seq(bbeta,dbeta,length=lt) 

      # Time-varying treatment effects
      simod=model(lt,m,bd,dd) 

      bmb=bm%*%betav

      # Model A
      fmucon=bmb+zm%*%simod$condeltav

      # Model B
      fmul=bmb+zm%*%simod$ldeltav 
      fmunl=bmb+zm%*%simod$nldeltav
      fmulag=bmb+zm%*%simod$lagdeltav

      # Model C
      rmu=bmb+xm%*%simod$cons

      tind=1:lt
      sda=NULL

      # Run 500 simulations for each cluster i
      for(i in 1:(m*ns)){

        indi=(i-1)*lt+tind 
        zmi=zm[indi,] 

        ## Cluster mean vector for 3 models
        
        # Model A
        fmuconi=fmucon[indi] 
        
        # Model B
        fmuli=fmul[indi]
        fmunli=fmunl[indi]
        fmulagi=fmulag[indi]
        
        # Model C
        rmui=rmu[indi]

        cvz=0
        for(kk in 1:m){
          indkk=((kk-1)*ns+1):(kk*ns) 
          zmikk=zmi[,indkk] 
          cvz=cvz+zmikk%*%mvrnorm(n = 1, mu = rep(0,ns), Sigma = diag(trsatv[kk],ns))
        }

        ## Sample outcomes for 3 models
        
        # Model A
        sfycon=t(mvrnorm(n*nsim,fmuconi,fcovy)) 
        fycon=c(sfycon) 

        # Model B
        sfyl=t(mvrnorm(n*nsim,fmuli,fcovy)) 
        fyl=c(sfyl)

        sfynl=t(mvrnorm(n*nsim,fmunli,fcovy)) 
        fynl=c(sfynl)

        sfylag=t(mvrnorm(n*nsim,fmulagi,fcovy)) 
        fylag=c(sfylag)

        # Model C
        sry=t(mvrnorm(n*nsim,rmui+cvz,fcovy)) 
        ry=c(sry)

        sdai=data.frame(fycon,fyl,fynl,fylag,ry,kronecker(onev,bm[indi,]),
                        kronecker(onev,zmi),kronecker(onev,xm[indi,]),
                        cid=rep(i,n*nsim*lt),time=rep(tind,n*nsim),
                        ind=rep(rep(1:n,each=lt),nsim),
                        simid=rep(1:nsim,each=n*lt))

        sda=rbind(sda,sdai)
    }

  names(sda)=c("fycon","fyl","fynl","fylag",
               "ry",paste0("b",1:lt),paste0("z",1:(m*ns)),paste0("x",1:m),
               "cid","time","id","simid")

  save(sda,file=paste0("concurrent_sda_time",lt,"_n",n,"_bd",bd,".RData"))
  
  }
 }
}
