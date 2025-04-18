source("simfun.R")

load("sobj.RData")
attach(sobj)

n = 30
es = c(2,4)

for (bd in c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)) {
  
  dd = bd + 0.5
  
  for (n in c(30, 100, 500)) {
    
    onev=rep(1,n*nsim)
    
    for(ts in 1:lts){

      lt=times[ts]
      ns=lt-1

      bm=bmf(lt,m)

      fo=fmatrix(lt,es[ts])
      zm=fo$zm
      nz=ncol(zm)
      xm=fo$xm
      nx=ncol(xm)
      
      fcovy=sigma2_a*matrix(1,lt,lt)+diag(sigma2_e,lt)

      betav=seq(bbeta,dbeta,length=lt)

      simod=model(lt,m,bd,dd)

      bmb=bm%*%betav

      fmucon=bmb+zm%*%simod$condeltav
      fmul=bmb+zm%*%simod$ldeltav
      fmunl=bmb+zm%*%simod$nldeltav
      fmulag=bmb+zm%*%simod$lagdeltav
      fmulag_my=bmb+zm%*%simod$lagdeltav_my

      rmu=bmb+xm%*%simod$cons
      
      tind=1:lt

      sda=NULL

      for(i in 1:(m*ns)){

        indi=(i-1)*lt+tind
        zmi=zm[indi,]

        fmuconi=fmucon[indi]
        fmuli=fmul[indi]
        fmunli=fmunl[indi]
        fmulagi=fmulag[indi]
        fmulagi_my=fmulag_my[indi]
        rmui=rmu[indi]

        cvz=0
        
        for(kk in 1:m){
          indkk=((kk-1)*ns+1):(kk*ns)
          zmikk=zmi[,indkk]
          cvz=cvz+zmikk%*%mvrnorm(1,rep(0,ns),diag(trsatv[kk],ns))
        }

        sfycon=t(mvrnorm(n*nsim,fmuconi,fcovy)) 
        fycon=c(sfycon)
        sfyl=t(mvrnorm(n*nsim,fmuli,fcovy)) 
        fyl=c(sfyl)
        sfynl=t(mvrnorm(n*nsim,fmunli,fcovy)) 
        fynl=c(sfynl)
        sfylag=t(mvrnorm(n*nsim,fmulagi,fcovy)) 
        fylag=c(sfylag)

        sfylag_my=t(mvrnorm(n*nsim,fmulagi_my,fcovy))
        fylag_my=c(sfylag_my)

        sry=t(mvrnorm(n*nsim,rmui+cvz,fcovy)) 
        ry=c(sry)

        sdai=data.frame(fycon,fyl,fynl,fylag,fylag_my,
                        ry,kronecker(onev,bm[indi,]),
                        kronecker(onev,zmi),kronecker(onev,xm[indi,]),cid=rep(i,n*nsim*lt),
                        time=rep(tind,n*nsim),ind=rep(rep(1:n,each=lt),nsim),simid=rep(1:nsim,each=n*lt))
        
        sda=rbind(sda,sdai)
        
      }

      names(sda)=c("fycon","fyl","fynl","fylag","fylag_my","ry",paste0("b",1:lt),
               paste0("z",1:(m*ns)),paste0("x",1:m),"cid","time","id","simid")

      save(sda,file=paste0("factorial_sda_time",lt,"_n",n,"_bd",bd,".RData"))

    }
  }
}

