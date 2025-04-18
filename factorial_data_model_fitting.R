source("simfun.R")

# Cluster TASK_ID
touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

n = 30
pb=c(0.025,0.975) 
yts=c("fycon","fyl","fynl","fylag","fylag_my") 

res = NULL

for (tt in c(1,3)) {
  
  load(paste0("sobj",tt,"_small_effect.RData"))
  attach(sobj)
  
  for(ts in 1:lts){
  
    lt=times[ts]
    ns=lt-1
  
    load(paste0("factorial_sda_time",lt,"_n",n,"_tt",tt,"_small_effect.RData"))
  
    for(ss in 1:length(yts)){
    
      sdas=subset(sda,simid==touse)
      yt=yts[ss]
    
      resA=fit_A(sdas,yt,lt,m) # constant effects
      resB=fit_B(sdas,yt,lt,m) # time-varying fixed effects
      resC=rfit3(sdas,yt,lt,m) # mean of time-varying random effects
    
    
      ## Bootstrap within cluster period
    
      bs_A = NULL
      bs_B = NULL
      bs_C = NULL
    
      for(bs in 1:(nbp)){
      
        bda=NULL
      
        for(ii in 1:(m*ns)){
          sdaii=subset(sdas,cid==ii)
          bi=NULL	
        
          for(jj in 1:lt){
            sdaij=subset(sdaii,time==jj)
            sind=sample(n,n,replace=T)
            bij=sdaij[sind,]
            bi=rbind(bi,bij) # sample n observations w. replacement for each cluster i at each time j 
          }
          
          bda=rbind(bda,bi)	
        }
      
      bs_A = cbind(bs_A, fit_A(bda,yt,lt,m))
      bs_B = cbind(bs_B, fit_B(bda,yt,lt,m))
      bs_C = cbind(bs_C, rfit3(bda,yt,lt,m))
     }
    
    # Bootstrap within cluster
    bs_A_c = NULL
    bs_B_c = NULL
    bs_C_c = NULL
    
    for(bs in 1:(nbp)){
      
      bda=NULL
      
      for(ii in 1:(m*ns)){
        sdaii=subset(sdas,cid==ii)
        n_c = nrow(sdaii)
        sind=sample(n_c,n_c,replace=T)
        bi=sdaii[sind,] # sample observations w. replacement for each cluster i 
        bda=rbind(bda,bi)	
      }
      
      bs_A_c = cbind(bs_A_c, fit_A(bda,yt,lt,m))
      bs_B_c = cbind(bs_B_c, fit_B(bda,yt,lt,m))
      bs_C_c = cbind(bs_C_c, rfit3(bda,yt,lt,m)) 
    }
    
    res_temp = rbind(c(resA,
                       apply(bs_A,1,quantile,prob=pb)[,1], apply(bs_A,1,quantile,prob=pb)[,2],
                       apply(bs_A, 1, sd),
                       apply(bs_A_c,1,quantile,prob=pb)[,1], apply(bs_A_c,1,quantile,prob=pb)[,2],
                       apply(bs_A_c, 1, sd),
                       n, 'A', yt, lt, touse,sigma2_a,sigma2_e), 
                     c(resB, 
                       apply(bs_B,1,quantile,prob=pb)[,1], apply(bs_B,1,quantile,prob=pb)[,2],
                       apply(bs_B, 1, sd),
                       apply(bs_B_c,1,quantile,prob=pb)[,1], apply(bs_B_c,1,quantile,prob=pb)[,2],
                       apply(bs_B_c, 1, sd),
                       n, 'B', yt, lt, touse,sigma2_a,sigma2_e), 
                     c(resC, 
                       apply(bs_C,1,quantile,prob=pb)[,1], apply(bs_C,1,quantile,prob=pb)[,2],
                       apply(bs_C, 1, sd),
                       apply(bs_C_c,1,quantile,prob=pb)[,1], apply(bs_C_c,1,quantile,prob=pb)[,2],
                       apply(bs_C_c, 1, sd),
                       n, 'C', yt, lt, touse,sigma2_a,sigma2_e))

      res = rbind(res, res_temp)
    }
  }
}

res = as.data.frame(res)

colnames(res) = c('X1_est','X2_est', 
                  'X1_LB_boots_clt_pd', 'X1_UB_boots_clt_pd', 'X2_LB_boots_clt_pd', 'X2_UB_boots_clt_pd',
                  "X1_SE_boots_clt_pd","X2_SE_boots_clt_pd",
                  'X1_LB_boots_clt', 'X1_UB_boots_clt', 'X2_LB_boots_clt', 'X2_UB_boots_clt',
                  "X1_SE_boots_clt","X2_SE_boots_clt",
                  "n", 'Model_fit', 'Outcome', 'Period', "simID","sigma2_a","sigma2_e")

detach(sobj)

