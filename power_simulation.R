source("helper_functions.R")

# Cluster TASK_ID
touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

pb=c(0.025,0.975) 
sigma2_a = 0.05
sigma2_e = 0.95
nsim=500
nbp=500
m=2
times = 5
lts=length(times)
bbeta=0.1
dbeta=0.5
trsatv=rep(0.2,2)

yts=c("fyl") 

res = NULL

for(design in c("concurrent", "factorial")){
 
  for (bd in c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)){

   for (n in c(30, 100, 500)) {
    
    for(ts in 1){
  
      lt=times[ts]
      ns=lt-1
  
      load(paste0(design,"_sda_time",lt,"_n",n,"_bd",bd,".RData"))
  
      for(ss in 1:length(yts)){
    
        sdas=subset(sda,simid==touse)
        yt=yts[ss]
  
        resB=fit_B(sdas,yt,lt,m)
   
        ## Bootstrap 
    
        bs_B = NULL
  
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
    
          bs_B = cbind(bs_B, fit_B(bda,yt,lt,m))
        }
    

        res_temp =  c(n, bd, resB,
                      apply(bs_B,1,quantile,prob=pb)[,1], 
                      apply(bs_B,1,quantile,prob=pb)[,2],
                      apply(bs_B,1,quantile,prob=pb_corr)[,1], 
                      apply(bs_B,1,quantile,prob=pb_corr)[,2],
                      'B', yt, lt, touse)
    
        res = rbind(res, res_temp)
        
     }
    }
   }
  }
  
  res$design = design
}

res = as.data.frame(res)

colnames(res) = c('n','bd', 'X1_est','X2_est', 
                  'X1_LB', 'X1_UB', 'X2_LB', 'X2_UB',
                  'X1_LB_adj', 'X1_UB_adj', 'X2_LB_adj', 'X2_UB_adj',
                  'Model_fit', 'Outcome', 'Period', "simID", "Design")
