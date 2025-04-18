library(dplyr)
library(xtable)

source("helper_functions.R")

load("sobj.RData")
attach(sobj)

for(ts in 1:lts){
  time=times[ts]
  ns=time-1
  simod=model(time,m,bd,dd)
  print(simod)
  print(paste0("constant: ", mean(simod$condeltav[1:ns]), " ", sd(simod$condeltav[1:ns])," ",
               mean(simod$condeltav[(1+ns):(m*ns)])," ", sd(simod$condeltav[(1+ns):(m*ns)])))
  print(paste0("linear:", mean(simod$ldeltav[1:ns]), " ", " ", sd(simod$ldeltav[1:ns])," ",
               mean(simod$ldeltav[(1+ns):(m*ns)])," ", sd(simod$ldeltav[(1+ns):(m*ns)])))
  print(paste0("lagged:", mean(simod$lagdeltav[1:ns]), " ", sd(simod$lagdeltav[1:ns])," ",
               mean(simod$lagdeltav[(1+ns):(m*ns)])," ", sd(simod$lagdeltav[(1+ns):(m*ns)])))
  print(paste0("nonlinear:", mean(simod$nldeltav[1:ns]), " ",sd(simod$nldeltav[1:ns])," ",
               mean(simod$nldeltav[(1+ns):(m*ns)])," ", sd(simod$nldeltav[(1+ns):(m*ns)])))
}

DF <- NULL
for (f in list.files()) {
  dat <- data.frame(Filename=f, read.csv(f, header = TRUE))
  DF <- rbind(DF, dat)
}

DF = DF %>%
  mutate(x1_truth = ifelse(Period == 5, 0.95, 
                           ifelse(Outcome == "fyl", mean(simod$ldeltav[1:ns]), 0.97)),
         x2_truth = ifelse(Period == 5, 1.35, 
                           ifelse(Outcome == "fyl", mean(simod$ldeltav[(1+ns):(m*ns)]), 1.33)))

# small effect size
DF = DF %>%
  mutate(x1_truth = ifelse((Period == 5)&(Outcome != "fyl"), 0.1, 
                            ifelse((Period == 5)&(Outcome == "fyl"), 0.095,
                                  ifelse((Period == 11)&(Outcome != "fyl"), 0.1, 0.097))),
         x2_truth = ifelse((Period == 5)&(Outcome != "fyl"), 0.14, 
                           ifelse((Period == 5)&(Outcome == "fyl"), 0.135,
                                  ifelse((Period == 11)&(Outcome != "fyl"), 0.13, 0.133))))
# large effect size
DF = DF %>%
  mutate(x1_truth = ifelse((Period == 5)&(Outcome != "fyl"), 0.28, 
                           ifelse((Period == 5)&(Outcome == "fyl"), 0.285,
                                  ifelse((Period == 11)&(Outcome != "fyl"), 0.29, 0.29))),
         x2_truth = ifelse((Period == 5)&(Outcome != "fyl"), 0.40, 
                           ifelse((Period == 5)&(Outcome == "fyl"), 0.405,
                                  ifelse((Period == 11)&(Outcome != "fyl"), 0.4, 0.4))))

summary = DF %>% 
  filter(sigma2_e==2.85 , sigma2_a==0.15, Period==5)%>% 
  group_by(Outcome, Model_fit)%>%
  summarise(bias_x1 = mean(X1_est - x1_truth),
            SD_x1 = sd(X1_est), 
            coverage_x1_clt_pd = 100*mean((X1_LB_boots_clt_pd <= x1_truth) & (X1_UB_boots_clt_pd >= x1_truth)),
            ci_length_x1_clt_pd = mean(X1_UB_boots_clt_pd-X1_LB_boots_clt_pd),
            mean_SE_x1_clt_pd = mean(X1_SE_boots_clt_pd), 
            bias_x2 = mean(X2_est - x2_truth),
            SD_x2 = sd(X2_est), 
            coverage_x2_clt_pd = 100*mean((X2_LB_boots_clt_pd <= x2_truth) & (X2_UB_boots_clt_pd >= x2_truth)),
            ci_length_x2_clt_pd = mean(X2_UB_boots_clt_pd-X2_LB_boots_clt_pd),
            mean_SE_x2_clt_pd = mean(X2_SE_boots_clt_pd))#, 
       

print(xtable(summary, digits=c(0,0,0,3,2,1,2,2,3,2,1,2,2)), 
      include.rownames=FALSE)

       