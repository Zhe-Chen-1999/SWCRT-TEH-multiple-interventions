library(dplyr)
library(xtable)

source("helper_functions.R")

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
ns = times-1

for (bd in c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)){
  dd = bd+0.5
  simod=model(times[1],m,bd,dd)
  print(paste0("linear:", mean(simod$ldeltav[1:ns]), " ", 
             mean(simod$ldeltav[(1+ns):(m*ns)])))
}

res = NULL

for(design in c("concurrent", "factorial")){
    
    setwd(paste0(design,'/res'))
    
    DF <- NULL
    for (f in list.files()) {
      dat <- data.frame(Filename=f, read.csv(f, header = TRUE))
      DF <- rbind(DF, dat)
    }
    
    summary = DF %>% 
        group_by(n, bd)%>% 
        summarise( design = design,
                   power_x1 = mean(X1_LB > 0),
                   power_x2 = mean(X2_LB > 0),
                   power_x1_adj = mean(X1_LB_adj > 0),
                   power_x2_adj = mean(X2_LB_adj > 0))
   
    res = rbind(res, summary)
}

############# Power Plot ##############

res_wider <- res %>%
  pivot_wider(
    names_from = design,
    values_from = c(power_x1, power_x2, power_x1_adj, power_x2_adj))  

plot_data <- res_wider %>%
  pivot_longer(cols = starts_with("power_"),
               names_to = "power",
               values_to = "value") %>%
  mutate(design = ifelse(power %in% c('power_x1_concurrent', 'power_x2_concurrent'), 
                         'Concurrent', 'Factorial'),
         trt = ifelse(power %in% c('power_x1_concurrent', 'power_x1_factorial'), 
                      'Treatment 1', 'Treatment 2'))%>%
  rowwise() %>% 
  mutate(x1_ave = model(5,2,bd,bd+0.5)$cons[1],
         x2_ave = model(5,2,bd,bd+0.5)$cons[2]) %>% 
  mutate(ave = ifelse(trt ==  'Treatment 1', x1_ave, x2_ave),
         n = as.factor(n))

ggplot(plot_data, aes(x = ave, y = value, color = n,  linetype = design)) +
  geom_line(linewidth = 1) +
  facet_wrap(~trt, scales = "free_x", nrow = 1) +
  labs(
    x = "Average Treatment Effect",
    y = "Power",
    color = "Design",
    linetype = "Model Fit") +
  scale_color_viridis_d()+
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 18),
    legend.position = "bottom",
    legend.title = element_text(size = 16),  # Adjust legend title size
    legend.text = element_text(size = 14),   # Adjust legend text size
    axis.title.x = element_text(size = 18),  # Increase x-axis title size
    axis.title.y = element_text(size = 18),  # Increase y-axis title size
    axis.text.x = element_text(size = 14),   # Increase x-axis label size
    axis.text.y = element_text(size = 14)    # Increase y-axis label size
  )


