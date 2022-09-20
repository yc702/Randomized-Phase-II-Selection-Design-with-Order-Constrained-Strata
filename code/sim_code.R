library(ggplot2)
library(survminer)
library(survival)
library(dplyr)
library(Hmisc)
library(foreach)
library(doParallel)
set.seed(123678)
options(digits=10)

## new simulation result 1
n <- c(20,30,40,50,60,70)
p1 <- 0.35
strata_diff <- 0.05
D1 <- 0.15
D2 <- 0.15
d1 <- c(0.02,0.05)
d2 <- c(0.02,0.05)
prop <- 0.4
rho = c(0,0.5)
output_old <- NULL
output_new <- NULL

source("bin_helper.R")

for (i in 1:length(n)){
  for (j in 1:length(d1)){
    result <- pickwin_original(n=n[i], p1=p1, strata_diff=strata_diff,
                               D1=D1, D2=D1,d1=d1[j], d2=d2[j],prop.strat=prop)
    output_old <- rbind(output_old,c(result$pcorr,result$pamb,n[i],p1,strata_diff,D1,D2,d1[j],d2[j],prop))
  }
  
}

for (i in 1:length(n)){
  for (j in 1:length(d1)){
    result <- pickwin_strat_sargent(n=n[i], p1=p1, strata_diff=strata_diff,  
                                    D1=D1, D2=D2,d1=d1[j], d2=d2[j],prop.strat=prop)
    output_new <- rbind(output_new,c(result$pcorr,result$pamb,n[i],p1,strata_diff,D1,D2,d1[j],d2[j],prop))
  }
  
}
output_old <- data.frame(output_old)
output_new <- data.frame(output_new)

colnames(output_old) <- c("Pcorr","Pamb","n","pi_b","strata_diff",
                          "D1","D2","d1","d2","prop")
colnames(output_new) <- c("Pcorr","Pamb","n","pi_b","strata_diff",
                          "D1","D2","d1","d2","prop")
output_old$Lambda_0 <- output_old$Pcorr
output_old$Lambda_05 <- output_old$Pcorr + 0.5*output_old$Pamb
output_new$Lambda_0 <- output_new$Pcorr
output_new$Lambda_05 <- output_new$Pcorr + 0.5*output_new$Pamb

plot_data <- data.frame(n = c(output_old$n,output_old$n,output_new$n,output_new$n),
                        Lambda = c(output_old$Lambda_0,output_old$Lambda_05,output_new$Lambda_0,output_new$Lambda_05),
                        rho = c(rep(rho,each=length(n)*2),rep(rho,each=length(n)*2)),
                        d1 = c(output_old$d1,output_old$d1,output_new$d1,output_new$d1),
                        d2 = c(output_old$d2,output_old$d2,output_new$d2,output_new$d2),
                        Group = c(rep("Original",4*length(n)),rep("Proposed",4*length(n))))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d1 <- factor(plot_data$d1)
plot_data$d2 <- factor(plot_data$d2)

ggplot(plot_data, aes(x=n, y=Lambda, group=interaction(Group,d1),color = Group,linetype=d1))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()

save(output_old,output_new,file = "sim1.RData")


## new simulation result 2
n <- 60
rho <- c(0,0.5)
p1 <- seq(0.25,0.55,0.05)
strata_diff <- 0.05
D1 <- 0.15
D2 <- 0.15
d1 <- 0.05
d2 <- 0.05
prop <- 0.4
output_old2 <- NULL
output_new2 <- NULL

for (i in 1:length(p1)){
  result <- pickwin_original(n=n, p1=p1[i], strata_diff=strata_diff,  
                                  D1=D1, D2=D2,d1=d1, d2=d2,prop.strat=prop)
  output_old2 <- rbind(output_old2,c(result$pcorr,result$pamb,n,p1[i],strata_diff,
                                     D1,D2,d1,d2,prop))
  
}

for (i in 1:length(p1)){
    result <- pickwin_strat_sargent(n=n, p1=p1[i], strata_diff=strata_diff,  
                                   D1=D1, D2=D2,d1=d1, d2=d2,prop.strat=prop)
    output_new2 <- rbind(output_new2,c(result$pcorr,result$pamb,n,p1[i],strata_diff,
                                       D1,D2,d1,d2,prop))
  
}

output_old2 <- data.frame(output_old2)
output_new2 <- data.frame(output_new2)

colnames(output_old2) <- c("Pcorr","Pamb","n","pi_b","strata_diff",
                           "D1","D2","d1","d2","prop")
colnames(output_new2) <- c("Pcorr","Pamb","n","pi_b","strata_diff",
                           "D1","D2","d1","d2","prop")
output_old2$Lambda_0 <- output_old2$Pcorr 
output_old2$Lambda_05 <- output_old2$Pcorr + 0.5*output_old2$Pamb
output_new2$Lambda_0 <- output_new2$Pcorr
output_new2$Lambda_05 <- output_new2$Pcorr + 0.5*output_new2$Pamb

plot_data <- data.frame(pi_b=c(output_old2$pi_b,output_old2$pi_b,
                               output_new2$pi_b,output_new2$pi_b),
                        Lambda = c(output_old2$Lambda_0,output_old2$Lambda_05,
                                   output_new2$Lambda_0,output_new2$Lambda_05),
                        rho = c(rep(rho,each=length(p1)),
                                rep(rho,each=length(p1))),
                        Group = c(rep("Original",2*length(p1)),
                                  rep("Proposed",2*length(p1))))

plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
ggplot(plot_data, aes(x=pi_b, y=Lambda, group=Group,color = Group))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous(expression(pi[b1]))+
  geom_line()+ylab(expression(paste(lambda)))+
  geom_point()

save(output_old2,output_new2,file = "sim2.RData")



### simulation for survival
## new simulation result 1

maxn <- c(20,30,40,50,60,70)
event_rate_A2 <- 0.03
strata_diff1 <- 0.03
strata_diff2 <- 0.03
trt_diff1 <- 0.08
trt_diff2 <- 0.08
d_diff1 <- c(0.02,0.05)
d_diff2 <- c(0.02,0.05)
rho = c(0,0.5)
prop=0.5
contr_result <- NULL
km_result <- NULL

cl <- makeCluster(6)
registerDoParallel(cl)
contr_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                        trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(maxn),.combine = 'rbind')%:%
    foreach(j=1:length(d_diff1),.combine = 'rbind') %dopar% {
      source("surv_helper.R")
     
      result <- sim_contr_fun(maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
                              trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      pcorr <- result[1]
      pamb <- result[2]
      
      c(pcorr,pamb,maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
        trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}

contr_result <- contr_combo (maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                             trt_diff1,trt_diff2,d_diff1,d_diff2)



cl <- makeCluster(6)
registerDoParallel(cl)
km_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                     trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(maxn),.combine = 'rbind')%:%
    foreach(j=1:length(d_diff1),.combine = 'rbind') %dopar% {
      source("surv_helper.R")
      result <- sim_km_fun(maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
                           trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      pcorr <- result[1]
      pamb <- result[2]
      
      c(pcorr,pamb,maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
        trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}

km_result <- km_combo(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                      trt_diff1,trt_diff2,d_diff1,d_diff2)




sim_data <- rbind(contr_result,km_result)
sim_data <- sim_data %>% data.frame() %>%
  mutate(Group=c(rep("Contr",nrow(km_result)),rep("KM",nrow(km_result))))
names(sim_data) <- c("Pcorr","Pamb","n","prop","event rate A2","strata diff1", "strata diff2",
                     "trt diff1","trt diff2","d diff1","d diff2","Group")
sim_data$Pcorr <- sim_data$Pcorr/5000
sim_data$Pamb <- sim_data$Pamb/5000

sim_data$Lambda_0 <- sim_data$Pcorr 
sim_data$Lambda_05 <- sim_data$Pcorr + 0.5*sim_data$Pamb

plot_data <- data.frame(n = c(sim_data$n,sim_data$n),
                        Lambda = c(sim_data$Lambda_0,sim_data$Lambda_05),
                        rho = c(rep(rho,each=length(sim_data$n))),
                        d1 = c(sim_data$`d diff1`,sim_data$`d diff1`),
                        d2 = c(sim_data$`d diff2`,sim_data$`d diff2`),
                        Group = c(ifelse(sim_data$Group=="Contr","Proposed","Original"),
                                  ifelse(sim_data$Group=="Contr","Proposed","Original")))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d1 <- factor(plot_data$d1)
plot_data$d2 <- factor(plot_data$d2)

ggplot(plot_data, aes(x=n, y=Lambda, group=interaction(Group,d1),color = Group,linetype=d1))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()

save(sim_data,file="sim_surv1.RData")


## new simulation result 2

maxn <- c(50)
prop=0.5
event_rate_A2 <- seq(0.03,0.13,0.02)
strata_diff1 <- 0.03
strata_diff2 <- 0.03
trt_diff1 <- 0.08
trt_diff2 <- 0.08
d_diff1 <- 0.05
d_diff2 <- 0.05
rho = c(0,0.5)
output_old <- NULL
output_new <- NULL
## x=3

cl <- makeCluster(6)
registerDoParallel(cl)
contr_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                        trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(event_rate_A2),.combine = 'rbind') %dopar% {
      source("surv_helper.R")
      result <- sim_contr_fun(maxn,prop,event_rate_A2[i],strata_diff1,strata_diff2,
                              trt_diff1,trt_diff2,d_diff1,d_diff2)
      pcorr <- result[1]
      pamb <- result[2]
      
      c(pcorr,pamb,maxn,prop,event_rate_A2[i],strata_diff1,strata_diff2,
        trt_diff1,trt_diff2,d_diff1,d_diff2)
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}

contr_result <- contr_combo (maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                             trt_diff1,trt_diff2,d_diff1,d_diff2)



cl <- makeCluster(6)
registerDoParallel(cl)
km_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                     trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(event_rate_A2),.combine = 'rbind') %dopar% {
    source("surv_helper.R")
    result <- sim_km_fun(maxn,prop,event_rate_A2[i],strata_diff1,strata_diff2,
                         trt_diff1,trt_diff2,d_diff1,d_diff2)
    pcorr <- result[1]
    pamb <- result[2]
    
    c(pcorr,pamb,maxn,prop,event_rate_A2[i],strata_diff1,strata_diff2,
      trt_diff1,trt_diff2,d_diff1,d_diff2)
    
  }
  return(start.p)
  on.exit(stopCluster(cl))
}

km_result <- km_combo(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                      trt_diff1,trt_diff2,d_diff1,d_diff2)



sim_data <- rbind(contr_result,km_result)
sim_data <- sim_data %>% data.frame() %>%
  mutate(Group=c(rep("Contr",nrow(km_result)),rep("KM",nrow(km_result))))
names(sim_data) <- c("Pcorr","Pamb","n","prop","event rate A2","strata diff1", "strata diff2",
                     "trt diff1","trt diff2","d diff2","d diff2","Group")
sim_data$Pcorr <- sim_data$Pcorr/5000
sim_data$Pamb <- sim_data$Pamb/5000


sim_data$Lambda_0 <- sim_data$Pcorr 
sim_data$Lambda_05 <- sim_data$Pcorr + 0.5*sim_data$Pamb

plot_data <- data.frame(h_a = c(sim_data$'event rate A2',sim_data$'event rate A2'),
                        Lambda = c(sim_data$Lambda_0,sim_data$Lambda_05),
                        rho = c(rep(rho,each=length(sim_data$'event rate A2'))),
                        Group = c(ifelse(sim_data$Group=="Contr","Proposed","Original"),
                                  ifelse(sim_data$Group=="Contr","Proposed","Original")))

plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
ggplot(plot_data, aes(x=h_a, y=Lambda, group=Group,color = Group))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+ylab(expression(paste(lambda)))+
  scale_x_continuous(expression(h[a2]))+
  geom_line()+
  geom_point()

save(sim_data,file = "sim_surv2.RData")



## case study:
## binary case
## new simulation result 1
n <- c(20,25,30,35,40,45,50)
p1 <- 0.4
strata_diff <- 0.1
D1 <- 0.2
D2 <- 0.2
d1 <- 0.05
d2 <- 0.05
prop <- 0.7
rho = c(0,0.5)
output_old <- NULL
output_new <- NULL
source("bin_helper.R")


for (i in 1:length(n)){
  for (j in 1:length(d1)){
    result <- pickwin_original(n=n[i], p1=p1, strata_diff=strata_diff,
                               D1=D1, D2=D1,d1=d1[j], d2=d2[j],prop.strat=prop)
    output_old <- rbind(output_old,c(result$pcorr,result$pamb,n[i],p1,strata_diff,D1,D2,d1[j],d2[j],prop))
  }
  
}

for (i in 1:length(n)){
  for (j in 1:length(d1)){
    result <- pickwin_strat_sargent(n=n[i], p1=p1, strata_diff=strata_diff,  
                                    D1=D1, D2=D2,d1=d1[j], d2=d2[j],prop.strat=prop)
    output_new <- rbind(output_new,c(result$pcorr,result$pamb,n[i],p1,strata_diff,D1,D2,d1[j],d2[j],prop))
  }
  
}
output_old <- data.frame(output_old)
output_new <- data.frame(output_new)

colnames(output_old) <- c("Pcorr","Pamb","n","pi_b","strata_diff",
                          "D1","D2","d1","d2","prop")
colnames(output_new) <- c("Pcorr","Pamb","n","pi_b","strata_diff",
                          "D1","D2","d1","d2","prop")
output_old$Lambda_0 <- output_old$Pcorr
output_old$Lambda_05 <- output_old$Pcorr + 0.5*output_old$Pamb
output_new$Lambda_0 <- output_new$Pcorr
output_new$Lambda_05 <- output_new$Pcorr + 0.5*output_new$Pamb

plot_data <- data.frame(n = c(output_old$n,output_old$n,output_new$n,output_new$n),
                        Lambda = c(output_old$Lambda_0,output_old$Lambda_05,output_new$Lambda_0,output_new$Lambda_05),
                        rho = c(rep(rho,each=length(n)),rep(rho,each=length(n))),
                        d1 = c(output_old$d1,output_old$d1,output_new$d1,output_new$d1),
                        d2 = c(output_old$d2,output_old$d2,output_new$d2,output_new$d2),
                        Group = c(rep("Original",2*length(n)),rep("Proposed",2*length(n))))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d1 <- factor(plot_data$d1)
plot_data$d2 <- factor(plot_data$d2)

ggplot(plot_data, aes(x=n, y=Lambda, group=Group,color = Group))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ 
  geom_hline(yintercept=0.8, linetype="dashed", 
                           color = "grey", size=0.5)
save(output_old2,output_new2,file="case_bin.RData")

## survival case
## new simulation result 1
maxn <- c(20,25,30,35,40,45)
event_rate_A2 <- -log(0.85)/2
strata_diff1 <- log(0.85)/2-log(0.75)/2
strata_diff2 <- log(0.85)/2-log(0.75)/2
trt_diff1 <- log(0.75)/2-log(0.6)/2
trt_diff2 <- log(0.85)/2-log(0.7)/2
d_diff1 <- 0.02
d_diff2 <- 0.02
rho = c(0,0.5)
prop=0.3
contr_result <- NULL
km_result <- NULL
## surv_helper_case.R has the same function as surv_helper.R
## except changing the accrual setting
source("surv_helper_case.R")


cl <- makeCluster(6)
registerDoParallel(cl)
contr_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                        trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(maxn),.combine = 'rbind')%:%
    foreach(j=1:length(d_diff1),.combine = 'rbind') %dopar% {
      source("surv_helper_case.R")
      
      result <- sim_contr_fun(maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
                              trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      pcorr <- result[1]
      pamb <- result[2]
      
      c(pcorr,pamb,maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
        trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}

contr_result <- contr_combo (maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                             trt_diff1,trt_diff2,d_diff1,d_diff2)



cl <- makeCluster(6)
registerDoParallel(cl)
km_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                     trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(maxn),.combine = 'rbind')%:%
    foreach(j=1:length(d_diff1),.combine = 'rbind') %dopar% {
      source("surv_helper_case.R")
      result <- sim_km_fun(maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
                           trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      pcorr <- result[1]
      pamb <- result[2]
      
      c(pcorr,pamb,maxn[i],prop,event_rate_A2,strata_diff1,strata_diff2,
        trt_diff1,trt_diff2,d_diff1[j],d_diff2[j])
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}

km_result <- km_combo(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                      trt_diff1,trt_diff2,d_diff1,d_diff2)




sim_data <- rbind(contr_result,km_result)
sim_data <- sim_data %>% data.frame() %>%
  mutate(Group=c(rep("Contr",nrow(km_result)),rep("KM",nrow(km_result))))
names(sim_data) <- c("Pcorr","Pamb","n","prop","event rate A1","strata diff1", "strata diff2",
                     "trt diff1","trt diff2","d diff1","d diff2","Group")
sim_data$Pcorr <- sim_data$Pcorr/5000
sim_data$Pamb <- sim_data$Pamb/5000

sim_data$Lambda_0 <- sim_data$Pcorr 
sim_data$Lambda_05 <- sim_data$Pcorr + 0.5*sim_data$Pamb

plot_data <- data.frame(n = c(sim_data$n,sim_data$n),
                        Lambda = c(sim_data$Lambda_0,sim_data$Lambda_05),
                        rho = c(rep(rho,each=length(sim_data$n))),
                        d1 = c(sim_data$`d diff1`,sim_data$`d diff1`),
                        d2 = c(sim_data$`d diff2`,sim_data$`d diff2`),
                        Group = c(ifelse(sim_data$Group=="Contr","Proposed","Original"),
                                  ifelse(sim_data$Group=="Contr","Proposed","Original")))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d1 <- factor(plot_data$d1)
plot_data$d2 <- factor(plot_data$d2)

ggplot(plot_data, aes(x=n, y=Lambda, group=Group,color = Group))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ 
  geom_hline(yintercept=0.8, linetype="dashed", 
             color = "grey", size=0.5)
save(sim_data,file = "case_surv.RData")
