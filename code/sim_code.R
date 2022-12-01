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
strata_diff <- 0.1
D1 <- 0.15
D2 <- 0.15
d1 <- c(0.03,0.07)
d2 <- c(0.03,0.07)
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
  geom_point()+ theme(text = element_text(size = 15)) 

save(output_old,output_new,file = "sim1.RData")




### simulation for survival
## new simulation result 1

maxn <- c(20,30,40,50,60,70)
event_rate_A2 <- -log(0.85)/3
strata_diff1 <- log(0.85)/3-log(0.75)/3
strata_diff2 <- log(0.85)/3-log(0.75)/3
trt_diff1 <- log(0.75)/3-log(0.55)/3
trt_diff2 <- log(0.85)/3-log(0.65)/3
d_diff1 <- c(0.03,0.07)
d_diff2 <- c(0.03,0.07)
rho = c(0,0.5)
prop=0.4
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
      ## set up for accrual
      x=3
      arrival_rate <- 1
      FUP = 3
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
      source("helper.R")
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
# sim_data$Lambda <- (sim_data$Pamb+sim_data$Pcorr)/5000
sim_data$Pcorr <- sim_data$Pcorr/5000
sim_data$Pamb <- sim_data$Pamb/5000

sim_data$Lambda_0 <- sim_data$Pcorr 
sim_data$Lambda_05 <- sim_data$Pcorr + 0.5*sim_data$Pamb

plot_data <- data.frame(n = c(sim_data$n,sim_data$n),
                        Lambda = c(sim_data$Lambda_0,sim_data$Lambda_05),
                        rho = c(rep(rho,each=length(sim_data$n))),
                        d1 = c(sim_data$`d diff1`,sim_data$`d diff1`),
                        d2 = c(sim_data$`d diff2`,sim_data$`d diff2`),
                        Group = c(ifelse(sim_data$Group=="Contr","With Constraints","Without Constraints"),
                                  ifelse(sim_data$Group=="Contr","With Constraints","Without Constraints")))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d1 <- factor(plot_data$d1)
plot_data$d2 <- factor(plot_data$d2)


ggplot(plot_data, aes(x=n, y=Lambda, group=interaction(Group,d1),color = Group,linetype=d1))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  scale_y_continuous(breaks = c(0.75,0.85,0.95))+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ theme(text = element_text(size = 15)) 

save(sim_data,file="sim_surv1.RData")



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
source("bin_function.R")


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
  geom_point()+ theme(text = element_text(size = 15)) + 
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
## changing the accrual setting


cl <- makeCluster(6)
registerDoParallel(cl)
contr_combo <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                        trt_diff1,trt_diff2,d_diff1,d_diff2){
  result <- NULL
  start.p <- foreach(i=1:length(maxn),.combine = 'rbind')%:%
    foreach(j=1:length(d_diff1),.combine = 'rbind') %dopar% {
      source("surv_helper.R")
      x=2
      arrival_rate <- 8
      FUP <- 2
      
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

ggplot(plot_data, aes(x=n, y=Lambda, group=Group,color = Group))+theme_classic() + 
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ theme(text = element_text(size = 15)) + 
  geom_hline(yintercept=0.8, linetype="dashed", 
             color = "grey", size=0.5)
save(sim_data,file = "case_surv.RData")



RNGkind("L'Ecuyer-CMRG") # the random generator to set parallel seed
set.seed(1236789)
options(digits=10)

### Simulation table 
## Binary
p1_seq <- c(0.05,0.2,0.35,0.5,0.65)
prop_seq <- c(0.4)
strata_diff <- 0.1
D1 <- 0.2
D2 <- 0.2
d1 <- 0.05
d2 <- 0.05
n_old1 <- NULL
n_new1 <- NULL

## explore different p1 seq
cl <- makeCluster(6)
registerDoParallel(cl)
contr_combo <- function(n, p1_seq, strata_diff,
                        D1, D2,d1, d2,prop_seq,sigma,rho,nstart){
  start.p <- foreach(i=1:length(p1_seq),.combine = 'rbind')%:%
    foreach(j=1:length(prop_seq),.combine = 'rbind') %dopar% {
      source("bin_helper.R")
      
      n <- nstart
      lambda <- 0
      while (lambda<sigma){
        n <- n+1
        result <- pickwin_strat_sargent(n=n, p1=p1_seq[i], strata_diff=strata_diff,
                                        D1=D1, D2=D1,d1=d1, d2=d2,prop.strat=prop_seq[j])
        lambda <- result$pcorr+result$pamb*rho
      }
      
      c(result$pcorr,result$pamb,n,p1_seq[i]+D2,prop_seq[j])
      
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}

n_new1 <- contr_combo (n, p1_seq, strata_diff,
                       D1, D2,d1, d2,prop_seq,sigma=0.8,rho = 0.5,nstart=14)


original_combo <- function(n, p1_seq, strata_diff,
                           D1, D2,d1, d2,prop_seq,sigma,rho,nstart){
  start.p <- foreach(i=1:length(p1_seq),.combine = 'rbind')%:%
    foreach(j=1:length(prop_seq),.combine = 'rbind') %dopar% {
      source("bin_helper.R")
      
      n <- nstart+i
      lambda <- 0
      while (lambda<sigma){
        n <- n+1
        result <- pickwin_original(n=n, p1=p1_seq[i], strata_diff=strata_diff,
                                   D1=D1, D2=D1,d1=d1, d2=d2,prop.strat=prop_seq[j])
        lambda <- result$pcorr+result$pamb*rho
      }
      
      c(result$pcorr,result$pamb,n,p1_seq[i]+D2,prop_seq[j])
      
      
    }
  return(start.p)
  on.exit(stopCluster(cl))
}


n_old1 <- original_combo (n, p1_seq, strata_diff,
                          D1, D2,d1, d2,prop_seq,sigma=0.8,rho = 0.5,nstart=15)



n_new1 <- data.frame(n_new1)
n_old1 <- data.frame(n_old1)
col_name_bin <- c("Pcorr","Pamb","Sample Size","Response Rate","Proportion")
colnames(n_new1) <- col_name_bin
colnames(n_old1) <- col_name_bin

n_new1$group <- "amb0.8_contr"
n_old1$group <- "amb0.8_orig"

combine_bin <- rbind(n_new1,n_old1)
save(combine_bin,file="case_bin_table_jsm.RData")


## make plots
combine_bin <- rbind(n_new1,n_old1)
ggplot(combine_bin, aes(x=`Response Rate`, y=`Sample Size`, group=group,color = group))+
  geom_line()+
  geom_point()+ theme(text = element_text(size = 15)) +theme_classic() + 
  facet_wrap(~Proportion,
             labeller = label_parsed)+ylab("Sample Size")



## survival outcome
## FUP = 4,5,6
sigma=0.8;rho = 0.5;nstart=30
FUP = 4
prop_seq <- 0.4
event_rate_A2_seq <- c(-log(0.95)/6,-log(0.8)/6,-log(0.65)/6,
                       -log(0.5)/6,-log(0.35)/6)
strata_diff_seq <- c(log(0.95)/6-log(0.85)/6,
                     log(0.8)/6-log(0.7)/6,
                     log(0.65)/6-log(0.55)/6,
                     log(0.5)/6-log(0.4)/6,
                     log(0.35)/6-log(0.25)/6)


trt_diff1_seq <- c(log(0.85)/6-log(0.65)/6,log(0.7)/6-log(0.5)/6,
                   log(0.55)/6-log(0.35)/6,
                   log(0.4)/6-log(0.2)/6,log(0.25)/6-log(0.05)/6)

trt_diff2_seq <- c(log(0.95)/6-log(0.75)/6,log(0.8)/6-log(0.6)/6,
                   log(0.65)/6-log(0.45)/6,log(0.5)/6-log(0.3)/6,
                   log(0.35)/6-log(0.15)/6)

d_diff1 <- 0.05
d_diff2 <- 0.05
n_old1 <- NULL
n_new1 <- NULL
nstart = c(15,27,27,27,16)
source("sim_table_helper.R")

## constrained survival
output <- NULL
result <- NULL
for(i in 1:length(event_rate_A2_seq)){
  for (j in 1:length(prop_seq)){
    maxn <- nstart[i]
    lambda <- 0
    
    while (lambda<sigma){
      maxn <- maxn+1
      corr <- 0
      err <- 0
      n1<- ceiling(prop_seq[j]*maxn)
      n2<- maxn-n1  
      event_rate_A2 <- event_rate_A2_seq[i]
      strata_diff1 <- strata_diff2 <- strata_diff_seq[i]
      trt_diff1 <- trt_diff1_seq[i]
      trt_diff2 <- trt_diff2_seq[i]
      
      result <- mclapply(1:8000, sim_contr_fun,mc.cores=6) %>% bind_rows()
      amb = (8000-sum(result$corr)-sum(result$err))/8000
      
      lambda <- mean(result$corr)+amb*rho
    }
    
    output <- rbind(output,c(mean(result$corr),amb,maxn,lambda,
                             event_rate_A2_seq[i],prop_seq[j],
                             strata_diff1,trt_diff1,FUP))
    
    
  }
}
output <- data.frame(output)


colnames(output) = c("Pcorr","Pamb","Sample Size","Lambda","Hazard2",
                     "Proportion","Strata Diff","Trt Diff","FUP")


n_new1 <- output


## km survival
output <- NULL
for(i in 1:length(event_rate_A2_seq)){
  for (j in 1:length(prop_seq)){
    maxn <- nstart[i]
    lambda <- 0
    
    while (lambda<sigma){
      maxn <- maxn+1
      corr <- 0
      err <- 0
      n1<- ceiling(prop_seq[j]*maxn)
      n2<- maxn-n1  
      event_rate_A2 <- event_rate_A2_seq[i]
      strata_diff1 <- strata_diff2 <- strata_diff_seq[i]
      trt_diff1 <- trt_diff1_seq[i]
      trt_diff2 <- trt_diff2_seq[i]
      
      result <- mclapply(1:5000, sim_km_fun,mc.cores=6) %>% 
        bind_rows()
      amb = (5000-sum(result$corr)-sum(result$err))/5000
      
      lambda <- mean(result$corr)+amb*rho
    }
    
    output <- rbind(output,c(mean(result$corr),amb,maxn,lambda,
                             event_rate_A2_seq[i],prop_seq[j],
                             strata_diff1,trt_diff1,FUP))
    
  }
}

output <- data.frame(output)
colnames(output) = c("Pcorr","Pamb","Sample Size","Lambda","Hazard2",
                     "Proportion","Strata Diff","Trt Diff","FUP")


n_old1 <- output
