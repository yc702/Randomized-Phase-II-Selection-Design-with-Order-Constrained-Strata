library(ggplot2)
library(survminer)
library(survival)
library(dplyr)
library(Hmisc)
library(foreach)
library(doParallel)
RNGkind("L'Ecuyer-CMRG") # the random generator to set parallel seed
set.seed(1236789)
options(digits=10)

source("function_table.R")

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

## explore different proportion, p1 seq
cl <- makeCluster(6)
registerDoParallel(cl)
contr_combo <- function(n, p1_seq, strata_diff,
                        D1, D2,d1, d2,prop_seq,sigma,rho,nstart){
  start.p <- foreach(i=1:length(p1_seq),.combine = 'rbind')%:%
    foreach(j=1:length(prop_seq),.combine = 'rbind') %dopar% {
      source("function_table.R")
      
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
      source("function_table.R")
      
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
source("helper_table.R")


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

