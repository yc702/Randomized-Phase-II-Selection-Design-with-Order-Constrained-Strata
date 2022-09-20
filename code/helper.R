library(survminer)
library(survival)
library(dplyr)
library(Hmisc)
RNGkind( "L'Ecuyer-CMRG")
set.seed(123678)
x=3
arrival_rate <- 1
FUP <- 3


Mg <- function(x,event_time,event_ind){
  sum(event_time[!duplicated(event_time)]<=x & event_ind[!duplicated(event_time)]==1)
}

Ng <- function(x,event_time){
  sum(event_time>=x)
}


Kg <- function(qn,nrisk,nevent,event_time,event_ind){
  mgx <- Mg(x,event_time,event_ind)
  
  if (mgx ==0 & abs(qn)>0.00001){
    
    return(-Ng(x,event_time))
    
  } else {
    ## k=inf if q=0
    if (abs(qn)<0.00001){
      k_hat <- Inf
      return(k_hat)
      
    } else if (qn == -Inf) {
      
      k_hat <- nevent[mgx]-nrisk[mgx]
      return(k_hat)
    } else{
      
      fr <- function(k){
        left <- 0
        for (j in 1:mgx){
          left <- left+log(1-nevent[j]/(nrisk[j]+k))
        }
        return(abs(left-qn))
      }
      k_hat = optimize(fr, c(-50, 50))$minimum
      
      k_hat <- max(k_hat,-Ng(x,event_time))
      return(k_hat)
      }
    
  } 
}

llkhd <- function(qq,nrisk,nevent,event_time,event_ind){
  mgx <- Mg(x,event_time,event_ind)
  
  if (mgx!=0){
    event_llkhd <- 0
    k_hat <- Kg(qq,nrisk,nevent,event_time,event_ind)
    
    ## if equal, there will be 0 in denominators which will lead to -inf values 
    ## n-Ng=d
    if(k_hat==-Ng(x,event_time)){
      return(-1000)
    } else{
      for (i in 1:mgx){
        event_llkhd <- event_llkhd+(nrisk[i]-nevent[i])*
          log(nrisk[i]+k_hat-nevent[i])-
          nrisk[i]*log(nrisk[i]+k_hat)
      }
      return(event_llkhd)
    }} else{
      return(Ng(x,event_time)*qq)
    }
  
}




sim_surv <- function(nmax,arrival_rate,event_rate,FUP){
  
  wait.t = rexp(nmax,rate = arrival_rate)
  arrival.t = cumsum(wait.t)
  event.t = rexp(nmax,rate=event_rate)
  tobs = arrival.t[nmax]
  tobs = tobs + FUP
  
  ## tobs being the total follow up time
  n.fail = sum(arrival.t+ event.t <= tobs)
  t.event = rep(0,nmax)
  t.ind = rep(0,nmax)
  for(j in 1:length(t.event)) {
    ## count events in the observed follow up time
    t.event[j] = ifelse(arrival.t[j]+event.t[j]<=tobs,event.t[j],tobs-arrival.t[j])
    t.ind[j] = ifelse(arrival.t[j]+event.t[j]<=tobs,1,0)
  }
  
  return(cbind(time = t.event,ind = t.ind))
  
}

sim_contr_fun <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                          trt_diff1,trt_diff2,d_diff1,d_diff2){
  corr <- 0
  err <- 0
  n1<- ceiling(prop*maxn)
  n2<- maxn-n1  

  # ## Assume only have two groups:
  ## helper function
  for(sim in 1:5000){
    
    trtA_1 <- sim_surv(n1,arrival_rate,event_rate=event_rate_A2+strata_diff1,FUP)
    trtA_2 <- sim_surv(n2,arrival_rate,event_rate=event_rate_A2,FUP)
    
    trtB_1 <- sim_surv(n1,arrival_rate,event_rate=event_rate_A2+trt_diff1+strata_diff2,FUP)
    trtB_2 <- sim_surv(n2,arrival_rate,event_rate=event_rate_A2+trt_diff2,FUP)
    
    testA <- data.frame(rbind(trtA_1,trtA_2),group=c(rep("Strata 1",n1),rep("Strata 2",n2)))
    testB <- data.frame(rbind(trtB_1,trtB_2),group=c(rep("Strata 1",n1),rep("Strata 2",n2)))
    
    testA <- testA %>%  arrange(group,time)
    testB <- testB %>%  arrange(group,time)
    
    fit <- survfit(Surv(time, ind)~group,data=testA)
    
    test1 <- testA[testA$group=="Strata 1",]
    test2 <- testA[testA$group=="Strata 2",]
    
    nrisk_1 <- summary(fit)$n.risk[1:nrow(test1[test1$ind==1,])]
    nrisk_2 <- summary(fit)$n.risk[(nrow(test1[test1$ind==1,])+1):(nrow(test1[test1$ind==1,])+nrow(test2[test2$ind==1,]))]
    nevent_1 <- summary(fit)$n.event[1:nrow(test1[test1$ind==1,])]
    nevent_2 <- summary(fit)$n.event[(nrow(test1[test1$ind==1,])+1):(nrow(test1[test1$ind==1,])+nrow(test2[test2$ind==1,]))]
    
    event_time_1 <- test1$time
    event_time_2 <- test2$time
    
    event_ind_1 <- test1$ind
    event_ind_2 <- test2$ind
    fn <- function(q){
      q1 <- q[1]
      
      llkhd1 <- -llkhd(q1,nrisk_1,nevent_1,event_time_1,event_ind_1)
      
      q2 <- q[2]
      
      llkhd2 <- -llkhd(q2,nrisk_2,nevent_2,event_time_2,event_ind_2)
      llkhd1+llkhd2
      
    }
    
    
    grr <- function(q){
      
      c(Kg(q[1],nrisk_1,nevent_1,event_time_1,event_ind_1),Kg(q[2],nrisk_2,nevent_2,event_time_2,event_ind_2))
    }
    q_optim <- constrOptim(c(-0.2,-0.1), fn, grad=grr, 
                           ui = t(rbind(c(-1,-1,0),c(1,0,-1))), ci = c(0,0,0))$par
    
    S_A1 <- exp(q_optim[1])
    S_A2 <- exp(q_optim[2])
    
    test1 <- testB[testB$group=="Strata 1",]
    test2 <- testB[testB$group=="Strata 2",]
    fit <- survfit(Surv(time, ind)~group,data=testB)
    
    nrisk_1 <- summary(fit)$n.risk[1:nrow(test1[test1$ind==1,])]
    nrisk_2 <- summary(fit)$n.risk[(nrow(test1[test1$ind==1,])+1):(nrow(test1[test1$ind==1,])+nrow(test2[test2$ind==1,]))]
    nevent_1 <- summary(fit)$n.event[1:nrow(test1[test1$ind==1,])]
    nevent_2 <- summary(fit)$n.event[(nrow(test1[test1$ind==1,])+1):(nrow(test1[test1$ind==1,])+nrow(test2[test2$ind==1,]))]
    
    event_time_1 <- test1$time
    event_time_2 <- test2$time
    
    event_ind_1 <- test1$ind
    event_ind_2 <- test2$ind
    
    q_optim <- constrOptim(c(-0.2,-0.1), fn, grad=grr, 
                           ui = t(rbind(c(-1,-1,0),c(1,0,-1))), ci = c(0,0,0))$par
    S_B1 <- exp(q_optim[1])
    S_B2 <- exp(q_optim[2])
    
    if(S_A1>S_B1+d_diff1 & S_A2>S_B2+d_diff2){
      corr <- corr+1
    }
    if(S_A1<S_B1+d_diff1 & S_A2<S_B2+d_diff2){
      err <- err+1
    }
  }
  amb = 5000-corr-err
  return(c(corr,amb))
}




sim_km_fun <- function(maxn,prop,event_rate_A2,strata_diff1,strata_diff2,
                       trt_diff1,trt_diff2,d_diff1,d_diff2){
  corr <- 0
  err <- 0
  n1<- ceiling(prop*maxn)
  n2<- maxn-n1  
  
  
  for(sim in 1:5000){
    trtA_1 <- sim_surv(n1,arrival_rate,event_rate=event_rate_A2+strata_diff1,FUP)
    trtA_2 <- sim_surv(n2,arrival_rate,event_rate=event_rate_A2,FUP)
    
    trtB_1 <- sim_surv(n1,arrival_rate,event_rate=event_rate_A2+trt_diff1+strata_diff2,FUP)
    trtB_2 <- sim_surv(n2,arrival_rate,event_rate=event_rate_A2+trt_diff2,FUP)
    
    testA <- data.frame(rbind(trtA_1,trtA_2),group=c(rep("Strata 1",n1),rep("Strata 2",n2)))
    testB <- data.frame(rbind(trtB_1,trtB_2),group=c(rep("Strata 1",n1),rep("Strata 2",n2)))
    
    testA <- testA %>%  arrange(group,time)
    testB <- testB %>%  arrange(group,time)
    
    test1 <- testA[testA$group=="Strata 1",]
    test2 <- testA[testA$group=="Strata 2",]
    fit1 <- survfit(Surv(time, ind)~1,data=test1)
    fit2 <- survfit(Surv(time, ind)~1,data=test2)
    
    time_1 <- summary(fit1)$time
    time_2 <- summary(fit2)$time
    
    if(length(time_1)==0){
      S_A1=1
    } else{
      if (max(time_1)<x){
        S_A1=summary(fit1,t=max(time_1))$surv
      } else{
        S_A1 <- summary(fit1,t=x)$surv
      }
    }
     
  
    if(length(time_2)==0){
      S_A2 <- 1
    } else{
      if (max(time_2)<x){
        S_A2 <- summary(fit2,t=max(time_2))$surv
      } else{
        S_A2 <- summary(fit2,t=x)$surv
      }
    }
    
    test1 <- testB[testB$group=="Strata 1",]
    test2 <- testB[testB$group=="Strata 2",]
    fit1 <- survfit(Surv(time, ind)~1,data=test1)
    fit2 <- survfit(Surv(time, ind)~1,data=test2)
    
    time_1 <- summary(fit1)$time
    time_2 <- summary(fit2)$time
    
    if(length(time_1)==0){
      S_B1 <- 1
    } else{
      if (max(time_1)<x){
        S_B1 <- summary(fit1,t=max(time_1))$surv
      } else{
        S_B1 <- summary(fit1,t=x)$surv
      }
    }
    
    
    if(length(time_2)==0){
      S_B2 <- 1
    } else{
      if (max(time_2)<x){
        S_B2 <- summary(fit2,t=max(time_2))$surv
      } else{
        S_B2 <- summary(fit2,t=x)$surv
      }
    }
    
    if(S_A1>S_B1+d_diff1 & S_A2>S_B2+d_diff2){
      corr <- corr+1
    }
    if(S_A1<S_B1+d_diff1 & S_A2<S_B2+d_diff2){
      err <- err+1
    }
    
  }
  amb = 5000-corr-err
  
  return(c(corr,amb))
  
}