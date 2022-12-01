## likelihood function for survival simulation

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



## simulation of survival times
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
