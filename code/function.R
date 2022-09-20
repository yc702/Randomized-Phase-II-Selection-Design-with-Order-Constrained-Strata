library(MASS)
library(mvtnorm)
library(quantreg)
library(quadprog)


QP<-function(p0,r,n){
  p_bar<-r/n
  S<-length(p0)
  Dmat<-matrix(0,S,S)
  diag(Dmat)<-n
  dvec<-r
  A_upper<-diag(S)
  A_lower<-matrix(0,nrow = (S-1),ncol = S)
  for(i in 1:(S-1)){
    A_lower[i,i:(i+1)]<-c(-1,1)
  }
  Amat<-rbind(A_upper,A_lower)
  # Amat <- A_lower
  Amat<-t(Amat)
  # b_vec<-c(p0,numeric(S-1))
  b_vec<-c(numeric(S-1+length(p0)))
  
  solve.QP(Dmat,dvec,Amat,b_vec)}


pickwin_strat_sargent<- function(n, p1, strata_diff,
                                 D1=0.15, D2=0.15,d1=0.05,d2=0.05, prop.strat=0.3)
{
  n1<- ceiling(prop.strat*n)
  n2<- n-n1
  Nk <- c(n1,n2)
  p <- c(0.1,0.1)
  pp_err <- 0
  pp_corr <- 0
  for (i1 in 0:n1){
    for (i2 in 0:n2){
      Rk1 <- c(i1,i2)
      p_1 <- QP(p,Rk1,Nk)$solution
      p11 <- p_1[1]
      p12 <- p_1[2]
      
      for (j1 in 0:n1) {
        for (j2 in 0:n2){
          Rk2 <- c(j1,j2)
          p_2 <- QP(p,Rk2,Nk)$solution
          p21 <- p_2[1]
          p22 <- p_2[2]
          
          if ((p21<p11+d1) & (p22<p12+d2)){
            
            pp_err <- pp_err+dbinom(i1, size=n1, prob=p1)*dbinom(i2, size=n2, prob=p1+strata_diff)*
                      dbinom(j1, size=n1, prob=p1+D1)*dbinom(j2, size=n2, prob=p1+strata_diff+D2)
          }
          if ((p21>p11+d1) & (p22>p12+d2)) {
            pp_corr <- pp_corr+dbinom(i1, size=n1, prob=p1)*dbinom(i2, size=n2, prob=p1+strata_diff)*
              dbinom(j1, size=n1, prob=p1+D1)*dbinom(j2, size=n2, prob=p1+strata_diff+D2)
          }

        }
        
      }
      
    }
    
  }
  p_amb = 1-pp_err-pp_corr
  return(data.frame(pcorr=pp_corr, pamb=p_amb))
}



## original idea
pickwin_original<- function(n, p1, strata_diff, D1=0.15, D2=0.15,d1=0.05,d2=0.05, prop.strat=0.3)
{
  n1<- ceiling(prop.strat*n)
  n2<- n-n1
  pp_err <- 0
  pp_corr <- 0
  for (i1 in 0:n1){
    for (i2 in 0:n2){
      p11 <- i1/n1
      p12 <- i2/n2
      
      for (j1 in 0:n1) {
        for (j2 in 0:n2){
          p21 <- j1/n1
          p22 <- j2/n2
          if ((p21<p11+d1) & (p22<p12+d2)){
            
            pp_err <- pp_err+dbinom(i1, size=n1, prob=p1)*dbinom(i2, size=n2, prob=p1+strata_diff)*
              dbinom(j1, size=n1, prob=p1+D1)*dbinom(j2, size=n2, prob=p1+strata_diff+D2)
          }
          if ((p21>p11+d1) & (p22>p12+d2)) {
            pp_corr <- pp_corr+dbinom(i1, size=n1, prob=p1)*dbinom(i2, size=n2, prob=p1+strata_diff)*
              dbinom(j1, size=n1, prob=p1+D1)*dbinom(j2, size=n2, prob=p1+strata_diff+D2)
          }

        }
        
      }
      
    }
    
  }
  p_amb = 1-pp_err-pp_corr
  return(data.frame(pcorr=pp_corr, pamb=p_amb))
}

