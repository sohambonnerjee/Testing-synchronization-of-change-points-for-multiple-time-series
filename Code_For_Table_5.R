## this yields the rows corresponding to "n=1000". 
## For the rows "n=500", change n in Line 153. 
row.sum<-function(x){
  if(length(nrow(x))==0){return(sum(x))}
  else{return(rowSums(x))}
}

kernel<- function(u){
  return(ifelse(abs(u)<=1, 1-u^2, 0))
}

SyncBootstage2<-function(X, n, d, nboot=5000){
  indiv_cusum<-matrix(nrow=n, ncol=d); common_tau<-0
  for(j in 1:d){
    indiv_cusum[,j]<- (cumsum(X[,j]-mean(X[,j])))^2
  }
  indiv_tau<- apply(indiv_cusum, 2, which.max) # estimate individual tau_j
  common_tau<- which.max(rowSums(indiv_cusum)) # estimate common tau
  
  T_n<- sum(sapply(1:d, function(j){max(indiv_cusum[,j])- indiv_cusum[common_tau, j]}))/n # test statistic
  
  ## estimate common means
  
  est_left_mean<- sapply(1:d, function(j){mean(X[(1:common_tau), j])})
  est_right_mean<- sapply(1:d, function(j){mean(X[((common_tau+1):n), j])})
  
  
  ## estimate covariance matrix under alternate
  
  est_true_left<- sapply(1:d, function(j){mean(X[(1:indiv_tau[j]), j])})
  est_true_right<- sapply(1:d, function(j){mean(X[((indiv_tau[j]+1):n), j])})
  est_true_mean<- matrix(nrow=n, ncol=d)
  for(j in 1:d){
    est_true_mean[,j]<-c(rep(est_true_left[j], indiv_tau[j]), rep(est_true_right[j], n-indiv_tau[j]))
  }
  
  B_n<- floor(n^{1/4})
  autocov<- function(k){
    cov_mat<- matrix(0,nrow=d, ncol=d)
    for(i in 1:(n-k)){
      cov_mat<-cov_mat+(1/n)* (X[i,] - est_true_mean[i,]) %*% t(X[i+k,] - est_true_mean[i+k,])
    }
    return(cov_mat)
  }
  
  
  Shat<- autocov(0)+ Reduce('+', lapply(1:B_n, function(k){kernel(k/B_n)*(autocov(k)+ t(autocov(k)))}))
  #Eg<-eigen(Shat) ; D<- Eg$values; V<-Eg$vectors
  # while(length(D[Re(D)<0])>0){
  #   B_n<-B_n-1
  #   Shat<- autocov(0)+ 2*Reduce('+', lapply(1:(B_n),function(k){kernel(k/(B_n))*autocov(k)}))
  #   D<- eigen(Shat)$values
  #   }
  
  ### two stage bootstrap
  est_boot_mean<- matrix(nrow=n, ncol=d)
  for(j in 1:d){
    est_boot_mean[,j]<-c(rep(est_left_mean[j], common_tau), rep(est_right_mean[j], n-common_tau))
  }
  
  T_n_boot<-c()
  Z<-list(); Z_indiv<-list()
  cusum_count<-matrix(nrow=nboot, ncol=d)
  
  #perform cusum for each dimension
  for ( B in 1:nboot) {
    Z[[B]]<- MASS::mvrnorm(n, mu=rep(0,d), Sigma = Shat)
    for (j in 1:d) {
      cusum_count[B,j]<- (max(indiv_cusum[,j]) >= max((cumsum((Z[[B]])[,j]-mean((Z[[B]])[,j])))^2))
    }
  }  
  
  pos_dim<- which(colMeans(cusum_count)>1-0.05) # which dimesions have jumps
  #if(length(pos_dim)==0){pval<- 1} # accept null if all dimensions have zero jumps
  # else{ # perform bootstrap on dimensions with significant jumps
  #T_n_trunc<- sum(sapply(pos_dim, function(j){max(indiv_cusum[,j])- indiv_cusum[common_tau, j]}))/n
  for(B in 1:nboot){
    for(j in pos_dim){
      Z[[B]][,j]<-Z[[B]][,j]+est_boot_mean[,j]
    }
    for(j in setdiff(1:d, pos_dim)){
      Z[[B]][,j]<-Z[[B]][,j]+mean(X[,j])
    }
    indiv_cusum_boot<-matrix(nrow=n, ncol=d); common_tau_boot<-0
    for(j in 1:d){
      indiv_cusum_boot[,j]<- (cumsum(Z[[B]][,j]-mean(Z[[B]][,j])))^2
    }
    
    #indiv_tau_boot<- apply(indiv_cusum_boot[, pos_dim], 2, which.max) # estimate individual tau_j
    
    common_tau_boot<- which.max(row.sum(indiv_cusum_boot[, 1:d])) # estimate common tau
    
    T_n_boot[B]<- sum(sapply(1:d, function(j){max(indiv_cusum_boot[,j])-
        indiv_cusum_boot[common_tau_boot, j]}))/n # test statistic
    
  }
  
  pval<- 1- mean(T_n>T_n_boot)  
  #}
  
  return(pval)
  #return(Shat)
}



# CovInnov: Cross-dimension covariance matrix
# coeff_mat: coefficients of the GJR garch model. size: d times 3.
mGJR_garch<- function(n, d, Cov_Innov, coeff_mat){
  if(nrow(Cov_Innov)!=d || ncol(Cov_Innov)!=d){print("Error in dimension; check input")}
  else{
    if(nrow(coeff_mat)!=d){print("Check coefficient matrix")}
    else{
      error<- MASS::mvrnorm(2*n, mu=rep(0,d), Sigma=Cov_Innov)
      
      X<- matrix(nrow=2*n, ncol=d)
      Sigma<- matrix(nrow=2*n, ncol=d)
      X[1,]<-rep(0,d); Sigma[1,]<-abs(error[1,])
      for(i in 2:(2*n)){
        Sigma[i,]<- sqrt( coeff_mat[,1]+ coeff_mat[,2] * (Sigma[i-1,])^2 + coeff_mat[,3]*(X[i-1,])^2 + coeff_mat[,4] * 
                            (X[i-1,])^2 * (X[i-1,]<0) )
        X[i,] <- Sigma[i,]* error[i,]
      }
      return(X[(n+1):(2*n),])
    }
  }
}
## helper function to induce cross-sectional dependence
RQK<-function(d,k,a){
  S<-matrix(nrow=d, ncol=d)
  for(i in 1:d){
    for(j in 1:d){
      S[i,j]<-(1+ (i-j)^2/(2*a*k^2))^(-a)
    }
  }
  return(S)
}

EXPK<- function(d, l){
  S<-matrix(nrow=d, ncol=d)
  for(i in 1:d){
    for(j in 1:d){
      S[i,j]<-exp(-(i-j)^2/(2*l^2))
    }
  }
  return(S)
}



### define your setting of simulation

n<-1000; d=4
Cov_Innov<-RQK(d, 1, 5);
coeff_mat<- t(diag(c(0.01,0.7, 0.1, 0.2)) %*% matrix(1, nrow=4, ncol=d))

## settings for jumps and changepoints
r1<- seq(0,0.1, by=0.01) ;
r2<- r1
pow_mat<-matrix(0,nrow=11, ncol=11)
for(w in 1:length(r1)){
  for (t in 1:length(r2)) {
    tau<- c(1/2, 1/2-r1[w], 1/2+r2[t], 1/2)
  delta<- c(1/log(n),1/(log(n)), -1/(log(n)), 0) 
  mean_mat<- matrix(nrow=n, ncol=d)
  for(j in 1:d){
    mean_mat[,j]<- c(rep(0, floor(n*tau[j])), rep(delta[j], n-floor(n*tau[j])))
  }
  niter<- 5000
  
  
  #source("twostage_modified.R")
  p_list<- c()
  for (iter in 1:niter) {
    e_innov<- mGJR_garch(n, d, Cov_Innov, coeff_mat)
    X<- mean_mat+e_innov
    p_list[iter]<- SyncBootstage2(X,n,d)
  }
  pow_mat[w,t]<-mean(p_list<0.05) # type-1-error or power
  }
}  
  
write.table(pow_mat, file = "table5.txt", row.names = FALSE, col.names = FALSE)
