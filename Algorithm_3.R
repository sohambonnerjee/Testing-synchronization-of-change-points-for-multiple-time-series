row.sum<-function(x){
  if(length(nrow(x))==0){return(sum(x))}
  else{return(rowSums(x))}
}

kernel<- function(u){
  return(ifelse(abs(u)<=1, 1-u^2, 0))
}

SyncBootstage2<-function(X, n, d, nboot=5000, B_n=floor(n^{1/4})){
  X<-as.matrix(X)
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
  
  
  autocov<- function(k){
    cov_mat<- matrix(0,nrow=d, ncol=d)
    for(i in 1:(n-k)){
      cov_mat<-cov_mat+(1/n)* (X[i,] - est_true_mean[i,]) %*% t(X[i+k,] - est_true_mean[i+k,])
    }
    return(cov_mat)
  }
  
  Shat<- autocov(0)+ Reduce('+', lapply(1:B_n, function(k){kernel(k/B_n)*(autocov(k)+ t(autocov(k)))}))
  
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
  
  return(list(pval=pval, common_tau=common_tau, indiv_tau=indiv_tau, T_n=T_n, boot=T_n_boot, Shat=Shat, pos_dim=pos_dim ))
}
