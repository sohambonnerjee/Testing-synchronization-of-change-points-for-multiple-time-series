## this yields the row corresponding to "n=1000". 
## For the row "n=500", change n in Line 155. 
row.sum<-function(x){
  if(length(nrow(x))==0){return(sum(x))}
  else{return(rowSums(x))}
}

kernel<- function(u){
  return(ifelse(abs(u)<=1, 1-u^2, 0))
}

SyncBootstage2<-function(X, n, d){
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
  
  ### two stage bootstrap
  est_boot_mean<- matrix(nrow=n, ncol=d)
  for(j in 1:d){
    est_boot_mean[,j]<-c(rep(est_left_mean[j], common_tau), rep(est_right_mean[j], n-common_tau))
  }
  
  T_n_boot<-c()
  Z<-list(); Z_indiv<-list()
  cusum_count<-matrix(nrow=1000, ncol=d)
  
  #perform cusum for each dimension
  for ( B in 1:1000) {
    Z[[B]]<- MASS::mvrnorm(n, mu=rep(0,d), Sigma = Shat)
    for (j in 1:d) {
      cusum_count[B,j]<- (max(indiv_cusum[,j]) >= max((cumsum((Z[[B]])[,j]-mean((Z[[B]])[,j])))^2))
    }
  }  
  
  pos_dim<- which(colMeans(cusum_count)>1-0.05) # which dimesions have jumps
  #if(length(pos_dim)==0){pval<- 1} # accept null if all dimensions have zero jumps
  # else{ # perform bootstrap on dimensions with significant jumps
  #T_n_trunc<- sum(sapply(pos_dim, function(j){max(indiv_cusum[,j])- indiv_cusum[common_tau, j]}))/n
  for(B in 1:1000){
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
}


# simulation of TAR 

# generate Multivariate Threshold Autoregressive Model

# thresh_vec: number of threshold parameters
# ar_mat: matrix of AR coefficients: rows: threshold+1, columns: dimension

mTAR<- function(n, d, thresh_vec, ar_mat, Innov){
  if(nrow(Innov)!=d || ncol(Innov)!=d){print("Error in dimension; check input")}
  else{
    nthresh<- length(thresh_vec)
    if(nrow(ar_mat)-nthresh!=1){print("Error; check the sizes of threshold/AR paramteres")}
    else{
      error<- MASS::mvrnorm(2*n, mu=rep(0,d), Sigma=Innov)
      thresh_vec<-sort(thresh_vec)
      thresh_fn<- function(x){
        ar_thresh<-c()
        for(j in 1:d){
          l_pos<-c();l_pos<- which( (x[j]-thresh_vec) >0)
          ar_thresh[j]<- ifelse(length(l_pos)==0, ar_mat[1,j], ar_mat[max(l_pos)+1, j])
        }
        
        return(ar_thresh)
      }
      X<- matrix(nrow=2*n, ncol=d)
      X[1,]<-rep(0,d)
      for (i in 2:(2*n)) {
        X[i,]<- thresh_fn(X[i-1,])* X[i-1,] + error[i,]
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
Cov_Innov<-RQK(d, 1, 5); Thresh_vec<-0; AR_mat<-matrix(rbind(rep(0.5, d), rep(-0.5, d)), nrow=2) ## settings to simulate innovations
X_dummy<- mTAR(5000, d, Thresh_vec, AR_mat, Cov_Innov); mu<- colMeans(X_dummy)
remove(X_dummy)

## settings for jumps and changepoints
r<- seq(0,0.1, by=0.01)
pow<-c()
for(w in 1:length(r)){
    tau<- c(1/2, 1/2-r[w], 1/2, 1/2)
    delta<- c(4/log(n), -4/log(n), 0, 0) 
    mean_mat<- matrix(nrow=n, ncol=d)
    for(j in 1:d){
      mean_mat[,j]<- c(rep(0, floor(n*tau[j])), rep(delta[j], n-floor(n*tau[j])))
    }
    
    niter<- 5000
    
    
    #source("twostage_modified.R")
    p_list<- c()
    for (iter in 1:niter) {
      e_innov<-sweep(mTAR(n, d, Thresh_vec, AR_mat, Cov_Innov), 2, mu) 
      X<- mean_mat+e_innov
      p_list[iter]<- SyncBootstage2(X,n,d)
    }
    pow[w]<-mean(p_list<0.05) # type-1-error or power
  }
  

write.table(pow, file = "table4.txt", row.names = FALSE, col.names = FALSE)



