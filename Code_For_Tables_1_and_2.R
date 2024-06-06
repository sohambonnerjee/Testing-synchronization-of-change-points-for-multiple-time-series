

d<-4; n=1000 #500
A<-matrix(nrow=d, ncol=d)
for(i in 1:d){
  for(j in 1:d){
    A[i,j]<- 0.3*exp(-abs(i-j))
  }
}

library(tsDyn)
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

Cov_Innov<-RQK(d, 1, 5)



S<- matrix(solve(diag(1, nrow=d*d) - kronecker(A, t(A)))  %*% 
             matrix(as.vector(Cov_Innov + A %*% solve(diag(1, nrow=d) - A ) %*% Cov_Innov +
                                t(A %*% solve(diag(1, nrow=d) - A ) %*% Cov_Innov)), ncol=1), nrow=d)

tau<- c(1/2, 1/2, 1/2, 1/2)
delta<-4
delta_1<- c(0,0,0,0) # model 1
delta_2 <- c(delta, 0,0,0) 
delta_3<-c(delta,delta,0,0)
delta_4<-c(delta,delta,delta,0)
delta_5<-c(delta,delta,delta,delta)

mean_mat1<- matrix(nrow=n, ncol=d) ; mean_mat2<- matrix(nrow=n, ncol=d) 
mean_mat3<- matrix(nrow=n, ncol=d) ;  mean_mat4<- matrix(nrow=n, ncol=d) 
mean_mat5<- matrix(nrow=n, ncol=d) 
for(j in 1:d){
  mean_mat1[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_1[j], n-floor(n*tau[j])))
  mean_mat2[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_2[j], n-floor(n*tau[j])))
  mean_mat3[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_3[j], n-floor(n*tau[j])))
  mean_mat4[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_4[j], n-floor(n*tau[j])))
  mean_mat5[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_5[j], n-floor(n*tau[j])))
}

niter<-5000

row.sum<-function(x){
  if(length(nrow(x))==0){return(sum(x))}
  else{return(rowSums(x))}
}

kernel1<- function(u, type="1"){
  if(type=="Parzen"){
    if(abs(u)<=1/2){
      return(1- 6 *u^2+ 6 * abs(u^3))
    }
    else{
      return(ifelse(abs(u)<=1, 2*(1- abs(u))^3, 0))
    }
  }
  else if(type=="Tukey-Hanning"){
    return(ifelse(abs(u)<=1, (1+cos(pi*u))/2, 0))
  }
  else if(type=="Block type"){
    if(abs(u)<=0.95){
      return(1)
    }
    else{
      return( ifelse(abs(u)<=1, (1+cos(20*(u-0.95)*pi))/2, 0) )
    }
  }
}


Cov_est_err<- function(X){
  X<-as.matrix(X); n<- nrow(X); d<-ncol(X)
  indiv_cusum<-matrix(nrow=n, ncol=d); common_tau<-0
  for(j in 1:d){
    indiv_cusum[,j]<- (cumsum(X[,j]-mean(X[,j])))^2
  }
  indiv_tau<- apply(indiv_cusum, 2, which.max) # estimate individual tau_j
  common_tau<- which.max(rowSums(indiv_cusum)) # estimate common tau
  
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
  
  B_n<- (floor(n^{1/5}): floor(n^{1/3}))
  autocov<- function(k){
    cov_mat<- matrix(0,nrow=d, ncol=d)
    for(i in 1:(n-k)){
      cov_mat<-cov_mat+(1/n)* (X[i,] - est_true_mean[i,]) %*% t(X[i+k,] - est_true_mean[i+k,])
    }
    return(cov_mat)
  }
  
  klist<-c("Parzen", "Tukey-Hanning", "Block type")
  
  err_mat<- matrix(nrow=length(klist), ncol=length(B_n))
  for (w in 1: length(klist)) {
    for (z in 1:length(B_n)) {
      Shat<- autocov(0)+ Reduce('+', lapply(1:B_n[z], function(k){kernel1(k/B_n[z], klist[w])*(autocov(k)+ t(autocov(k)))}))
      err_mat[w,z]<- max(abs(eigen(Shat -S)$values)) 
    }
  }
  
  
  return(err_mat)
  
}

t1<- list(); t2<- list(); t3<- list(); t4<- list(); t5<-list()

for (i in 1: niter) {
  e<-VAR.sim(A, n=n, include = "none", varcov = Cov_Innov)
  X1<- mean_mat1+e; t1[[i]]<- Cov_est_err(X1)
  X2<- mean_mat2+e; t2[[i]]<- Cov_est_err(X2)
   X3<- mean_mat3+e; t3[[i]]<- Cov_est_err(X3)
  X4<- mean_mat4+e; t4[[i]]<- Cov_est_err(X4)
 X5<- mean_mat5+e; t5[[i]]<- Cov_est_err(X5)
}

write_csv(data.frame(sqrt(Reduce(`+`, lapply(t1, function(x) (x - mean(x))^2))/niter)), "t1sd.csv")

write_csv(data.frame(Reduce(`+`, t1)/niter), "t1.csv")

write_csv(data.frame(sqrt(Reduce(`+`, lapply(t2, function(x) (x - mean(x))^2))/niter)), "t2sd.csv")

write_csv(data.frame(Reduce(`+`, t2)/niter), "t2.csv")

write_csv(data.frame(sqrt(Reduce(`+`, lapply(t3, function(x) (x - mean(x))^2))/niter)), "t3sd.csv")

write_csv(data.frame(Reduce(`+`, t3)/niter), "t3.csv")

write_csv(data.frame(sqrt(Reduce(`+`, lapply(t4, function(x) (x - mean(x))^2))/niter)), "t4sd.csv")

write_csv(data.frame(Reduce(`+`, t4)/niter), "t4.csv")

write_csv(data.frame(sqrt(Reduce(`+`, lapply(t5, function(x) (x - mean(x))^2))/niter)), "t5sd.csv")

write_csv(data.frame(Reduce(`+`, t5)/niter), "t5.csv")


data_list<-list()

for( i in 1:5){
  # Create the file name
  file_name <- paste0("t", i, ".csv")
  file_name_sd <- paste0("t", i, "sd.csv")
  
  # Create the file path
  file_path <- paste0("", file_name)
  file_path_sd <- paste0("", file_name_sd)
  
  
  # Load the CSV file
  data <- round(read.csv(file_path),3)
  data_sd<- round(read.csv(file_path_sd), 3)
  
  
  # Add data to the list
  data_list[[i]] <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))
  
  for (k in 1:nrow(data)) {
    for (l in 1:ncol(data)) {
      data_list[[i]][k, l] <- paste0(data[k, l], "(", data_sd[k, l], ")")
    }
  }
  
}

combined_df<- do.call(rbind, data_list)
write_csv(data.frame(combined_df), "comb1000.csv") #for n=500, write comb500.csv

matrix_to_latex <- function(mat) {
  nrows <- nrow(mat)
  ncols <- ncol(mat)
  
  latex_code <- "\\begin{pmatrix}\n"
  
  for (i in 1:nrows) {
    for (j in 1:ncols) {
      latex_code <- paste(latex_code, format(mat[i, j], digits = 3))
      if (j < ncols) {
        latex_code <- paste(latex_code, " & ")
      }
    }
    if (i < nrows) {
      latex_code <- paste(latex_code, " \\\\ \n")
    }
  }
  
  latex_code <- paste(latex_code, "\n\\end{pmatrix}")
  
  return(latex_code)
}

# Example usage:

latex_code <- matrix_to_latex(S)
cat(latex_code)

