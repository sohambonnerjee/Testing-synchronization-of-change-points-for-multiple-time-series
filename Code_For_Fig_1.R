library(ggplot2)
library(gridExtra)
library(tsDyn)

d<-4
A<-matrix(nrow=d, ncol=d)
for(i in 1:d){
  for(j in 1:d){
    A[i,j]<- 0.3*exp(-abs(i-j))
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

Cov_Innov<-RQK(d, 1, 5)

tau<- rep(1/2, d)
delta<-1.5
delta_c <- list()
for( j in 1:(d+1)){
  delta_c[[j]] <- c(rep(delta, j-1), rep(0, d-j+1))
}

n=500
mean_mat<- list()
for(j in 1:(d+1)){
  mean_mat[[j]]<- matrix(nrow=n, ncol=d)
  for(k in 1:d){
    mean_mat[[j]][,k]<-c(rep(0, floor(n*tau[k])), rep(delta_c[[j]][k], n-floor(n*tau[k])))
  }
}
# for(j in 1:d){
#   mean_mat1[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_1[j], n-floor(n*tau[j])))
#   mean_mat2[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_2[j], n-floor(n*tau[j])))
#   mean_mat3[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_3[j], n-floor(n*tau[j])))
#   mean_mat4[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_4[j], n-floor(n*tau[j])))
#   mean_mat5[,j]<- c(rep(0, floor(n*tau[j])), rep(delta_5[j], n-floor(n*tau[j])))
# }

test_stat<- function(X){
  X<-as.matrix(X); n<- nrow(X); d<-ncol(X)
  indiv_cusum<-matrix(nrow=n, ncol=d); common_tau<-0
  for(j in 1:d){
    indiv_cusum[,j]<- (cumsum(X[,j]-mean(X[,j])))^2
  }
  indiv_tau<- apply(indiv_cusum, 2, which.max) # estimate individual tau_j
  common_tau<- which.max(rowSums(indiv_cusum)) # estimate common tau
  
  T_n<- sum(sapply(1:d, function(j){max(indiv_cusum[,j])- indiv_cusum[common_tau, j]}))/n # test statistic
  return(T_n)
}

niter<-5000; ts<-list()

for(j in 1:(d+1)){
  ts[[j]]<-rep(0, niter)
}

for (i in 1: niter) {
  e<-VAR.sim(A, n=n, include = "none", varcov = Cov_Innov)
  for(j in 1:(d+1)){
    ts[[j]][i]<-test_stat(mean_mat[[j]]+e)
  }
  
  # X1<- mean_mat1+e; t1[i]<- test_stat(X1)
  # X2<- mean_mat2+e; t2[i]<- test_stat(X2)
  # X3<- mean_mat3+e; t3[i]<- test_stat(X3)
  # X4<- mean_mat4+e; t4[i]<- test_stat(X4)
  # X5<- mean_mat5+e; t5[i]<- test_stat(X5)
}

data <- data.frame(
  value = unlist(ts),
  group = rep(paste0("Model ", 1:(d+1)), each = niter)
)

# Plot the density of the five vectors
p500<-ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(title = "", x = "Value", y = "Density", fill = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


######

n=1000
mean_mat<- list()
for(j in 1:(d+1)){
  mean_mat[[j]]<- matrix(nrow=n, ncol=d)
  for(k in 1:d){
    mean_mat[[j]][,k]<-c(rep(0, floor(n*tau[k])), rep(delta_c[[j]][k], n-floor(n*tau[k])))
  }
}

niter<-5000; ts<-list()

for(j in 1:(d+1)){
  ts[[j]]<-rep(0, niter)
}

for (i in 1: niter){
  e<-VAR.sim(A, n=n, include = "none", varcov = Cov_Innov)
  for(j in 1:(d+1)){
    ts[[j]][i]<-test_stat(mean_mat[[j]]+e)
  }
}

data <- data.frame(
  value = unlist(ts),
  group = rep(paste0("Model ", 1:(d+1)), each = niter)
)

# Plot the density of the five vectors
p1000<-ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(title = "", x = "Value", y = "Density", fill = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


grid.arrange(p500, p1000, nrow=1)
