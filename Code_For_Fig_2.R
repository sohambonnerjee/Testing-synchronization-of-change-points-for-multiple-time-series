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

niter<-5000; ts<-list(); ts_boot<-list()

for(j in 1:(d+1)){
  ts[[j]]<-rep(0, niter)
  ts_boot[[j]]<- rep(0, niter)
}

S<- matrix(solve(diag(1, nrow=d*d) - kronecker(A, t(A)))  %*% 
             matrix(as.vector(Cov_Innov + A %*% solve(diag(1, nrow=d) - A ) %*% Cov_Innov +
                                t(A %*% solve(diag(1, nrow=d) - A ) %*% Cov_Innov)), ncol=1), nrow=d)

for (i in 1: niter) {
  e<-VAR.sim(A, n=n, include = "none", varcov = Cov_Innov)
  Z<- MASS::mvrnorm(n, mu=rep(0,d), Sigma=S)
  for(j in 1:(d+1)){
    ts[[j]][i]<-test_stat(mean_mat[[j]]+e)
    ts_boot[[j]][i]<-test_stat(mean_mat[[j]]+Z)
  }
}

data<- data.frame(
  value = unlist(ts), value_boot=unlist(ts_boot),
  group = rep(paste0("Model ", 1:(d+1)), each = niter)
)

p1<- ggplot(data[1:5000,], aes(x=sort(value), y=sort(value_boot)))+geom_point()+theme(legend.title = 
                                                                   element_blank(), legend.text =  element_text(size = 18))+
  geom_abline(slope=1, intercept = 0, col="red")+labs(x=expression(Quantiles~of~T[n]) , y=expression(Quantiles~of~T[n]^Z))+
  ggtitle(label="") + theme_minimal()

p2<- ggplot(data[5001:10000,], aes(x=sort(value), y=sort(value_boot)))+geom_point()+theme(legend.title = 
                                                                                        element_blank(), legend.text =  element_text(size = 18))+
  geom_abline(slope=1, intercept = 0, col="red")+labs(x=expression(Quantiles~of~T[n]) , y=expression(Quantiles~of~T[n]^Z))+
  ggtitle(label="") + theme_minimal()

p3<- ggplot(data[10001:15000,], aes(x=sort(value), y=sort(value_boot)))+geom_point()+theme(legend.title = 
                                                                                            element_blank(), legend.text =  element_text(size = 18))+
  geom_abline(slope=1, intercept = 0, col="red")+labs(x=expression(Quantiles~of~T[n]) , y=expression(Quantiles~of~T[n]^Z))+
  ggtitle(label="") + theme_minimal()


p4<- ggplot(data[15001:20000,], aes(x=sort(value), y=sort(value_boot)))+geom_point()+theme(legend.title = 
                                                                                             element_blank(), legend.text =  element_text(size = 18))+
  geom_abline(slope=1, intercept = 0, col="red")+labs(x=expression(Quantiles~of~T[n]) , y=expression(Quantiles~of~T[n]^Z))+
  ggtitle(label="") + theme_minimal()

p5<- ggplot(data[20001:25000,], aes(x=sort(value), y=sort(value_boot)))+geom_point()+theme(legend.title = 
                                                                                             element_blank(), legend.text =  element_text(size = 18))+
  geom_abline(slope=1, intercept = 0, col="red")+labs(x=expression(Quantiles~of~T[n]) , y=expression(Quantiles~of~T[n]^Z))+
  ggtitle(label="") + theme_minimal()


grid.arrange(p1,p2,p3,p4,p5, nrow=1)  
