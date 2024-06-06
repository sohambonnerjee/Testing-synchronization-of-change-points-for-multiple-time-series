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

tau<- c(1/2, 1/2, 1/2, 1/2)
delta<- c(1/log(n),1/(log(n)), -1/(log(n)), 0) 
mean_mat<- matrix(nrow=n, ncol=d)
for(j in 1:d){
  mean_mat[,j]<- c(rep(0, floor(n*tau[j])), rep(delta[j], n-floor(n*tau[j])))
}

e_innov<- mGJR_garch(n, d, Cov_Innov, coeff_mat)
X<- mean_mat+e_innov

library(ggplot2)
df <- data.frame(Time = (1:n)/n, 
                 Series1 = X[,1],
                 Series2 = X[,2],
                 Series3 = X[,3],
                 Series4 = X[,4], 
                 Vline = tau)
library(tidyr)
df_long <- gather(df, key = "Series", value = "Value", -Time, -Vline)

# Plotting
vertical_lines <- data.frame(
  Series = c("Series1", "Series2", "Series3", "Series4"),
  Vline = tau
)

# Original colors
series_colors <- c("Series1" = "#E54E61", "Series2" = "#ACD85F", "Series3" = "#5FD8CF", "Series4" = "#CD5FD8")



# Plotting
p1<-ggplot(df_long, aes(x = Time, y = Value, color = Series)) +
  geom_line() +
  geom_vline(data = vertical_lines, aes(xintercept = Vline), linetype = "dashed") +
  facet_wrap(~ Series, scales = "free_y", ncol = 1) +
  scale_color_manual(values = series_colors) +  # Restore original colors
  theme_minimal() +
  theme(panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_text(size = 15),  # Adjust size of x-axis label
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),  # Adjust size of x-axis tick labels
        axis.text.y = element_text(size = 8)) +
  labs(x = "", y = "X") +
  guides(color = FALSE)

tau<- c(1/2, 1/2-0.1, 1/2+0.1, 1/2)
delta<- c(1/log(n),1/(log(n)), -1/(log(n)), 0) 
mean_mat<- matrix(nrow=n, ncol=d)
for(j in 1:d){
  mean_mat[,j]<- c(rep(0, floor(n*tau[j])), rep(delta[j], n-floor(n*tau[j])))
}

e_innov<- mGJR_garch(n, d, Cov_Innov, coeff_mat)
X<- mean_mat+e_innov

df <- data.frame(Time = (1:n)/n, 
                 Series1 = X[,1],
                 Series2 = X[,2],
                 Series3 = X[,3],
                 Series4 = X[,4], 
                 Vline = tau)
library(tidyr)
df_long <- gather(df, key = "Series", value = "Value", -Time, -Vline)

# Plotting
vertical_lines <- data.frame(
  Series = c("Series1", "Series2", "Series3", "Series4"),
  Vline = tau
)

# Original colors
series_colors <- c("Series1" = "#E54E61", "Series2" = "#ACD85F", "Series3" = "#5FD8CF", "Series4" = "#CD5FD8")



# Plotting
p2<-ggplot(df_long, aes(x = Time, y = Value, color = Series)) +
  geom_line() +
  geom_vline(data = vertical_lines, aes(xintercept = Vline), linetype = "dashed") +
  facet_wrap(~ Series, scales = "free_y", ncol = 1) +
  scale_color_manual(values = series_colors) +  # Restore original colors
  theme_minimal() +
  theme(panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_text(size = 15),  # Adjust size of x-axis label
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),  # Adjust size of x-axis tick labels
        axis.text.y = element_text(size = 8)) +
  labs(x = "", y = "X") +
  guides(color = FALSE)


