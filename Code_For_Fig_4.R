
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


tau<- c(1/2, 1/2, 1/2, 1/2)
delta<- c(4/log(n), -4/log(n), 4/log(n), 0) 
mean_mat<- matrix(nrow=n, ncol=d)
for(j in 1:d){
  mean_mat[,j]<- c(rep(0, floor(n*tau[j])), rep(delta[j], n-floor(n*tau[j])))
}

e_innov<-sweep(mTAR(n, d, Thresh_vec, AR_mat, Cov_Innov), 2, mu) 
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

tau<- c(1/2, 1/2-0.1, 1/2+0.01, 1/2)
delta<- c(4/log(n), -4/log(n), 4/log(n), 0) 
mean_mat<- matrix(nrow=n, ncol=d)
for(j in 1:d){
  mean_mat[,j]<- c(rep(0, floor(n*tau[j])), rep(delta[j], n-floor(n*tau[j])))
}

e_innov<-sweep(mTAR(n, d, Thresh_vec, AR_mat, Cov_Innov), 2, mu) 
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


