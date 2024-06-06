library(readr)
X_1 <- read_csv("missipi-river/Vicksberg2324 - Sheet1.csv")
X_2 <- read_csv("missipi-river/Tennesi2324 - Sheet1.csv")

shorten<-function(data){
  data<-data[,c(3,4,6)]
  colnames(data)<-c("Date", "Time", "Level")
  dat1<- data %>% group_by(Date) %>% summarise(glevel=mean(Level), n=n())
  dat1$Date<- format(dat1$Date, "%m-%d")
  ind<-which(dat1$Date %in% c("09-30","03-12", "02-29"))
  if(length(ind)>0){dat1<-dat1[-ind, ]}
  return(dat1)
}

dates<- shorten(X_1)$Date
tot_dat<-as.data.frame(matrix(nrow=241, ncol=2))
for(i in 1:2){
  tot_dat[,i]<- (shorten(get(paste0("X_", i))))[,2]
}
#colnames(tot_dat)<-c("2017-18", "2018-19", "2019-20", "2020-21", "2021-22", "2022-23", "2023-24")
source("Algorithm2.R")

res<-SyncBootstage2(X=as.matrix(tot_dat), n=241, d=2)



date_range <- seq(as.Date("2022-09-01"), as.Date("2023-05-01"), by = "day")[-c(30,193)]

colnames(tot_dat)<-c("Vicksburg, MS", "Memphis, TN")

my_data_long <- cbind( date_range, tot_dat[,c(1:2)]) %>%
  gather(key = "variable", value = "value", -1)
#colnames(my_data_long)<-c("X1", "variable", "value")

vertical_dates <- date_range[res$indiv_tau]

# Plot the data using ggplot2
ggplot(my_data_long, aes(x =date_range, 
                         y = value, color = variable,  group = variable)) +
  geom_line() +
  labs(title = "",
       x = "Dates",
       y = "Discharge, cubic feet per second",
       color = "Legend Title")+
  facet_wrap(~variable, nrow = 2, ncol = 1, scales = "fixed") +
  scale_x_date(date_breaks = "1 month", date_labels = "%m/%d") +
  theme(legend.position = "none",  # Remove legend box
        panel.background = element_rect(fill = "white"),  # Set background color to white
        panel.grid.major = element_line(color = "gray", linetype = "dotted"),  # Add major grid lines
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  geom_vline(data = data.frame(variable = unique(my_data_long$variable), 
                               X1 = vertical_dates), 
             aes(xintercept = X1), linetype = "dashed", color = "red")

