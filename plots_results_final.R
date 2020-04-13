library(rstan)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
load("~/Dropbox/contagion_studies_cloud/results_contagion_poisson_hurdle_new_rate.RData")
typesim="hurdle"
appendix ="hurdle_new_rate"


today <- as.Date("2020-03-23")
end_fit <- as.Date("2020-03-15")
data <- read.csv2("~/Dropbox/COVID-19/archived_data/archived_time_series/time_series_19-covid-Confirmed_archived_0325.csv", sep=",", header=TRUE)
colnames(data)[5:ncol(data)] <- sapply(seq(as.Date("2020/1/22"), by = "day", to=today), toString)
data_recovered <- read.csv2("~/Dropbox/COVID-19/archived_data/archived_time_series/time_series_19-covid-Recovered_archived_0325.csv",sep=",", header=TRUE)
colnames(data_recovered)[5:ncol(data)] <- sapply(seq(as.Date("2020/1/22"), by = "day", to=today), toString)
data_death <- read.csv2("~/Dropbox/COVID-19/archived_data/archived_time_series/time_series_19-covid-Deaths_archived_0325.csv", sep=",", header=TRUE)
colnames(data_death)[5:ncol(data)] <- sapply(seq(as.Date("2020/1/22"), by = "day", to=today), toString)
T = ncol(data) - 4 
names <- sapply(1:141, FUN =function(x){
  return(ifelse(toString(data$Province.State[x]) !="",toString(data$Province.State[x]), toString(data$Country.Region[x])  ))
})



Time = dim(data)[2] -4
N = t(apply(data[,5:ncol(data)], 1,FUN=function(x){ diff(as.numeric(x))}))
D = t(apply(data_death[,5:ncol(data)], 1,FUN=function(x){ diff(as.numeric(x))}))
R = t(apply(data_recovered[,5:ncol(data)], 1,FUN=function(x){ diff(as.numeric(x))}))
print(c(min(N), min(D),min(R)))
D[which(D<0)] = 0
R[which(R<0)] = 0
N = N + D + R
print(which(N<0, arr.ind = TRUE))
print(data[248:259,c(1:2,(47+4):(50+4))])

N[which(N<0)] = 0





index_fr = which(data$Province.State == "France")
index_it = which(data$Country.Region == "Italy")
index_germany = which(data$Country.Region == "Germany")
index_spain = which(data$Country.Region == "Spain")
index_sw = which(data$Country.Region == "Switzerland")
index_uk = which(data$Country.Region == "United Kingdom")[3]
index_japan = which(data$Country.Region == "Japan")
index_iran = which(data$Country.Region == "Iran")
index_thai = which(data$Country.Region == "Thailand")
index_sing = which(data$Country.Region == "Singapore")
index_korea = which(data$Country.Region == "Korea, South")
index_dp = which(data$Province.State == "Diamond Princess")[2]
index_cal = which(data$Province.State == "California")
index_wash = which(data$Province.State == "Washington")
index_hubei = which(data$Province.State == "Hubei")
index_hk = which(data$Province.State == "Hong Kong")
index_guizhou = which(data$Province.State == "Guizhou")
index_ny = which(data$Province.State == "New York")
index_us = which(data$Country.Region == "US")[1:52] ### Not the  granular stuff

indexes <- c(index_fr, index_it, index_germany, index_spain, 
             index_sw, index_uk, index_iran, index_japan,
             index_thai,index_sing, index_korea,index_dp,index_cal,
             index_wash,index_ny, index_hubei,index_hk, index_guizhou )



xx<- N[c(index_hubei), 1:36]
xx <- rbind(xx, N[c(index_hubei,index_hk, index_guizhou,index_fr, index_it, index_germany, index_spain,
                    index_sw, index_uk,index_iran,
                    index_japan, index_thai, index_sing, index_korea, index_cal,
                    index_wash, index_ny), 18:53])
xx<- rbind(xx,apply(N[index_us, 18:53],2, sum) )
names= c("Hubei0","Hubei1", "Hong Kong",  "Guizhou",
         "France", "Italy",  "Germany", "Spain", "Switzerland", "UK","Iran",
         "Japan", "Thailand",
         "Singapore", "Korea, South","California",
         "Washington", "New York",   "US" )
indexes <- c(index_hubei, index_hubei,index_hk, index_guizhou,index_fr, index_it, index_germany, index_spain,
             index_sw, index_uk,index_iran,
             index_japan, index_thai, index_sing, index_korea, index_cal,
             index_wash, index_ny, "NA")

population <- data.frame(row.names = names,
                         pop = c(58.5 , 58.5 , 7.4, 34.75, 67, 60.6, 82.9,47,8.57,66.87,
                                 82.9, 126.9, 69.6, 5.7, 51.2, 40,  7.55, 19.5, 329.5))

alphas = test$data.w


library(R0)
library(projections)
set.seed(1)
library(distcrete)
library(epitrix)
mu <- 3.96
sigma <- 4.75
cv <- sigma / mu
params <- gamma_mucv2shapescale(mu, cv)

## Outbreak during 1918 influenza pandemic in Germany)
G=19
CI = matrix(0,G+1,2*3)
for (g in 1:G){
  mGT<-generation.time("gamma", c(mu, sigma))
  t<-estimate.R( xx[g,], mGT, begin=1, end=36, methods=c("EG", "AR"),
                 pop.size=population$pop[g] * 1e6, nsim=1000)
  CI[g,1] =t$estimates$EG$R
  CI[g,2:3] =t$estimates$EG$conf.int
  CI[g,4] =t$estimates$AR$R
  CI[g,5:6] =t$estimates$AR$conf.int
  
}
g=20
t<-estimate.R(apply(xx[c(2:15,19),], 2, sum), mGT, begin=1, end=36, methods=c("EG", "AR"),
              pop.size= sum(population$pop[c(2:15,19)])* 1e6, nsim=1000)
CI[g,1] =t$estimates$EG$R
CI[g,2:3] =t$estimates$EG$conf.int
CI[g,4] =t$estimates$AR$R
CI[g,5:6] =t$estimates$AR$conf.int


#### Trace plot
#stan_trace(test$fit,pars = "tau")
#stan_trace(test$fit,pars = "R0")

#### Chains
df_of_draws <- as.data.frame(test$fit)
chainname = c()
itername = c()
for (i in 1:8){
  chainname<- c(chainname, rep(i,1000))
  itername<- c(itername, 1:1000)
}
df_of_draws$chains = as.factor(chainname)
df_of_draws$iter = itername

Time_vec = seq(as.Date("2020/02/09"), by = "day",length.out = 36)
Time_vec1 = seq(as.Date("2020/01/23"), by = "day",length.out = 36)
si <- distcrete("gamma", shape = params$shape,
                scale = params$scale,
                interval = 1, w = 0)


interval_R0 = 420:438
interval_theta= 21:39


B=1000
G = 19
K = 20
TimeToPredict = 60
Di= 20
Predictions <- array(dim=c(B, G, 36 + TimeToPredict))

for ( b in 1:B){
  #### sample tau
  print(b)
  Predictions[b,,1:36] = xx
  tau= df_of_draws$tau[sample(1:4000,1)]
  for (g in 1:G){
    cbar= df_of_draws[sample(1:4000,1), 1+g]
    R0 = cbar  * tau * Di
    for (t in 37:TimeToPredict){
      if(typesim=="hurdle"){
        p = 0
        z = 0
        while(  p * (1-z) ==0){
          p= rpois(1, 5 + R0 * sum(alphas * (Predictions[b, g, (t-K):(t-1)])));
          z = rbinom(1,1,df_of_draws[sample(1:8000,1),interval_theta[g]]);        
        }
        #p = rnorm(1, 2 * sqrt(3.0/8. + R0 * sum(alphas * (Predictions[b, g, (t-K):(t-1)]))), df_of_draws[sample(1:8000,1), interval_eps[g]]);
        Predictions[b, g, t]  = p * (1-z)  #(p/2)^2-3.0/8 #rpois(1, 1 + R0 * sum(alphas_s * (Predictions[b, g, (t-K):(t-1)])))
      }else{
        Predictions[b, g, t]  = rpois(1, 5 + R0 * sum(alphas * (Predictions[b, g, (t-K):(t-1)])))
      }
      
    }
    
  }
}



interval=18:61
perf <- data.frame(matrix(0,19,15))


g = 1
interval1 = 1:44
dat2 <- c()
for (i in 1:36){
  dat2 <- c(dat2, rep(toString(Time_vec1[i]), xx[g, i]))
}
pred <- project(incidence(dat2), R = CI[g,1], si = si, n_days = 34, n_sim = 1000)
pred_h <- project(incidence(dat2), R = CI[1,1], si = si, n_days = 34, n_sim = 1000)
pred_w <- project(incidence(dat2), R = CI[20,1], si = si, n_days = 34, n_sim = 1000)
pred_us <- project(incidence(dat2), R = CI[19,1], si = si, n_days = 34, n_sim = 1000)
temp_df = data.frame(Y = apply(Predictions[,g,1:length(interval1)],2, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}),
                     time=seq(as.Date("2020/01/23"), by = "day", length.out = length(interval)),
                     type=c(rep('Train',36),rep('Test',length(interval)-36)),
                     L = apply(Predictions[,g,1:length(interval1)],2,FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}),
                     U = apply(Predictions[,g,1:length(interval1)],2,FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}),
                     N= sapply(N[as.numeric(indexes[g]),interval1],
                               FUN=function(x){ifelse(x==0,0,log(x))})
)
temp_df$Y_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_h[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_h[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_h[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                     apply(pred_us[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                     apply(pred_us[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                     apply(pred_us[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
perf[1,] = c(mean(abs(temp_df$Y[37:44] - temp_df$N[37:44])),
              mean(abs(temp_df$Y_EG[37:44] - temp_df$N[37:44])),
              mean(abs(temp_df$Y_EG_hubei[37:44] - temp_df$N[37:44])),
              mean(abs(temp_df$Y_EG_world[37:44] - temp_df$N[37:44])),
              mean(abs(temp_df$Y_EG_us[37:44] - temp_df$N[37:44])),
              mean(temp_df$U[37:44] - temp_df$L[37:44]),
              mean(temp_df$U_EG[37:44] - temp_df$L_EG[37:44]),
              mean(temp_df$U_EG_hubei[37:44] - temp_df$L_EG_hubei[37:44]),
              mean(temp_df$U_EG_world[37:44] - temp_df$L_EG_world[37:44]),
              mean(temp_df$U_EG_us[37:44] - temp_df$L_EG_us[37:44]),
              mean((temp_df$N[37:44]>temp_df$L[37:44])* (temp_df$N[37:44])<temp_df$U[37:44]),
              mean((temp_df$N[37:44]>temp_df$L_EG[37:44])* (temp_df$N[37:44])<temp_df$U_EG[37:44]),
              mean((temp_df$N[37:44]>temp_df$L_EG_hubei[37:44])* (temp_df$N[37:44])<temp_df$U_EG_hubei[37:44]),
              mean((temp_df$N[37:44]>temp_df$L_EG_world[37:44])* (temp_df$N[37:44])<temp_df$U_EG_world[37:44]),
              mean((temp_df$N[37:44]>temp_df$L_EG_us[37:44])* (temp_df$N[37:44])<temp_df$U_EG_us[37:44])
              )

pd <- position_dodge(0.1)  
ggplot(temp_df, aes(x=time, y= Y))+
  geom_point(aes(x=time, y=Y_EG), fill="orchid1", size=1, shape=16)+  
  geom_errorbar(data=temp_df,aes(x = time, ymin=L_EG, ymax=U_EG),
                width=1, size=1, col="plum1",alpha=0.8) +
  geom_errorbar(data=temp_df,aes(x = time, ymin=L, ymax=U),
                width=1, size=1, col="lightgreen",position = pd)  + 
  geom_line() +
  geom_point(size=1, shape=21, fill="white")+
  theme_bw()+
  geom_point(aes(x=time, y=N, col=type),size=3, shape=17)
ggsave(paste("~/Dropbox/contagion_studies/", appendix,"hubei0.pdf",sep=''))



for (g in 2:18){
dat2 <- c()
for (i in 1:36){
  dat2 <- c(dat2, rep(toString(Time_vec[i]), xx[g, i]))
}
pred <- project(incidence(dat2), R = CI[g,1], si = si, n_days = 34, n_sim = 1000)
pred_h <- project(incidence(dat2), R = CI[1,1], si = si, n_days = 34, n_sim = 1000)
pred_w <- project(incidence(dat2), R = CI[20,1], si = si, n_days = 34, n_sim = 1000)
pred_us <- project(incidence(dat2), R = CI[19,1], si = si, n_days = 34, n_sim = 1000)
temp_df = data.frame(Y = apply(Predictions[,g,1:length(interval)],2, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}),
                     time=seq(as.Date("2020/02/09"), by = "day", length.out = length(interval)),
                     type=c(rep('Train',36),rep('Test',length(interval)-36)),
                     L = apply(Predictions[,g,1:length(interval)],2,FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}),
                     U = apply(Predictions[,g,1:length(interval)],2,FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}),
                     N= sapply(N[as.numeric(indexes[g]),interval],
                               FUN=function(x){ifelse(x==0,0,log(x))})
)
temp_df$Y_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred_h[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred_h[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred_h[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_us[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_us[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_us[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))


pd <- position_dodge(0.1)  
ggplot(temp_df, aes(x=time, y= Y))+
  geom_point(aes(x=time, y=Y_EG), fill="orchid1", size=1, shape=16)+  
  geom_errorbar(data=temp_df,aes(x = time, ymin=L_EG, ymax=U_EG),
                width=1, size=1, col="plum1",alpha=0.8) +
  geom_errorbar(data=temp_df,aes(x = time, ymin=L, ymax=U),
                width=1, size=1, col="lightgreen",position = pd)  + 
  geom_line() +
  geom_point(size=1, shape=21, fill="white")+
  theme_bw()+
  geom_point(aes(x=time, y=N, col=type),size=3, shape=17)
  ggsave(paste("~/Dropbox/contagion_studies/",appendix, names[g], ".pdf",sep=''))

  
  perf[g,] = c(mean(abs(temp_df$Y[37:44] - temp_df$N[37:44])),
               mean(abs(temp_df$Y_EG[37:44] - temp_df$N[37:44])),
               mean(abs(temp_df$Y_EG_hubei[37:44] - temp_df$N[37:44])),
               mean(abs(temp_df$Y_EG_world[37:44] - temp_df$N[37:44])),
               mean(abs(temp_df$Y_EG_us[37:44] - temp_df$N[37:44])),
               mean(temp_df$U[37:44] - temp_df$L[37:44]),
               mean(temp_df$U_EG[37:44] - temp_df$L_EG[37:44]),
               mean(temp_df$U_EG_hubei[37:44] - temp_df$L_EG_hubei[37:44]),
               mean(temp_df$U_EG_world[37:44] - temp_df$L_EG_world[37:44]),
               mean(temp_df$U_EG_us[37:44] - temp_df$L_EG_us[37:44]),
               mean((temp_df$N[37:44]>temp_df$L[37:44])* (temp_df$N[37:44])<temp_df$U[37:44]),
               mean((temp_df$N[37:44]>temp_df$L_EG[37:44])* (temp_df$N[37:44])<temp_df$U_EG[37:44]),
               mean((temp_df$N[37:44]>temp_df$L_EG_hubei[37:44])* (temp_df$N[37:44])<temp_df$U_EG_hubei[37:44]),
               mean((temp_df$N[37:44]>temp_df$L_EG_world[37:44])* (temp_df$N[37:44])<temp_df$U_EG_world[37:44]),
               mean((temp_df$N[37:44]>temp_df$L_EG_us[37:44])* (temp_df$N[37:44])<temp_df$U_EG_us[37:44])
  )
  
  
}

g=19
interval=18:ncol(N)
dat2 <- c()
for (i in 1:36){
  dat2 <- c(dat2, rep(toString(Time_vec[i]), xx[g, i]))
}
pred <- project(incidence(dat2), R = CI[g,1], si = si, n_days = 34, n_sim = 1000)
pred_h <- project(incidence(dat2), R = CI[1,1], si = si, n_days = 34, n_sim = 1000)
pred_w <- project(incidence(dat2), R = CI[20,1], si = si, n_days = 34, n_sim = 1000)
pred_us <- project(incidence(dat2), R = CI[19,1], si = si, n_days = 34, n_sim = 1000)

temp_df = data.frame(Y = apply(Predictions[,g,1:length(interval)],2, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}),
                     time =seq(as.Date("2020/02/09"), by = "day", length.out = length(interval)),
                     type=c(rep('Train',36),rep('Test',5)),
                     L = apply(Predictions[,g,1:length(interval)],2,FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025, na.rm = TRUE)))}),
                     U = apply(Predictions[,g,1:length(interval)],2,FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975, na.rm = TRUE)))}),
                     N= sapply(apply(N[index_us,interval],2,sum),
                               FUN=function(x){ifelse(x==0,0,log(x))})
)
temp_df$Y_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                  apply(pred[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_h[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_h[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_hubei  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_h[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_world  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                        apply(pred_w[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))
temp_df$Y_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                     apply(pred_us[1:8], 1, FUN=function(x){ifelse(mean(x)==0, 0, log(mean(x)))}))
temp_df$U_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                     apply(pred_us[1:8], 1, FUN=function(x){ifelse(quantile(x,0.975)==0, 0, log(quantile(x,0.975)))}))
temp_df$L_EG_us  = c(sapply(xx[g,], FUN=function(x){ifelse(x==0, 0, log(x))}), 
                     apply(pred_us[1:8], 1, FUN=function(x){ifelse(quantile(x,0.025)==0, 0, log(quantile(x,0.025)))}))

ggplot(temp_df, aes(x=time, y= Y)) +  
  geom_point(aes(x=time, y=Y_EG), fill="orchid1", size=1, shape=16)+  
  geom_errorbar(data=temp_df,aes(x = time, ymin=L_EG, ymax=U_EG),
                width=1, size=1, col="plum1",alpha=0.8) +
  geom_errorbar(data=temp_df,aes(x = time, ymin=L, ymax=U),
                width=1, size=1, col="lightgreen",position = pd)  + 
  geom_line() +
  geom_point(size=1, shape=21, fill="white")+
  theme_bw()+
  geom_point(aes(x=time, y=N, col=type),size=3, shape=17)
ggsave(paste("~/Dropbox/contagion_studies/", appendix, "USA.pdf",sep=''))


perf[g,] = c(mean(abs(temp_df$Y[37:44] - temp_df$N[37:44])),
             mean(abs(temp_df$Y_EG[37:44] - temp_df$N[37:44])),
             mean(abs(temp_df$Y_EG_hubei[37:44] - temp_df$N[37:44])),
             mean(abs(temp_df$Y_EG_world[37:44] - temp_df$N[37:44])),
             mean(abs(temp_df$Y_EG_us[37:44] - temp_df$N[37:44])),
             mean(temp_df$U[37:44] - temp_df$L[37:44]),
             mean(temp_df$U_EG[37:44] - temp_df$L_EG[37:44]),
             mean(temp_df$U_EG_hubei[37:44] - temp_df$L_EG_hubei[37:44]),
             mean(temp_df$U_EG_world[37:44] - temp_df$L_EG_world[37:44]),
             mean(temp_df$U_EG_us[37:44] - temp_df$L_EG_us[37:44]),
             mean((temp_df$N[37:44]>temp_df$L[37:44])* (temp_df$N[37:44])<temp_df$U[37:44]),
             mean((temp_df$N[37:44]>temp_df$L_EG[37:44])* (temp_df$N[37:44])<temp_df$U_EG[37:44]),
             mean((temp_df$N[37:44]>temp_df$L_EG_hubei[37:44])* (temp_df$N[37:44])<temp_df$U_EG_hubei[37:44]),
             mean((temp_df$N[37:44]>temp_df$L_EG_world[37:44])* (temp_df$N[37:44])<temp_df$U_EG_world[37:44]),
             mean((temp_df$N[37:44]>temp_df$L_EG_us[37:44])* (temp_df$N[37:44])<temp_df$U_EG_us[37:44])
)

### Quantify performance



##### Plot quatieis
library(reshape2)
dummy2 <- data.frame(variable = sapply(1:19,
                                       FUN=function(x){paste("cbar[",toString(x), "]", sep="")}),
                     Z = apply(df_of_draws[,2:20],2,FUN=function(x){quantile(x,0.025)}),
                     ZZ = apply(df_of_draws[,2:20],2,FUN=function(x){quantile(x,0.975)}))

dummy3 <- data.frame(variable = sapply(1:19,
                                       FUN=function(x){paste("R0[",toString(x), "]", sep="")}),
                     Z = apply(df_of_draws[,interval_R0],2,FUN=function(x){quantile(x,0.025)}),
                     ZZ = apply(df_of_draws[,interval_R0],2,FUN=function(x){quantile(x,0.975)}))



melt(df_of_draws[,interval_R0]) %>%
  ggplot( aes(x=value, col=variable)) +  theme_bw()+
  geom_density(fill="skyblue", color="black", alpha=0.8) + 
  facet_wrap(~ variable, , scales = "free", nrow=4) +
  geom_vline(data = dummy3, aes(xintercept = Z), color="gray")+
  geom_vline(data = dummy3, aes(xintercept =ZZ), color="gray")
ggsave(paste("~/Dropbox/contagion_studies/", appendix, "R0.pdf",sep=''))


melt(df_of_draws[,1:20]) %>%
  ggplot( aes(x=value, col=variable)) +  theme_bw()+
  geom_density(fill="skyblue", color="black", alpha=0.8) + 
  facet_wrap(~ variable, , scales = "free", nrow=4) +
  geom_vline(data = dummy2, aes(xintercept = Z), color="gray")+
  geom_vline(data = dummy2, aes(xintercept =ZZ), color="gray")
ggsave(paste("~/Dropbox/contagion_studies/", appendix, "cbar.pdf",sep=''))

