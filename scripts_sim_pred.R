#library(tidyverse, lib="~/R_libs")
print(R.version)
library(StanHeaders)
library(ggplot2)
library(dplyr)
library(rstan)
#load("/scratch/users/cdonnat/contagion/CI.RData")
#setwd("/scratch/users/cdonnat/contagion")
args = commandArgs(trailingOnly=TRUE)  ### Pass the seed + name of saved file where we want to save the results.
j = as.numeric(args[1])
session = as.numeric(args[2])
r0 = as.numeric(args[3])
level_cont = 0.01 * as.numeric(args[4])
typerd = args[5]
print(typerd)
if (typerd == "cauchy"){
print("cauchy")
load("/scratch/users/cdonnat/contagion/params_Cauchy.RData")
}else{
load("results_contagion_switch_normal_gammaR0_rndom.RData")
}
appendix= paste(typerd, "_predictions_scenario_", args[1],"r0_", toString(r0),"_levelcont_", as.numeric(args[4]), "_session", toString(session), sep="")

df_of_draws = as.data.frame(test$fit)
alphas = test$data.w
print(alphas)

chainname = c()
itername = c()
for (i in 1:8){
  chainname<- c(chainname, rep(i,1000))
  itername<- c(itername, 1:1000)
}
df_of_draws$chains = as.factor(chainname)
df_of_draws$iter = itername
print(df_of_draws[1:3,1:4])
df_of_draws = df_of_draws[df_of_draws$chains %in% c(1,3,4,5,7,8),]
print(dim(df_of_draws))
interval_cbar = 1104:1122
interval_R0 = 1085:1103


today <- as.Date("2020-03-18")
end_fit <- as.Date("2020-03-15")
data <- read.csv2("time_series_19-covid-Confirmed.csv", sep=",", header=TRUE)
colnames(data)[5:ncol(data)] <- sapply(seq(as.Date("2020/1/22"), by = "day", to=today), toString)
data_recovered <- read.csv2("time_series_19-covid-Recovered.csv",sep=",", header=TRUE)
colnames(data_recovered)[5:ncol(data)] <- sapply(seq(as.Date("2020/1/22"), by = "day", to=today), toString)
data_death <- read.csv2("time_series_19-covid-Deaths.csv",
			                        sep=",", header=TRUE)
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



interval_fit1=1:36
interval_fit=18:53
xx<- N[c(index_hubei), interval_fit1]
xx <- rbind(xx, N[c(index_hubei,index_hk, index_guizhou,index_fr, index_it, index_germany, index_spain,
                    index_sw, index_uk,index_iran,
                    index_japan, index_thai, index_sing, index_korea, index_cal,
                    index_wash, index_ny), interval_fit])
xx<- rbind(xx,apply(N[index_us, interval_fit],2, sum) )
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

population$pop = population$pop * 10^6
population$`0-24` = rep(0,19)
population$`25-54` = rep(0,19)
population$`55-64` = rep(0,19)
population$`65+` = rep(0,19)
population[1:4,2:ncol(population)]=c(17.26 + 11.48,46.81,12.08,12.34)
population[3,2:ncol(population)]=c(12.81 + 8.81, 42.66, 17.24,18.48)
population[5,2:ncol(population)]=c(18.36 + 11.88, 36.83, 12.47,20.46)
population[6,2:ncol(population)]=c(13.45 + 9.61, 40.86, 14, 22.08)
population[7,2:ncol(population)]=c(12.89 + 9.81, 38.58, 15.74, 22.99)
population[8,2:ncol(population)]=c(15.02 + 9.9, 43.61, 12.99, 18.49)
population[9,2:ncol(population)]=c(15.34 + 10.39, 42.05, 13.48, 18.73)
population[10,2:ncol(population)]=c(17.63 + 11.49, 39.67, 12.73, 18.48)
population[11,2:ncol(population)]=c(24.11 + 13.36, 48.94, 7.72, 5.87)
population[12,2:ncol(population)]=c(12.49 + 9.47, 36.8, 12.06, 29.18)
population[13,2:ncol(population)]=c(16.45 + 13.02, 45.69, 13.01, 11.82)
population[14,2:ncol(population)]=c(12.8 + 15.01, 50.73, 10.58, 10.89)
population[15,2:ncol(population)]=c(12.77 + 11.18, 44.66, 15.47, 15.92)
population[16,2:ncol(population)]=c(22.9+ 14.2, 44.4, 7.7, 10.8)
population[17,2:ncol(population)]=c(21.3+ 13.9, 45.2, 8.4, 11.2 )
population[18,2:ncol(population)]=c(20.39+ 13.87, 38.92, 12.86, 16.85)
population[19,2:ncol(population)]=c(18.46+ 12.91, 38.92, 12.86, 16.85)





row.names(population) =names
hosp_rate_l = 0.01* c(1.6, 25/35 * 14.3 + 10 /35 * 21.2, 20.5, 1/3 * 28.6 + 1/3 * 30.5 + 1/3 * 31.3)
hosp_rate_u = 0.01* c(2.5, 25/35 * 20.8 + 10 /35 * 28.3, 30.1, 1/3 * 43.5 + 1/3 * 58.7 + 1/3 * 70.3)
icu_rate_l =  0.01* c(100/2500, 25/35 * 2.0 + 10 /35 * 5.4,   4.7, 1/3 * 8.1  + 1/3 * 10.5 + 1/3 * 6.3)
#icu_rate_u = c(0.01,  25/35 * 4.2 + 10 /35 * 10.4, 11.2, 1/3 * 18.2, 18.8+ 1/3* 31.0 + 1/3* 29.0)
death_rate_l = 0.01* c(100/2500, 25/35 * 0.1 + 10 /35 * 0.5, 1.4,   1/3 * 2.7+ 1/3* 4.3 + 1/3* 10.4)
death_rate_u = 0.01* c(100/2500,25/35 * 0.2 + 10 /35 * 0.8, 2.6,    1/3 * 4.9+ 1/3* 10.5+ 1/3* 27.3)

# 

### Lowering cbar to 1 everywhere
B= 500
G = 19
TimeToPredict = 200
Di= 20
K=20

mu_hosp = 22 - 10
IQR_hosp = 7
sigma_hosp = IQR_hosp/(27/20)

mu_icu = 14.5  - 8
IQR_icu = (19-12)
sigma_icu = IQR_icu/(27/20)

mu_deaths = 18.5 - 8
IQR_death = (22-15)
sigma_deaths = IQR_death /(27/20)

level = c(c(1,0.2,0.4,0.5,0.75, 0.85,0.9,1.0) , rep(1.0, 9))
Predictions2 <- array(0,dim=c(B,  G, TimeToPredict))
Predictions3 <- array(0,dim=c(B,  G, TimeToPredict))
Deaths <- array(0,dim=c(B, G, TimeToPredict))
Hospitals <- array(0,dim=c(B, G, TimeToPredict))
ICU <- array(0,dim=c(B, G, TimeToPredict))

mortality <- data_death[indexes,66]/data[indexes,66]
mortality[19] = sum(data_death[index_us,66])/sum(data[index_us,66])
#print(data_death[indexes,c(1:2,64:66)])  
#print(data[indexes,c(1:2,64:66)])
#print(data_death[indexes,66])
mean_R0 = apply(df_of_draws[,interval_R0], 2, mean) 
print(colnames(df_of_draws)[interval_R0])
print((df_of_draws)[1:10,interval_R0])

intervals_r <- list()
for (g in 1:G){
  intervals_r[[g]] <- seq(from=400 +g, to= 400 + G*36, by =G)
}

print(colnames(df_of_draws)[interval_cbar])
print(colnames(df_of_draws)[interval_R0])
print(colnames(df_of_draws)[intervals_r[[5]]])

for ( b in 1:B){
    #### sample tau
    print(b)
    Predictions2[b,,1:36] = xx
    Predictions3[b,,1:36] = xx
   # Deaths[b,,1:36] = rbind(D[index_hubei, interval_fit1],D[indexes[1:18],interval_fit],
   #                         apply(D[index_us, interval_fit],2, sum) )
    tau= df_of_draws$tau[sample(1:nrow(df_of_draws),1)]
    if (j== length(level)-1){tau = tau*0.5}
    if (j== length(level)){tau = tau*0.75}
    alphas_s = alphas
    for (g in c(5:10,13,16:19)){
      print(paste("g = ", toString(g)))
      if (j==1){
        cbar =1
      }else{
        cbar= df_of_draws[sample(1:nrow(df_of_draws),1), 1+g]  * level[j]
        #print(cbar)
      }
      if(r0>20){
          R0 = cbar  * tau * Di
      }else{
          R0 = CI[r0,1] * level[j]
      }
      inc_days=1
      t = 37
      nb_hospitals = 0
      while ( (t<TimeToPredict)&& (nb_hospitals< 0.01 * population$pop[g] )){
          #print(paste("t = ", toString(t)))          
          if(r0>20){
              R0 =  df_of_draws[sample(1:nrow(df_of_draws),1), sample(intervals_r[[g]],1)]  * level[j]  #* df_of_draws$tau[sample(1:nrow(df_of_draws),1)] * Di
              R02 =  mean_R0[g] * level[j]
           }else{
              R0 = ifelse(r0==20,CI[r0,1] * level[j], df_of_draws[sample(1:nrow(df_of_draws),1), 1+r0]  * level[j]* df_of_draws$tau[sample(1:nrow(df_of_draws),1)] * Di)
              R02 =  ifelse(r0==20, mean(mean_R0) * level[j], mean_R0[r0] * level[j])
          }
          if(inc_days > 1 && j == 9){ R0 = R0 * level_cont; R02 = R02 * level_cont;} # ifelse(r0>20,   level_cont, level_cont)} 
          if(inc_days > 2  && j == 10){ R0 = R0 * level_cont; R02 = R02 * level_cont;} #ifelse(r0>20,level_cont, CI[r0,1]  * level_cont)}
	  if(inc_days > 3 && j == 11){ R0 = R0 * level_cont; R02 = R02 * level_cont;} #ifelse(r0>20,level_cont, CI[r0,1]  * level_cont)}
	  if(inc_days > 4 && j == 12){ R0 = R0 * level_cont; R02 = R02 * level_cont; } #ifelse(r0>20,level_cont, CI[r0,1]  * level_cont)}#  * tau * Di *level_cont}
	  if(inc_days > 5 && j == 13){ R0 = R0 * level_cont; R02 = R02 * level_cont;} #ifelse(r0>20,level_cont, CI[r0,1]  * level_cont)}#*level_cont}
	  if(inc_days > 6 && j == 14){ R0 = R0 * level_cont; R02 = R02 * level_cont;}# ifelse(r0>20,level_cont, CI[r0,1]  * level_cont)}
	  if(inc_days > 7 && j == 15){ R0 = R0 * level_cont; R02 = R02 * level_cont;} #ifelse(r0>20, level_cont, CI[r0,1]  * level_cont)}
          #print(paste("t = ", toString(t), " R0 = ", toString(R0)))          
          p = rnorm(1, 2 * sqrt(3.0/8. + R0 * sum(alphas* (Predictions2[b, g, (t-K):(t-1)]))), 1);
          #print(paste("t = ", toString(t)))          
          p3 = rnorm(1, 2 * sqrt(3.0/8. + R02 * sum(alphas* (Predictions3[b, g, (t-K):(t-1)]))), 1);
          #print((Predictions2[b, g, (t-K):(t-1)]))
	  #print(alphas_s)
          #print(c(R02, R0,"yo",sum(alphas_s* (Predictions2[b, g, (t-K):(t-1)])),p))
          Predictions2[b, g, t]  = max(c(0,floor((p/2)^2-3.0/8.0)))
          Predictions3[b, g, t]  = max(c(0,floor((p3/2)^2-3.0/8.0)))
          #### Each of these gets asigned an hospitalization rate
          #### Upon hospitalization a time of death or recover
          p_host = 0.01* rgamma(1,sum(population[g,2:ncol(population)]* hosp_rate_l) * 10, 10)
          p_icu = min(c(1.0,0.01* rgamma(1,sum(population[g,2:ncol(population)]* icu_rate_l)*10,10)/p_host))  # p(ICU and hos  )/ p(hosp)
          p_mort = min(c(0.01* rgamma(1,sum(population[g,2:ncol(population)]* death_rate_l)*10,10)/p_icu, 1.0))
          #print(c(Predictions2[b, g, t], p_host, p_icu, p_mort))          
          nb_hospitals = ifelse( Predictions2[b,g, t]>0, rbinom(1, Predictions2[b,g, t], p_host), 0)

          nb_icu = ifelse(nb_hospitals==0, 0, rbinom(1, nb_hospitals, p_icu))
          nb_case_fatality =ifelse(nb_icu==0, 0, rbinom(1, nb_icu, p_mort))
          #print(c(t, Predictions2[b,g, t],nb_icu, nb_case_fatality, nb_hospitals,0.01 * population$pop[g], p_host))
          if((nb_hospitals>0) && (nb_hospitals-nb_icu>0)){
            when_discharged = data.frame(patients = 1:(nb_hospitals-nb_icu),
                                         day= sapply(rnorm(nb_hospitals-nb_icu, mu_hosp, sigma_hosp), FUN=function(x){ifelse(x<0, 0, floor(x))}))%>%count(day)
            for (i in 1:nrow(when_discharged)){
              Hospitals[b,g, t:min(c(TimeToPredict, t+ when_discharged[i,]$day))] = Hospitals[b, g, t:min(c(TimeToPredict, t+ when_discharged[i,]$day))] + when_discharged[i,]$n  ;
            }
          }
          if((nb_icu>0) && (nb_icu-nb_case_fatality>0)){
            when_discharge_icu = data.frame(patients = 1:(nb_icu-nb_case_fatality),
                                            day= sapply(rnorm(nb_icu-nb_case_fatality, mu_icu, sigma_icu), FUN=function(x){ifelse(x<0, 0, floor(x))}))%>%count(day)
            for (i in 1:nrow(when_discharge_icu)){
              Hospitals[b,g, t:min(c(TimeToPredict,(t+ when_discharge_icu[i,]$day)))] = Hospitals[b, g, t:min(c(TimeToPredict,(t+ when_discharge_icu[i,]$day)))] + when_discharge_icu[i,]$n  ;
              ICU[b,g, t:min(c(TimeToPredict,(t+ when_discharge_icu[i,]$day)))] = ICU[b, g, t:min(c(TimeToPredict,(t+ when_discharge_icu[i,]$day)))] + when_discharge_icu[i,]$n  ;
            }
            
          }
          if(nb_case_fatality>0){
            when_dies = data.frame(patients = 1:(nb_case_fatality),
                                   day= sapply(rnorm(nb_case_fatality, mu_deaths, sigma_deaths), FUN=function(x){ifelse(x<0, 0, floor(x))}))%>%count(day)         
            for (i in 1:nrow(when_dies)){
              Hospitals[b, g, t:min(c(TimeToPredict,(t+ when_dies[i,]$day)))] = Hospitals[b, g, t:min(c(TimeToPredict,(t+ when_dies[i,]$day)))] + when_dies[i,]$n  ;
              ICU[b, g, t:min(c(TimeToPredict,(t+ when_dies[i,]$day)))] = ICU[b, g, t:min(c(TimeToPredict,(t+ when_dies[i,]$day)))] + when_dies[i,]$n  ;
              Deaths[b, g, min(c(TimeToPredict,(t+ when_dies[i,]$day)))] = Deaths[b, g, min(c(TimeToPredict,(t+ when_dies[i,]$day)))] + when_dies[i,]$n  ;
            }          
          }
          inc_days = ifelse(inc_days ==7, 1, inc_days+1)
          t = t + 1 
          #print(paste("t = ", toString(t), "inc= ", toString(inc_days)))    
    }
    if (t<TimeToPredict){
             Hospitals[b,g, t] = NA
             ICU[b,g, t] = NA
             Deaths[b,g, t] = NA
             Predictions2[b,g, t] = NA
    }
  }
  
}
test= list(Deaths = Deaths, ICU = ICU, Hospitals = Hospitals, Predictions = Predictions2, Predictions_static = Predictions3)
save(test,file= paste(appendix, ".RData", sep=""))








