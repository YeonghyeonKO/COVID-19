#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)
library(lubridate)

# Analysis Period( ~ 9/18)
max_date <- set_date("2020/9/18") 

#### Data import ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

continent <- as.data.frame(unique(df[,c("countriesAndTerritories","continentExp")]))
colnames(continent)<-c("country","continent")
continent$country<-as.character(continent$country)

population <- as.data.frame(unique(df[,c("countriesAndTerritories","popData2019")]))
colnames(population)<-c("country","popData2019")
population$country<-as.character(population$country)

# preprocess data
df_sum = preprocessing_data()

#### Analysis ####
main_country = c("South_Korea","Japan","United_Kingdom","France",
                 "Germany","Italy","Iran","China","United_States_of_America",
                 "Canada","Spain","Brazil","Russia","India")

df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)

##### segmented poisson ####
par(mfrow=c(1,1))

# 한국 피팅
segPoisson("South_Korea",break_point = c(50,75,130,220))


# 주요 국가 - 중국과 같이 종식된 나라는 끝맺음 구간이 필요
break_sum = df_result[,13:17]
rownames(break_sum) = df_result$country
for(i in main_country){
  temp = which(rownames(break_sum)==i)
  break_point = df_result$break1[temp]
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  
  segPoisson(Country=i,break_point = break_point,save_image = FALSE)
  print(break_point)
}


#### segmented logistic ####

## 주요 국가 - 중국과 같이 종식된 나라는 끝맺음 구간이 필요
break_sum = df_result[,13:17]
rownames(break_sum) = df_result$country
for(i in main_country){
  temp = which(rownames(break_sum)==i)
  break_point = df_result$break1[temp]
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  segLogistic(Country=i, break_point = break_point)
}

# segLogistic(Country="Iran", break_point = c(57, 122, 249))
# segLogistic(Country="Spain", break_point = c(61, 164))
segLogistic(Country = "France", break_point = c(61, 170))
