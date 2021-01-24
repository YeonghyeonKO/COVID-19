#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

# Analysis Period( ~ 8/31)
max_date <- set_date("2020/9/1") 

#### Data import ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

continent <- as.data.frame(unique(df[,c("countriesAndTerritories","continentExp")]))
colnames(continent)<-c("country","continent")
continent$country<-as.character(continent$country)

# preprocess data
df_sum = preprocessing_data()

#### Analysis ####
main_country = c("South_Korea","Japan","United_Kingdom","France",
                  "Germany","Italy","Iran","China","United_States_of_America",
                  "Canada","Spain","Brazil","Russia","India")

df_result = derivative_analysis(Country = country,gkf_bandwidth = 14)

######### segmented poisson ###############

# 한국 피팅
segPoisson("South_Korea",break_point = c(80,130,200))


# 주요 국가
break_sum = df_result[,13:17]
rownames(break_sum) = df_result$country
for(i in main_country){
  temp = which(rownames(break_sum)==i)
  break_point = c()
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  segPoisson(Country=i,break_point = break_point)
  print(break_point)
}


