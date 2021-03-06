#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

# Analysis Period( ~ 9/22)
max_date <- set_date("2020/9/22") 

#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data
df_sum = preprocessing_data()

# Peak Detectiong Analysis
df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)

derivative_analysis(Country = "Germany",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)


#### segmented logistic ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #

par(mfrow=c(1,1))
break_sum = df_result[,13:17]
rownames(break_sum) = df_result$country

for(i in main_country){
  temp = which(rownames(break_sum)==i)
  break_point = break_sum$break1[temp]
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  segLogistic(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
              max_iter = 1000)
}


#### segmented logistic - daily cases ####
for(i in main_country){
  temp = which(rownames(break_sum)==i)
  break_point = break_sum$break1[temp]
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  segLogistic_daily(Country=i,break_point = break_point,
                    save_image = TRUE,prediction_plot = TRUE,
                    max_iter = 1000)
}


#### separately fit logistic models #### 
for(i in main_country){
  temp = which(rownames(break_sum)==i)
  break_point = break_sum$break1[temp]
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  segLogistic_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
              max_iter = 1000)
}



#### country classification by # of waves ####
# number of waves is counted by # of breakpoint.
sum(table(df_result$break_number))
table(df_result$break_number)



