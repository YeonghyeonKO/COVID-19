#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)


# Analysis Period( ~ 11/30)
max_date <- set_date("2020/11/30") 
getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201226_Validation_of_models_for_each_country")


#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
# df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
# write.csv(df, "covid_data.csv")
df <- read.csv("covid_data.csv")


## generate vector of country names
# country <- as.character(unique(df$countriesAndTerritories))
# country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]
pop <- na.omit(df[,c(8,11)] %>% distinct())
country <- pop %>% filter(popData2019 >= 1000000) %>% select(countriesAndTerritories)
country <- as.vector(country[[1]])


## preprocess data (until 20/11/30)
df_sum = preprocessing_data()[1:336,]


## Simple Moving average (window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}


## Peak Detectiong Analysis
df_result = derivative_analysis(Country=country, gkf_bandwidth=14, 
                                first_break=0, save_image=TRUE, save_excel=TRUE)



## Segmented Logistic Model
par(mfrow=c(1,1))
break_sum = cbind(df_result$first_day,df_result[,11:15])
rownames(break_sum) = df_result$country
country = df_result$country

coef_sep_Logi <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Logi) <- country
colnames(coef_sep_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", "a2_Logi", "b2_Logi", "c2_Logi",
                             "startpoint", "breakpoint")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = break_sum$break1[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_sep_Logi[i, 1:3] <- segLogistic_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                                    max_iter = 1000)[1,]
      coef_sep_Logi[i, 4:6] <- segLogistic_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                                    max_iter = 1000)[2,]  
      coef_sep_Logi[i, 7:8] <- break_sum[i, 1:2]
      
    },
    error = function(e){
      next
    },
    finally = next
  )
}



A <- coef_sep_Logi

write.csv(coef_sep_Logi, "201130_Separated_Logistic_coefficients.csv")


