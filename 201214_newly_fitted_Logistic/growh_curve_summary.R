#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)
library(gridExtra)
library(tidyverse)

getwd()
setwd("~/Documents/GitHub/COVID-19/201214_fitting_Korea")
coef_result <- read_csv("coef_result.csv")


# Analysis Period( ~ 8/31)
max_date <- set_date("2020/8/31") 
# Analysis Period( ~ 9/15)
max_date <- set_date("2020/9/15") 
# Analysis Period( ~ 9/30)
max_date <- set_date("2020/9/30") 
# Analysis Period( ~ 10/15)
max_date <- set_date("2020/10/15")
# Analysis Period( ~ 10/31)
max_date <- set_date("2020/10/31")

# For Taiwan

# derivative_analysis(Country = "Taiwan",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)
# derivative_analysis(Country = "South_Korea",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)
# derivative_analysis(Country = "India",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)
# derivative_analysis(Country = "UK",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)
# derivative_analysis(Country = "France",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)
ggplot() +
  geom_line(aes(df_sum$Date, df_sum$Taiwan), col=2) +
  geom_line(aes(df_sum$Date, df_sum$South_Korea), col=3)


#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data (until 20/08/31)
df_sum = preprocessing_data()[1:245,]
# preprocess data (until 20/9/15)
df_sum = preprocessing_data()[1:260,]
# preprocess data (until 20/9/30)
df_sum = preprocessing_data()[1:275,]
# preprocess data (until 20/10/15)
df_sum = preprocessing_data()[1:290,]
# preprocess data (until 20/10/31)
df_sum = preprocessing_data()[1:306,]

ts.plot(df_sum["Paraguay"])
# Simple Moving average (window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}

# Peak Detectiong Analysis
df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_break = 50,save_image = TRUE,save_excel = TRUE)

derivative_analysis(Country = "South_Korea",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)



#### segmented logistic ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #

par(mfrow=c(1,1))
break_sum = df_result[,13:17]
rownames(break_sum) = df_result$country
country = df_result$country


#### separately fit logistic models #### 
coef_sep_Logi <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Logi) <- country
colnames(coef_sep_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", "a2_Logi", "b2_Logi", "c2_Logi",
                             "startpoint", "breakpoint")

country = c("South_Korea", "India", "France", "UK", "Taiwan")

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
      coef_temp = data.frame()
      coef_temp = segLogistic_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
                                       max_iter = 1000)
      if(!is.null(coef_temp)){
        coef_sep_Logi[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_sep_Logi[i, 4:6] <- coef_temp[2,]
        }
      }
      coef_sep_Logi[i, 7:8] <- break_sum[i, 1:2]
    }
  )
}


for(i in 1:length(country)){
  if(!is.na(coef_sep_Logi$breakpoint[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_sep_Logi$a2_Logi[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_sep_Logi$a1_Logi[i])){
      m = 1
    }else{
      m = 0
    }
  }

  # Retry fitting 
  i_country = country[i]
  if(M-m>0){
    tryCatch(
      expr={
        temp = which(rownames(break_sum)==i_country)
        break_point = break_sum$break1[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        coef_temp = segLogistic_separate(Country=i_country,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                         max_iter = 1000)
      }
    )
    
    # Case by Case  
    if(M==2&m==0){
      if(!is.null(coef_temp)){
        coef_sep_Logi[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_sep_Logi[i, 4:6] <- coef_temp[2,]
        }
      }
    }
    
    if(M==2&m==1){
      if(!is.null(coef_temp)){
        if(length(coef_temp[,1])==2){
          coef_sep_Logi[i, 4:6] <- coef_temp[2,]
        }
      }
    }
    
    if(M==1&m==0){
      if(!is.null(coef_temp)){
        coef_sep_Logi[i, 1:3] <- coef_temp[1,]
      }
    }
  }
}

for(i in 1:length(country)){
  
}


#### segmented Bertalanffy ####
coef_sep_Ber <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Ber) <- country
colnames(coef_sep_Ber) <- c("a1_Ber", "b1_Ber", "c1_Ber", "a2_Ber", "b2_Ber", "c2_Ber",
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
      coef_temp = data.frame()
      coef_temp = segBertalanffy_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                       max_iter = 1000)
      if(!is.null(coef_temp)){
        coef_sep_Ber[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_sep_Ber[i, 4:6] <- coef_temp[2,]
        }
      }
      coef_sep_Ber[i, 7:8] <- break_sum[i, 1:2]
    }
  )
}


#### segmented Gompertz ####
coef_sep_Gom <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Gom) <- country
colnames(coef_sep_Gom) <- c("a1_Gom", "b1_Gom", "c1_Gom", "a2_Gom", "b2_Gom", "c2_Gom",
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
      coef_temp = data.frame()
      coef_temp = segGompertz_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                       max_iter = 1000)
      if(!is.null(coef_temp)){
        coef_sep_Gom[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_sep_Gom[i, 4:6] <- coef_temp[2,]
        }
      }
      coef_sep_Gom[i, 7:8] <- break_sum[i, 1:2]
    }
  )
}


for(i in 1:length(country)){
  if(!is.na(coef_sep_Gom$breakpoint[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_sep_Gom$a2_Gom[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_sep_Gom$a1_Gom[i])){
      m = 1
    }else{
      m = 0
    }
  }
  
  # Retry fitting 
  i_country = country[i]
  if(M-m>0){
    tryCatch(
      expr={
        temp = which(rownames(break_sum)==i_country)
        break_point = break_sum$break1[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        coef_temp = segGompertz_separate(Country=i_country,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                         max_iter = 1000)
      }
    )
    
    # Case by Case  
    if(M==2&m==0){
      if(!is.null(coef_temp)){
        coef_sep_Gom[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_sep_Gom[i, 4:6] <- coef_temp[2,]
        }
      }
    }
    
    if(M==2&m==1){
      if(!is.null(coef_temp)){
        if(length(coef_temp[,1])==2){
          coef_sep_Gom[i, 4:6] <- coef_temp[2,]
        }
      }
    }
    
    if(M==1&m==0){
      if(!is.null(coef_temp)){
        coef_sep_Gom[i, 1:3] <- coef_temp[1,]
      }
    }
  }
}

#### Result ####
coefficient_result = cbind(coef_sep_Logi[,1:6],coef_sep_Gom)
write.csv(coefficient_result,"coef_result.csv")

# number of fitted parameter pairs
length(which(!is.na(coef_sep_Logi$a1_Logi)))
length(which(!is.na(coef_sep_Logi$a2_Logi)))

length(which(!is.na(coef_sep_Gom$a1_Gom)))
length(which(!is.na(coef_sep_Gom$a2_Gom)))

# scatter plot : logistic vs gompertz
par(mfrow=c(1,1))
plot(log(coefficient_result$a1_Logi),log(coefficient_result$a1_Gom))
plot(log(coefficient_result$a2_Logi),log(coefficient_result$a2_Gom))

plot(log(coefficient_result$b1_Logi),log(coefficient_result$b1_Gom))
plot(log(coefficient_result$b2_Logi),log(coefficient_result$b2_Gom))

plot(log(coefficient_result$c1_Logi),log(coefficient_result$c1_Gom))
plot(log(coefficient_result$c2_Logi),log(coefficient_result$c2_Gom))



