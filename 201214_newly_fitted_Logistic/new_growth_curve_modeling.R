#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)


# Analysis Period( ~ 8/31)
max_date <- set_date("2020/8/31") 


#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
# df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
# write.csv(df, "covid_data.csv")
df <- read.csv("covid_data.csv")

setwd("~/Documents/GitHub/COVID-19/201214_newly_fitted_Logistic")
getwd()

# generate vector of country names
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data (until 20/08/31)
df_sum = preprocessing_data()[1:245,]

# Peak Detectiong Analysis
df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_day_criteria = 50,save_image = TRUE,save_excel = FALSE)

# Simple Moving average (window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}


#   <Segmented Growth Curve model fitting>  
# - use break_sum data.frame which is union of first_break and breakpoint
# - 1st break_point is not real breakpoint but starting point of model fitting            
# - from 2nd point, real breakpoint is counted                                           
# - (ex) break_point=c(50, 70, 120) => 50 is the date at which model fitting is started    
# -      so the # of real breakpoint is 2 (70, 120)                                             
par(mfrow=c(1,1))
break_sum = cbind(df_result$first_day,df_result[,11:15])
rownames(break_sum) = df_result$country
country = df_result$country

#### separately fit logistic models ####
coef_sep_Logi <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Logi) <- country
colnames(coef_sep_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", "a2_Logi", "b2_Logi", "c2_Logi",
                             "startpoint", "breakpoint")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = first_day[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_temp = data.frame()
      coef_temp = segLogistic(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
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

# retry fitting for unfitted countries
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
        break_point = first_day[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        coef_temp = segLogistic(Country=i_country,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
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

#### segmented Bertalanffy ####
coef_sep_Ber <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Ber) <- country
colnames(coef_sep_Ber) <- c("a1_Ber", "b1_Ber", "c1_Ber", "a2_Ber", "b2_Ber", "c2_Ber",
                            "startpoint", "breakpoint")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = first_day[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_temp = data.frame()
      coef_temp = segBertalanffy(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
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

# retry fitting for unfitted countries
for(i in 1:length(country)){
  if(!is.na(coef_sep_Ber$breakpoint[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_sep_Ber$a1_Ber[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_sep_Ber$a2_Ber[i])){
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
        break_point = first_day[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        coef_temp = segBertalanffy(Country=i_country,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                   max_iter = 1000)
      }
    )
    
    # Case by Case  
    if(M==2&m==0){
      if(!is.null(coef_temp)){
        coef_sep_Ber[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_sep_Ber[i, 4:6] <- coef_temp[2,]
        }
      }
    }
    
    if(M==2&m==1){
      if(!is.null(coef_temp)){
        if(length(coef_temp[,1])==2){
          coef_sep_Ber[i, 4:6] <- coef_temp[2,]
        }
      }
    }
    
    if(M==1&m==0){
      if(!is.null(coef_temp)){
        coef_sep_Ber[i, 1:3] <- coef_temp[1,]
      }
    }
  }
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
      break_point = first_day[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_temp = data.frame()
      coef_temp = segGompertz(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
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

# retry fitting for unfitted countries
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
        break_point = first_day[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        coef_temp = segGompertz(Country=i_country,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
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

coefficient_result = cbind(coef_sep_Logi[,1:6],coef_sep_Gom)
#write.csv(coefficient_result,"coef_result.csv")




