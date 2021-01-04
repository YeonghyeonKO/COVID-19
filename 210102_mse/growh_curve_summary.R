#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

#### Study Period Setting( ~ 11/30) ####
max_date <- set_date("2020/11/30") 

#### Importing Data #### 
# Importing daily confirmed cases data
getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/210102_mse/data")
df <- read.csv("COVID-19-geographic-disbtribution-worldwide-2020-12-14.csv")

#### Data Preprocessing & Segmentation ####
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# Preprocess data by country
df_sum = preprocessing_data()[1:(max_date),]

# Segmentation
df_result = segmentation(Country = country,gkf_bandwidth = 14,first_day_criteria = 50,save_image = TRUE,save_excel = FALSE)
# df_result = read.csv("C:/Users/김학용/Desktop/코로나 연구 인턴/주차별 작업사항/최종 정리(날짜 수정)/result/segmentation_result.csv")

# Simple Moving Average(window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}

#### Segmented Growth Curve model fitting ####  
# - use break_sum data.frame which is union of first_break and breakpoint
# - 1st break_point is not real breakpoint but starting point of model fitting            
# - from 2nd point, real breakpoint is counted                                           
# - (ex) break_point=c(50, 70, 120) => 50 is the date at which model fitting is started    
# -      so the # of real breakpoint is 2 (70, 120)                                             
par(mfrow=c(1,1))
break_sum = cbind(df_result[,c(3,11:15)])
rownames(break_sum) = df_result$country
country = as.character(df_result$country)

#### 1) Segmented Logistic model Fitting ####
coef_seg_Logi <- data.frame(matrix(NA, ncol = 9, nrow = length(country)))
rownames(coef_seg_Logi) <- country
colnames(coef_seg_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", "a2_Logi", "b2_Logi", "c2_Logi",
                             "a3_Logi", "b3_Logi", "c3_Logi")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = break_sum$first_day[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_temp = data.frame()
      coef_temp = segLogistic(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                       max_iter = 1000)
      if(!is.null(coef_temp)){
        coef_seg_Logi[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])>=2){
          coef_seg_Logi[i, 4:6] <- coef_temp[2,]
          if(length(coef_temp[,1])==3){
            coef_seg_Logi[i, 7:9] <- coef_temp[3,]
          }
        }
      }
    }
  )
}

# retry fitting for unfitted countries
for(i in 1:length(country)){
  M = ifelse(df_result$break_number[i] <= 2, df_result$break_number[i]+1, 3)
  
  if(!is.na(coef_seg_Logi$a3_Logi[i])){
    m = 3
  }else{
    if(!is.na(coef_seg_Logi$a2_Logi[i])){
      m = 2          # number of fitted parameter pairs
    }else{
      if(!is.na(coef_seg_Logi$a1_Logi[i])){
        m = 1
      }else{
        m = 0
      }
    }
  }
    
  # Retry fitting 
  i_country = as.character(country[i])
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
        coef_temp = segLogistic(Country=i_country,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                         max_iter = 1000)
        
        if(!is.null(coef_temp)){
          coef_seg_Logi[i, 1:3] <- coef_temp[1,]
          if(length(coef_temp[,1])>=2){
            coef_seg_Logi[i, 4:6] <- coef_temp[2,]
            if(length(coef_temp[,1])==3){
              coef_seg_Logi[i, 7:9] <- coef_temp[3,]
            }
          }
        }
        
      }
    )
  }
}

write.csv(coef_seg_Logi,"result/segmented_Logistic_result.csv")

#### 2) Segmented Bertalanffy model Fitting ####
coef_seg_Ber <- data.frame(matrix(NA, ncol = 9, nrow = length(country)))
rownames(coef_seg_Ber) <- country
colnames(coef_seg_Ber) <- c("a1_Ber", "b1_Ber", "c1_Ber", "a2_Ber", "b2_Ber", "c2_Ber",
                            "a3_Ber", "b3_Ber", "c3_Ber")

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
      coef_temp = segBertalanffy(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                       max_iter = 1000)
      if(!is.null(coef_temp)){
        coef_seg_Ber[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])>=2){
          coef_seg_Ber[i, 4:6] <- coef_temp[2,]
          if(length(coef_temp[,1])==3){
            coef_seg_Ber[i, 7:9] <- coef_temp[3,]
          }
        }
      }
      
    }
  )
}

# retry fitting for unfitted countries
for(i in 1:length(country)){
  M = ifelse(df_result$break_number[i] <= 2, df_result$break_number[i]+1, 3)
  
  if(!is.na(coef_seg_Ber$a3_Ber[i])){
    m = 3
  }else{
    if(!is.na(coef_seg_Ber$a2_Ber[i])){
      m = 2          # number of fitted parameter pairs
    }else{
      if(!is.na(coef_seg_Ber$a1_Ber[i])){
        m = 1
      }else{
        m = 0
      }
    }
  }
  
  # Retry fitting 
  i_country = as.character(country[i])
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
        coef_temp = segBertalanffy(Country=i_country,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                max_iter = 1000)
        
        if(!is.null(coef_temp)){
          coef_seg_Ber[i, 1:3] <- coef_temp[1,]
          if(length(coef_temp[,1])>=2){
            coef_seg_Ber[i, 4:6] <- coef_temp[2,]
            if(length(coef_temp[,1])==3){
              coef_seg_Ber[i, 7:9] <- coef_temp[3,]
            }
          }
        }
        
      }
    )
  }
}

write.csv(coef_seg_Ber,"result/segmented_bertalanffy_result.csv")

#### 3) Segmented Gompertz model Fitting ####
coef_seg_Gom <- data.frame(matrix(NA, ncol = 9, nrow = length(country)))
rownames(coef_seg_Gom) <- country
colnames(coef_seg_Gom) <- c("a1_Gom", "b1_Gom", "c1_Gom", "a2_Gom", "b2_Gom", "c2_Gom",
                            "a3_Gom", "b3_Gom", "c3_Gom")

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
      coef_temp = segGompertz(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                       max_iter = 1000)
      if(!is.null(coef_temp)){
        coef_seg_Gom[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])>=2){
          coef_seg_Gom[i, 4:6] <- coef_temp[2,]
          if(length(coef_temp[,1])==3){
            coef_seg_Gom[i, 7:9] <- coef_temp[3,]
          }
        }
      }
      
    }
  )
}


# retry fitting for unfitted countries
for(i in 1:length(country)){
  M = ifelse(df_result$break_number[i] <= 2, df_result$break_number[i]+1, 3)
  
  if(!is.na(coef_seg_Gom$a3_Gom[i])){
    m = 3
  }else{
    if(!is.na(coef_seg_Gom$a2_Gom[i])){
      m = 2          # number of fitted parameter pairs
    }else{
      if(!is.na(coef_seg_Gom$a1_Gom[i])){
        m = 1
      }else{
        m = 0
      }
    }
  }
  
  # Retry fitting 
  i_country = as.character(country[i])
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
        coef_temp = segGompertz(Country=i_country,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                   max_iter = 1000)
        
        if(!is.null(coef_temp)){
          coef_seg_Gom[i, 1:3] <- coef_temp[1,]
          if(length(coef_temp[,1])>=2){
            coef_seg_Gom[i, 4:6] <- coef_temp[2,]
            if(length(coef_temp[,1])==3){
              coef_seg_Gom[i, 7:9] <- coef_temp[3,]
            }
          }
        }
        
      }
    )
  }
}

write.csv(coef_seg_Gom,"result/segmented_Gompertz_result.csv")
