#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

#### Study Period Setting( ~ 08/31) ####
max_date <- set_date("2020/08/31") 

#### Importing Data #### 
# Importing daily confirmed cases data
getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/210104_MSE_MSSE/Summary based Analysis 최종본")
df <- read.csv("data/COVID-19-geographic-disbtribution-worldwide-2020-12-14.csv")
# df <- read.csv("C:/Users/김학용/Desktop/코로나 연구 인턴/주차별 작업사항/최종 정리(8.31)/data/COVID-19-geographic-disbtribution-worldwide-2020-12-14.csv")

#### Data Preprocessing & Segmentation ####
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# Preprocess data by country
df_sum = preprocessing_data()[2:(max_date+1),]

# Segmentation
df_result = segmentation(Country = country,gkf_bandwidth = 14,first_day_criteria = 50,save_image = FALSE,save_excel = FALSE)
# df_result = read.csv("C:/Users/김학용/Desktop/코로나 연구 인턴/주차별 작업사항/최종 정리(날짜 수정)/result/segmentation_result.csv")

# Simple Moving Average(window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}

#### Segmented Growth Curve model fitting ####  
# - use break_sum data.frame which is union of first_break and breakpoint
par(mfrow=c(1,1))
break_sum = cbind(df_result[,c(3,11:15)])
rownames(break_sum) = df_result$country
country = as.character(df_result$country)

#### 1) Segmented Logistic model Fitting ####
coef_seg_Logi = data.frame(matrix(NA, ncol = 7, nrow = length(country)))
rownames(coef_seg_Logi) <- country
colnames(coef_seg_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", "a2_Logi", "b2_Logi", "c2_Logi", 
                             "MSE", "MAPE", "MSSE")

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
      result_temp = segLogistic(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                       max_iter = 1000)
      coef_temp = result_temp[[1]]
      
      if(!is.null(coef_temp)){
        coef_seg_Logi[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])>=2){
          coef_seg_Logi[i, 4:6] <- coef_temp[2,]
        }
        coef_seg_Logi[i, 7] <- result_temp[[2]]
        coef_seg_Logi[i, 8] <- result_temp[[3]]
        coef_seg_Logi[i, 9] <- result_temp[[4]]
      }
    }
  )
}

# retry fitting for unfitted countries
for(i in 1:length(country)){
  if(!is.na(break_sum$break1[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_seg_Logi$a2_Logi[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_seg_Logi$a1_Logi[i])){
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
        break_point = break_sum$first_day[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        result_temp = segLogistic(Country=i_country,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                max_iter = 1000)
        coef_temp = result_temp[[1]]
      }
    )
    
    # Case by Case  
    if(M==2&m==0){
      if(!is.null(coef_temp)){
        coef_seg_Logi[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_seg_Logi[i, 4:6] <- coef_temp[2,]
        }
        coef_seg_Logi[i, 7] <- result_temp[[2]]
        coef_seg_Logi[i, 8] <- result_temp[[3]]
        coef_seg_Logi[i, 9] <- result_temp[[4]]
      }
    }
    
    if(M==2&m==1){
      if(!is.null(coef_temp)){
        if(length(coef_temp[,1])==2){
          coef_seg_Logi[i, 4:6] <- coef_temp[2,]
        }
        coef_seg_Logi[i, 7] <- result_temp[[2]]
        coef_seg_Logi[i, 8] <- result_temp[[3]]
        coef_seg_Logi[i, 9] <- result_temp[[4]]
      }
    }
    
    if(M==1&m==0){
      if(!is.null(coef_temp)){
        coef_seg_Logi[i, 1:3] <- coef_temp[1,]
        coef_seg_Logi[i, 7] <- result_temp[[2]]
        coef_seg_Logi[i, 8] <- result_temp[[3]]
        coef_seg_Logi[i, 9] <- result_temp[[4]]
      }
    }
  }
}

write.csv(coef_seg_Logi,"result/segmented_Logistic_result.csv")
country_Logi = rownames(coef_seg_Logi)[which(coef_seg_Logi$MSSE<=0.4)]
# among 134 countries, 124 countries are fitted and 5 countries are excluded for high msse
# Final result is 119 countries


#### 2) Segmented Gompertz model Fitting ####
coef_seg_Gom <- data.frame(matrix(NA, ncol = 7, nrow = length(country)))
rownames(coef_seg_Gom) <- country
colnames(coef_seg_Gom) <- c("a1_Gom", "b1_Gom", "c1_Gom", "a2_Gom", "b2_Gom", "c2_Gom", 
                            "MSE", "MAPE", "MSSE")

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
      result_temp = segGompertz(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                       max_iter = 1000)
      coef_temp = result_temp[[1]]
      
      if(!is.null(coef_temp)){
        coef_seg_Gom[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])>=2){
          coef_seg_Gom[i, 4:6] <- coef_temp[2,]
        }
        coef_seg_Gom[i, 7] <- result_temp[[2]]
        coef_seg_Gom[i, 8] <- result_temp[[3]]
        coef_seg_Gom[i, 9] <- result_temp[[4]]
      }
    }
  )
}


# retry fitting for unfitted countries
for(i in 1:length(country)){
  if(!is.na(break_sum$break1[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_seg_Gom$a2_Gom[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_seg_Gom$a1_Gom[i])){
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
        break_point = break_sum$first_day[temp]
        for(j in 2:5){
          if(!is.na(break_sum[temp,j])){
            break_point = c(break_point,break_sum[temp,j])
          }
        }
        coef_temp = data.frame()
        result_temp = segGompertz(Country=i_country,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                max_iter = 1000)
        coef_temp = result_temp[[1]]
      }
    )
    
    # Case by Case  
    if(M==2&m==0){
      if(!is.null(coef_temp)){
        coef_seg_Gom[i, 1:3] <- coef_temp[1,]
        if(length(coef_temp[,1])==2){
          coef_seg_Gom[i, 4:6] <- coef_temp[2,]
        }
        coef_seg_Gom[i, 7] <- result_temp[[2]]
        coef_seg_Gom[i, 8] <- result_temp[[3]]
        coef_seg_Gom[i, 9] <- result_temp[[4]]
      }
    }
    
    if(M==2&m==1){
      if(!is.null(coef_temp)){
        if(length(coef_temp[,1])==2){
          coef_seg_Gom[i, 4:6] <- coef_temp[2,]
        }
        coef_seg_Gom[i, 7] <- result_temp[[2]]
        coef_seg_Gom[i, 8] <- result_temp[[3]]
        coef_seg_Gom[i, 9] <- result_temp[[4]]
      }
    }
    
    if(M==1&m==0){
      if(!is.null(coef_temp)){
        coef_seg_Gom[i, 1:3] <- coef_temp[1,]
        coef_seg_Gom[i, 7] <- result_temp[[2]]
        coef_seg_Gom[i, 8] <- result_temp[[3]]
        coef_seg_Gom[i, 9] <- result_temp[[4]]
      }
    }
  }
}

write.csv(coef_seg_Gom,"result/segmented_Gompertz_result.csv")
country_Gom = rownames(coef_seg_Gom)[which(coef_seg_Gom$MSSE<=0.4)]
# among 134 countries, 119 countries are fitted and 5 countries are excluded for high msse
# Final result is 114 countries

