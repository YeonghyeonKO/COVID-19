#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

#### Study Period Setting( ~ 21/03/20) ####
max_date <- set_date("2021/03/20") 

#### Importing Data #### 
# Importing daily confirmed cases data
getwd()
setwd("/Users/yeonghyeonko/Documents/GitHub/COVID-19/210329_final_update/")
df <- read.csv("data/OxCGRT_210320.csv")

#### Data Preprocessing & Segmentation ####
country <- as.character(unique(df$CountryCode))
# country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# Preprocess data by country
df_sum = preprocessing_data()[0:max_date+1,]



dfdf <- data.frame(1:(ncol(df_sum)-1), apply(df_sum[,-1], 2, sum))
colnames(dfdf) <- c("country", "sum")

dfdf[,1] <- rownames(dfdf)

left_join(df_result, 
          dfdf,
          by="country")


write.csv(left_join(df_result, 
                    dfdf,
                    by="country"),
          "asdf.csv")

# Segmentation
df_result = segmentation(Country = country,
                         gkf_bandwidth = 14,
                         first_day_criteria = 50,
                         save_image = TRUE,
                         save_excel = TRUE)


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
coef_seg_Logi = data.frame(matrix(NA, ncol = 12, nrow = length(country)))
rownames(coef_seg_Logi) <- country
colnames(coef_seg_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", 
                             "a2_Logi", "b2_Logi", "c2_Logi", 
                             "a3_Logi", "b3_Logi", "c3_Logi", 
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
        if(length(coef_temp[,1])>=3){
          coef_seg_Logi[i, 7:9] <- coef_temp[3,]
        }
        coef_seg_Logi[i, 10] <- result_temp[[2]]
        coef_seg_Logi[i, 11] <- result_temp[[3]]
        coef_seg_Logi[i, 12] <- result_temp[[4]]
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
sum(!is.na(coef_seg_Logi$a1_Logi)) # the number of fitted countries (134)
# ( ~ 20/08/31) among 134 countries, 124 countries are fitted and 119 countries are remained by msse ( <= 0.4 )
# ( ~ 20/12/14) among 154 countries, 118 countries are fitted and 101 countries are remained by msse ( <= 0.4 )
# ( ~ 21/03/20) among 156 countries, 134 countries are fitted and 98 countries are remained by msse ( <= 0.4 )

# Final result is 119 (~ 20/08/31), 101 (~ 20/12/14), 98 (~ 21/03/20) countries respectively.



#### 2) Segmented Gompertz model Fitting ####
coef_seg_Gom <- data.frame(matrix(NA, ncol = 12, nrow = length(country)))
rownames(coef_seg_Gom) <- country
colnames(coef_seg_Gom) <- c("a1_Gom", "b1_Gom", "c1_Gom", 
                            "a2_Gom", "b2_Gom", "c2_Gom", 
                            "a3_Gom", "b3_Gom", "c3_Gom", 
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
        if(length(coef_temp[,1])>=3){
          coef_seg_Gom[i, 7:9] <- coef_temp[3,]
        }
        coef_seg_Gom[i, 10] <- result_temp[[2]]
        coef_seg_Gom[i, 11] <- result_temp[[3]]
        coef_seg_Gom[i, 12] <- result_temp[[4]]
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
sum(!is.na(coef_seg_Gom$a1_Gom)) # the number of fitted countries (123)
# ( ~ 20/08/31) among 134 countries, 119 countries are fitted and 114 countries are remained by msse ( <= 0.4 )
# ( ~ 20/12/14) among 154 countries, 111 countries are fitted and 101 countries are remained by msse ( <= 0.4 )
# ( ~ 21/03/20) among 156 countries, 123 countries are fitted and 84 countries are remained by msse ( <= 0.4 )

# Final result is 114 (~ 20/08/31), 101 (~ 20/12/14), 84 (~ 21/03/20) countries respectively.



KOR_data <- read.csv("result/segmented_Logistic_result.csv") %>%
  filter(X == "KOR")
?dplyr

b = df_result %>% filter(country=="KOR") %>% select(6:8)
y = cumsum(df_sum["KOR"][-c(1:b$first_day),])
x = 1:length(y)

plot(x,y,cex=0.5,
     xlab="Days since 2020-02-21",
     ylab="Cumulative Cases", 
     main=" KOR / Logistic",
     sub="Simple Moving Avarage (window size = 7)") 
abline(v=b,lty=2)
legend("topleft",c("data","fit1","fit2","fit3"),
       pch=c(1,-1,-1,-1),
       lty=c(0,1,1,1),
       col=c(1,2,3,4),
       lwd=c(1,2,2,2))

initial_list = data.frame(a1=KOR_data$a1_Logi, b1=KOR_data$b1_Logi, c1=KOR_data$c1_Logi)
fit1 <- nls2(y[1:b[[1]]] ~ a1/(1+exp(b1-c1*x[1:b[[1]]])), 
     start = initial_list,
     algorithm = "plinear-random",
     control = nls.control(maxiter=500))
pred_y = predict(fit1,newdata=data.frame(x))
lines(1:b$peak1, pred_y, col=2,lwd=3)

initial_list = data.frame(a2=KOR_data$a2_Logi, b2=KOR_data$b2_Logi, c2=KOR_data$c2_Logi)
fit2 <- nls2(y[b[[1]]:b[[2]]] ~ y[b[[1]]] + a2/(1+exp(b2-c2*x[ b[[1]]:b[[2]] ])), 
             start = initial_list,
             algorithm = "plinear-random",
             control = nls.control(maxiter=500))
pred_y = predict(fit2)

lines(b$peak1:b$peak2, pred_y + (y[b[[1]]]-pred_y[1]), col=3,lwd=3)

lines(pred_y, col=3,lwd=3)

initial_list = data.frame(a3=KOR_data$a3_Logi, b3=KOR_data$b3_Logi, c3=KOR_data$c3_Logi)
fit3 <- nls2(y[b[[2]]:b[[3]]] ~ a3/(1+exp(b3-c3*x[b[[2]]:b[[3]]])), 
             start = initial_list,
             algorithm = "plinear-random",
             control = nls.control(maxiter=500))
pred_y = predict(fit3, new)
lines(x3+(b[3]-1), pred_y, col=4,lwd=3)


lines(1:177,  a2/(1+exp(b2-c2*x[ b[[1]]:b[[2]] ])))
