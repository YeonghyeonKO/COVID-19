#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

# Analysis Period( ~ 8/31)
max_date <- set_date("2020/8/31") 

#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data (until 20/08/31)
df_sum = preprocessing_data()[1:245,]

# Simple Moving average (window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}

# Peak Detectiong Analysis
df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)

derivative_analysis(Country = "Germany",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)


###############################################################################

#### segmented logistic ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #

par(mfrow=c(1,1))
break_sum = df_result[,13:17]
rownames(break_sum) = df_result$country

main_country <- c("South_Korea", "Japan", "Italy", "Iran", "United_States_of_America",
                  "Canada", "Spain", "Russia")

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
      segLogistic(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
                  max_iter = 1000)
      
    },
    error = function(e){
      next
    },
    finally = next
  )
}

country = c("Austria", "United_States_of_America", "South_Korea")

#### segmented logistic - daily cases ####
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
      segLogistic_daily(Country=i,break_point = break_point,
                        save_image = TRUE,prediction_plot = TRUE,
                        max_iter = 1000)
    },
    error = function(e){
      next
    },
    finally = next
  )
}


#### separately fit logistic models #### 
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
      segLogistic_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = TRUE,
                  max_iter = 1000)
    },
    error = function(e){
    next
  },
  finally = next
  )
}



###############################################################################

#### segmented Bertalanffy ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #


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
      segBertalanffy(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
                  max_iter = 1000)
      
    },
    error = function(e){
      next
    },
    finally = next
  )
}


#### segmented Bertalanffy - daily cases ####
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
      segBertalanffy_daily(Country=i,break_point = break_point,
                        save_image = TRUE,prediction_plot = TRUE,
                        max_iter = 1000)
    },
    error = function(e){
      next
    },
    finally = next
  )
}


#### separately fit Bertalanffy models #### 
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
      segBertalanffy_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
                           max_iter = 1000)
    },
    error = function(e){
      next
    },
    finally = next
  )
}


###############################################################################

#### segmented Gompertz ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #

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
      segGompertz(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = TRUE,
                     max_iter = 1000)
      
    },
    error = function(e){
      next
    },
    finally = next
  )
}

#### segmented Gompertz - daily cases ####
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
      segGompertz_daily(Country=i,break_point = break_point,
                           save_image = TRUE,prediction_plot = TRUE,
                           max_iter = 1000)
    },
    error = function(e){
      next
    },
    finally = next
  )
}


#### separately fit Gompertz models #### 
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
      segGompertz_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = TRUE,
                           max_iter = 1000)
    },
    error = function(e){
      next
    },
    finally = next
  )
}



#### country classification by # of waves ####
# number of waves is counted by # of breakpoint.
sum(table(df_result$break_number))
table(df_result$break_number)

pp1 <- ggplot()+
  geom_point(aes(1:length(df_sum$Malawi), df_sum$Malawi), alpha=0.5) +
  labs(title = "Malawi") +
  theme_bw()

pp2 <- ggplot()+
  geom_point(aes(1:length(df_sum$Benin), df_sum$Benin), alpha=0.5) +
  labs(title = "Benin") +
  theme_bw()

gridExtra::grid.arrange(pp1, pp2)


setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201116_Permutation_Test")
coef <- read.csv("coef_result.csv")

which(coef$b1_Logi - coef$b2_Logi>0)
which(coef$b1_Gom - coef$b2_Gom>0)



# Andorra (b2_Logi=401) / Kyrgyzstan (b1_Gom=5915) / Benin (b1_Gom=3454) / Kuwait (b2_Gom)
par(mfrow=c(1,2))
boxplot(coef[-c(4,18,81,82),]$b1_Logi, main="b1_Logi"); boxplot(coef[-c(4,18,81,82),]$b2_Logi,  main="b2_Logi")
boxplot(coef[-c(4,18,81,82),]$b1_Gom, main="b1_Gom"); boxplot(coef[-c(4,18,81,82),]$b2_Gom,  main="b2_Gom")

summary(coef[-c(4,18,81,82),]$b1_Gom)
summary(coef[-c(4,18,81,82),]$b2_Gom)
var(coef[-c(4,18,81,82),]$b2_Gom)

sd(coef[-c(4,18,81,82),]$b1_Logi[!is.na(coef[-c(4,18,81,82),]$b1_Logi)])
sd(coef[-c(4,18,81,82),]$b2_Logi[!is.na(coef[-c(4,18,81,82),]$b2_Gom)])
sd(coef[-c(4,18,81,82),]$b1_Gom[!is.na(coef[-c(4,18,81,82),]$b1_Gom)])
sd(coef[-c(4,18,81,82),]$b2_Gom[!is.na(coef[-c(4,18,81,82),]$b2_Gom)])
