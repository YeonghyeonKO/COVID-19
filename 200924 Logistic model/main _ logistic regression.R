#rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(nls2)
library(zoo)

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



country <- c("South_Korea", "Japan", "Italy", "Iran", "United_States_of_America",
             "Canada", "Spain", "Russia")

a2_sep <- c(63366, 62543, 94862, 302481, 4712386, 15732, 641896, 155988)
b2_sep <- c(7.31, 14.4, 13.41, 5.4, 9.11, 19.53, 12.07, 29.54)
c2_sep <- c(0.027, 0.081, 0.065, 0.04, 0.06, 0.138, 0.062, 0.164)
a2_seg <- c(281861, 62670, 334783, 292288, 4756758, 19728, 829985, 183613)
b2_seg <- c(7.91, 14.34, 9.1, 5.71, 8.96, 13.27, 10.51, 20.64)
c2_seg <- c(0.023, 0.08, 0.036, 0.043, 0.059, 0.095, 0.051, 0.116)


coef_2 <- as.data.frame(
  cbind(country, a2_seg, b2_seg, c2_seg,
        a2_sep, b2_sep, c2_sep)) %>%
    
  # alphabetically ordered
  arrange(country) %>%
  
  bind_cols(df_result[which(df_result$country %in% country),][13]) %>%

  # start
  bind_cols(df_result[which(df_result$country %in% country),][14]) %>%

  # end
  bind_cols(as.vector(na.fill(
    df_result[which(df_result$country %in% country),][15], 267)))

colnames(coef_2)[8:10] <- c("initial", "start", "end")
coef_2[, c(2:7)] <- sapply(coef_2[, c(2:7)], as.numeric)



## 아래 두 문장의 문제를 해결하여 ggplot() 2,3번째 줄 "df_sum$Canada"만 고치면 됨.
# noquote 해봤는데 안됨...
get(paste0("df_sum$", "Canada"))
df_sum$Canada

for (i in 1:8){
  # country <- noquote(coef_2[i, 1])
  x <- coef_2[i,]$start:coef_2[i,]$end
  y <- coef_2[i,]$a2_seg / (1 + exp(coef_2[i,]$b2_seg - coef_2[i,]$c2_seg * x))
  
  ggplot() +
    geom_point(aes(1:length(df_sum$Canada), cumsum(df_sum$Canada))) +
    geom_point(aes(x, cumsum(df_sum$Canada)[x]), col = 2) +
    geom_vline(aes(xintercept = c(coef_2[i,]$start, coef_2[i,]$end)), col = 2, linetype = 4) +
    theme_bw() + 
    labs(title = paste0("2nd Segmented Logistic coefficient in ", coef_2[i, 1]),
         x = "Days", y = "Cumulative Confirmed Cases")
}


eq <- a / (1 + exp(b - c*x))
