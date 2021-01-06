#rm(list=ls())
library(dplyr)
library(PerformanceAnalytics)

#### Importing Data #### 
# Importing daily confirmed cases data
getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/210104_MSE_MSSE/Summary based Analysis 최종본/")

coef_Logi <- read.csv("result/segmented_Logistic_result.csv")[1:7]
coef_Gomp <- read.csv("result/segmented_Gompertz_result.csv")[1:7]



#### Correlation Analysis ####
# correlation coefficients of Logistic Models
chart.Correlation(log(coef_Logi[, c(2:7)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(2:4)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(5:7)]), histogram = TRUE, pch=19)

chart.Correlation(log(coef_Logi[, c(2,5)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(3,6)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(4,7)]), histogram = TRUE, pch=19)


# correlation coefficients of Gompertz Models
chart.Correlation(log(coef_Gomp[, c(2:7)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(2:4)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(5:7)]), histogram = TRUE, pch=19)

chart.Correlation(log(coef_Gomp[, c(2,5)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(3,6)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(4,7)]), histogram = TRUE, pch=19)


# correlation coefficients of both models in the same segment.
chart.Correlation(data.frame(log(coef_Logi[,c(2:4)]),log(coef_Gomp[, c(2:4)])), histogram = TRUE, pch=19)
chart.Correlation(data.frame(log(coef_Logi[,c(5:7)]),log(coef_Gomp[, c(5:7)])), histogram = TRUE, pch=19)

