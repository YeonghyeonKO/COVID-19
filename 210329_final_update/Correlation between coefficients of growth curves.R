#rm(list=ls())
library(dplyr)
library(PerformanceAnalytics)

#### Importing Data #### 
# Importing daily confirmed cases data
getwd()
setwd("/Users/yeonghyeonko/Documents/GitHub/COVID-19/210324_update_date/")

coef_Logi <- read.csv("result/segmented_Logistic_result.csv")[1:10]
colnames(coef_Logi) <- c("Country", 
                         "α1(Logistic)", "β1(Logistic)", "γ1(Logistic)",
                         "α2(Logistic)", "β2(Logistic)", "γ2(Logistic)",
                         "α3(Logistic)", "β3(Logistic)", "γ3(Logistic)")

coef_Gomp <- read.csv("result/segmented_Gompertz_result.csv")[1:10]
colnames(coef_Gomp) <- c("Country", 
                         "α1(Gompertz)", "β1(Gompertz)", "γ1(Gompertz)",
                         "α2(Gompertz)", "β2(Gompertz)", "γ2(Gompertz)",
                         "α3(Gompertz)", "β3(Gompertz)", "γ3(Gompertz)")

coef_total <- cbind(coef_Logi, coef_Gomp)[,-11]

#### Correlation Analysis ####
# correlation coefficients of Logistic Models
chart.Correlation(log(coef_Logi[, c(2:10)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(2:4)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(5:7)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(8:10)]), histogram = TRUE, pch=19)


chart.Correlation(log(coef_Logi[, c(2,5,8)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(3,6,9)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Logi[, c(4,7,10)]), histogram = TRUE, pch=19)


# correlation coefficients of Gompertz Models
chart.Correlation(log(coef_Gomp[, c(2:10)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(2:4)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(5:7)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(8:10)]), histogram = TRUE, pch=19)

chart.Correlation(log(coef_Gomp[, c(2,5,8)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(3,6,9)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_Gomp[, c(4,7,10)]), histogram = TRUE, pch=19)


# correlation coefficients of both models in the same segment.
chart.Correlation(log(coef_total[,c(2:4, 11:13)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_total[,c(5:7, 14:16)]), histogram = TRUE, pch=19)
chart.Correlation(log(coef_total[,c(8:10, 17:19)]), histogram = TRUE, pch=19)

