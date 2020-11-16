#### 1. Import Data
library(GGally)
library(tidyverse)
library(PerformanceAnalytics)
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201109_Permutation_Test")
coef <- read.csv("coef_result.csv")

outlier_countries <- c("Andorra", "Aruba", "Benin", "French_Polynesia", 
                           "Kyrgyzstan", "Kuwait", "Mongolia", "Niger",
                           "Sao_Tome_and_Principe", "Seychelles")

exclude_countries <- c("Albania", "Argentina", "Cameroon", "Colombia",
                       "Dominican_Republic", "India", "Indonesia", "Kosovo",
                       "Moldova", "Mozambique", "Namibia", "Puero_Rico",
                       "uzbekistan", "Venezuela", "Zambia")

coef <- filter(read.csv("coef_result.csv"), !X %in% small_cases_countries)[, -1]
coef <- filter(read.csv("coef_result.csv"), !X %in% exclude_countries)[, -1]
log_coef <- log(coef)
log_coef_st <- scale(log(coef))


chart.Correlation(log_coef[, c(1:6)], histogram = TRUE, pch=19)

#### 2. Correlation between Segments / Models
## segment 1에서 계수별 비교
chart.Correlation(coef[, c(2,3,8,9)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(2,3,8,9)], histogram = TRUE, pch=19)
chart.Correlation(log_coef_st[, c(2,3,8,9)], histogram = TRUE, pch=19)
# b1_Logi와 b1_Gom, c1_Logi와 c1_Gom이 각각 상관계수가 높음 -> permutation test

## segment 2에서 계수별 비교
chart.Correlation(coef[, c(5,6,11,12)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(5,6,11,12)], histogram = TRUE, pch=19)
chart.Correlation(log_coef_st[, c(5,6,11,12)], histogram = TRUE, pch=19)
# segment 1에서와 마찬가지로 Logistic, Gompertz 모델의 계수 b와 c가 각각 상관계수가 높음.
# 다만 Logistic에서의 계수 b와 c가 상당히 상관계수가 높게 나와 의아함 -> 

## segment간 Logistic 계수별 비교
chart.Correlation(coef[, c(2,5,3,6)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(2,5,3,6)], histogram = TRUE, pch=19)
chart.Correlation(log_coef_st[, c(2,5,3,6)], histogram = TRUE, pch=19)
# Logistic의 계수 비교에서 segment 2일 때 계수 b와 c의 상관계수가 높음.

## segment간 Gompertz 계수별 비교
chart.Correlation(coef[, c(8,11,9,12)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(8,11,9,12)], histogram = TRUE, pch=19)
chart.Correlation(log_coef_st[, c(8,11,9,12)], histogram = TRUE, pch=19)
# Gompertz의 계수 비교에서 segment 2일 때 계수 b와 c의 상관계수가 높음.


#### 3. Permutation Test
## segment별 계수 비교
# b1_Logi vs b2_Logi
set.seed(1)
chart.Correlation(log_coef[, c(2,5)], histogram = TRUE, pch=19)
perm_data <- log_coef[, c(2, 5)] %>% filter(!is.na(b2_Logi), !is.nan(b1_Logi))
t_stat <- t.test(perm_data[,1], perm_data[,2])$statistic
t_perm <- paired.perm.test(perm_data[,1], perm_data[,2], n.perm=1000)
hist(t_perm, xlim=c(-20,20),
     main="Permutation Test for b1_Logi and b2_Logi")
abline(v=abs(t_stat),lty=2,col=2)
pvalue=mean(abs(t_perm)>=abs(t_stat)); pvalue

# c1_Logi vs c2_Logi
set.seed(1)
chart.Correlation(log_coef[, c(3,6)], histogram = TRUE, pch=19)
perm_data <- log_coef[, c(3, 6)] %>% filter(!is.na(c2_Logi), !is.nan(c1_Logi))
t_stat <- t.test(perm_data[,1], perm_data[,2])$statistic
t_perm <- perm.test(perm_data[,1], perm_data[,2], n.perm=1000)
hist(t_perm, main="Permutation Test for c1_Logi and c2_Logi")
abline(v=abs(t_stat),lty=2,col=2)
pvalue=mean(abs(tp)>=abs(real_test))
pvalue=mean(abs(t_perm)>=abs(t_stat)); pvalue

# b1_Gom vs b2_Gom
set.seed(1)
chart.Correlation(log_coef[, c(8,11)], histogram = TRUE, pch=19)
perm_data <- log_coef[, c(8, 11)] %>% filter(!is.na(b2_Gom), !is.nan(b1_Gom))
t_stat_ <- t.test(perm_data[,1], perm_data[,2])$statistic
t_perm <- perm.test(perm_data[,1], perm_data[,2], n.perm=1000)
hist(t_perm, main="Permutation Test for b1_Gom and b2_Gom")
abline(v=abs(t_stat),lty=2,col=2)
pvalue=mean(abs(t_perm)>=abs(t_stat)); pvalue

# c1_Gom vs c2_Gom
set.seed(1)
chart.Correlation(log_coef[, c(9,12)], histogram = TRUE, pch=19)
perm_data <- log_coef[, c(9, 12)] %>% filter(!is.na(c2_Gom), !is.nan(c1_Gom))
t_stat_ <- t.test(perm_data[,1], perm_data[,2])$statistic
t_perm <- perm.test(perm_data[,1], perm_data[,2], n.perm=1000)
hist(t_perm, main="Permutation Test for c1_Gom and c2_Gom")
abline(v=abs(t_stat),lty=2,col=2)
pvalue=mean(abs(t_perm)>=abs(t_stat)); pvalue

