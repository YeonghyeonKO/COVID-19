#### 1. Import Data
library(GGally)
library(tidyverse)
library(PerformanceAnalytics)
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201116_Permutation_Test")
coef <- read.csv("coef_result.csv")

which(coef$b1_Logi - coef$b2_Logi>0)
which(coef$b1_Gom - coef$b2_Gom>0)

outlier_countries <- c("Andorra", "Aruba", "Benin", "French_Polynesia", 
                           "Kyrgyzstan", "Kuwait", "Mongolia", "Niger",
                           "Sao_Tome_and_Principe", "Seychelles")

exclude_countries <- c("Albania", "Argentina", "Cameroon", "Colombia",
                       "Dominican_Republic", "India", "Indonesia", "Kosovo",
                       "Moldova", "Mozambique", "Namibia", "Puero_Rico",
                       "uzbekistan", "Venezuela", "Zambia")

coef <- filter(read.csv("coef_result.csv"), !X %in% outlier_countries)[, -1]
coef <- filter(read.csv("coef_result.csv"), !X %in% exclude_countries)[, -1]
log_coef <- log(coef)
log_coef_st <- scale(log(coef))

chart.Correlation(coef[, c(1:6)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(1:6)], histogram = TRUE, pch=19)

chart.Correlation(log_coef[, c(1,4)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(2,5)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(3,6)], histogram = TRUE, pch=19)

chart.Correlation(log_coef[, c(7,10)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(8,11)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(9,12)], histogram = TRUE, pch=19)

chart.Correlation(log_coef[, c(1:3)], histogram = TRUE, pch=19)
coef$b2_Logi - coef$b1_Logi

summary(lm(a1_Logi ~ c1_Logi, data=log_coef))
summary(lm(a2_Logi ~ c2_Logi, data=log_coef))
summary(lm(a1_Gom ~ c1_Gom, data=log_coef))
summary(lm(a2_Gom ~ c2_Gom, data=log_coef))

chart.Correlation(log_coef[, c(7:9)], histogram = TRUE, pch=19)

chart.Correlation(log_coef[, c(4:6)], histogram = TRUE, pch=19)
chart.Correlation(log_coef[, c(10:12)], histogram = TRUE, pch=19)


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
### corr coef 비교
chart.Correlation(log_coef[, c(1,4,2,5,3,6)])


### segment별 계수 비교
## a1_Logi vs a2_Logi
set.seed(1)
chart.Correlation(log_coef[, c(1,4)], histogram = TRUE, pch=19)
a_logi <- coef[, c(1, 4)] %>% filter(!is.na(a2_Logi), !is.nan(a1_Logi))
a_logi <- log_coef[, c(1, 4)] %>% filter(!is.na(a2_Logi), !is.nan(a1_Logi))
attach(a_logi)
boxplot(a_logi, main="Boxplots for Logistic coef. a between two segments")
paired_t_stat <- t.test(a1_Logi, a2_Logi, paired=TRUE)$statistic
# Paired t-test
paired_t_perm <- paired.perm.test(a_logi, n.perm=1000)
hist(paired_t_perm, xlim=c(-20,20),
     main="Paired Permutation Test for a1_Logi & a2_Logi")
abline(v=abs(paired_t_stat),lty=2,col=2)
pvalue=mean(abs(paired_t_perm)>=abs(paired_t_stat)); pvalue


## b1_Logi vs b2_Logi
set.seed(1)
chart.Correlation(log_coef[, c(2,5)], histogram = TRUE, pch=19)
b_logi <- log_coef[, c(2, 5)] %>% filter(!is.na(b2_Logi), !is.nan(b1_Logi))
attach(b_logi)
boxplot(b_logi, main="Boxplots for Logistic coef. b between two segments")
paired_t_stat <- t.test(b1_Logi, b2_Logi, paired=TRUE)$statistic
# Paired t-test
paired_t_perm <- paired.perm.test(b_logi, n.perm=1000)
hist(paired_t_perm, xlim=c(-20,20),
     main="Paired Permutation Test for b1_Logi & b2_Logi")
abline(v=abs(paired_t_stat),lty=2,col=2)
pvalue=mean(abs(paired_t_perm)>=abs(paired_t_stat)); pvalue


## c1_Logi vs c2_Logi
set.seed(1)
chart.Correlation(log_coef[, c(3,6)], histogram = TRUE, pch=19)
c_logi <- log_coef[, c(3, 6)] %>% filter(!is.na(c2_Logi), !is.nan(c1_Logi))
attach(c_logi)
boxplot(c_logi, main="Boxplots for Logistic coef. c between two segments")
paired_t_stat <- t.test(c1_Logi, c2_Logi, paired=TRUE)$statistic
# Paired t-test
paired_t_perm <- paired.perm.test(c_logi, n.perm=1000)
hist(paired_t_perm, xlim=c(-20,20),
     main="Paired Permutation Test for c1_Logi & c2_Logi")
abline(v=abs(paired_t_stat),lty=2,col=2)
pvalue=mean(abs(paired_t_perm)>=abs(paired_t_stat)); pvalue


### a1_Gom vs a2_Gom
set.seed(1)
chart.Correlation(log_coef[, c(7,10)], histogram = TRUE, pch=19)
a_Gom <- log_coef[, c(7, 10)] %>% filter(!is.na(a2_Gom), !is.nan(a1_Gom))
attach(a_Gom)
boxplot(a_Gom, main="Boxplots for Gompertz coef. a between two segments")
paired_t_stat <- t.test(a1_Gom, a2_Gom, paired=TRUE)$statistic
# Paired t-test
paired_t_perm <- paired.perm.test(a_Gom, n.perm=1000)
hist(paired_t_perm, xlim=c(-20,20),
     main="Paired Permutation Test for a1_Gom & a2_Gom")
abline(v=abs(paired_t_stat),lty=2,col=2, lwd=2)
pvalue=mean(abs(paired_t_perm)>=abs(paired_t_stat)); pvalue


### b1_Gom vs b2_Gom
set.seed(1)
chart.Correlation(log_coef[, c(8,11)], histogram = TRUE, pch=19)
b_Gom <- log_coef[, c(8, 11)] %>% filter(!is.na(b2_Gom), !is.nan(b1_Gom))
attach(b_Gom)
boxplot(b_Gom, main="Boxplots for Gompertz coef. b between two segments")
paired_t_stat <- t.test(b1_Gom, b2_Gom, paired=TRUE)$statistic
# Paired t-test
paired_t_perm <- paired.perm.test(b_Gom, n.perm=1000)
hist(paired_t_perm, xlim=c(-20,20),
     main="Paired Permutation Test for b1_Gom & b2_Gom")
abline(v=abs(paired_t_stat),lty=2,col=2, lwd=2)
pvalue=mean(abs(paired_t_perm)>=abs(paired_t_stat)); pvalue


### c1_Gom vs c2_Gom
set.seed(1)
chart.Correlation(log_coef[, c(9,12)], histogram = TRUE, pch=19)
c_Gom <- log_coef[, c(9, 12)] %>% filter(!is.na(c2_Gom), !is.nan(c1_Gom))
attach(c_Gom)
boxplot(c_Gom, main="Boxplots for Gompertz coef. c between two segments")
paired_t_stat <- t.test(c1_Gom, c2_Gom, paired=TRUE)$statistic
# Paired t-test
paired_t_perm <- paired.perm.test(c_Gom, n.perm=1000)
hist(paired_t_perm, xlim=c(-20,20),
     main="Paired Permutation Test for c1_Gom & c2_Gom")
abline(v=abs(paired_t_stat),lty=2,col=2, lwd=2)
pvalue=mean(abs(paired_t_perm)>=abs(paired_t_stat)); pvalue
