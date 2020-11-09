#### 1. Import Data
library(GGally)
library(tidyverse)
library(PerformanceAnalytics)
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201102 Segmented Growth Curve coefficients Comparison")
coef <- read.csv("coef_result.csv")

colnames(coef); coef$X; summary(coef)

ggpairs(coef[2:7])
chart.Correlation(log(coef[2:7]), histogram = FALSE, pch=19)
chart.Correlation(log(coef[c(-4, -157), 2:7]), histogram = FALSE, pch=19)

chart.Correlation(log(coef[8:13]), histogram = FALSE, pch=19)
chart.Correlation(log(coef[c(-18, -82), 8:13]), histogram = FALSE, pch=19)


chart.Correlation(log(coef[,c(4,9,10)]), histogram = FALSE, pch=19)
chart.Correlation(log(coef[,c(7,12,13)]), histogram = FALSE, pch=19)


log_coef_st <- scale(log(coef[, c(2:13)]))

chart.Correlation(log_coef_st[, c(3,8,9)], histogram = FALSE, pch=19)
chart.Correlation(log_coef_st[, c(6,11,12)], histogram = FALSE, pch=19)


t.test(log_coef_st[,3], log_coef_st[,8])

t_stat_1 <- t.test(log(coef[,4]), log(coef[,9]))$statistic

t_stat <- t.test(log_coef_st[,3], log_coef_st[,8])$statistic

set.seed(1)
tperm = perm.test(log_coef_st[,3], log_coef_st[,8], n.perm=1000) # 1000번 permutation 하여, 1000 개의 t-statistc 계산

hist(tperm)

abline(v=abs(t_stat),lty=2,col=2) # 실제 t-value 가 분포의 극단치에서 보임(reject H0으로 보임)

# Empirical p-value
pvalue=mean(abs(tperm) >= abs(t_stat)) # 위 그래프에서 Red (real t-value) 오른쪽에 있는 개수들의 평균을 구함
pvalue
# original p-value=1과 같은 값

# Conlcusion: 저체중아 여부에 따라 산모의 체중의 모평균은 통계적으로 유의미한 차이가 있다.





# 189 obs, 10 variables
str(birthwt)



# 2. 정상군과 실험군 분류 
normal=birthwt[birthwt[,"low"]==0,"lwt"]
normal
str(normal)

case=birthwt[birthwt[,"low"]==1,"lwt"]
case
str(case)

# original p-value=0.01308
t.test(normal,case)

# 3. 두 그룹의 산모 체중에 대한 t 검정 값
real_test=t.test(normal,case)$statistic
real_test

# 4. 두 그룹간의 permutation test
getwd()
setwd("/Users/yeonghyeon/Documents/고영현/수업/데이터마이닝방법 및 실습/실습/6. permutation Test")
source("permfunc.R") # 첨부파일 다운하여 working directory에 놓고 실행

set.seed(1)
tperm = perm.test(normal, case, n.perm=1000) # 1000번 permutation 하여, 1000 개의 t-statistc 계산

hist(tperm)
abline(v=abs(real_test),lty=2,col=2) # 실제 t-value 가 분포의 극단치에서 보임(reject H0으로 보임)

# 5. Empirical p-value
pvalue=mean(abs(tperm)>=abs(real_test)) # 위 그래프에서 Red (real t-value) 오른쪽에 있는 개수들의 평균을 구함
# 0.017 : significant ; original p-value=0.01308 와 유사한 값이 나옴 (seed 값에 따라 조금씩은 바뀜)
# Conlcusion: 저체중아 여부에 따라 산모의 체중의 모평균은 통계적으로 유의미한 차이가 있다.
pvalue


# age 검정
str(birthwt)
# 2. 정상군과 실험군 분류 
normal=birthwt[birthwt[,"low"]==0,"age"]
normal
str(normal)

case=birthwt[birthwt[,"low"]==1,"age"]
case
str(case)

# original p-value=0.07834
t.test(normal,case)

# 3. 두 그룹의 산모 age에 대한 t 검정 값=1.773672 (threshold for permutation test)
real_test=t.test(normal,case)$statistic
real_test

# 4. 두 그룹간의 permutation test
getwd()
setwd("C:/Users/keonvin/Desktop/데이터마이닝 방법 및 실습 (박태성 교수님;TA 과목)/KV Park/9주차 실습(permutation test)")
source("permfunc.R") # 첨부파일 다운하여 working directory에 놓고 실행

set.seed(1)
tperm=perm.test(normal,case,n.perm=1000) # 1000번 permutation 하여, 1000 개의 t-statistc 계산

hist(tperm)
abline(v=abs(real_test),lty=2,col=2) # 실제 t-value 가 분포의 극단치에서 보임(reject H0으로 보임)

# 5. Empirical p-value
pvalue=mean(abs(tperm)>=abs(real_test)) # 위 그래프에서 Red (real t-value) 오른쪽에 있는 개수들의 평균을 구함
# 0.086 : significant ; original p-value=0.07834 와 유사한 값이 나옴 (seed 값에 따라 조금씩은 바뀜)
# Conlcusion: 저체중아 여부에 따라 산모의 나이의 모평균은 통계적으로 유의미한 차이가 없다.
pvalue