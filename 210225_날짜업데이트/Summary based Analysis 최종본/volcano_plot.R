library(dplyr)
library(ggplot2)


getwd()
setwd("/Users/yeonghyeonko/Documents/GitHub/COVID-19/210225_날짜업데이트/Summary based Analysis 최종본/")

slr1 <- read.csv("result/Linear Regression(1st Segment, standardized).csv")
slr2 <- read.csv("result/Linear Regression(2nd Segment, standardized).csv")

Logi_a <- cbind(slr1[,c(2,3,4)], a2_Logi=slr2[,4])
Logi_b <- cbind(slr1[,c(2,3,5)], b2_Logi=slr2[,5])
Logi_c <- cbind(slr1[,c(2,3,6)], c2_Logi=slr2[,6])

Gomp_a <- cbind(slr1[,c(2,3,7)], a2_Gomp=slr2[,7])
Gomp_b <- cbind(slr1[,c(2,3,8)], b2_Gomp=slr2[,8])
Gomp_c <- cbind(slr1[,c(2,3,9)], c2_Gomp=slr2[,9])


df_a <- data.frame(
  Explanatory = factor(unique(Logi_a$Explanatory)),
  coefficient = c(Logi_a$a1_Logi[which((1:nrow(Logi_a))%%3==1)], Logi_a$a2_Logi[which((1:nrow(Logi_a))%%3==1)]),
  p.value = c(Logi_a$a1_Logi[which((1:nrow(Logi_a))%%3==0)], Logi_a$a2_Logi[which((1:nrow(Logi_a))%%3==0)]),
  segment = rep(1:2, each=length(unique(Logi_a$Explanatory)))
)

# non-log scaled
ggplot(df_a, aes(x=coefficient, y=p.value, 
                 group=Explanatory, color=Explanatory)) +
  geom_point(shape=df_a$segment) +
  geom_line() +
  labs(x="coefficient", y="p-value")

# log scaled (only p.value)
ggplot(df_a, aes(x=coefficient, y=log(p.value), 
                 group=Explanatory, color=Explanatory)) +
  geom_point(shape=df_a$segment) +
  geom_line() +
  labs(x="coefficient", y="log(p-value)")

# log scaled (only coefficient)
ggplot(df_a, aes(x=log(coefficient), y=p.value, 
                 group=Explanatory, color=Explanatory)) +
  geom_point(shape=df_a$segment) +
  geom_line() +
  labs(x="log(coefficient)", y="p-value")

# log scaled (both)
ggplot(df_a, aes(x=log(coefficient), y=log(p.value), 
                 group=Explanatory, color=Explanatory)) +
  geom_point(shape=df_a$segment) +
  geom_line() +
  labs(x="log(coefficient)", y="log(p-value)")
