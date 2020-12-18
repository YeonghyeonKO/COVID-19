#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)


# Analysis Period( ~ 8/31)
max_date <- set_date("2020/8/31") 


#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
# df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
# write.csv(df, "covid_data.csv")
df <- read.csv("covid_data.csv")
coef <- read.csv("coef_result.csv")

setwd("~/Documents/GitHub/COVID-19/201218_coefficient_b")
getwd()

# generate vector of country names
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data (until 20/08/31)
df_sum = preprocessing_data()[1:245,]

Malawi <- df_sum['Malawi']
Malawi_st <- coef[which(coef$X=="Malawi"),]$startpoint
Malawi_bp <- coef[which(coef$X=="Malawi"),]$breakpoint

Andorra <- df_sum['Andorra']
Andorra_st <- coef[which(coef$X=="Andorra"),]$startpoint
Andorra_bp <- coef[which(coef$X=="Andorra"),]$breakpoint

Benin <- df_sum['Benin']
Benin_st <- coef[which(coef$X=="Benin"),]$startpoint
Benin_bp <- coef[which(coef$X=="Benin"),]$breakpoint

ggplot() +
  geom_point(aes(x=1:nrow(Malawi), y=malawi$Malawi)) +
  geom_vline(aes(xintercept=Malawi_st), color=2, linetype=2) +
  geom_vline(aes(xintercept=Malawi_bp), color=2) +
  theme_bw() + 
  labs(title="Malawi", x="Days", y="Daily Cases")

ggplot() +
  geom_point(aes(x=1:nrow(Andorra), y=Andorra$Andorra)) +
  geom_vline(aes(xintercept=Andorra_st), color=2, linetype=2) +
  geom_vline(aes(xintercept=Andorra_bp), color=2) +
  theme_bw() + 
  labs(title="Andorra", x="Days", y="Daily Cases")

ggplot() +
  geom_point(aes(x=1:nrow(Benin), y=benin$Benin)) +
  geom_vline(aes(xintercept=Benin_st), color=2, linetype=2) +
  geom_vline(aes(xintercept=Benin_bp), color=2) +
  theme_bw() + 
  labs(title="Benin", x="Days", y="Daily Cases")


coef[which(coef$b1_Logi - coef$b2_Logi >0 ),]
coef[which(coef$b1_Gom - coef$b2_Gom >0 ),]


# Benin, Malawi, Andorra 제거
coef <- coef[-which(coef$X%in%c("Benin", "Andorra", "Malawi")),]

boxplot(coef$b1_Logi, coef$b2_Logi, main="b_Logistic", xlab="Segment")
boxplot(coef$b1_Gom, coef$b2_Gom, main="b_Gompertz", xlab="Segment")
# outlier : Kuwait(b2_Gom), Kyrgyzstan(b1_Gom)

hist(coef$b1_Logi - coef$b2_Logi, main="b1 - b2 for Logistic Model")
hist(coef$b1_Gom - coef$b2_Gom, main="b1 - b2 for Gompertz Model", breaks = 100)
coef <- coef[-which(coef$X%in%c("Kuwait", "Kyrgyzstan")),]
hist(coef$b1_Gom - coef$b2_Gom, main="b1 - b2 for Gompertz Model", breaks = 100)



#### coefficient b와 c와의 관계
ggplot(coef) +
  geom_point(aes(coef$b1_Logi, coef$c1_Logi)) +
  theme_bw()

ggplot(coef) +
  geom_point(aes(coef$b2_Logi- coef$b1_Logi, coef$c2_Logi)) +
  theme_bw()



ggplot(coef) +
  geom_point(aes(coef$b1_Gom, coef$c1_Gom)) +
  theme_bw()

ggplot(coef) +
  geom_point(aes(coef$b2_Gom, coef$c2_Gom)) +
  theme_bw()

ggplot(coef) +
  geom_point(aes(coef$b2_Gom - coef$b1_Logi, coef$c2_Gom)) +
  theme_bw()
