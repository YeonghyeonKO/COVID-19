library(gridExtra)
library(locfit)
library(Metrics)
library(mgcv)
library(segmented)
library(splines)
library(tidyverse)


getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201219_Non_Linear_Model_Korea")
covid <- read_csv("KOR_additive_measure_cluster_201218.csv",
                  locale=locale(encoding="euc-kr"))[-c(1:10),c(3,19:26)]
colnames(covid)


#### 1. Plotting ####
## 전국 일일 확진자
p1_a <- ggplot(covid) +
  geom_point(aes(Days, Domestic_DailyConfirmed), alpha=0.5) +
  theme_classic(base_family="NanumGothic") +
  labs(x="From 2020/01/20", y="일일 확진자", title="전국 신규 확진자"); p1_a

## 수도권 일일 확진자
p1_b <- ggplot(covid) +
  geom_point(aes(Date, Capital_DailyConfirmed), alpha=0.5) +
  theme_classic(base_family="NanumGothic") +
  labs(x="From 2020/01/20", y="일일 확진자", title="수도권 신규 확진자"); p1_b

## 비수도권 일일 확진자
p1_c <- ggplot(covid) +
  geom_point(aes(Date, NonCapital_DailyConfirmed), alpha=0.5) +
  theme_classic(base_family="NanumGothic") +
  labs(x="From 2020/01/20", y="일일 확진자", title="비수도권 신규 확진자"); p1_c

grid.arrange(p1_a, p1_b, p1_c)



#### 2. 전국 일일 확진자 모델링 ####
train_idx <- which(covid$Date <= "2020-12-11")
train <- covid[train_idx,]
train_gov <- train$GovernmentIndex
train_gov_lag <- lag(train_gov, 12, default=0)
# lag = 12부터 Government Index의 효과가 발동.

test_idx <- which(covid$Date > "2020-12-11")
test <- covid[test_idx,]
test_gov <- test$GovernmentIndex

total_idx <- which(covid$Date <= "2020-12-18")
total <- covid[total_idx,]
total_gov <- total$GovernmentIndex
total_gov_lag <- lag(total_gov, 12, default=0)

pred_idx <- 335:354

y_train <- train$Domestic_DailyConfirmed
y_test <- test$Domestic_DailyConfirmed
y_total <- total$Domestic_DailyConfirmed


#### 2.1. Polynomial Regression ####
fit_1_dom <- lm(y_train ~ poly(train_idx, degree=6) + poly(train_gov_lag))
fit_1_dom_pred <- lm(y_total ~ poly(total_idx, degree=6) + poly(total_gov_lag))
summary(fit_1_dom); summary(fit_1_dom_pred)
predict_1_dom <- predict(fit_1_dom, data.frame(train_idx=test_idx, train_gov_lag=test_gov))
predict_1_dom_pred_0 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0))
predict_1_dom_pred_0.5 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0.5))
predict_1_dom_pred_1 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1))
predict_1_dom_pred_1.5 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1.5))
predict_1_dom_pred_2 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2))
predict_1_dom_pred_2.5 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2.5))
predict_1_dom_pred_3 <- predict(fit_1_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=3))

# train / test
p2_1_dom <- ggplot() +
  geom_point(aes(train_idx, y_train), alpha=0.5) +
  geom_line(aes(train_idx, fitted(fit_1_dom)), color=2, alpha=1) +
  
  geom_line(aes(train_idx, train_gov*500), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./500, name="사회적 거리두기 (단계)"))+
  
  geom_point(aes(test_idx, y_test), alpha=0.5) +
  geom_line(aes(test_idx, predict_1_dom), color=2, size=1, linetype=3) +
  
  geom_vline(aes(xintercept=tail(train_idx,1)), color=8, linetype=6, size=0.5) + 
  
  annotate("text", x=330, y=100, label="12/11", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Polynomial Regression (Degree = 6)\nLag for Government Index = 12"); p2_1_dom
training_rmse_1_dom <- rmse(y_train, fitted(fit_1_dom)); training_rmse_1_dom
testing_rmse_1_dom <- rmse(y_test, predict_1_dom); testing_rmse_1_dom

# predict
p2_1_dom_pred <- ggplot() +
  geom_point(aes(total_idx, y_total), alpha=0.5) +
  geom_line(aes(total_idx, fitted(fit_1_dom_pred)), color=2) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)"))+
  
  geom_line(aes(pred_idx, predict_1_dom_pred_0), color=3) +
  geom_line(aes(pred_idx, predict_1_dom_pred_0.5), color=4) +
  geom_line(aes(pred_idx, predict_1_dom_pred_1), color=5) +
  geom_line(aes(pred_idx, predict_1_dom_pred_1.5), color=6) +
  geom_line(aes(pred_idx, predict_1_dom_pred_2), color=7) +
  geom_line(aes(pred_idx, predict_1_dom_pred_2.5), color=8) +
  geom_line(aes(pred_idx, predict_1_dom_pred_3), color=9) +
  
  geom_vline(aes(xintercept=tail(total_idx,1)), color=8, linetype=6, size=0.5) + 
  
  annotate("text", x=330, y=100, label="12/18", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Polynomial Regression (Degree = 6)\nLag for Government Index = 12"); p2_1_dom_pred
# 향후 20일간 (12/18일 이후) 일일 확진자 예측값
cbind(predict_1_dom_pred_0, predict_1_dom_pred_0.5, predict_1_dom_pred_1, 
      predict_1_dom_pred_1.5, predict_1_dom_pred_2, predict_1_dom_pred_2.5, predict_1_dom_pred_3)





#### 2.2. Cubic Spline ####
train_gov_lag <- lag(train_gov, 12, default=0)
total_gov_lag <- lag(total_gov, 12, default=0)

fit_2_dom <- lm(y_train ~ bs(train_idx, knots=c(15, 190, 267)) + train_gov_lag)
fit_2_dom_pred <- lm(y_total ~ bs(total_idx,  knots=c(15, 190, 267)) + total_gov_lag)
summary(fit_2_dom); summary(fit_2_dom_pred)
predict_2_dom <- predict(fit_2_dom, data.frame(train_idx=test_idx, train_gov_lag=test_gov))
predict_2_dom_pred_0 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0))
predict_2_dom_pred_0.5 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0.5))
predict_2_dom_pred_1 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1))
predict_2_dom_pred_1.5 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1.5))
predict_2_dom_pred_2 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2))
predict_2_dom_pred_2.5 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2.5))
predict_2_dom_pred_3 <- predict(fit_2_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=3))

# train / test
p2_2_dom <- ggplot() +
  geom_point(aes(train_idx, y_train), alpha=0.5) +
  geom_line(aes(train_idx, fitted(fit_2_dom)), color=2, alpha=0.5) +
  
  geom_point(aes(test_idx, y_test), alpha=0.5) +
  geom_line(aes(test_idx, predict_2_dom), color=2, size=1, linetype=3) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)"))+
  
  geom_vline(aes(xintercept=tail(train_idx,1)), color=8, linetype=6, size=0.5) + 
  geom_vline(aes(xintercept=c(15, 190, 267)), color=8, linetype=3) + 
  
  annotate("text", x=330, y=100, label="12/11", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Cubic Spline (Knots = 15, 190, 267)\nLag for Government Index = 12"); p2_2_dom
training_rmse_2_dom <- rmse(y_train, fitted(fit_2_dom)); training_rmse_2_dom
testing_rmse_2_dom <- rmse(y_test, predict_2_dom); testing_rmse_2_dom

# predict
p2_2_dom_pred <- ggplot() +
  geom_point(aes(total_idx, y_total), alpha=0.5) +
  geom_line(aes(total_idx, fitted(fit_2_dom_pred)), color=2, alpha=0.5) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)")) +
  
  geom_line(aes(pred_idx, predict_1_dom_pred_0), color=3) +
  geom_line(aes(pred_idx, predict_1_dom_pred_0.5), color=4) +
  geom_line(aes(pred_idx, predict_1_dom_pred_1), color=5) +
  geom_line(aes(pred_idx, predict_1_dom_pred_1.5), color=6) +
  geom_line(aes(pred_idx, predict_1_dom_pred_2), color=7) +
  geom_line(aes(pred_idx, predict_1_dom_pred_2.5), color=8) +
  geom_line(aes(pred_idx, predict_1_dom_pred_3), color=9) +
  
  geom_vline(aes(xintercept=tail(total_idx,1)), color=8, linetype=6, size=0.5) + 
  geom_vline(aes(xintercept=c(15, 190, 267)), color=8, linetype=3) + 
  
  annotate("text", x=330, y=100, label="12/18", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Cubic Spline (Knots = 15, 190, 267)\nLag for Government Index = 12"); p2_2_dom_pred
# 향후 20일간 (12/18일 이후) 일일 확진자 예측값
cbind(predict_2_dom_pred_0, predict_2_dom_pred_0.5, predict_2_dom_pred_1, 
      predict_2_dom_pred_1.5, predict_2_dom_pred_2, predict_2_dom_pred_2.5, predict_2_dom_pred_3)




# #### 2.3. Smoothing Spline ####
# train_gov_lag <- lag(train_gov, 12, default=0)
# total_gov_lag <- lag(total_gov, 12, default=0)
# fit_3_dom <- gam(y_train ~ s())
# 
# 
# fit_3_dom <- smooth.spline(x=train_idx, y=y_train, cv = TRUE, df=10)
# fit_3_dom_pred <- smooth.spline(x=total_idx, y=y_total, cv = TRUE, df=10)
# fit_3_dom; fit_3_dom_pred
# predict_3_dom <- predict(fit_3_dom, data.frame(train_idx=test_idx))$y[[1]]
# predict_3_dom_pred <- predict(fit_3_dom_pred, data.frame(total_idx=pred_idx))$y[[1]]
# 
# # train / test
# p2_3_dom <- ggplot() +
#   geom_point(aes(train_idx, y_train), alpha=0.5) +
#   geom_line(aes(train_idx, fitted(fit_3_dom)), color=2, alpha=0.5) +
#   
#   geom_point(aes(test_idx, y_test), alpha=0.5) +
#   geom_line(aes(test_idx, predict_3_dom), color=2, size=1, linetype=3) +
#   
#   geom_vline(aes(xintercept=tail(train_idx,1)), color=8, linetype=6, size=0.5) + 
#   
#   annotate("text", x=300, y=1000, label="12/11", size=5, color=8) +
#   theme_classic(base_family="NanumGothic") +
#   labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
#        subtitle="Smoothing Spline (Df = 10)"); p2_3_dom
# training_rmse_3_dom <- rmse(y_train, fitted(fit_3_dom)); training_rmse_3_dom
# testing_rmse_3_dom <- rmse(y_test, predict_3_dom); testing_rmse_3_dom
# 
# # predict
# p2_3_dom_pred <- ggplot() +
#   geom_point(aes(total_idx, y_total), alpha=0.5) +
#   geom_line(aes(total_idx, fitted(fit_3_dom_pred)), color=2, alpha=0.5) +
#   
#   geom_line(aes(pred_idx, predict_3_dom_pred), color=2, size=1, linetype=3) +
#   
#   geom_vline(aes(xintercept=tail(total_idx,1)), color=8, linetype=6, size=0.5) + 
#   
#   annotate("text", x=310, y=1200, label="12/18", size=5, color=8) +
#   theme_classic(base_family="NanumGothic") +
#   labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
#        subtitle="Smoothing Spline (Df = 10)"); p2_3_dom_pred
# predict_3_dom_pred
# # 향후 10일간 (12/18일 이후) 일일 확진자 예측값



#### 2.4 loess function ####
train_gov_lag <- lag(train_gov, 12, default=0)
total_gov_lag <- lag(total_gov, 12, default=0)

fit_4_dom <- locfit(y_train ~ train_idx + train_gov_lag, alpha=0.2)
fit_4_dom_pred <- locfit(y_total ~ total_idx + total_gov_lag, alpha=0.2)
summary(fit_4_dom); summary(fit_4_dom_pred)

predict_4_dom <- predict(fit_4_dom, data.frame(train_idx=test_idx, train_gov_lag=test_gov))
predict_4_dom_pred_0 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0))
predict_4_dom_pred_0.5 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0.5))
predict_4_dom_pred_1 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1))
predict_4_dom_pred_1.5 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1.5))
predict_4_dom_pred_2 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2))
predict_4_dom_pred_2.5 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2.5))
predict_4_dom_pred_3 <- predict(fit_4_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=3))

# train / test
p2_4_dom <- ggplot() +
  geom_point(aes(train_idx, y_train), alpha=0.5) +
  geom_line(aes(train_idx, fitted(fit_4_dom)), color=2) +
  
  geom_point(aes(test_idx, y_test), alpha=0.5) +
  geom_line(aes(test_idx, predict_4_dom), color=2, size=1, linetype=3) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)"))+
  
  geom_vline(aes(xintercept=tail(train_idx,1)), color=8, linetype=6, size=0.5) + 
  
  annotate("text", x=330, y=100, label="12/11", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Local Regression (alpha = 0.2))\nLag for Government Index = 12"); p2_4_dom
training_rmse_4_dom <- rmse(y_train, fitted(fit_4_dom)); training_rmse_4_dom
testing_rmse_4_dom <- rmse(y_test, predict_4_dom); testing_rmse_4_dom

# predict
p2_4_dom_pred <- ggplot() +
  geom_point(aes(total_idx, y_total), alpha=0.5) +
  geom_line(aes(total_idx, fitted(fit_4_dom_pred)), color=2) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)")) +
  
  geom_line(aes(pred_idx, predict_4_dom_pred_0), color=3) +
  geom_line(aes(pred_idx, predict_4_dom_pred_0.5), color=4) +
  geom_line(aes(pred_idx, predict_4_dom_pred_1), color=5) +
  geom_line(aes(pred_idx, predict_4_dom_pred_1.5), color=6) +
  geom_line(aes(pred_idx, predict_4_dom_pred_2), color=7) +
  geom_line(aes(pred_idx, predict_4_dom_pred_2.5), color=8) +
  geom_line(aes(pred_idx, predict_4_dom_pred_3), color=9) +
  
  geom_vline(aes(xintercept=tail(total_idx,1)), color=8, linetype=6, size=0.5) + 
  geom_vline(aes(xintercept=c(15, 190, 267)), color=8, linetype=3) + 
  
  annotate("text", x=330, y=100, label="12/18", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Local Regression (alpha = 0.2))\nLag for Government Index = 12"); p2_4_dom_pred
# 향후 20일간 (12/18일 이후) 일일 확진자 예측값
cbind(predict_4_dom_pred_0, predict_4_dom_pred_0.5, predict_4_dom_pred_1, 
      predict_4_dom_pred_1.5, predict_4_dom_pred_2, predict_4_dom_pred_2.5, predict_4_dom_pred_3)




#### 2.5. GAM ####
train_gov_lag <- lag(train_gov, 12, default=0)
total_gov_lag <- lag(total_gov, 12, default=0)

fit_5_dom <- gam(y_train ~ s(train_idx)  + ti(train_gov_lag), method="REML")
fit_5_dom_pred <- gam(y_total ~ s(total_idx) + ti(total_gov_lag), method="REML")
summary(fit_5_dom); summary(fit_5_dom_pred)

predict_5_dom <- predict(fit_5_dom, data.frame(train_idx=test_idx, train_gov_lag=test_gov))
predict_5_dom_pred_0 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0))
predict_5_dom_pred_0.5 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=0.5))
predict_5_dom_pred_1 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1))
predict_5_dom_pred_1.5 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=1.5))
predict_5_dom_pred_2 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2))
predict_5_dom_pred_2.5 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=2.5))
predict_5_dom_pred_3 <- predict(fit_5_dom_pred, data.frame(total_idx=pred_idx, total_gov_lag=3))

# train / test
p2_5_dom <- ggplot() +
  geom_point(aes(train_idx, y_train), alpha=0.5) +
  geom_line(aes(train_idx, fitted(fit_5_dom)), color=2) +
  
  geom_point(aes(test_idx, y_test), alpha=0.5) +
  geom_line(aes(test_idx, predict_5_dom), color=2, size=1, linetype=3) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)"))+
  
  geom_vline(aes(xintercept=tail(train_idx,1)), color=8, linetype=6, size=0.5) + 
  
  annotate("text", x=330, y=100, label="12/11", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="GAM (method = `RELM`)\nLag for Government Index = 12"); p2_5_dom
training_rmse_5_dom <- rmse(y_train, fitted(fit_5_dom)); training_rmse_5_dom
testing_rmse_5_dom <- rmse(y_test, predict_5_dom); testing_rmse_5_dom

# predict
p2_5_dom_pred <- ggplot() +
  geom_point(aes(total_idx, y_total), alpha=0.5) +
  geom_line(aes(total_idx, fitted(fit_5_dom_pred)), color=2) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)")) +
  
  geom_line(aes(pred_idx, predict_5_dom_pred_0), color=3) +
  geom_line(aes(pred_idx, predict_5_dom_pred_0.5), color=4) +
  geom_line(aes(pred_idx, predict_5_dom_pred_1), color=5) +
  geom_line(aes(pred_idx, predict_5_dom_pred_1.5), color=6) +
  geom_line(aes(pred_idx, predict_5_dom_pred_2), color=7) +
  geom_line(aes(pred_idx, predict_5_dom_pred_2.5), color=8) +
  geom_line(aes(pred_idx, predict_5_dom_pred_3), color=9) +
  
  geom_vline(aes(xintercept=tail(total_idx,1)), color=8, linetype=6, size=0.5) + 
  geom_vline(aes(xintercept=c(15, 190, 267)), color=8, linetype=3) + 
  
  annotate("text", x=330, y=100, label="12/18", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="GAM (method = `RELM`)\nLag for Government Index = 12"); p2_5_dom_pred
# 향후 20일간 (12/18일 이후) 일일 확진자 예측값
cbind(predict_5_dom_pred_0, predict_5_dom_pred_0.5, predict_5_dom_pred_1, 
      predict_5_dom_pred_1.5, predict_5_dom_pred_2, predict_5_dom_pred_2.5, predict_5_dom_pred_3)




#### 2.6. Segmented Poisson Regression ####
train_gov_lag <- lag(train_gov, 12, default=0)
total_gov_lag <- lag(total_gov, 12, default=0)

fit_6_dom <- glm(y_train ~ log(train_idx) + train_idx + train_gov_lag, family=poisson)
fit_6_dom_pred <- glm(y_total ~ log(total_idx) + total_gov_lag + total_gov_lag, family=poisson)
summary(fit_6_dom); summary(fit_6_dom_pred)

breakpoints <- c(15, 190, 267)
y1 <- y_train[breakpoints[1]:breakpoints[2]]; y2 <- y_train[breakpoints[2]:breakpoints[3]]; y3 <- y_train[breakpoints[3]:length(y_train)]
train_idx1 <- train_idx[breakpoints[1]:breakpoints[2]]; train_idx2 <- train_idx[breakpoints[2]:breakpoints[3]]; train_idx3 <- train_idx[breakpoints[3]:length(y_train)]
train_gov_lag1 <- train_gov_lag[breakpoints[1]:breakpoints[2]]; train_gov_lag2 <- train_gov_lag[breakpoints[2]:breakpoints[3]]; train_gov_lag3 <- train_gov_lag[breakpoints[3]:length(y_train)]
fit_6_dom1 <- glm(y1 ~ log(train_idx1) + train_idx1 + train_gov_lag1, family=poisson)
fit_6_dom2 <- glm(y2 ~ log(train_idx2) + train_idx2 + train_gov_lag2, family=poisson)
fit_6_dom3 <- glm(y3 ~ log(train_idx3) + train_idx3 + train_gov_lag3, family=poisson)

y_total1 <- y_total[breakpoints[1]:breakpoints[2]]; y_total2 <- y_total[breakpoints[2]:breakpoints[3]]; y_total3 <- y_total[breakpoints[3]:length(y_train)]
total_idx1 <- total_idx[breakpoints[1]:breakpoints[2]]; total_idx2 <- total_idx[breakpoints[2]:breakpoints[3]]; total_idx3 <- total_idx[breakpoints[3]:length(y_train)]
total_gov_lag1 <- total_gov_lag[breakpoints[1]:breakpoints[2]]; total_gov_lag2 <- total_gov_lag[breakpoints[2]:breakpoints[3]]; total_gov_lag3 <- total_gov_lag[breakpoints[3]:length(y_train)]
fit_6_dom_pred1 <- glm(y_total1 ~ log(total_idx1) + total_idx1 + total_gov_lag1, family=poisson)
fit_6_dom_pred2 <- glm(y_total2 ~ log(total_idx2) + total_idx2 + total_gov_lag2, family=poisson)
fit_6_dom_pred3 <- glm(y_total3 ~ log(total_idx3) + total_idx3 + total_gov_lag3, family=poisson)

predict_6_dom <- exp(predict(fit_6_dom3, data.frame(train_idx3=test_idx, train_gov_lag3=test_gov)))
predict_6_dom_pred_0 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=0)))
predict_6_dom_pred_0.5 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=0.5)))
predict_6_dom_pred_1 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=1)))
predict_6_dom_pred_1.5 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=1.5)))
predict_6_dom_pred_2 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=2)))
predict_6_dom_pred_2.5 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=2.5)))
predict_6_dom_pred_3 <- exp(predict(fit_6_dom_pred3, data.frame(total_idx3=pred_idx, total_gov_lag3=3)))

# train / test
p2_6_dom <- ggplot() +
  geom_point(aes(train_idx, y_train), alpha=0.5) +
  geom_line(aes(train_idx1, fitted(fit_6_dom1)), color=2) +
  geom_line(aes(train_idx2, fitted(fit_6_dom2)), color=2) +
  geom_line(aes(train_idx3, fitted(fit_6_dom3)), color=2) +
  
  geom_point(aes(test_idx, y_test), alpha=0.5) +
  geom_line(aes(test_idx, predict_6_dom), color=2, size=1, linetype=3) +
  
  geom_line(aes(train_idx, train_gov*800), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./800, name="사회적 거리두기 (단계)"))+
  
  geom_vline(aes(xintercept=tail(train_idx,1)), color=8, linetype=6, size=0.5) + 
  geom_vline(aes(xintercept=c(15, 190, 267)), color=8, linetype=3) + 
  
  annotate("text", x=330, y=100, label="12/11", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Segmented Poisson (Breakpoint = 15, 190, 267)\nLag for Government Index = 12"); p2_6_dom
training_rmse_6_dom <- rmse(y_train, fitted(fit_6_dom)); training_rmse_6_dom
testing_rmse_6_dom <- rmse(y_test, predict_6_dom); testing_rmse_6_dom

# predict
p2_6_dom_pred <- ggplot() +
  geom_point(aes(total_idx, y_total), alpha=0.5) +
  geom_line(aes(total_idx1, fitted(fit_6_dom_pred1)), color=2) +
  geom_line(aes(total_idx2, fitted(fit_6_dom_pred2)), color=2) +
  geom_line(aes(total_idx3, fitted(fit_6_dom_pred3)), color=2) +
  
  geom_line(aes(train_idx, train_gov*7000), color=4) +
  scale_y_continuous(sec.axis = sec_axis(~./7000, name="사회적 거리두기 (단계)")) +
  
  geom_line(aes(pred_idx, predict_6_dom_pred_0), color=3) +
  geom_line(aes(pred_idx, predict_6_dom_pred_0.5), color=4) +
  geom_line(aes(pred_idx, predict_6_dom_pred_1), color=5) +
  geom_line(aes(pred_idx, predict_6_dom_pred_1.5), color=6) +
  geom_line(aes(pred_idx, predict_6_dom_pred_2), color=7) +
  geom_line(aes(pred_idx, predict_6_dom_pred_2.5), color=8) +
  geom_line(aes(pred_idx, predict_6_dom_pred_3), color=9) +
  
  geom_vline(aes(xintercept=tail(total_idx,1)), color=8, linetype=6, size=0.5) + 
  geom_vline(aes(xintercept=c(15, 190, 267)), color=8, linetype=3) + 
  
  annotate("text", x=330, y=100, label="12/18", size=5, color=8) +
  theme_classic(base_family="NanumGothic") +
  labs(x="날짜 (월)", y="일일 확진자", title="COVID-19 전국 신규 확진자",
       subtitle="Segmented Poisson (Breakpoint = 15, 190, 267)\nLag for Government Index = 12"); p2_6_dom_pred
# 향후 20일간 (12/18일 이후) 일일 확진자 예측값
cbind(predict_6_dom_pred_0, predict_6_dom_pred_0.5, predict_6_dom_pred_1, 
      predict_6_dom_pred_1.5, predict_6_dom_pred_2, predict_6_dom_pred_2.5, predict_6_dom_pred_3)

