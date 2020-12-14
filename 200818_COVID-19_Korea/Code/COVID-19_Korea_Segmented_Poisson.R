#install.packages("segmented")
library(segmented)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggrepel)

path <- "C:/Users/김학용/Desktop/COVID-19_Korea_200818[277]/"
data_path <- paste0(path, "Data/")
save_path <- paste0(path, "Plot/")

# COVID_19_Confirmed_Korea_Merged_fin: From Kaggle & Public Data Portal
df <- read.csv(paste0(data_path, "COVID_19_Confirmed_Korea_Merged_fin.csv"))

# COVID_19_Confirmed_Korea_Stat: From KOSTAT
#df <- read.csv(paste0(data_path, "COVID_19_Confirmed_Korea_Stat.csv"))

############################################ 0. EDA ############################################

p_1_1 <- ggplot(df, aes(x = Days_after_First_Cases, y = Cases, col = Country_Region)) +
  geom_line() +
  geom_text_repel(
    data = subset(df, Days_after_First_Cases == max(Days_after_First_Cases)),
    aes(x = max(Days_after_First_Cases), label = Country_Region),
    size = 3,
    nudge_x = 45,
    segment.color = NA
  ) +
  labs(title = "COVID-19 Confirmed Cases", subtitle = paste0("Korea, ", as.character(min(as.Date(df$Date))), "~", as.character(max(as.Date(df$Date)))), x = "Days from First Case", y = "Confirmed Cases")+
  theme_bw() +
  theme(
    plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
    plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    axis.text.x = element_text(hjust = 1, size = 12),
    legend.position = "none")
p_1_1
#ggsave(p_1_1, filename = "Confirmed_case_in_Korea_by_Region_total.png", path = save_path, width = 12, height = 7, dpi = 300)

p_1_2 <- ggplot(filter(df, Cases > 1), aes(x = Days_after_First_Cases, y = log10(Cases), col = Country_Region)) +
  geom_line() +
  geom_text_repel(
    data = subset(df, Days_after_First_Cases == max(Days_after_First_Cases)),
    aes(x = max(Days_after_First_Cases), label = Country_Region),
    size = 3,
    nudge_x = 45,
    segment.color = NA
  ) +
  labs(title = "COVID-19 Confirmed Cases", subtitle = paste0("Korea, ", as.character(min(as.Date(df$Date))), "~", as.character(max(as.Date(df$Date)))), x = "Days from First Case", y = "Log10(Confirmed Cases)")+
  theme_bw() +
  theme(
    plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
    plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    axis.text.x = element_text(hjust = 1, size = 12),
    legend.position = "none")
p_1_2
#ggsave(p_1_2, filename = "Confirmed_case_in_Korea_log_by_Region_total.png", path = save_path, width = 12, height = 7, dpi = 300)


p_2_1 <- ggplot(df, aes(x = Days_after_First_Cases, y = Difference, col = Country_Region)) +
  geom_line() +
  geom_text_repel(
    data = subset(df, Days_after_First_Cases == max(Days_after_First_Cases)),
    aes(x = max(Days_after_First_Cases), label = Country_Region),
    size = 3,
    nudge_x = 45,
    segment.color = NA
  ) +
  labs(title = "COVID-19 Daily Cases", subtitle = paste0("Korea, ", as.character(min(as.Date(df$Date))), "~", as.character(max(as.Date(df$Date)))), x = "Days from First Case", y = "Confirmed Cases")+
  theme_bw() +
  theme(
    plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
    plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    axis.text.x = element_text(hjust = 1, size = 12),
    legend.position = "none")
p_2_1
#ggsave(p_2_1, filename = "Daily_case_in_Korea_by_Region_total.png", path = save_path, width = 12, height = 7, dpi = 300)






province <- unique(as.character(df$Country_Region))

{
  n <- 0
  coef_Poi <- data.frame(); coef_seg_Poi <- data.frame()
  current_max <- c()
  pred_week_4 <- c()
  pred_week_8 <- c()
}

for(p in province){

  n <- n+1
  df_sub <- filter(df, Country_Region == p)

  ### Poisson regression
  fit <- glm(Difference ~ log(Days_after_First_Cases) + Days_after_First_Cases, data = df_sub, family = poisson)
  
  SSE_Poi <- sum((predict(fit, data.frame(Days_after_First_Cases = df_sub$Days_after_First_Cases), type = "response") - df_sub$Difference)^2)
  SST_Poi <- sum((df_sub$Difference - mean(df_sub$Difference))^2)
  R_Poi <- 1 - SSE_Poi / SST_Poi
  
  for (k in 1:3){coef_Poi[k, n] <- summary(fit)$coefficients[k, 1]}
  coef_Poi[4, n] <- R_Poi
  colnames(coef_Poi)[n] = p
  rownames(coef_Poi) <- c(rownames(summary(fit)$coefficients)[1:3], "R-Square (Daily)")
  
  
  ### Segmented Poisson regression coefficients & R-Sqaure  
  seg_fit <- segmented(fit, seg.Z = ~ log(Days_after_First_Cases) + Days_after_First_Cases, npsi = 1) # npsi: break point의 개수
  psi <- round(seg_fit$psi[2])
  
  fit1 <- glm(Difference[1:psi] ~ log(Days_after_First_Cases[1:psi]) + Days_after_First_Cases[1:psi], 
              data = df_sub, family = poisson)
  fit2 <- glm(Difference[psi+1:nrow(df_sub)] ~ log(Days_after_First_Cases[psi+1:nrow(df_sub)]) + Days_after_First_Cases[psi+1:nrow(df_sub)], 
              data = df_sub, family = poisson)
  
  SSE_seg_Poi <-sum((predict(seg_fit, data.frame(Days_after_First_Cases = df$Days_after_First_Cases), type = "response") - df_sub$Difference)^2)
  SST_seg_Poi <- sum((df_sub$Difference - mean(df_sub$Difference))^2)
  R_seg_Poi <- 1 - SSE_seg_Poi / SST_seg_Poi
  
  for (k in 1:3){coef_seg_Poi[k, n] <- summary(fit1)$coefficients[k, 1]}
  tryCatch(
    expr = {
      for (k in 4:6){coef_seg_Poi[k, n] <- summary(fit2)$coefficients[k-3, 1]}
    },
    error = function(e){
      for (k in 4:6){coef_seg_Poi[k, n] <- NA}
    },
    finally = NULL
  )
  coef_seg_Poi[7, n] <- psi
  coef_seg_Poi[8, n] <- R_seg_Poi
  colnames(coef_seg_Poi)[n] = p
  
  
  max_date <- max(as.Date(df_sub$Date))
  week_4_date <- max(as.Date(df_sub$Date))+28
  week_8_date <- max(as.Date(df_sub$Date))+56
  
  df_sub$Date <- as.character(df_sub$Date)
  df_sub_new <- data.frame(Date = as.character(as.Date(origin = max(as.Date(df_sub$Date)), 1:56)),
                           Case_Type = "Confirmed",
                           Cases = NA,
                           Difference = NA,
                           Country_Region = as.character(df_sub$Country_Region[1]),
                           Days_after_First_Cases = (max(df_sub$Days_after_First_Cases)+1):max(df_sub$Days_after_First_Cases+56))
  df_sub_ch <- rbind(df_sub,df_sub_new)
  
  
  pred_pois <- data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases,
                          pred = exp(predict(fit, df_sub_ch)))
  
  pred_seg_pois <- data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases,
                          pred = exp(predict(seg_fit, df_sub_ch)))
  
  current_max <- c(current_max, max(df_sub$Cases))
  pred_week_4 <- c(pred_week_4, pred_seg_pois$pred[pred_seg_pois$Days_after_First_Cases == df_sub_ch$Days_after_First_Cases[df_sub_ch$Date == week_4_date]])
  pred_week_8 <- c(pred_week_8, pred_seg_pois$pred[pred_seg_pois$Days_after_First_Cases == df_sub_ch$Days_after_First_Cases[df_sub_ch$Date == week_8_date]])
  
  ### 일일 확진자 plot
  Daily <- ggplot(df_sub_ch, aes(x = Days_after_First_Cases, y = Difference)) +
    geom_point(alpha = 0.8, size = 2, color = "grey") +
    scale_y_continuous(limits = c(0, 1.5 * max(df_sub$Difference)))+
    geom_line(pred_pois, mapping = aes(x = Days_after_First_Cases, y = pred), col = "midnightblue") +
    geom_line(pred_seg_pois, mapping = aes(x = Days_after_First_Cases, y = pred), col = "tomato") +
    geom_vline(xintercept = df_sub_ch$Days_after_First_Cases[df_sub_ch$Date == max_date], 
               col = "grey", size = 1, linetype = 2, alpha = 0.8) +
    labs(x = NULL, y = "Daily Cases") +
    theme_bw()+
    theme(axis.title.y = element_text(size = 15))
  Daily
  
  ### 누적 확진자 plot
  Cumulative <- ggplot(df_sub_ch, aes(x = Days_after_First_Cases, y = Cases)) +
    geom_point(alpha = 0.8, size = 2, color = "grey") +
    scale_y_continuous(limits = c(0, 1.5 * max(df_sub$Cases)))+
    geom_line(pred_pois, mapping = aes(x = Days_after_First_Cases, y = cumsum(pred)), col = "midnightblue") +
    geom_line(pred_seg_pois, mapping = aes(x = Days_after_First_Cases, cumsum(pred)), col = "tomato") +
    geom_vline(xintercept = df_sub_ch$Days_after_First_Cases[df_sub_ch$Date == max_date], 
               col = "grey", size = 1, linetype = 2, alpha = 0.8) +
    labs(x = NULL, y = "Cumulative Cases") +
    theme_bw()+
    theme(axis.title.y = element_text(size = 15))
  Cumulative
  
  ### 최종 plot
  p_merged <- grid.arrange(Daily, Cumulative,
                      top = textGrob(p, rot=0, gp = gpar(fontsize=20, bold=T)),
                      bottom = textGrob(paste0("Days From ",as.character(min(as.Date(df_sub$Date)))), rot=0, gp = gpar(fontsize=15)), nrow = 2)
  p_merged
 # ggsave(paste0(save_path,"Segmented Poisson in ", p, ".png"), p_merged, dpi = 300, width = 8, height = 6)
}


# df_pred: Predictive Daily Cases
df_pred <- data.frame(Max_Cases = current_max,
                     Predict_Week_4 = pred_week_4,
                     Predict_Week_8 = pred_week_8)
row.names(df_pred) <- province

row.names(coef_Poi) <- c("(Intercept)", "log(Days)", "Days", "R-square")

row.names(coef_seg_Poi) <- c("(Intercept)_1", "log(Days)_1", "Days_1",
                             "(Intercept)_2", "log(Days)_2", "Days_2",
                             "Break_Point", "R-square")

write.csv(coef_Poi, paste0(save_path, "Poisson_Korea.csv"), row.names = T) # Poisson Coefficients
write.csv(coef_seg_Poi, paste0(save_path, "Segmented_Poisson_Korea.csv"), row.names = T) # Segmented Poisson Coefficients
write.csv(df_pred, paste0(save_path, "Segmented_Poisson_Korea_Predict.csv"), row.names = T) # Predicted Daily Cases



unique(df$Country_Region)
library(forecast)

# 수도권
y1<-filter(df,Country_Region == "수도권(서울+경기+인천)")$Difference
gkf1<-ksmooth(1:length(y1),y1,kernel="normal",bandwidth = 14)$y
ma1<-ma(y1,7)
plot(y1)
lines(gkf1,col=2,lwd=2)
lines(ma1,col=4)

plot(cumsum(y1))
lines(cumsum(gkf1),col=2)
lines(4:(3+length(na.omit(ma1))),cumsum(na.omit(ma1)),col=4)

#poisson regression
x<-1:length(gkf1)
fit1_gkf1 <- glm(gkf1 ~ log(x) + x,  family = poisson)
seg1_gkf1<-segmented(fit1_gkf1, seg.Z = ~ log(x) + x, npsi = 3)

y_test<-gkf1[195:211]
x_test<-1:length(y_test)
log_x_test<-log(x_test)
fit_test <- glm(y_test ~ log_x_test + x_test,  family = poisson)

plot(y_test,xlim=c(0,1000),ylim=c(0,500))
lines(fitted(fit_test),col=2)
lines(exp(predict.glm(fit_test,data.frame(log_x_test=log(1:1000),x_test=1:1000))),col=4,lty=2)
P<-function(p,t) exp(p[1])*(t^(p[2]))*exp(p[3]*t) 
P(fit_test$coefficients,1:1000)



plot(gkf1,xlim=c(0,300),ylim=c(0,1000))
lines(fitted(fit1_gkf1),col=2)
lines(fitted(seg1_gkf1),col=4,lwd=2)
lines(exp(predict.segmented(seg1_gkf1,data.frame(x=c(1:300)))),col=4,lty=2)
predict(seg1_gkf1,)

# 비수도권
y2<-filter(df,Country_Region == "비수도권")$Difference
gkf2<-ksmooth(1:length(y2),y2,kernel="normal",bandwidth = 14)
ma2<-ma(y2,7)
plot(y2)
lines(gkf2,col=2,lwd=2)
lines(ma2,col=4,lwd=2,lty=2)


plot(cumsum(y2))
lines(cumsum(gkf2$y),col=2,lwd=2)
lines(4:(3+length(na.omit(ma2))),cumsum(na.omit(ma2)),col=4,lwd=2,lty=2)

# 전국

y3<-filter(df,Country_Region == "전국(지역)")$Difference
gkf3<-ksmooth(1:length(y3),y3,kernel="normal",bandwidth = 14)$y
ma3<-ma(y3,7)
plot(y3)
lines(gkf3,col=2,lwd=2)
lines(ma3,col=4)

plot(cumsum(y1))
lines(cumsum(gkf1),col=2)
lines(4:(3+length(na.omit(ma1))),cumsum(na.omit(ma1)),col=4)









library(nls2)

  y<-y_test
  x<-1:length(y_test)
  ## grid setting for initial value
  st = data.frame(a = c(0, 10), b = c(0, 10), c = c(-10, 0))
  fit11 = nls2(y ~ exp(a)*(x^(b))*exp(c*x), start = st,
               algorithm = "plinear-random",control = nls.control(maxiter = 3000))
  initial1 = coef(fit11)
  initial1 = coef(fit11)
  ## fitting
  nlsfit12 = try(nls(y ~ exp(a)*(x^(b))*exp(c*x), 
                  start = list(a = initial1[1],b = initial1[2],c = initial1[3]),
                  upper=c(Inf,2,0),
                  algorithm = "port"))

  plot(x,y,type="p",xlim=c(0,100),ylim=c(0,500))
  lines(x,(fitted(fit12)),col=2)
  lines(predict(fit12,data.frame(x=c(1:100))),col=4,lty=2)
  
  
  
  
  logistic_country <- c(logistic_country,i)
  if(is.null(logistic_result)){
    logistic_result <- coef(fit12)
  }else{
    logistic_result <- rbind(logistic_result,coef(fit12))
  }
}




























