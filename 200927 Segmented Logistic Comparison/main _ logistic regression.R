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
# i = "Canada"
# get(paste0("df_sum$", "Canada"))
# df_sum$i
# df_sum$(noquote(i))
# as.expression("Canada")
# df_sum$as.name("Canada")
# 
# df_sum$Canada


save_path = "C:/Users/kyh/GitHub/COVID-19/200927 Segmented Logistic Comparison/Plot"


for (i in 1:8){
  ind <- which(colnames(df_sum)==coef_2$country[i])
  x <- coef_2$start[i] : coef_2$end[i]
  y_seg <- coef_2$a2_seg[i] / (1 + exp(coef_2$b2_seg[i] - coef_2$c2_seg[i] * (x- coef_2$initial[i])))
  y_sep <- coef_2$a2_sep[i] / (1 + exp(coef_2$b2_sep[i] - coef_2$c2_sep[i] * (x- coef_2$initial[i])))
  
  
  p <- ggplot() +
    geom_point(aes(coef_2$initial[i]:length(df_sum[,ind]), cumsum(df_sum[,ind])[coef_2$initial[i]:length(df_sum[,ind])]),
               alpha = 0.1) +
    geom_point(aes(x, cumsum(df_sum[,ind])[x]), col = 1, alpha = 0.1) +
    geom_vline(aes(xintercept = c(coef_2[i,]$start, coef_2[i,]$end)), col = 2, linetype = 4) +
    
    # segmented Logistic model
    geom_line(aes(x, y_seg - y_seg[1] + cumsum(df_sum[,ind])[coef_2$start[i]]), col = 2, size = 1.7) +
    
    # separated Logistic model
    geom_line(aes(x, y_sep - y_sep[1] + cumsum(df_sum[,ind])[coef_2$start[i]]), col = 4, size = 1.7) +
    
    labs(title = paste0("2nd Segmented vs Separated Logistic Model"),
         subtitle = paste0(coef_2[i,1], "  / breakpoint : ", coef_2[i,]$start, " ~ ", coef_2[i,]$end),
         x = paste0(df_sum[coef_2$initial[i], 1], "~", df_sum[length(df_sum[,ind]), 1]),
         y = "Cumulative Confirmed Cases") +
    annotate("text", x=coef_2$initial[i]*1.32, y = max(cumsum(df_sum[,ind])), label="— Segmented", col = 2, size=5) +
    annotate("text", x=coef_2$initial[i]*1.3, y = max(cumsum(df_sum[,ind]))*0.9, label="— Separated", col = 4, size=5) +
    # geom_segment(aes(x=coef_2$initial[i], y=max(cumsum(df_sum[,ind])), xend=coef_2$initial[i]*1.2, yend=max(cumsum(df_sum[,ind]))),color = 2, size=1.5) +
    # geom_segment(aes(x=coef_2$initial[i], y=max(cumsum(df_sum[,ind]))*0.9, xend=coef_2$initial[i]*1.2, yend=max(cumsum(df_sum[,ind]))*0.9),color = 4, size=1.5) +
    
    annotate("rect", xmin = coef_2[i,]$start, xmax = coef_2[i,]$end,
             ymin = cumsum(df_sum[,ind])[coef_2[i,]$start], ymax = cumsum(df_sum[,ind])[coef_2[i,]$end], 
             alpha=0.05, fill="red") +
  
    theme_bw() + 
    theme(
      plot.title=element_text(size=20, hjust=0.5, face="bold", colour="black", vjust=2),
      plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
      axis.text=element_text(size=14),
      axis.title=element_text(size=14),
      axis.text.x = element_text(hjust = 1, size = 12))
  
  ggsave(p, filename = paste0("2nd Segmented vs Separated Logistic Model in ", coef_2$country[i], ".png"), 
         path = save_path, width = 12, height = 7, dpi = 300)
  
  
  
  y_seg - y_seg[1] + cumsum(df_sum[,ind])[coef_2$start[i]]
  y_sep - y_sep[1] + cumsum(df_sum[,ind])[coef_2$start[i]]
  
}





ggplot() +
  geom_histogram(aes(coef_2[,2] - coef_2[, 5], color = country))

ggplot() +
  

ggplot() +
  geom_histogram(aes(coef_2[,3] - coef_2[, 6]))

ggplot() +
  geom_histogram(aes(coef_2[,4] - coef_2[, 7]))


eq <- a / (1 + exp(b - c*x))

