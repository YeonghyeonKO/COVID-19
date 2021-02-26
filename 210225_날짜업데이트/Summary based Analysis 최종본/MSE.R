library(dplyr)
library(ggplot2)
library(nls2)

getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/210102_mse/")

# Preprocessing data (2020/08/31)
df <- read.csv("최종 정리(0831)/COVID-19-geographic-disbtribution-worldwide.csv")
max_date <- set_date("2020/8/31") 
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]
df_sum = preprocessing_data()[1:245,]
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}

break_result <- read.csv("최종 정리(0831)/categorizaioin_result.csv")
coef_result <- read.csv("최종 정리(0831)/coef_result.csv")

eq_log <- Y ~ a/(1 + exp(b - c*X))
mse = data.frame(country=country, mse=numeric(213), scaled_mse=numeric(213))


for (i in country){
  i = "Sri_Lanka"
  coef <- coef_result[which(coef_result$X==i),]
  
  # 1st Segment
  if (length(coef$a1_Logi)==0 || is.na(coef$a1_Logi)){
    mse[which(mse$country==i),"mse"] = NA
    mse[which(mse$country==i),"scaled_mse"] = NA
    next
  }
  Date <- 1:(max_date+1)
  Cases <- cumsum(df_sum[,i])
  startpoint <- coef$startpoint
  breakpoint <- coef$breakpoint
  a1 <- coef$a1_Logi; b1 <- coef$b1_Logi; c1 <- coef$c1_Logi
  fit_Cases1 <- a1 / (1 + exp(b1 - c1*(Date - startpoint)))

  plot(Date, Cases, main=i)
  lines(Date, fit_Cases1, col=2)
  
  
  
  # 2nd Segment
  if (length(coef$a2_Logi)==0 || is.na(coef$a2_Logi)){
    next
  }
  intercept <- Cases[breakpoint]
  a2 <- coef$a2_Logi; b2 <- coef$b2_Logi; c2 <- coef$c2_Logi
  fit_Cases2 <- a2 / (1 + exp(b2 - c2*(Date - breakpoint))) + intercept
  lines(Date, fit_Cases2)
  
  # mse[which(mse$country==i),"mse"] = mean((Cases - fit_Cases)^2)
  # mse[which(mse$country==i),"scaled_mse"] = mean((scale(Cases) - scale(fit_Cases))^2)
  # cat("\n#### ", i, " ####\n")
  # cat("MSE : ", mse[which(mse$country==i),"mse"], "\n")
  # cat("Scaled MSE : ", mse[which(mse$country==i),"scaled_mse"], "\n")
}

mse
