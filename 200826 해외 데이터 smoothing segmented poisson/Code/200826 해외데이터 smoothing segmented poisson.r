#################
#### library ####
#################

library(segmented)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggrepel)



############################
#### data preprocessing ####
############################

path <- "C:/Users/kyh/Desktop/고영현/인턴십/해외 코로나/200826 해외 데이터 smoothing segmented poisson/"
data_path <- paste0(path, "Data/")
save_path <- paste0(path, "Plot/")
result_path <- paste0(path, "Result/")

df <- read.csv(paste0(data_path, "corona_ecdc.csv")) %>%
    mutate(date = make_datetime(year, month, day),
           country = countriesAndTerritories) %>%
    filter(date < "2020-07-11",) %>%
    select(country, date, cases, deaths) %>%
    arrange(country, date)

for (i in which(df$cases < 0)){df$cases[i] = 0} # Case < 0 처리 (0으로 대체)



##################################################################
#### Simple Moving Average (window size = 3, 5, 7, 9 ,11, 13) ####
##################################################################

for (j in seq(from = 3, by = 2, length.out = 6)){
    assign(paste0("cases_smoothing_", j), pracma::movavg(df$cases, j, type = "s"))
    df <- cbind(df, get(paste0("cases_smoothing_", j)))
    colnames(df)[ncol(df)] = paste0("cases_smoothing_", j)
}

# # segmented poisson에서 적합되지 않는 국가들 list
# country2 <- c("China", "Japan", "Nepal", "Puerto_Rico", "Tunisia", "New_Zealand", "Burkina_Faso")
# 
# `%notin%` <- purrr::negate(`%in%`)

num <- df %>%
    group_by(country) %>%
    mutate(total = cumsum(cases)) %>%
    filter(total > 0) %>%
    # filter(total > 0, country %notin% country2) %>%
    summarise(count = n())


n <- 0; coef_Poi <- data.frame(matrix(NA, ncol = 6, nrow = nrow(num))); coef_Poi_smoothing <- data.frame()
colnames(coef_Poi) <- c(names(fit$coefficients), paste0(names(fit$coefficients), "_smoothing"))
rownames(coef_Poi) <- num$country



###########################################################
#### Poisson Model과 Smoothing Poisson Model 계수 비교 ####
###########################################################

for (k in num$country){
    n <- n + 1
    days <- num$count[n]
    
    data_0710 <- df %>%
        filter(country == k) %>%
        mutate(total = cumsum(cases)) %>%
        filter(total >= 1) %>%
        mutate(Days_after_first_Case = 1:days) %>%
        filter(date < "2020-07-11")
    
    if (data_0710$total[days] < 300){
        next
    }
    
    fit <- glm(cases ~ log(Days_after_first_Case) + Days_after_first_Case, data = data_0710, family = poisson)
    fit_s <- glm(cases_smoothing_7 ~ log(Days_after_first_Case) + Days_after_first_Case, data = data_0710, family = poisson)
    
    coef_Poi[k, ] <- c(fit$coefficients, fit_s$coefficients)
}



########################
#### .csv 파일 저장 ####
########################

write.csv(coef_Poi, paste0(result_path, "Coef_Poi.csv"), row.names = T)
