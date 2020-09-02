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

path <- "C:/Users/kyh/Desktop/고영현/인턴십/해외 코로나/200902 해외 데이터 smoothing segmented poisson/"
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

# for (j in seq(from = 3, by = 2, length.out = 6)){
#     assign(paste0("cases_smoothing_", j), pracma::movavg(df$cases, j, type = "s"))
#     df <- cbind(df, get(paste0("cases_smoothing_", j)))
#     colnames(df)[ncol(df)] = paste0("cases_smoothing_", j)
# }



# segmented poisson에서 적합되지 않는 국가들 list
country2 <- c("China", "Japan", "Nepal", "Puerto_Rico", "Tunisia", "New_Zealand", 
              "Burkina_Faso", "Comoros", "Equatorial_Guinea", "Botswana", "Sri_Lanka")
            

`%notin%` <- purrr::negate(`%in%`)

num <- df %>%
    group_by(country) %>%
    mutate(total = cumsum(cases)) %>%
    filter(total > 0, country %notin% country2) %>%
    summarise(count = n())



###########################################################
#### Poisson Model과 Smoothing Poisson Model 계수 비교 ####
###########################################################
set.seed(1234)

for (j in seq(from = 3, by = 2, length.out = 6)){
    n <- 0
    coef_Poi <- data.frame(matrix(NA, ncol = 6, nrow = nrow(num)))
    coef_seg_Poi <- data.frame(matrix(NA, ncol = 14, nrow = nrow(num)))
    rownames(coef_Poi) <- num$country
    rownames(coef_seg_Poi) <- num$country

    assign(paste0("cases_smoothing"), pracma::movavg(df$cases, j, type = "s"))
    df[5] <- cases_smoothing
    colnames(df)[5] = paste0("cases_smoothing")
    
    for (k in num$country){
        n <- n + 1
        days <- num$count[n]
        
        data_s <- df %>%
            filter(country == k) %>%
            mutate(total = cumsum(cases)) %>%
            filter(total >= 1) %>%
            mutate(Days_after_first_Case = 1:days) %>%
            filter(date < "2020-07-11")
        
        if (data_s$total[days] < 300){
            next
        }
        
        # Poisson Model
        fit <- glm(cases ~ log(Days_after_first_Case) + Days_after_first_Case, data = data_s, family = poisson)
        fit_s <- glm(cases_smoothing ~ log(Days_after_first_Case) + Days_after_first_Case, data = data_s, family = poisson)
        #fit_s <- glm(as.numeric(unlist(data_s[(j+7)/2])) ~ log(Days_after_first_Case) + Days_after_first_Case, data = data_s, family = poisson)
        
        coef_Poi[k, ] <- c(fit$coefficients, fit_s$coefficients)
        
        # Segmented Model
        tryCatch(
            expr = {
                seg_fit <- segmented(fit, seg.Z = ~ log(Days_after_first_Case) + Days_after_first_Case, npsi = 1)
                psi <- round(seg_fit$psi[2])
                fit1 <- glm(cases[1:psi] ~ log(Days_after_first_Case[1:psi]) + Days_after_first_Case[1:psi], 
                            data = data_s, family = poisson)
                fit2 <- glm(cases[psi+1:num$count[n]] ~ log(Days_after_first_Case[psi+1:num$count[n]]) + Days_after_first_Case[psi+1:num$count[n]], 
                            data = data_s, family = poisson)
                coef_seg_Poi[k, 1:7] <- c(summary(fit1)$coefficients[,1], summary(fit2)$coefficients[,1], psi)
            },
            
            error = function(e){
                next
            },
            
            finally = {
                NULL
            }
        )
        
        
        tryCatch(
            expr = {
                seg_fit_s <- segmented(fit_s, seg.Z = ~ log(Days_after_first_Case) + Days_after_first_Case, npsi = 1)
                psi_s <- round(seg_fit_s$psi[2])
                fit1_s <- glm(cases[1:psi_s] ~ log(Days_after_first_Case[1:psi_s]) + Days_after_first_Case[1:psi_s], 
                              data = data_s, family = poisson)
                fit2_s <- glm(cases[psi_s+1:num$count[n]] ~ log(Days_after_first_Case[psi_s+1:num$count[n]]) + Days_after_first_Case[psi_s+1:num$count[n]], 
                              data = data_s, family = poisson)
                coef_seg_Poi[k, 8:14] <- c(summary(fit1_s)$coefficients[,1], summary(fit2_s)$coefficients[,1], psi_s)
            },
            
            error = function(e){
                next
            },
            
            finally = {
                NULL
            }
        )
    }
    
    colnames(coef_Poi) <- c(names(fit$coefficients), paste0(names(fit$coefficients), "_smoothing"))
    colnames(coef_seg_Poi) <- c("(Intercept)_1", "log(Days_after_first_Case)_1", "Days_after_first_Case_1",
                                "(Intercept)_2", "log(Days_after_first_Case)_2", "Days_after_first_Case_2",	"psi",
                                "(Intercept)_1_s", "log(Days_after_first_Case)_1_s", "Days_after_first_Case_1_s", 
                                "(Intercept)_2_s", "log(Days_after_first_Case)_2_s", "Days_after_first_Case_2_s", "psi_s")
    
    
    ########################
    #### .csv 파일 저장 ####
    ########################
    
    write.csv(coef_Poi, paste0(result_path, "Coef_Poi_", j, ".csv"), row.names = T)
    write.csv(coef_seg_Poi, paste0(result_path, "Coef_seg_Poi_", j, ".csv"), row.names = T)
}
