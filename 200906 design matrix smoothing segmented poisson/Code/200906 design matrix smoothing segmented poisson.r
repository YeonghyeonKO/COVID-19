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

path <- "C:/Users/kyh/Desktop/고영현/인턴십/해외 코로나/200906 design matrix smoothing segmented poisson/"
data_path <- paste0(path, "Data/")
save_path <- paste0(path, "Plot/")
result_path <- paste0(path, "Result/")

df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv",
               na.strings = "", fileEncoding = "UTF-8-BOM") %>%
    mutate(date = make_datetime(year, month, day),
           country = countriesAndTerritories) %>%
    filter(date < "2020-09-01",) %>%
    select(country, date, cases, deaths) %>%
    arrange(country, date)

for (i in which(df$cases < 0)){df$cases[i] = 0} # Case < 0 처리 (0으로 대체)



##################################################################
#### Simple Moving Average (window size = 3, 5, 7, 9 ,11, 13) ####
##################################################################

# segmented poisson에서 적합되지 않는 국가들 list
country2 <- c("China", "Japan", "Nepal", "Puerto_Rico", "Tunisia", 
              "New_Zealand", "Burkina_Faso", "Comoros", "Equatorial_Guinea", 
              "Botswana", "Sri_Lanka", "Western_Sahara")
            

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
    coef_seg_Poi <- data.frame(matrix(NA, ncol = 12, nrow = nrow(num)))
    rownames(coef_Poi) <- num$country
    rownames(coef_seg_Poi) <- num$country
    
    for (k in num$country){
        n <- n + 1
        days <- num$count[n]
        
        data_s <- df %>%
            filter(country == k) %>%
            mutate(total = cumsum(cases),
                   cases_smoothing = pracma::movavg(cases, j, type = "s")) %>%
            filter(total >= 1) %>%
            mutate(Day = 1:days)
        
        if (data_s$total[days] < 1000){
            next
        }
        
        # Poisson Model (original & smoothing)
        fit <- glm(cases ~ Day + log(Day), data = data_s, family = poisson)
        fit_s <- glm(cases_smoothing ~ Day + log(Day), data = data_s, family = poisson)
        coef_Poi[k, ] <- c(fit$coefficients, fit_s$coefficients)
        
        
        # Segmented Poisson Model (original & smoothing)
        tryCatch(
            expr = {
                # design matrix for original data
                seg_fit <- segmented(fit, seg.Z = ~ log(Day) + Day, npsi = 1)
                psi <- round(seg_fit$psi[2])
                
                X <- data_s$Day
                X_design = data.frame(cbind(X, log(X + 1), 
                                            c(rep(0, psi-1), X[psi:days] - X[psi-1]), 
                                            c(rep(0, psi-1), log(X[psi:days] - X[psi-1] + 1))))
                
                
                # design matrix for smoothing data
                seg_fit_s <- segmented(fit_s, seg.Z = ~ log(Day) + Day, npsi = 1)
                psi_s <- round(seg_fit_s$psi[2])
                
                X_design_s = data.frame(cbind(X, log(X + 1), 
                                            c(rep(0, psi_s-1), X[psi_s:days] - X[psi_s-1]), 
                                            c(rep(0, psi_s-1), log(X[psi_s:days] - X[psi_s-1] + 1))))

                                
                # Segmented Poisson Model (original & smoothing)
                seg_fit <- try(glm(data_s$cases ~  . , family = poisson(), data = X_design))
                seg_fit_s <- try(glm(data_s$cases_smoothing ~  . , family = poisson(), data = X_design_s))
                
                coef_seg_Poi[k, 1:6] <- c(summary(seg_fit)$coefficients[,1], psi)
                coef_seg_Poi[k, 7:12] <- c(summary(seg_fit_s)$coefficients[,1], psi_s)
            },
            
            error = function(e){
                next
            },
            
            finally = {
                NULL
            }
        )
        
        p <- ggplot() +
                # Segmented Poisson plot
                geom_point(aes(data_s$Day, data_s$cases), alpha = 0.2) +
                geom_line(aes(data_s$Day, fitted(seg_fit)), size = 1) +
                geom_vline(aes(xintercept = psi), linetype = 4) +
                
                # Smoothing Segmented Poisson plot
                geom_point(aes(data_s$Day, data_s$cases_s), col = 2, alpha = 0.4) +
                geom_line(aes(data_s$Day, fitted(seg_fit_s)), col = 2, size = 1) +
                geom_vline(aes(xintercept = psi_s), linetype = 4, col = 2) +
                
                labs(title = paste0("Segmented Poisson Model of ", k),
                     subtitle = paste0("Simple Moving Average (Window size = ", j , ")",
                                       "\nBreakpoint : ", psi, " (original)",
                                       "\nBreakpoint : ", psi_s, " (smoothing)"),
                     x = "Days after the first case ( ~ 2020/08/31)", y = "Daily Cases") +
                theme_bw() +
                theme(plot.title = element_text(size = 15),
                      plot.subtitle = element_text(size = 12))
        
        ggsave(p, filename = paste0("Segmented_Poisson_Model_of_", k, "(", j, ").png"), 
               path = save_path, width = 12, height = 7, dpi = 300)    
        
    
    colnames(coef_Poi) <- c(names(fit$coefficients), paste0(names(fit$coefficients), "_smoothing"))
    colnames(coef_seg_Poi) <- c("(Intercept)", "Day_1", "log(Day)_1",
                                "Day_2", "log(Day)_2", "psi",
                                
                                "(Intercept)_s", "Day_1", "log(Day)_1_s", 
                                "Day_2_s", "log(Day)_2_s", "psi_s")
    }
    
    
    ########################
    #### .csv 파일 저장 ####
    ########################
    
    write.csv(coef_Poi, paste0(result_path, "Coef_Poi_", j, ".csv"), row.names = T)
    write.csv(coef_seg_Poi, paste0(result_path, "Coef_seg_Poi_", j, ".csv"), row.names = T)
}


