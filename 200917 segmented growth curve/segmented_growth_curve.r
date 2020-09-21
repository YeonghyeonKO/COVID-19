#################
#### library ####
#################

# library(minpack.lm)
library(segmented)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggrepel)
library(knitr)
library(nls2)


############################
#### data preprocessing ####
############################

path <- "C:/Users/kyh/Desktop/ê³ ì˜?˜„/?¸?„´?‹­/?•´?™¸ ì½”ë¡œ?‚˜/200906 design matrix smoothing segmented poisson/"
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

for (i in which(df$cases < 0)){df$cases[i] = 0} # Case < 0 ì²˜ë¦¬ (0?œ¼ë¡? ??€ì²?)

num <- df %>%
    group_by(country) %>%
    mutate(total = cumsum(cases)) %>%
    filter(total > 0) %>%
    summarise(count = n())


#####################
#### nls fitting ####
#####################

## Logistic  model
eq_log <- Y ~ a/(1 + exp(b + c*X))

fit_log <- nls2(eq_log, start = grid_Brazil, algorithm = "plinear-random", control = nls.control(maxiter = 3000))
# initial_log_B <- coef(fit_log_B)

n = 11; k = "Australia"

for (k in num$country){
    n <- n + 1
    days <- num$count[n]
    
    data_s <- df %>%
        filter(country == k) %>%
        mutate(total = cumsum(cases)) %>%
        filter(total >= 1) %>%
        mutate(Day = 1:days)
    
    if (data_s$total[days] < 1000){
        next
    }
    
    breakpoint <- 130
    
    #### NLS2 Model (for initialization)
    X <- data_s$Day; Y <- data_s$total
    g <- cut(X, c(-Inf, breakpoint, Inf), label = FALSE)
    
    grid_nls <- data.frame(a = c(max(Y) - 50000, max(Y) + 50000), b = c(-100, 100), c = c(-100, 100))
    fit_log_nls2 <- nls2(eq_log, start = grid_nls, algorithm = "plinear-random", control = nls.control(maxiter = 3000))
    
    initial_log <- coef(fit_log_nls2)
    
    # ì´ˆê¸°ê°’ì„ ë°”íƒ•?œ¼ë¡? ? ?•©
    fit_log_nls <- nls(eq_log, start = list(a = initial_log[1],
                                          b = initial_log[2],
                                          c = initial_log[3]))
    
    ggplot() +
        geom_point(aes(X, Y), shape = 1) +
        geom_line(aes(X, coef(fit_log_nls)[1] / (1 + exp(coef(fit_log_nls)[2] + coef(fit_log_nls)[3] * X))),
                  color = "coral", lwd = 1) +
        labs(title = paste0("COVID-19 Cumulative cases in ", k, " : Logistic Model"),
             x = paste0("From ", data_s$date[1], " To 08/31 in 2020"), y = "Cumulative Cases")
}



X[breakpoint:days]


#Use Segmented to get estimates for parameters with 2 breakpoints
breakpoint <- 130

g <- cut(X, c(-Inf, breakpoint, Inf), label = FALSE)

fit_log_nls2_1 <- nls2(Y[1:breakpoint] ~ a/(1 + exp(b + c * X[1:breakpoint])), 
                     start = grid_nls, 
                     algorithm = "plinear-random", 
                     control = nls.control(maxiter = 1000))

fit_log_nls2_2 <- nls2(Y[(breakpoint+1) : days] - Y[breakpoint] ~ a/(1 + exp(b + c * X[1:(days-breakpoint)])), 
                       start = grid_nls, 
                       algorithm = "plinear-random", 
                       control = nls.control(maxiter = 1000))

initial_log <- cbind(coef(fit_log_nls2_1), coef(fit_log_nls2_2))
initial_log[1, 2] <- initial_log[1,2] + Y[breakpoint]

initial_log_1 <- coef(fit_log_nls2_1)

ggplot() +
    geom_point(aes(X, Y), alpha = 0.3) +
    geom_point(aes(X[(breakpoint+1):days], Y[(breakpoint+1) : days] - Y[breakpoint]), col = 2) + 
    geom_line(aes(X[1:breakpoint], fitted(fit_log_nls2_1))) +
    geom_line(aes(X[(breakpoint+1):days], fitted(fit_log_nls2_2)), col = 2)


for (i in 1:1000){
    print(i)
    tryCatch(
        expr = {
            set.seed(i)
            fm <- nls(Y ~ a/(1 + exp(b + c*X)), 
                      start = list(a = initial_log_1, 
                                   b = initial_log_1, 
                                   c = initial_log_1))
            },
            
            # fm <- nls(Y ~ a[g]/(1 + exp(b[g] + c[g]*X)), start = list(a = initial_log[1,],
            #                                                       b = initial_log[2,],
            #                                                       c = initial_log[3,]))
        
        error = function(e) {
            print(e)
            break
            },
        finally = next
    )
}

?nls


plot(Y ~ X, df)
lines(fitted(fm) ~ X, df, col = "red")
