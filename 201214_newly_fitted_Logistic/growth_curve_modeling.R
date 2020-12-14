#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)
library(gridExtra)

setwd("~/Documents/GitHub/COVID-19/201214_fitting_Korea")
getwd()

# Analysis Period( ~ 8/31)
max_date <- set_date("2020/8/31") 
# Analysis Period( ~ 9/15)
max_date <- set_date("2020/9/15") 
# Analysis Period( ~ 9/20)
max_date <- set_date("2020/9/20") 
# Analysis Period( ~ 9/25)
max_date <- set_date("2020/9/25") 
# Analysis Period( ~ 9/30)
max_date <- set_date("2020/9/30") 
# Analysis Period( ~ 10/15)
max_date <- set_date("2020/10/15")
# Analysis Period( ~ 10/31)
max_date <- set_date("2020/10/31")

#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]


# preprocess data (until 20/08/31)
df_sum = preprocessing_data()[1:245,]
# preprocess data (until 20/9/15)
df_sum = preprocessing_data()[1:260,]
# preprocess data (until 20/9/20)
df_sum = preprocessing_data()[1:265,]
# preprocess data (until 20/9/25)
df_sum = preprocessing_data()[1:270,]
# preprocess data (until 20/9/30)
df_sum = preprocessing_data()[1:275,]
# preprocess data (until 20/10/15)
df_sum = preprocessing_data()[1:290,]
# preprocess data (until 20/10/31)
df_sum = preprocessing_data()[1:306,]

par(mfrow=c(1,1))
ts.plot(df_sum["South_Korea"])
ts.plot(df_sum["Taiwan"])
# Simple Moving average (window size = 7)
for(i in country){
  df_sum[i][,1] <- pracma::movavg(df_sum[i][,1], 7, type = "s")  
}

# Peak Detectiong Analysis
df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_break = 50,save_image = TRUE,save_excel = TRUE)

derivative_analysis(Country = "South_Sudan",gkf_bandwidth = 14,first_break = 50,save_image = FALSE,save_excel = FALSE)


###############################################################################

#### segmented logistic ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #

#### separately fit logistic models #### 
coef_sep_Logi <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Logi) <- country
colnames(coef_sep_Logi) <- c("a1_Logi", "b1_Logi", "c1_Logi", "a2_Logi", "b2_Logi", "c2_Logi",
                            "startpoint", "breakpoint")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = break_sum$break1[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_sep_Logi[i, 1:3] <- segLogistic_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                  max_iter = 1000)[1,]
      coef_sep_Logi[i, 4:6] <- segLogistic_separate(Country=i,break_point = break_point,save_image = TRUE,prediction_plot = FALSE,
                                                    max_iter = 1000)[2,]  
      coef_sep_Logi[i, 7:8] <- break_sum[i, 1:2]
      
    },
    error = function(e){
    next
  },
  finally = next
  )
}


A <- coef_sep_Logi

write.csv(coef_sep_Logi, "200925_Separated_Logistic_coefficients.csv")


coef_200831 <- read_csv("200831_Separated_Logistic_coefficients.csv")
na.omit(coef_200831)

coef_200925 <- read_csv("200925_Separated_Logistic_coefficients.csv")
na.omit(coef_200925)


###############################################################################

#### segmented Bertalanffy ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #


#### separately fit Bertalanffy models #### 

coef_sep_Ber <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Ber) <- country
colnames(coef_sep_Ber) <- c("a1_Ber", "b1_Ber", "c1_Ber", "a2_Ber", "b2_Ber", "c2_Ber",
                            "startpoint", "breakpoint")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = break_sum$break1[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      
      coef_sep_Ber[i, 1:3] <- segBertalanffy_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                                    max_iter = 1000)[1,]
      coef_sep_Ber[i, 4:6] <- segBertalanffy_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = FALSE,
                                                    max_iter = 1000)[2,]  
      coef_sep_Ber[i, 7:8] <- break_sum[i, 1:2]
    },
    error = function(e){
      next
    },
    finally = next
  )
}

B <- coef_sep_Ber

write.csv(coef_sep_Ber, "Separated_Bertalanffy_coefficients.csv")

###############################################################################

#### segmented Gompertz ####
# - 1st breakpoint is not real breakpoint but starting point of logistic regression            #
# - from 2nd point, real breakpoint is counted                                                 #
# - ex) breakpoint=c(50, 70, 120) => 50 is the first day since cumulative # is bigger than 50  # 
# -     so the # of real breakpoint is 2(70, 120)                                              #

#### separately fit Gompertz models #### 

coef_sep_Gom <- data.frame(matrix(NA, ncol = 8, nrow = length(country)))
rownames(coef_sep_Gom) <- country
colnames(coef_sep_Gom) <- c("a1_Gom", "b1_Gom", "c1_Gom", "a2_Gom", "b2_Gom", "c2_Gom",
                            "startpoint", "breakpoint")

for(i in country){
  tryCatch(
    expr={
      temp = which(rownames(break_sum)==i)
      break_point = break_sum$break1[temp]
      for(j in 2:5){
        if(!is.na(break_sum[temp,j])){
          break_point = c(break_point,break_sum[temp,j])
        }
      }
      coef_sep_Gom[i, 1:3] <- segGompertz_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot =FALSE,
                           max_iter = 1000)[1,]
      coef_sep_Gom[i, 4:6] <- segGompertz_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot =FALSE,
                           max_iter = 1000)[2,]
      coef_sep_Gom[i, 7:8] <- break_sum[i, 1:2]
    },
    error = function(e){
      next
    },
    finally = next
  )
}

C <- coef_sep_Gom

write.csv(coef_sep_Gom, "Separated_Gompertz_coefficients.csv")


#### Comparison of the coefficients of Separated Growth Curve Model ####
## First Segmentation

# All 3 models
length(country[!is.na(coef_sep_Logi[1]) & !is.na(coef_sep_Ber[1]) & !is.na(coef_sep_Gom[1])])
fitted_country_Logi_Ber_Gom <- country[!is.na(coef_sep_Logi[1]) & !is.na(coef_sep_Ber[1]) & !is.na(coef_sep_Gom[1])]

common_coef_all_1 <- cbind(fitted_country_Logi_Ber_Gom, coef_sep_Logi[fitted_country_Logi_Ber_Gom, 7:8]) %>%
  cbind(coef_sep_Logi[fitted_country_Logi_Ber_Gom, 1:3]) %>%
  cbind(coef_sep_Ber[fitted_country_Logi_Ber_Gom, 1:3]) %>%
  cbind(coef_sep_Gom[fitted_country_Logi_Ber_Gom, 1:3]) %>%
  rename(country = fitted_country_Logi_Ber_Gom)

rownames(common_coef_all_1) <- 1:nrow(common_coef_all_1)
write.csv(common_coef_all_1, "common_coef_all_1.csv")


# Logistic & Gompertz models
length(country[!is.na(coef_sep_Logi[1]) & !is.na(coef_sep_Gom[1])])
fitted_country_Logi_Gom <- country[!is.na(coef_sep_Logi[1]) & !is.na(coef_sep_Gom[1])]

common_coef_1 <- cbind(fitted_country_Logi_Gom, coef_sep_Logi[fitted_country_Logi_Gom, 7:8]) %>%
  cbind(coef_sep_Logi[fitted_country_Logi_Gom, 1:3]) %>%
  cbind(coef_sep_Gom[fitted_country_Logi_Gom, 1:3]) %>%
  rename(country = fitted_country_Logi_Gom)

rownames(common_coef_1) <- 1:nrow(common_coef_1)
write.csv(common_coef_1, "common_coef_1.csv")



## Second Segmentation

# All 3 models
length(country[!is.na(coef_sep_Logi[4]) & !is.na(coef_sep_Ber[4]) & !is.na(coef_sep_Gom[4])])
fitted_country_Logi_Ber_Gom <- country[!is.na(coef_sep_Logi[4]) & !is.na(coef_sep_Ber[4]) & !is.na(coef_sep_Gom[4])]

common_coef_all_2 <- cbind(fitted_country_Logi_Ber_Gom, coef_sep_Logi[fitted_country_Logi_Ber_Gom, 7:8]) %>%
  cbind(coef_sep_Logi[fitted_country_Logi_Ber_Gom, 4:6]) %>%
  cbind(coef_sep_Ber[fitted_country_Logi_Ber_Gom, 4:6]) %>%
  cbind(coef_sep_Gom[fitted_country_Logi_Ber_Gom, 4:6]) %>%
  rename(country = fitted_country_Logi_Ber_Gom)

# rownames(common_coef_all_2) <- 1:nrow(common_coef_all_2)
write.csv(common_coef_all_2, "common_coef_all_2.csv")


# Logistic & Gompertz models
length(country[!is.na(coef_sep_Logi[4]) & !is.na(coef_sep_Gom[4])])
fitted_country_Logi_Gom <- country[!is.na(coef_sep_Logi[4]) & !is.na(coef_sep_Gom[4])]

common_coef_2 <- cbind(fitted_country_Logi_Gom, coef_sep_Logi[fitted_country_Logi_Gom, 7:8]) %>%
  cbind(coef_sep_Logi[fitted_country_Logi_Gom, 4:6]) %>%
  cbind(coef_sep_Gom[fitted_country_Logi_Gom, 4:6]) %>%
  rename(country = fitted_country_Logi_Gom)

rownames(common_coef_2) <- 1:nrow(common_coef_2)
write.csv(common_coef_2, "common_coef_2.csv")



#### Analysis of common_coef_1 & common_coef_2 ####
par( mfrow = c(2,1))

## raw data
summary(common_coef_1[-1:-3])
hist(common_coef_1[,4], breaks=50, main="a1_Logi"); hist(common_coef_1[,7], breaks=50, main="a1_Gom"); 
hist(common_coef_1[,5], breaks=50, main="b1_Logi"); hist(common_coef_1[,8], breaks=50, main="b1_Gom");
hist(common_coef_1[,6], breaks=50, main="c1_Logi"); hist(common_coef_1[,9], breaks=50, main="c1_Gom");

summary(common_coef_2[-1:-3])
hist(common_coef_2[,4], breaks=50, main="a2_Logi"); hist(common_coef_2[,7], breaks=50, main="a2_Gom"); 
hist(common_coef_2[,5], breaks=50, main="b2_Logi"); hist(common_coef_2[,8], breaks=50, main="b2_Gom");
hist(common_coef_2[,6], breaks=50, main="c2_Logi"); hist(common_coef_2[,9], breaks=50, main="c2_Gom");


## standardized data
scaled_common_coef_1 <- data.frame(scale(common_coef_1[-1:-3]))
summary(scaled_common_coef_1)
hist(scaled_common_coef_1[,1], breaks=50, main="a1_Logi"); hist(scaled_common_coef_1[,4], breaks=50, main="a1_Gom"); 
hist(scaled_common_coef_1[,2], breaks=50, main="b1_Logi"); hist(scaled_common_coef_1[,5], breaks=50, main="b1_Gom"); 
# except Kyrgyzstan (outlier)
hist(scale(common_coef_1[-68,2]), breaks=50, main="b1_Logi"); hist(scale(common_coef_1[-68,5]), breaks=50, main="b1_Gom");
hist(scale(common_coef_1[,3]), breaks=50, main="c1_Logi"); hist(scale(common_coef_1[,6]), breaks=50, main="c1_Gom");
par( mfrow = c(1,1))

scaled_common_coef_2 <- data.frame(scale(common_coef_2[-1:-3]))
summary(scaled_common_coef_2)
hist(scaled_common_coef_2[,1], breaks=50, main="a2_Logi"); hist(scaled_common_coef_2[,4], breaks=50, main="a2_Gom"); 
hist(scaled_common_coef_2[,2], breaks=50, main="b2_Logi"); hist(scaled_common_coef_2[,5], breaks=50, main="b2_Gom"); 
# except Kyrgyzstan (outlier)
hist(scale(common_coef_2[-68,2]), breaks=50, main="b2_Logi"); hist(scale(common_coef_2[-68,5]), breaks=50, main="b2_Gom");
hist(scale(common_coef_2[,3]), breaks=50, main="c2_Logi"); hist(scale(common_coef_2[,6]), breaks=50, main="c2_Gom");
par( mfrow = c(1,1))


## Plotting
a1 <- ggplot(scaled_common_coef_1) + 
  geom_point(aes(common_coef_1$country, a1_Logi), col=2) +
  geom_point(aes(common_coef_1$country, a1_Gom), col=4) +
  labs(title="Comparison between Logistic vs Gompertz", x="", y="a1 coefficient") +
  annotate("text", x="Spain", y=9, label="• Logistic", col = 2, size=5) +
  annotate("text", x="Spain", y=7, label="• Gompertz", col = 4, size=5) +
  theme_bw()

b1 <- ggplot(scaled_common_coef_1) + 
  geom_point(aes(common_coef_1$country, b1_Logi), col=2) +
  geom_point(aes(common_coef_1$country, b1_Gom), col=4) +
  labs(x="", y="b1 coefficient") +
  theme_bw()

c1 <- ggplot(scaled_common_coef_1) + 
  geom_point(aes(common_coef_1$country, c1_Logi), col=2) +
  geom_point(aes(common_coef_1$country, c1_Gom), col=4) +
  labs(x="countries", y="c1 coefficient") +
  theme_bw()

p1 <- grid.arrange(a1, b1, c1, nrow = 3)
ggsave(p1, filename = "Logistic vs Gompertz coef. in the 1st segmentation.png", 
       width = 11, height = 7, dpi = 300)

a2 <- ggplot(scaled_common_coef_2) + 
  geom_point(aes(common_coef_2$country, a2_Logi), col=2) +
  geom_point(aes(common_coef_2$country, a2_Gom), col=4) +
  labs(title="Comparison between Logistic vs Gompertz", x="", y="a2 coefficient") +
  annotate("text", x="Peru", y=6, label="• Logistic", col = 2, size=5) +
  annotate("text", x="Peru", y=4.5, label="• Gompertz", col = 4, size=5) +
  theme_bw()

b2 <- ggplot(scaled_common_coef_2) + 
  geom_point(aes(common_coef_2$country, b2_Logi), col=2) +
  geom_point(aes(common_coef_2$country, b2_Gom), col=4) +
  labs(x="", y="b2 coefficient") +
  theme_bw()

c2 <- ggplot(scaled_common_coef_2) + 
  geom_point(aes(common_coef_2$country, c2_Logi), col=2) +
  geom_point(aes(common_coef_2$country, c2_Gom), col=4) +
  labs(x="countries", y="c2 coefficient") +
  theme_bw()

p2 <-grid.arrange(a2, b2, c2, nrow = 3)
ggsave(p2, filename = "Logistic vs Gompertz coef. in the 2nd segmentation.png", 
       width = 11, height = 7, dpi = 300)
