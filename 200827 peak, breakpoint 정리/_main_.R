#rm(list=ls())
library(dplyr)
library(ggplot2)
library(nls2)

# Analysis Period( ~ 8/25)
max_date <- set_date("2020/8/26")

#### Data import ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data
df_sum = preprocessing_data()


#### Analysis ####
main_country <- c("South_Korea","Japan","United_Kingdom","France",
                  "Germany","Italy","Iran","China","United_States_of_America",
                  "Canada","Spain","Brazil","Russia","India")

df_result<-derivative_analysis(Country = country,gkf_bandwidth = 14,save_image = TRUE,save_excel = TRUE)

#### Segmented Logistic regression
cumulative_analysis(Country = main_country)


