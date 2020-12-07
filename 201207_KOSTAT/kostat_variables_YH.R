library(tidyverse)

getwd()
setwd("/Users/yeonghyeon/Documents/GitHub/COVID-19/201207_KOSTAT")

coefs <- read_csv("covariate(YH,DE,KH).csv")
kor_eng <- read_csv("Kor_Eng.csv", locale=locale('ko',encoding='euc-kr'))

child_vaccination <- read_csv("child_vaccination.csv", locale=locale('ko',encoding='euc-kr')) %>%
  rename(Korean=1, child_vaccination=2)

malnutrition <- read_csv("malnutrition.csv", locale=locale('ko',encoding='euc-kr'))[,-1] %>%
  rename(Korean=1, malnutrition=2)

migration_rate <- read_csv("migration_rate.csv", locale=locale('ko',encoding='euc-kr')) %>%
  rename(Korean=1, migration_rate=2)

prop_urban_pop <- read_csv("prop_urban_pop.csv", locale=locale('ko',encoding='euc-kr')) %>%
  rename(Korean=1, prop_urban_pop=2)


new_coefs <- kor_eng %>% 
  left_join(child_vaccination, by="Korean") %>%
  left_join(malnutrition, by="Korean") %>%
  left_join(migration_rate, by="Korean") %>%
  left_join(prop_urban_pop, by="Korean")

coefs_result <- inner_join(coefs, new_coefs, by="iso3")[,-c(21:23)][,c(21:24, 1:20)]
write_csv(coefs_result, "kostat_vars.csv")
