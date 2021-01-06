#rm(list=ls())
library(dplyr)
library(ggplot2)

# Option
log_option = TRUE # If True, all parameters log-scaled.

# Path
path <- "C:/Users/김학용/Desktop/코로나 연구 인턴/주차별 작업사항/최종 정리(8.31)/"
# path <- "/Users/yeonghyeon/Documents/GitHub/COVID-19/210104_MSE_MSSE/Summary based Analysis 최종본/"
data_path <- paste0(path, "data/")
result_path <- paste0(path, "result/")
plot_path <- paste0(path, "plot/SLR/")

#### Data Load ####
# Import coefficient data
df_logistic = read.csv(paste0(result_path,"segmented_Logistic_result.csv"))
df_logistic = df_logistic[which(coef_seg_Logi$MSSE<=0.4),1:7]

df_gompertz = read.csv(paste0(result_path,"segmented_Gompertz_result.csv"))
df_gompertz = df_gompertz[which(coef_seg_Gom$MSSE<=0.4),1:7]

df_coef = full_join(df_logistic,df_gompertz,by="X")
rownames(df_coef) = df_coef$X
df_coef = df_coef[,-1]

country = rownames(df_coef)


## KOSTAT data preprocessing - Travel/Trade
KOSTAT = read.csv(paste0(data_path,"kostat_vars.csv"))
Country = data.frame(Country=country)
KOSTAT = inner_join(Country,KOSTAT,by="Country")
df_KOSTAT = KOSTAT[,c(1:19)]
for(i in 12:16){
  df_KOSTAT[,i] = as.numeric(df_KOSTAT[,i])
}


## ourworld data preprocessing
ourworld = read.csv(paste0(data_path,"OurWorld_20201101.csv"))
country_ourworld = as.character(unique(ourworld$location))
df_ourworld = NULL
for(j in country_ourworld){
  inx = min(which(ourworld$location==j))
  temp_row = ourworld[inx,c(1:2,35:49)]   
  df_ourworld = rbind(df_ourworld,temp_row)  
}

# country_name matching - only "Eswatini" is mismatched
country_ourworld = gsub(" ","_",country_ourworld)
index_country = which(country_ourworld %in% c("Cote_d'Ivoire","Czech_Republic",
                                              "Democratic_Republic_of_Congo", "Faeroe_Islands", 
                                              "Guinea-Bissau", "Macedonia", "Tanzania", "United_States")) 
country_ourworld[index_country] = c("Cote_dIvoire", "Czechia", "Democratic_Republic_of_the_Congo",                        
                                    "Faroe_Islands","Guinea_Bissau","North_Macedonia",
                                    "United_Republic_of_Tanzania", "United_States_of_America")
rownames(df_ourworld) = country_ourworld
df_ourworld$Country = country_ourworld
df_ourworld = df_ourworld[,3:18]
# intersect(country,country_ourworld) 
# setdiff(country,country_ourworld)

Covariate = full_join(df_KOSTAT,df_ourworld,by="Country")
colnames(Covariate)

#### Linear Regression ####
# Time-independent variable sources are (1)ourworld.csv, (2)KOSTAT.csv
# First, we analized respectively and then we merge the two results.

## List for X and Y variables
Y_list1 <- c("a1_Logi", "b1_Logi", "c1_Logi", "a1_Gom", "b1_Gom", "c1_Gom")
Y_list2 <- c("a2_Logi", "b2_Logi", "c2_Logi", "a2_Gom", "b2_Gom", "c2_Gom")

X_list <- colnames(Covariate)[-1]
df = inner_join(cbind(df_coef,Country=rownames(df_coef)), Covariate, by="Country")

## SLR
## df_sum: final result
## df_res: result table for each Y
df_sum1 <- NULL
df_sum2 <- NULL

loop_count = 0 
for(Y_list in list(Y_list1,Y_list2)){
  df_sum <- NULL
  loop_count = loop_count + 1    
  for(Y in Y_list){
    # initialize
    {
      X_box <- c()
      Y_box <- c()
      Res_box <-c()
      Res_type_box <- c()
      n_box <- c()
    }
    
    for(X in X_list){
      # Make dataframe with one X variable, one Y variable
      df_sub <- data.frame(X = df[[X]], Y = df[[Y]])
      fitted <- lm(Y~X, df_sub)
      res <- summary(fitted)
      
      # 3 result values(coefficient, R square, P value) will be added to "box" variable
      Res_box <- c(Res_box,res$coefficients[2], res$r.squared, anova(fitted)$'Pr(>F)'[1])
      Res_type_box <- c(Res_type_box, c("Coefficient","R-square", "P-value"))
      X_box <- c(X_box, rep(X,3))
      Y_box <- c(Y_box, rep(Y,3))
      n_box <- c(n_box, rep(length(which(!is.na(df_sub$X) & !is.na(df_sub$Y))), 3))
      
      # # scatter plot
      # p <- ggplot(df_sub)+
      #   geom_point(aes(x=X, y=Y))+
      #   labs(title = paste0("SLR model for parmeter ", Y),
      #        subtitle = paste0(Y, " ~ ", X), x = X, y = Y)+
      #   theme_bw() +
      #   theme(
      #     plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
      #     plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
      #     axis.text=element_text(size=14, face = "bold", colour = "black"),
      #     axis.text.x = element_text(size = 14, hjust = 0.5),
      #     axis.title=element_text(size=16, colour = "black"))
      # 
      # ggsave(p, filename = paste0(Y, " to ", X, " SLR Scatter Plot.png"), path = plot_path, dpi = 300, width = 8, height = 8.5)
      
    }
    
    df_res <- data.frame("#Samples" = n_box, Explanatory = X_box, Res_Type = Res_type_box, Res = Res_box)
    
    # Adding result table for each Y (df_res) to final result table(df_sum)
    if(is.null(df_sum)){
      df_sum <- df_res
    }else{
      df_sum <- cbind(df_sum, df_res$Res)
    }
    
R version 3.5.3 (2019-03-11) -- "Great Truth"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from C:/Users/源?숈슜/Desktop/肄붾줈???곌뎄 ?명꽩/二쇱감蹂??묒뾽?ы빆/理쒖쥌 ?뺣━(8.31)/.RData]

> View(df_sub)
> print(j)
[1] "International"
> for(j in list(1,2)){
+   print(j)
+ }
[1] 1
[1] 2
> for(j in list(Y_list1,Y_list2)){
+   print(j)
+ }
Error in list(Y_list1, Y_list2) : object 'Y_list1' not found
> ## List for X and Y variables
> Y_list1 <- c("a1_Logi", "b1_Logi", "c1_Logi", "a1_Gom", "b1_Gom", "c1_Gom")
> Y_list2 <- c("a2_Logi", "b2_Logi", "c2_Logi", "a2_Gom", "b2_Gom", "c2_Gom")
> for(j in list(Y_list1,Y_list2)){
+   print(j)
+ }
[1] "a1_Logi" "b1_Logi" "c1_Logi" "a1_Gom"  "b1_Gom"  "c1_Gom" 
[1] "a2_Logi" "b2_Logi" "c2_Logi" "a2_Gom"  "b2_Gom"  "c2_Gom" 
> print(k)
Error in print(k) : object 'k' not found
> for(j in list(Y_list1,Y_list2)){
+   for(k in j){
+     print(k)
+   }
+ }
[1] "a1_Logi"
[1] "b1_Logi"
[1] "c1_Logi"
[1] "a1_Gom"
[1] "b1_Gom"
[1] "c1_Gom"
[1] "a2_Logi"
[1] "b2_Logi"
[1] "c2_Logi"
[1] "a2_Gom"
[1] "b2_Gom"
[1] "c2_Gom"
> a==Y_list1
       awe
[1,] FALSE
[2,] FALSE
[3,] FALSE
[4,] FALSE
> ## SLR
> ## df_sum: final result
> ## df_res: result table for each Y
> df_sum1 <- NULL
> df_sum2 <- NULL
> loop_count = 0
> ## SLR
> ## df_sum: final result
> ## df_res: result table for each Y
> df_sum1 <- NULL
> df_sum2 <- NULL
> loop_count = 0
> for(Y_list in list(Y_list1,Y_list2)){
+   df_sum <- NULL
+   loop_count = loop_count + 1    
+   for(Y in Y_list){
+     # initialize
+     {
+       X_box <- c()
+       Y_box <- c()
+       Res_box <-c()
+       Res_type_box <- c()
+       n_box <- c()
+     }
+     
+     for(X in X_list){
+       # Make dataframe with one X variable, one Y variable
+       df_sub <- data.frame(X = df[[X]], Y = df[[Y]])
+       fitted <- lm(Y~X, df_sub)
+       res <- summary(fitted)
+       
+       # 3 result values(coefficient, R square, P value) will be added to "box" variable
+       Res_box <- c(Res_box,res$coefficients[2], res$r.squared, anova(fitted)$'Pr(>F)'[1])
+       Res_type_box <- c(Res_type_box, c("Coefficient","R-square", "P-value"))
+       X_box <- c(X_box, rep(X,3))
+       Y_box <- c(Y_box, rep(Y,3))
+       n_box <- c(n_box, rep(length(which(!is.na(df_sub$X) & !is.na(df_sub$Y))), 3))
+       
+       # # scatter plot
+       # p <- ggplot(df_sub)+
+       #   geom_point(aes(x=X, y=Y))+
+       #   labs(title = paste0("SLR model for parmeter ", Y),
+       #        subtitle = paste0(Y, " ~ ", X), x = X, y = Y)+
+       #   theme_bw() +
+       #   theme(
+       #     plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
+       #     plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
+       #     axis.text=element_text(size=14, face = "bold", colour = "black"),
+       #     axis.text.x = element_text(size = 14, hjust = 0.5),
+       #     axis.title=element_text(size=16, colour = "black"))
+       # 
+       # ggsave(p, filename = paste0(Y, " to ", X, " SLR Scatter Plot.png"), path = plot_path, dpi = 300, width = 8, height = 8.5)
+       
+     }
+     
+     df_res <- data.frame("#Samples" = n_box, Explanatory = X_box, Res_Type = Res_type_box, Res = Res_box)
+     
+     # Adding result table for each Y (df_res) to final result table(df_sum)
+     if(is.null(df_sum)){
+       df_sum <- df_res
+     }else{
+       df_sum <- cbind(df_sum, df_res$Res)
+     }
+     
+     names(df_sum)[ncol(df_sum)] <- Y
+   }
+   
    names(df_sum)[ncol(df_sum)] <- Y
  }
  
  if(loop_count==1){
    df_sum1 = df_sum
  }else{
    df_sum2 = df_sum
  } 
}



#### Write result #### 
write.csv(df_sum1, paste0(result_path, "Linear Regressioin(1st Segment).csv"), row.names = F)
write.csv(df_sum2, paste0(result_path, "Linear Regressioin(2nd Segment).csv"), row.names = F)






