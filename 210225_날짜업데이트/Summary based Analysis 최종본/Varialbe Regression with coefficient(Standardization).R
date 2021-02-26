#rm(list=ls())
library(dplyr)
library(ggplot2)
library(forcats)

# Path
getwd()
path <- "/Users/yeonghyeonko/Documents/GitHub/COVID-19/210225_날짜업데이트/Summary based Analysis 최종본/"
data_path <- paste0(path, "data/")
result_path <- paste0(path, "result/")
plot_path <- paste0(path, "plot/SLR/")

#### Data Load ####
# Import coefficient data
df_logistic = read.csv(paste0(result_path,"segmented_Logistic_result.csv"))
df_logistic = df_logistic[which(df_logistic$MSSE<=0.4),1:7]

df_gompertz = read.csv(paste0(result_path,"segmented_Gompertz_result.csv"))
df_gompertz = df_gompertz[which(df_gompertz$MSSE<=0.4),1:7]

df_coef = full_join(df_logistic,df_gompertz,by="X")
rownames(df_coef) = df_coef$X
df_coef = df_coef[,-1]

country = rownames(df_coef)


## KOSIS data preprocessing - Travel/Trade
KOSIS = read.csv(paste0(data_path,"KOSIS.csv"))
Country = data.frame(Country=country)
KOSIS = inner_join(Country,KOSIS,by="Country")
df_KOSIS = KOSIS[,c(1:17)]
for(i in 2:17){
  df_KOSIS[,i] = as.numeric(df_KOSIS[,i])
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

Covariate = full_join(df_KOSIS,df_ourworld,by="Country")
colnames(Covariate)

# Distribution of Variables
for(i in 2:length(Covariate)){
  hist(Covariate[,i], main=colnames(Covariate)[i])
}

# Standardization
Covariate_standardized = Covariate
for(i in 2:length(Covariate)){
  Covariate_standardized[,i] = (Covariate[,i] - mean(Covariate[,i],na.rm = TRUE))/sd(Covariate[,i],na.rm = TRUE)
}
for(i in 2:length(Covariate)){
  hist(Covariate_standardized[,i])
}


#### Linear Regression ####
# Time-independent variable sources are (1)ourworld.csv, (2)KOSIS.csv
# First, we analized respectively and then we merge the two results.

## List for X and Y variables
Y_list1 <- c("a1_Logi", "b1_Logi", "c1_Logi", "a1_Gom", "b1_Gom", "c1_Gom")
Y_list2 <- c("a2_Logi", "b2_Logi", "c2_Logi", "a2_Gom", "b2_Gom", "c2_Gom")

X_list <- colnames(Covariate_standardized)[-1]
df = inner_join(cbind(df_coef,Country=rownames(df_coef)), Covariate_standardized, by="Country")

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
    names(df_sum)[ncol(df_sum)] <- Y
  } 
  
  if(loop_count==1){
    df_sum1 = df_sum
  }else{
    df_sum2 = df_sum
  } 
}



#### Write result #### 
write.csv(df_sum1, paste0(result_path, "Linear Regression(1st Segment, standardized).csv"), row.names = F)
write.csv(df_sum2, paste0(result_path, "Linear Regression(2nd Segment, standardized).csv"), row.names = F)



#### Bar Plot for standardized SLR coefficients between Growth Curve Parameter and Time-independent variable ####
## Logistic Model
for (i in 4:6){
  df_Logi <- data.frame(
    Explanatory = factor(unique(df_sum1$Explanatory)),
    coefficient = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==1)], df_sum2[,i][which((1:nrow(df_sum2))%%3==1)]),
    p.value = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==0)], df_sum2[,i][which((1:nrow(df_sum2))%%3==0)]),
    segment = factor(rep(1:2, each=length(unique(df_sum1$Explanatory)))),
    category = factor(c(7,7,1,1,7,5,4,4,7,3,3,6,5,5,2,1,1,1,2,2,2,1,7,7,7,7,7,7,7,2,5))
    # Territory&population : 1 / Age : 2 / Trade : 3 / Environment : 4 / 
    # Education&Culture&Science : 5 / National Accounts : 6 / 
    # Health&Sociecy&Welfare : 7
  ) %>% arrange(category)
  
  
  # volcano plot
  ggplot(df_Logi, aes(x=coefficient, y=p.value,
                        group=Explanatory, color=category)) +
    geom_point(aes(shape=segment)) +
    geom_line() +
    labs(x="coefficient", y="p-value", 
         title=colnames(df_sum1)[i], subtitle="volcano plot") +
    theme_classic() +
    theme(plot.title=element_text(size=15, hjust=0.5, face="bold", vjust=2))

  # bar plot
  ggplot(df_Logi %>%
           mutate(Explanatory = fct_reorder(Explanatory, desc(category))),
         aes(x=Explanatory, y=p.value, fill=segment, label=category)) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(position=position_dodge(width = 0), aes(x=Explanatory, y=1)) +
    geom_hline(yintercept=0.05, color = 2, size=1, linetype=8) +
    labs(x="", title=paste0("SLR for ", colnames(df_sum1)[i])) +
    coord_flip() +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    # scale_fill_hue(l=70, c=100) +
    theme_classic() +
    theme(plot.title=element_text(size=15, hjust=0.5, face="bold", vjust=2))
}



## Gompertz Model
for (i in 7:9){
  df_Gomp <- data.frame(
    Explanatory = factor(unique(df_sum1$Explanatory)),
    coefficient = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==1)], df_sum2[,i][which((1:nrow(df_sum2))%%3==1)]),
    p.value = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==0)], df_sum2[,i][which((1:nrow(df_sum2))%%3==0)]),
    segment = factor(rep(1:2, each=length(unique(df_sum1$Explanatory)))),
    category = factor(c(7,7,1,1,7,5,4,4,7,3,3,6,5,5,2,1,1,1,2,2,2,1,7,7,7,7,7,7,7,2,5))
    # Territory&population : 1 / Age : 2 / Trade : 3 / Environment : 4 / 
    # Education&Culture&Science : 5 / National Accounts : 6 / 
    # Health&Sociecy&Welfare : 7
  ) %>% arrange(category)
  
  
  # volcano plot
  ggplot(df_Gomp, aes(x=coefficient, y=p.value,
                      group=Explanatory, color=category)) +
    geom_point(aes(shape=segment)) +
    geom_line() +
    labs(x="coefficient", y="p-value", 
         title=colnames(df_sum1)[i], subtitle="volcano plot") +
    theme_classic() +
    theme(plot.title=element_text(size=15, hjust=0.5, face="bold", vjust=2))
  
  # bar plot
  ggplot(df_Logi %>%
           mutate(Explanatory = fct_reorder(Explanatory, desc(category))),
         aes(x=Explanatory, y=p.value, fill=segment, label=category)) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(position=position_dodge(width = 0), aes(x=Explanatory, y=1)) +
    geom_hline(yintercept=0.05, color = 2, size=1, linetype=8) +
    labs(x="", title=paste0("SLR for ", colnames(df_sum1)[i])) +
    coord_flip() +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    # scale_fill_hue(l=70, c=100) +
    theme_classic() +
    theme(plot.title=element_text(size=15, hjust=0.5, face="bold", vjust=2))
}
