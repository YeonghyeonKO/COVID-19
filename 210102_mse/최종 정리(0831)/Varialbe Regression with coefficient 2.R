rm(list=ls())
library(dplyr)
library(ggplot2)

# Option
log_option = TRUE # If True, all parameters log-scaled.

# Path
path <- "C:/Users/김학용/Desktop/코로나 연구 인턴/주차별 작업사항/최종 정리/"
data_path <- paste0(path, "prep/")
result_path <- ifelse(log_option,paste0(path, "regression analysis/regression result(log scaled)/"),
                      paste0(path, "regression analysis/regression result/"))
plot_path <- ifelse(log_option,paste0(path, "regression analysis/regression result(log scaled)/plot/"),
                    paste0(path, "regression analysis/regression result/plot/"))

#### Data Load ####
## segmented growthcurve model parameter
if(log_option){
  df_coef = log(read.csv(paste0(path,"coef_result.csv"),row.names = 1)[,1:12])
}else{
  df_coef = read.csv(paste0(path,"coef_result.csv"),row.names = 1)[,1:12]
}
country = rownames(df_coef)

## KOSTAT data preprocessing - Travel/Trade
covariate = read.csv(paste0(data_path,"kostat_vars.csv"))
Country = data.frame(Country=country)
KOSTAT = inner_join(Country,covariate,by="Country")
df_KOSTAT = KOSTAT[,c(1:19)]

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


#### Linear Regression ####
# Time-independent variable sources are (1)ourworld.csv, (2)KOSTAT.csv
# First, we analized respectively and then we merge the two results.

## List for X and Y variables
Y_list <- c("a1_Logi", "b1_Logi", "c1_Logi", "a1_Gom", "b1_Gom", "c1_Gom",
            "a2_Logi", "b2_Logi", "c2_Logi", "a2_Gom", "b2_Gom", "c2_Gom")

#X_list <- colnames(df_ourworld)[-c(1,2)]
#df = inner_join(cbind(df_coef,Country=rownames(df_coef)),cbind(df_ourworld,Country=rownames(df_ourworld)),by="Country")

X_list <- colnames(df_KOSTAT)[-1]
df = inner_join(cbind(df_coef,Country=rownames(df_coef)),df_KOSTAT,by="Country") 

## SLR
## df_sum: final result
## df_res: result table for each Y
df_sum <- NULL

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
    
      # scatter plot
      p <- ggplot(df_sub)+
        geom_point(aes(x=X, y=Y))+
        labs(title = paste0("SLR model for parmeter ", Y),
             subtitle = paste0(Y, " ~ ", X), x = X, y = Y)+
        theme_bw() +
        theme(
          plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
          plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
          axis.text=element_text(size=14, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 14, hjust = 0.5),
          axis.title=element_text(size=16, colour = "black"))

      ggsave(p, filename = paste0(Y, " to ", X, " SLR Scatter Plot.png"), path = plot_path, dpi = 300, width = 8, height = 8.5)

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


#### Write result ####

# write.csv(df_sum, paste0(result_path, "regression_ourworld.csv"), row.names = F)
# write.csv(df_sum, paste0(result_path, "regression_ourworld(log scaled).csv"), row.names = F)
# write.csv(df_sum, paste0(result_path, "regression_KOSTAT.csv"), row.names = F)
 # write.csv(df_sum, paste0(result_path, "regression_KOSTAT(log scaled).csv"), row.names = F)

