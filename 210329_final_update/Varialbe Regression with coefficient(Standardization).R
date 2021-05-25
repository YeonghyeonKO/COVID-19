#rm(list=ls())
library(dplyr)
library(ggplot2)
library(forcats)
library(ggtext)

# Path
getwd()
path <- "/Users/yeonghyeonko/Documents/GitHub/COVID-19/210329_final_update/"
data_path <- paste0(path, "data/")
result_path <- paste0(path, "result/")
plot_path <- paste0(path, "plot/SLR/")

#### Data Load ####
# Import coefficient data
df_logistic = read.csv(paste0(result_path,"segmented_Logistic_result.csv"))
df_logistic = df_logistic[which(df_logistic$MSSE<=0.4),1:10]

df_gompertz = read.csv(paste0(result_path,"segmented_Gompertz_result.csv"))
df_gompertz = df_gompertz[which(df_gompertz$MSSE<=0.4),1:10]

df_coef = full_join(df_logistic,df_gompertz,by="X")
rownames(df_coef) = df_coef$X
df_coef = df_coef[,-1]

country = rownames(df_coef)


## KOSIS data preprocessing - Travel/Trade
KOSIS = read.csv(paste0(data_path,"KOSIS.csv"))[,-20] %>%
  mutate(Country=iso3)
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


rownames(df_ourworld) = df_ourworld$iso_code
df_ourworld$Country = df_ourworld$iso_code
df_ourworld = df_ourworld[-c(214:215),3:18]
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
Y_list3 <- c("a3_Logi", "b3_Logi", "c3_Logi", "a3_Gom", "b3_Gom", "c3_Gom")

X_list <- colnames(Covariate_standardized)[-1]
df = inner_join(cbind(df_coef,Country=rownames(df_coef)), Covariate_standardized, by="Country")

## SLR
## df_sum: final result
## df_res: result table for each Y
df_sum1 <- NULL
df_sum2 <- NULL
df_sum3 <- NULL

loop_count = 0  
for(Y_list in list(Y_list1,Y_list2,Y_list3)){
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
  }
  if(loop_count==2){
    df_sum2 = df_sum
  }
  if(loop_count==3){
    df_sum3 = df_sum
  }
}

df_sum1_ <- df_sum1[-seq(2, 93, 3),]
a <- data.frame(matrix(NA, nrow=31, ncol=6))
for (i in 1:31) {
  for (j in 4:9){
    a[i,j-3] <- paste0(substr(df_sum1_[2*i - 1, j], 1,13), " (", substr(df_sum1_[2*i,j], 1, 6), ")")
  }
}
a <- cbind(size=df_sum1$X.Samples[seq(2, 93, 3)], factor=df_sum1$Explanatory[seq(2, 93, 3)], a)
a <- a %>% select(1,2,3,6,4,7,5,8)
colnames(a) <- c("size", "factor", "a1_logi", "a1_gom","b1_logi", "b1_gom","c1_logi", "c1_gom" )


df_sum2_ <- df_sum2[-seq(2, 93, 3),]
b <- data.frame(matrix(NA, nrow=31, ncol=6))
for (i in 1:31) {
  for (j in 4:9){
    b[i,j-3] <- paste0(substr(df_sum2_[2*i - 1, j], 1,13), " (", substr(df_sum2_[2*i,j], 1, 6), ")")
  }
}
b <- cbind(size=df_sum2$X.Samples[seq(2, 93, 3)], factor=df_sum2$Explanatory[seq(2, 93, 3)], b)
b <- b %>% select(1,2,3,6,4,7,5,8)
colnames(b) <- c("size", "factor", "a2_logi", "a2_gom","b2_logi", "b2_gom","c2_logi", "c2_gom" )

df_sum3_ <- df_sum3[-seq(2, 93, 3),]
c <- data.frame(matrix(NA, nrow=31, ncol=6))
for (i in 1:31) {
  for (j in 4:9){
    c[i,j-3] <- paste0(substr(df_sum3_[2*i - 1, j], 1,13), " (", substr(df_sum3_[2*i,j], 1, 6), ")")
  }
}
c <- cbind(size=df_sum3$X.Samples[seq(2, 93, 3)], factor=df_sum3$Explanatory[seq(2, 93, 3)], c)
c <- c %>% select(1,2,3,6,4,7,5,8)
colnames(c) <- c("size", "factor", "a3_logi", "a3_gom","b3_logi", "b3_gom","c3_logi", "c3_gom" )

write.csv(a, "seg1.csv")
write.csv(b, "seg2.csv")
write.csv(c, "seg3.csv")


#### Write result #### 
write.csv(df_sum1, paste0(result_path, "Linear Regression(1st Segment, standardized).csv"), row.names = F)
write.csv(df_sum2, paste0(result_path, "Linear Regression(2nd Segment, standardized).csv"), row.names = F)
write.csv(df_sum3, paste0(result_path, "Linear Regression(3rd Segment, standardized).csv"), row.names = F)

SLR_seg1 <- read.csv(paste0(result_path, "Linear Regression(1st Segment, standardized).csv"))[seq(3, 93, 3),]
SLR_seg2 <- read.csv(paste0(result_path, "Linear Regression(2nd Segment, standardized).csv"))[seq(3, 93, 3),]
SLR_seg3 <- read.csv(paste0(result_path, "Linear Regression(3rd Segment, standardized).csv"))[seq(3, 93, 3),]

significant_SLR_seg <- cbind(SLR_seg1$Explanatory, 
                             SLR_seg1[,4:9] < 0.05,
                             SLR_seg2[,4:9] < 0.05,
                             SLR_seg3[,4:9] < 0.05)

significant_SLR_a <- na.omit(significant_SLR_seg[,c(1,seq(2,19,3))])
significant_SLR_b <- na.omit(significant_SLR_seg[,c(1,seq(3,19,3))])
significant_SLR_c <- na.omit(significant_SLR_seg[,c(1,seq(4,19,3))])



#### Bar Plot for standardized SLR coefficients between Growth Curve Parameter and Time-independent variable ####
## Logistic Model
coef_name <- c("a_Logi", "b_Logi", "c_Logi",
               "a_Gom", "b_Gom", "c_Gom")
for (i in 4:6){
  i=6
  df_Logi <- data.frame(
    Explanatory = factor(unique(df_sum1$Explanatory)),
    coefficient = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==1)], 
                    df_sum2[,i][which((1:nrow(df_sum2))%%3==1)],
                    df_sum3[,i][which((1:nrow(df_sum3))%%3==1)]),
    p.value = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==0)], 
                df_sum2[,i][which((1:nrow(df_sum2))%%3==0)],
                df_sum3[,i][which((1:nrow(df_sum3))%%3==0)]),
    segment = factor(rep(1:3, each=length(unique(df_sum1$Explanatory)))),
    category = factor(c(7,7,1,1,7,5,4,4,7,3,3,6,5,5,2,1,1,1,2,2,2,1,7,7,7,7,7,7,7,2,5))
    # Territory&population : 1 / Age : 2 / Trade : 3 / Environment : 4 / 
    # Education&Culture&Science : 5 / National Accounts : 6 / 
    # Health&Sociecy&Welfare : 7
  ) %>% arrange(category)
  
  
  # volcano plot
  vol_plot <- ggplot(df_Logi, aes(x=coefficient, y=-log(p.value),
                        group=Explanatory, color=category)) +
    geom_point(aes(shape=segment)) +
    geom_line() +
    labs(x="coefficient", y="-log(p.value)", 
         title=paste0("SLR Volcano Plot for ", coef_name[i-3])) +
    theme_classic() +
    theme(plot.title=element_text(size=15, hjust=0.5, face="bold", vjust=2))
  ggsave(paste0(plot_path, coef_name[i-3],"_volcano_plot.png"), plot = vol_plot)

  # bar plot
  bar_plot <- ggplot(df_Logi %>%
           mutate(Explanatory = fct_reorder(Explanatory, desc(category))),
         aes(x=Explanatory, y=-log(p.value), fill=segment, label=category)) +
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept=-log(0.05), color = 2, size=0.5, linetype=3) +
    labs(x="National Factors", title = "P-value of *&theta;*<sub>1</sub> of *&gamma;* ~ National Factors(Logistic)") +
    coord_flip() +
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#DD07E0")) +
    # scale_fill_hue(l=70, c=100) +
    theme_classic() +
    theme(
      plot.title= element_markdown(size=10, hjust=0.5, vjust=2, face="bold"),
      legend.title = element_text(face="bold", size=10),
      legend.text = element_text(size=8)
      )
  ggsave(paste0(plot_path, coef_name[i-3],"_bar_plot.png"), plot = bar_plot)
}

set.seed(20)
df <- data.frame(d13C = rnorm(20, -23, 5),
                 DIC = rnorm(20, 4, 0.2),
                 d13CDIC = rnorm(20, -8, 2))

ggplot(df, aes(x = d13C, y = d13CDIC)) +
  geom_point(aes(fill = DIC), pch = 21, cex = 5) + 
  labs(
    x = "*&delta;*<sup>13</sup>C (&permil; VPDB)",
    y = "*&delta;*<sup>13</sup>C<sub>DIC</sub> (&permil; VPDB)",
    title = "*&theta;*<sub>1</sub>"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    plot.title = element_markdown()
  )


## Gompertz Model
for (i in 7:9){
  i=7
  df_Gomp <- data.frame(
    Explanatory = factor(unique(df_sum1$Explanatory)),
    coefficient = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==1)], 
                    df_sum2[,i][which((1:nrow(df_sum2))%%3==1)],
                    df_sum3[,i][which((1:nrow(df_sum3))%%3==1)]),
    p.value = c(df_sum1[,i][which((1:nrow(df_sum1))%%3==0)], 
                df_sum2[,i][which((1:nrow(df_sum2))%%3==0)],
                df_sum3[,i][which((1:nrow(df_sum3))%%3==0)]),
    segment = factor(rep(1:3, each=length(unique(df_sum1$Explanatory)))),
    category = factor(c(7,7,1,1,7,5,4,4,7,3,3,6,5,5,2,1,1,1,2,2,2,1,7,7,7,7,7,7,7,2,5))
    # Territory&population : 1 / Age : 2 / Trade : 3 / Environment : 4 / 
    # Education&Culture&Science : 5 / National Accounts : 6 / 
    # Health&Sociecy&Welfare : 7
  ) %>% arrange(category)
  
  
  # volcano plot
  vol_plot <- ggplot(df_Gomp, aes(x=coefficient, y=-log(p.value),
                                  group=Explanatory, color=category)) +
    geom_point(aes(shape=segment)) +
    geom_line() +
    labs(x="coefficient", y="p-value", 
         title=paste0("SLR Volcano Plot for ", coef_name[i-3])) +
    theme_classic() +
    theme(plot.title=element_text(size=15, hjust=0.5, face="bold", vjust=2))
  ggsave(paste0(plot_path, coef_name[i-3],"_volcano_plot.png"), plot = vol_plot)
  
  # bar plot
  bar_plot <- ggplot(df_Gomp %>%
                       mutate(Explanatory = fct_reorder(Explanatory, desc(category))),
                     aes(x=Explanatory, y=-log(p.value), fill=segment, label=category)) +
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept=-log(0.05), color = 2, size=0.5, linetype=3) +
    labs(x="National Factors", title = "P-value of *&theta;*<sub>1</sub> of *&alpha;* ~ National Factors(Gompertz)") +
    coord_flip() +
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#DD07E0")) +
    # scale_fill_hue(l=70, c=100) +
    theme_classic() +
    theme(
      plot.title= element_markdown(size=10, hjust=0.5, vjust=2, face="bold"),
      legend.title = element_text(face="bold", size=10),
      legend.text = element_text(size=8)
    )
  ggsave(paste0(plot_path, coef_name[i-3],"_bar_plot.png"), plot = bar_plot)
}
