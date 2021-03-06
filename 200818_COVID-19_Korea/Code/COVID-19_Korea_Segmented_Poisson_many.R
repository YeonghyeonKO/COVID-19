#install.packages("segmented")
library(segmented)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggrepel)

path <- path <- "C:/Users/김학용/Desktop/COVID-19_Korea_200818[277]/"
data_path <- paste0(path, "Data/")
save_path <- paste0(path, "Plot/")
result_path <- paste0(path, "Result/")

# COVID_19_Confirmed_Korea_Merged_fin: From Kaggle & Public Data Portal
df <- read.csv(paste0(data_path, "COVID_19_Confirmed_Korea_Merged_fin.csv"))

# COVID_19_Confirmed_Korea_Stat: From KOSTAT
#df <- read.csv(paste0(data_path, "COVID_19_Confirmed_Korea_Stat.csv"))

############################################ 0. EDA ############################################

p_1_1 <- ggplot(df, aes(x = Days_after_First_Cases, y = Cases, col = Country_Region)) +
  geom_line() +
  geom_text_repel(
    data = subset(df, Days_after_First_Cases == max(Days_after_First_Cases)),
    aes(x = max(Days_after_First_Cases), label = Country_Region),
    size = 3,
    nudge_x = 45,
    segment.color = NA
  ) +
  labs(title = "COVID-19 Confirmed Cases", subtitle = paste0("Korea, ", as.character(min(as.Date(df$Date))), "~", as.character(max(as.Date(df$Date)))), x = "Days from First Case", y = "Confirmed Cases")+
  theme_bw() +
  theme(
    plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
    plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    axis.text.x = element_text(hjust = 1, size = 12),
    legend.position = "none")
p_1_1
ggsave(p_1_1, filename = "Confirmed_case_in_Korea_by_Region_total.png", path = save_path, width = 12, height = 7, dpi = 300)

p_1_2 <- ggplot(filter(df, Cases > 1), aes(x = Days_after_First_Cases, y = log10(Cases), col = Country_Region)) +
  geom_line() +
  geom_text_repel(
    data = subset(df, Days_after_First_Cases == max(Days_after_First_Cases)),
    aes(x = max(Days_after_First_Cases), label = Country_Region),
    size = 3,
    nudge_x = 45,
    segment.color = NA
  ) +
  labs(title = "COVID-19 Confirmed Cases", subtitle = paste0("Korea, ", as.character(min(as.Date(df$Date))), "~", as.character(max(as.Date(df$Date)))), x = "Days from First Case", y = "Log10(Confirmed Cases)")+
  theme_bw() +
  theme(
    plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
    plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    axis.text.x = element_text(hjust = 1, size = 12),
    legend.position = "none")
p_1_2
ggsave(p_1_2, filename = "Confirmed_case_in_Korea_log_by_Region_total.png", path = save_path, width = 12, height = 7, dpi = 300)


p_2_1 <- ggplot(df, aes(x = Days_after_First_Cases, y = Difference, col = Country_Region)) +
  geom_line() +
  geom_text_repel(
    data = subset(df, Days_after_First_Cases == max(Days_after_First_Cases)),
    aes(x = max(Days_after_First_Cases), label = Country_Region),
    size = 3,
    nudge_x = 45,
    segment.color = NA
  ) +
  labs(title = "COVID-19 Daily Cases", subtitle = paste0("Korea, ", as.character(min(as.Date(df$Date))), "~", as.character(max(as.Date(df$Date)))), x = "Days from First Case", y = "Confirmed Cases")+
  theme_bw() +
  theme(
    plot.title=element_text(size=25, hjust=0.5, face="bold", colour="black", vjust=2),
    plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    axis.text.x = element_text(hjust = 1, size = 12),
    legend.position = "none")
p_2_1
ggsave(p_2_1, filename = "Daily_case_in_Korea_by_Region_total.png", path = save_path, width = 12, height = 7, dpi = 300)

province <- unique(as.character(df$Country_Region))

{
  n <- 0
  coef_Poi <- data.frame(); coef_seg_Poi <- data.frame()
  R_square_list <- c()
  MSE_list <- c()
  psi_list <- c()
}

for(p in province){
  
  n <- n+1
  df_sub <- filter(df, Country_Region == p)

  # Date
  max_date <- max(as.Date(df_sub$Date))
  week_4_date <- max(as.Date(df_sub$Date))+28
  week_8_date <- max(as.Date(df_sub$Date))+56
  
  # For prediction
  df_sub$Date <- as.character(df_sub$Date)
  df_sub_new <- data.frame(Date = as.character(as.Date(origin = max(as.Date(df_sub$Date)), 1:56)),
                           Case_Type = "Confirmed",
                           Cases = NA,
                           Difference = NA,
                           Country_Region = as.character(df_sub$Country_Region[1]),
                           Days_after_First_Cases = (max(df_sub$Days_after_First_Cases)+1):max(df_sub$Days_after_First_Cases+56))
  df_sub_ch <- rbind(df_sub,df_sub_new)
  
  ### Poisson regression
  fit <- glm(Difference ~ log(Days_after_First_Cases) + Days_after_First_Cases, data = df_sub, family = poisson)
  
  SSE_Poi <- sum((predict(fit, data.frame(Days_after_First_Cases = df_sub$Days_after_First_Cases), type = "response") - df_sub$Difference)^2)
  SST_Poi <- sum((df_sub$Difference - mean(df_sub$Difference))^2)
  R_Poi <- 1 - SSE_Poi / SST_Poi
  
  pred_pois <- data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases,
                          pred = exp(predict(fit, df_sub_ch)), Type = "Poisson")
  
  for (k in 1:3){coef_Poi[k, n] <- summary(fit)$coefficients[k, 1]}
  coef_Poi[4, n] <- R_Poi
  colnames(coef_Poi)[n] = p
  rownames(coef_Poi) <- c(rownames(summary(fit)$coefficients)[1:3], "R-Square (Daily)")
  
  pred_seg_pois <- NULL
  fitted_idx <- c(1)
  
  ### Segmented Poisson regression coefficients & R-Sqaure  
  psi_var <- c()
  for(n_seg in 1:3){
    #n_seg <- 1
    seg_fit <- NULL
    set.seed(777)
    try(seg_fit <- segmented(fit, seg.Z = ~ log(Days_after_First_Cases) + Days_after_First_Cases, npsi = n_seg)) # npsi: break point의 개수
    
    if(is.null(seg_fit)){
      psi_var <- c(psi_var, NA)
      pred_seg_pois <- data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases,
                                  pred = NA, Type = paste0("Segmented_",as.character(n_seg+1)))
    }else{
      
      if(n_seg == 1){
        psi_var <- c(psi_var, round(seg_fit$psi[2]))
      }else{
        psi_var <- c(psi_var, paste(as.character(round(seg_fit$psi[,2])), collapse = ","))
      }
      r <- 1
      psi_0 <- 0
      pred_seg_pois <- NULL
      while(r <= n_seg+1){
        
        if(r == 1){
          psi <- round(seg_fit$psi[2])
        }else if(r == n_seg+1){
          psi <- nrow(df_sub)
        }else{
          psi <- round(seg_fit$psi[r,2])
        }
        
        seg_fit_part <- glm(Difference[psi_0+1:psi] ~ log(Days_after_First_Cases[psi_0+1:psi]) + Days_after_First_Cases[psi_0+1:psi],
                   data = df_sub, family = poisson)
        if(n_seg == 2){
          if(is.null(coef_Poi)){
            coef_Poi <- data.frame(coefficients(seg_fit_part))
          }else{
            coef_Poi <- cbind(coef_Poi, data.frame(coefficients(seg_fit_part)))
          }
          
        }
        
        if(is.null(pred_seg_pois)){
          pred_seg_pois <- data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases[psi_0+1:psi],
                                      pred = exp(predict(seg_fit_part, df_sub_ch)),
                                      Type = paste0("Segmented_",as.character(n_seg+1)))
        }else{
          pred_seg_pois <- rbind(pred_seg_pois, data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases[psi_0+1:psi],
                                                           pred = NA,
                                                           Type = paste0("Segmented_",as.character(n_seg+1))))
        }
        psi_0 <- psi
        r <- r+1
      }
      pred_seg_pois <- data.frame(Days_after_First_Cases = df_sub_ch$Days_after_First_Cases,
                                    pred = exp(predict(seg_fit, df_sub_ch)),
                                    Type = paste0("Segmented_",as.character(n_seg+1)))
      fitted_idx <- c(fitted_idx, n_seg+1)
      
    }
    pred_pois <- rbind(pred_pois, pred_seg_pois)
  }

  
  ### Save the results: R-square, MSE, Break Point

  R_square <- c()
  MSE <- c()
  for(n_seg in 0:5){
    seg_type <- levels(pred_pois$Type)[n_seg+1]
    pred_pois_sub <- pred_pois %>% filter(Type == seg_type)
    R_square <- c(R_square, 1 - (sum((df_sub$Difference - pred_pois_sub$pred[1:nrow(df_sub)])^2)) / sum((df_sub$Difference- mean(df_sub$Difference))^2))
    MSE <- c(MSE, sum((df_sub$Difference - pred_pois_sub$pred[1:nrow(df_sub)])^2) / nrow(df_sub))
  }
  
  R_square <- list(R_square)
  names(R_square) <- p
  R_square_list <- c(R_square_list, R_square)

  MSE <- list(MSE)
  names(MSE) <- p
  MSE_list <- c(MSE_list, MSE)
  
  psi_var <- list(psi_var)
  names(psi_var) <- p
  psi_list <- c(psi_list, psi_var)
    
  
  ### 일일 확진자 plot
  Daily <- ggplot(df_sub_ch, aes(x = Days_after_First_Cases, y = Difference)) +
    geom_point(alpha = 0.8, size = 2, color = "grey") +
    scale_y_continuous(limits = c(0, 1.5 * max(df_sub$Difference)))+
    geom_line(pred_pois, mapping = aes(x = Days_after_First_Cases, y = pred, group = Type), col = "midnightblue", size =1.2) +
    geom_vline(xintercept = df_sub_ch$Days_after_First_Cases[df_sub_ch$Date == max_date], 
               col = "grey", size = 1, linetype = 2, alpha = 0.8) +
    labs(x = NULL, y = "Daily Cases") +
    theme_bw()+
    theme(plot.title=element_text(size=20, hjust=0.5, colour="black", vjust=2),
          plot.subtitle=element_text(size=16, hjust=0.5, face="italic", color="maroon", vjust=2),
          axis.title.y = element_text(size = 15),
          legend.position = "none") +
    facet_wrap( ~ Type, ncol = 2) 
  Daily

  
  ### 누적 확진자 plot
  pred_pois_sum <- pred_pois %>% group_by(Type) %>% mutate(pred = cumsum(pred))
  Cumulative <- ggplot(df_sub_ch, aes(x = Days_after_First_Cases, y = Cases)) +
    geom_point(alpha = 0.8, size = 2, color = "grey") +
    scale_y_continuous(limits = c(0, 1.5 * max(df_sub$Cases)))+
    geom_line(pred_pois_sum, mapping = aes(x = Days_after_First_Cases, y = pred, group = Type, col = "midnightblue"), size =1.2) +
    geom_vline(xintercept = df_sub_ch$Days_after_First_Cases[df_sub_ch$Date == max_date], 
               col = "grey", size = 1, linetype = 2, alpha = 0.8) +
    labs(x = NULL, y = "Cumulative Cases") +
    theme_bw()+
    theme(axis.title.y = element_text(size = 15),
          legend.position = "none") +
    facet_wrap( ~ Type, ncol = 2) 
  Cumulative
  
  ### 최종 plot
  p_merged <- grid.arrange(Daily, Cumulative,
                      top = textGrob(p, rot=0, gp = gpar(fontsize=20, bold=T)),
                      bottom = textGrob(paste0("Days From ",as.character(min(as.Date(df_sub$Date)))), rot=0, gp = gpar(fontsize=15)), ncol = 2)

  ggsave(paste0(save_path,"Segmented Poisson in ", p, ".png"), p_merged, dpi = 300, width = 10, height = 6)
}


# # df_pred: Predictive Daily Cases
# df_pred <- data.frame(Max_Cases = current_max,
#                      Predict_Week_4 = pred_week_4,
#                      Predict_Week_8 = pred_week_8)
# row.names(df_pred) <- province
# 
# row.names(coef_Poi) <- c("(Intercept)", "log(Days)", "Days", "R-square")
# 
# row.names(coef_seg_Poi) <- c("(Intercept)_1", "log(Days)_1", "Days_1",
#                              "(Intercept)_2", "log(Days)_2", "Days_2",
#                              "Break_Point", "R-square")
# 
# write.csv(coef_Poi, paste0(result_path, "Poisson_Korea.csv"), row.names = T) # Poisson Coefficients
# write.csv(coef_seg_Poi, paste0(result_path, "Segmented_Poisson_Korea.csv"), row.names = T) # Segmented Poisson Coefficients
# write.csv(df_pred, paste0(result_path, "Segmented_Poisson_Korea_Predict.csv"), row.names = T) # Predicted Daily Cases

df_R_sqr <- data.frame(R_square_list)
df_MSE <- data.frame(MSE_list)
df_psi <- data.frame(psi_list)

row.names(df_R_sqr) <- row.names(df_MSE) <- levels(pred_pois$Type)
row.names(df_psi) <- levels(pred_pois$Type)[-1]

write.csv(df_R_sqr, paste0(result_path, "R_square_Pois.csv"), row.names = T)
write.csv(df_MSE, paste0(result_path, "MSE_Pois.csv"), row.names = T)
write.csv(df_psi, paste0(result_path, "Break_Point_Pois.csv"), row.names = T)
