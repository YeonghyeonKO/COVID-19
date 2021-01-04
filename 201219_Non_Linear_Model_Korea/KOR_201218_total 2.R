rm(list = ls()); gc()

library(car)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# set date for data
data_date <- "201218"
today <- Sys.Date()
today <- format(today, format="%y%m%d")


# set directories for data and saving
path <- "/home/gjhuh/covid_19/"
setwd(path)
data_path <- paste0(path, "Korea_Prep/Data/")
save_path <- paste0(path, "segmented_poisson/", today, "_report/")
dir.create(save_path)

# Loading Data
bp_df <- read.csv(paste0(data_path,"segmented_point_manual_201218.csv"), stringsAsFactors = F, encoding = "UTF-8")
data_add_orig = read.csv(paste0(data_path,"KOR_additive_measure_cluster_", data_date, ".csv"), stringsAsFactors = F, encoding = "UTF-8", row.names = 1)

formula_variables = c("ClosingIndex", "RestrictionIndex", "EconomicIndex", "HealthIndex", "StringencyIndex", "GovernmentIndex")
multiple_variables <- c("ClosingIndex", "RestrictionIndex", "EconomicIndex", "HealthIndex", "GovernmentIndex")
province_list <- list(Domestic = "전국(지역)", Capital = "수도권(서울+경기+인천)", NonCapital = "비수도권")

merged_data <- data_add_orig[, c("Date", "Days", "Domestic_ConfirmedCases", "Domestic_DailyConfirmed", 
                                 "Capital_ConfirmedCases", "Capital_DailyConfirmed", 
                                 "NonCapital_ConfirmedCases", "NonCapital_DailyConfirmed", 
                                 formula_variables)]

merged_data$GovernmentIndex <- merged_data$GovernmentIndex * (100/3)

# setting range of lagging effect
lag_days <- seq(0, 10)

{
  # correlation checking stage
  corr_mat <- matrix(0, nrow = length(formula_variables), ncol = length(formula_variables))
  colnames(corr_mat) <- rownames(corr_mat) <- formula_variables
  for (policy_1 in formula_variables) {
    for (policy_2 in formula_variables) {
      corr_mat[policy_1, policy_2] = cor(merged_data[policy_1], merged_data[policy_2], method = "pearson")
    }
  }
  # par(mar = c(5, 4, 4, 4) + 0.3)
  # heatmap(corr_mat, scale = "none")
}


single_model_idx <- c("Without Index", formula_variables)
multiple_model_idx <- "Multiple Indices"

index_list <- list()
for (each_policy in setdiff(formula_variables, c("GovernmentIndex"))) {
  temp_idx <- c()
  for (each_idx in seq(100, 0, length.out = 6)) {
    if (each_idx >= min(merged_data[(nrow(merged_data)-90):nrow(merged_data), each_policy])) {
      temp_idx <- c(temp_idx, each_idx)
    } else {
      temp_idx <- c(temp_idx, each_idx)
      break
    }
  }
  index_list[[each_policy]] <- temp_idx
}
index_list[["GovernmentIndex"]] <- rev(seq(2, 3, by = 0.5))*(100/3)


data_names <- colnames(merged_data)
colors <- brewer.pal(n = 8, name = "Dark2")

for(sel_prov in names(province_list)){
  # sel_prov <- names(province_list)[1]
  # sel_prov <- names(province_list)[2]
  target_variable <- colnames(merged_data)[which(colnames(merged_data) == paste(sel_prov, "DailyConfirmed", sep = "_"))]
  save_path_region <- paste0(save_path, sel_prov, "/")
  dir.create(save_path_region)
  
  mult_coef_total <- mult_pvalue_total <- mult_min_max_total <- min_max_total <- seg4_coef_total <- seg4_pvalue_total <- c()
  
  plots <- list()
  legend_list <- list()
  
  for(lag_day in lag_days){
    # lag_day <- 0
    # lag_day <- 5
    # lag_day <- 8
    
    data = merged_data
    data_last <- NA
    residual_list <- list()
    
    if(lag_day > 0){
      # 확진자 수 정보는 없지만 정책 정보는 있는 후반부 데이터
      data_last <- data[(nrow(data) - lag_day + 1):nrow(data),]
      data_last[, grepl("Confirmed", colnames(data_last))] <- NA
      
      data <- cbind(data[-(1:lag_day),!(colnames(data) %in% formula_variables)], 
                    data[1:(nrow(data) - lag_day),formula_variables])
      data <- subset(data, select = data_names)
    }
    
    data <- data %>% filter(Days > 0)
    model_idx <- nrow(data)
    
    # assigning breakpoints
    df_break_sel <- bp_df %>% filter(province == province_list[[sel_prov]])
    idx_seg4 <- which(!is.na(df_break_sel$psi_3))
    
    
    ########################################### breakpoint 3 ###########################################
    
    
    # making design matrix for seg4 segmented poisson regression
    {
      seg4_coef_colnames <- c("b0", "b11", "b21", "b12", "b22", "b13", "b23", "b14", "b24")
      seg4_default_colnames <- c(seg4_coef_colnames, target_variable)
      
      seg4_bp <- c(1, floor(c(df_break_sel$psi_1[idx_seg4], df_break_sel$psi_2[idx_seg4], df_break_sel$psi_3[idx_seg4]))) # breakpoints
      seg4_dsg_mat = matrix(0, nrow = nrow(data), ncol = (1+length(seg4_bp)*2+1)) # 10 = 1 + 4*2 + 1
      
      # days, log(days)
      for (i in 1:length(seg4_bp)) {
        each_bp <- seg4_bp[i]
        seg4_dsg_mat[each_bp:nrow(seg4_dsg_mat), (1:2)+(2*i)-1] <- 
          matrix(c(((1:(nrow(seg4_dsg_mat)-each_bp+1))-1), log(1:(nrow(seg4_dsg_mat)-each_bp+1))), byrow = F, ncol = 2)
      }
      
      # add confirmed case data
      seg4_dsg_mat[, ncol(seg4_dsg_mat)] = data[, target_variable]
      
      seg4_dsg_mat = data.frame(seg4_dsg_mat)
      colnames(seg4_dsg_mat) = seg4_default_colnames
      
      seg4_dsg_mat <- cbind(seg4_dsg_mat, data[, formula_variables])
      
      # # add remained data from lagging day
      # if(lag_day != 0){
      #   seg4_dsg_mat[(nrow(data)+1):(nrow(data)+lag_day), formula_variables] <- data_last[,formula_variables]
      # }
      
      # add rows for prediction
      seg4_dsg_mat[(nrow(data)+1):(nrow(data)+56), target_variable] = NA
      for (i in 1:length(seg4_bp)) {
        seg4_dsg_mat[nrow(data) + (1:56), i*2] <- seg4_dsg_mat[nrow(data), i*2] + (1:56)
        seg4_dsg_mat[nrow(data) + (1:56), i*2+1] <- log(seg4_dsg_mat[nrow(data), i*2] + (1:56)+1)
      }
      
      seg4_dsg_mat[, 1] = 1
    }
    
    seg4_model_formula_list <- as.list(c(
      # Days + log Days
      paste0(target_variable, " ~ ", 
             paste(seg4_coef_colnames, collapse = "+"), 
             " -1"), 
      # Days + log Days + Policy
      paste0(target_variable, " ~ ", 
             paste(paste(seg4_coef_colnames, collapse = "+"), 
                   formula_variables, sep = "+"),
             " -1"), 
      paste0(target_variable, " ~ ", 
             paste(c(seg4_coef_colnames, multiple_variables), collapse = "+"), 
             " -1")))
    
    names(seg4_model_formula_list) <- c(single_model_idx, multiple_model_idx)
    
    # seg4 model list
    seg4_model_list = list()
    
    for(m in 1:length(seg4_model_formula_list)){
      set.seed(777)
      glm_result <- NULL
      try({
        glm_result = glm(as.formula(seg4_model_formula_list[[m]]), family = 'poisson', data = seg4_dsg_mat[1:model_idx, ], maxit = 100)
      })
      
      if (!is.null(glm_result)) {
        seg4_model_list[[names(seg4_model_formula_list)[m]]] = glm_result
      }
    }
    
    # result saving single
    {
      seg4_coef = matrix(nrow = length(single_model_idx), ncol = length(seg4_coef_colnames) + 3) # policy + 2 MSE
      rownames(seg4_coef) = single_model_idx
      colnames(seg4_coef) <- c(seg4_coef_colnames, "Policy","train_MSE", "test_MSE")
      seg4_pvalue = seg4_coef
      for (each_policy in single_model_idx) {
        temp_coef <- seg4_model_list[[each_policy]]$coefficients
        temp_p <- summary(seg4_model_list[[each_policy]])$coefficients[,4]
        
        seg4_coef[each_policy, seg4_coef_colnames] <- temp_coef[which(names(temp_coef) %in% seg4_coef_colnames)]
        seg4_coef[each_policy, "Policy"] <- temp_coef[strsplit(each_policy, split = "_")[[1]][1]]
        seg4_pvalue[each_policy, seg4_coef_colnames] <- temp_p[which(names(temp_p) %in% seg4_coef_colnames)]
        seg4_pvalue[each_policy, "Policy"] <- temp_p[strsplit(each_policy, split = "_")[[1]][1]]
        
        y_hat <- exp(predict(seg4_model_list[[each_policy]], newdata = seg4_dsg_mat[1:model_idx, ]))
        y_act <- seg4_dsg_mat[1:model_idx, target_variable]
        seg4_coef[each_policy, "train_MSE"] <- mean((y_hat - y_act)^2)
        
      }
      
      seg4_coef_total <- rbind(seg4_coef_total, cbind(seg4_coef, rep(lag_day, nrow(seg4_coef))))
      seg4_pvalue_total <- rbind(seg4_pvalue_total, cbind(seg4_pvalue, rep(lag_day, nrow(seg4_pvalue))))
    }
    
    pred_mat <- seg4_dsg_mat[model_idx:nrow(seg4_dsg_mat), seg4_coef_colnames]
    
    min_max_mat = matrix(nrow = length(names(seg4_model_formula_list)), ncol = 8) # policy + 2 MSE
    rownames(min_max_mat) = names(seg4_model_formula_list)
    colnames(min_max_mat) <- c(1:7, "diff")
    
    # plots
    for (each_policy in single_model_idx) {
      p <- ggplot(data.frame(x = 1:model_idx, y = data[1:model_idx, target_variable]), aes(x = x, y = y)) +
        geom_line(color = "gray50") +
        labs(title = each_policy, x = "Days from First Case (1/20)", y = "Daily Cases", subtitle = paste0("# lag = ", as.character(lag_day))) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = 2),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14), 
              legend.position="right", 
              legend.title=element_text(size=10), 
              legend.text=element_text(size=8), 
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
              legend.key.width = unit(1, "cm")) + 
        coord_cartesian(ylim = c(0, 3050))
      
      if (each_policy == "Without Index") {
        p <- p + geom_line(data = data.frame(x = (1:nrow(seg4_dsg_mat)), y = exp(predict(seg4_model_list[[each_policy]], newdata = seg4_dsg_mat))), 
                           mapping = aes(x = x, y = y), color = "red")
      } else {
        p <- p + geom_line(data = data.frame(x = (1:model_idx), y = exp(predict(seg4_model_list[[each_policy]], newdata = seg4_dsg_mat[1:model_idx, ]))), 
                           mapping = aes(x = x, y = y), color = "red")
        
        index_range <- index_list[[each_policy]]
        pred_list = c()
        plot_df <- c()
        plot_default_idx = model_idx # + length(test_idx)
        for (idx in length(index_range):1) {
          policy_index <- c(seg4_dsg_mat[model_idx, each_policy], rep(index_range[idx], nrow(pred_mat)-1))
          temp_data <- cbind(pred_mat, policy_index)
          colnames(temp_data) <- c(seg4_coef_colnames, each_policy)
          
          temp_plot_df <- data.frame(x = (plot_default_idx + 1:nrow(pred_mat) -1), 
                                     y = exp(predict(seg4_model_list[[each_policy]], newdata = temp_data)), 
                                     index = as.character(index_range[idx]))
          plot_df <- rbind(plot_df, temp_plot_df)
          pred_list <- cbind(pred_list, exp(predict(seg4_model_list[[each_policy]], newdata = temp_data)))
        }
        plot_df$index <- factor(plot_df$index, levels = index_range)
        max_row <- pred_list[which.max(as.numeric(apply(pred_list, 1, function(vec){return(max(vec) - min(vec))}))),]
        min_max_mat[each_policy, 1:length(max_row)] =  max_row
        min_max_mat[each_policy, "diff"] = max(max_row)-min(max_row)
        
        if (grepl(pattern = "Government", each_policy)) {
          model_labels <- factor(as.character(index_range/(100/3)))
        } else {
          model_labels <- factor(as.character(index_range))
        }
        
        col_list <- colors[1:length(index_range)]
        p <- p + geom_line(data=plot_df, 
                           aes(x=x, y=y, colour = index, group = index), show.legend = T) + 
          scale_color_manual(name = "Index",  labels = model_labels, values = col_list)
        
        if (is.null(legend_list[[each_policy]])) {
          legend_list[[each_policy]] <- g_legend(p)
        }
      }
      
      # saving plots
      plots[[paste(each_policy, lag_day, sep = "_")]] <- p
    }
    
    min_max_total <- rbind(min_max_total, cbind(min_max_mat, rep(lag_day, nrow(min_max_mat))))
  }
  
  for (each_plot in names(plots)) {
    g <- plots[[each_plot]]
    ggsave(filename = paste0(save_path_region, paste(target_variable, "seg4", each_plot, "lag", sep = "_"), ".png"),
           plot = g, units = 'in', width = 8, height = 6, dpi = 300)
  }
  
  
  # save results
  write.csv(seg4_coef_total, paste0(save_path_region, paste(target_variable, "seg4", sep = "_"), "_coef.csv"))
  write.csv(seg4_pvalue_total, paste0(save_path_region, paste(target_variable, "seg4", sep = "_"), "_pvalue.csv"))
  write.csv(min_max_total, paste0(save_path_region, paste(target_variable, "seg4", sep = "_"), "_min_max_total.csv"))
  write.csv(mult_coef_total, paste0(save_path_region, paste(target_variable, "seg4", sep = "_"), "_mult_coef.csv"))
  write.csv(mult_pvalue_total, paste0(save_path_region, paste(target_variable, "seg4", sep = "_"), "_mult_pvalue.csv"))
  write.csv(mult_min_max_total, paste0(save_path_region, paste(target_variable, "seg4", sep = "_"), "_mult_min_max_total.csv"))
}
