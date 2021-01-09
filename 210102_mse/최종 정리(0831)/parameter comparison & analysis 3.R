#### Parameter Comparison
library(GGally)
coefficient_result=read.csv("coef_result.csv",row.names = 1)

# number of fitted parameter pairs
length(which(!is.na(coefficient_result$a1_Logi)))
length(which(!is.na(coefficient_result$a2_Logi)))

length(which(!is.na(coefficient_result$a1_Gom)))
length(which(!is.na(coefficient_result$a2_Gom)))

# scatter plot : logistic seg1 vs logistic seg1
ggpairs(log(coefficient_result[,1:6]))

# scatter plot : gompertz seg1 vs gompertz seg2
ggpairs(log(coefficient_result[,7:12]))

# scatter plot : logistic seg1 vs gompertz seg1
ggpairs(log(coefficient_result[,c(1:3,7:9)]))





