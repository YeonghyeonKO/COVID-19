#### Parameter Comparison
library(GGally)
coefficient_result=read.csv("coef_result.csv",row.names = 1)

# number of fitted parameter pairs
length(which(!is.na(coef_seg_Logi$a1_Logi)))
length(which(!is.na(coef_seg_Logi$a2_Logi)))
length(which(!is.na(coef_seg_Logi$a3_Logi)))

length(which(!is.na(coef_seg_Gom$a1_Gom)))
length(which(!is.na(coef_seg_Gom$a2_Gom)))
length(which(!is.na(coef_seg_Gom$a3_Gom)))


# scatter plot : logistic seg1 vs logistic seg1
ggpairs(log(coef_seg_Logi[,1:9]))

# scatter plot : gompertz seg1 vs gompertz seg2
ggpairs(log(coef_seg_Gom[,1:9]))

# scatter plot : logistic seg1 vs gompertz seg1
ggpairs(cbind(log(coef_seg_Logi[,1:3]),log(coef_seg_Gom[,1:3])))





