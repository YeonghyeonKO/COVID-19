coefficient_result = read.csv("coef_result.csv",row.names = 1)
coef_sep_Logi = coefficient_result[,c(1:6,14)]
coef_sep_Gom = coefficient_result[,c(7:12,14)]

# Logistic model counting
country = rownames(coef_sep_Logi)
count1 = 0
count2 = 0
count3 = 0
segment2_len = length(which(!is.na(coef_sep_Logi$breakpoint)))
segment1_len = 165 - segment2_len

for(i in 1:length(country)){
  if(!is.na(coef_sep_Logi$breakpoint[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_sep_Logi$a2_Logi[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_sep_Logi$a1_Logi[i])){
      m = 1
    }else{
      m = 0
    }
  }

  # Case by Case  
  if(M==2&m==0){
    count1 = count1 + 1
  }
  
  if(M==2&m==1){
    count2 = count2 + 1
  }
  
  if(M==1&m==0){
    count3 = count3 + 1
  }
}
for(i in 1){
  cat("-------fitted model counting-------\n")
  
  cat("\n")
  cat("num of countries whose segment num is 2 :", segment2_len,"\n")
  cat("segment num is 2, fitted segment num is 2 :", segment2_len-(count1+count2),"\n")
  cat("segment num is 2, fitted segment num is 1 :", segment2_len-count1,"\n")
  
  cat("\n")
  cat("num of countries whose segment num is 1 :", segment1_len,"\n")
  cat("segment num is 1, fitted segment num is 1 :", segment1_len-count3) 
}
count1;count2;count3



# Gompertz model counting
count1 = 0
count2 = 0
count3 = 0
country = rownames(coef_sep_Gom)
segment2_len = length(which(!is.na(coef_sep_Gom$breakpoint)))
segment1_len = 165 - segment2_len

for(i in 1:length(country)){
  if(!is.na(coef_sep_Gom$breakpoint[i])){
    M = 2          # expected number of parameter pairs
  }else{
    M = 1
  }
  
  if(!is.na(coef_sep_Gom$a2_Gom[i])){
    m = 2          # number of fitted parameter pairs
  }else{
    if(!is.na(coef_sep_Gom$a1_Gom[i])){
      m = 1
    }else{
      m = 0
    }
  }
  
  # Case by Case  
  if(M==2&m==0){
    count1 = count1 + 1
  }
  
  if(M==2&m==1){
    count2 = count2 + 1
  }
  
  if(M==1&m==0){
    count3 = count3 + 1
  }
}
for(i in 1){
  cat("-------fitted model counting-------\n")
  
  cat("\n")
  cat("num of countries whose segment num is 2 :", segment2_len,"\n")
  cat("segment num is 2, fitted segment num is 2 :", segment2_len-(count1+count2),"\n")
  cat("segment num is 2, fitted segment num is 1 :", segment2_len-count1,"\n")
  
  cat("\n")
  cat("num of countries whose segment num is 1 :", segment1_len,"\n")
  cat("segment num is 1, fitted segment num is 1 :", segment1_len-count3) 
}
count1;count2;count3










