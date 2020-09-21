current_state_analysis = function(Country, 
                                   exclusion_criteria=1000,
                                   gkf_bandwidth=14,
                                start_date=1,end_date=max_date){
  
  # objects for storing results 
  current_state_result <- NULL
  country_list <- c()

  # Categorization Algorithm
  for(i in Country){
    
    ### Initialization
    temp <- which(colnames(df_sum)==i)
    
    # since 2020/01/01
    t=start_date
    x = t:end_date - (t-1)
    y <-df_sum[t:end_date,temp]
    if(max(y)<exclusion_criteria) next
    
    # gaussian kernel estimation & 1st, 2nd derivative
    gkf <- ksmooth(x,y,kernel="normal",bandwidth = gkf_bandwidth)$y  # gaussian kernel 
    deriv_gkf <- diff(gkf)    
    
    ind = 0
    date = c()
    current_state = c()
    
    for(j in seq(from=100,to=end_date,by=1)){
      ind = ind + 1
      date[ind] = j
      current_state[ind] = gkf[j]/max(deriv_gkf[1:j-1])
    }
    
    if(is.null(current_state_result)){
      current_state_result = data.frame(current_state)
    }else{
      current_state_result = cbind(current_state_result,current_state)
    }
    
    country_list = c(country_list,i)
  }
  
  colnames(current_state_result) = country_list
  rownames(current_state_result) = date
  return(current_state_result)
}

data.frame(c(1,2,3,4))
is.null(a)

current = current_state_analysis("Germany")
current
plot(current$Germany)
lines(df_sum$Germany[101:248]*0.005)



# Global
par(mfrow=c(1,1))
hist((df_result$current_state),main="Histogram of Current State( ~ 8/31)",xlab=expression(S[t_max]),breaks = 40,xlim=c(0,80),freq=FALSE,ylim=c(0,0.1))
barplot(table(df_result$peak_number),main="Histogram of number of Peaks( ~ 8/31)",xlab=" # of peak")
hist(df_result$peak1, main="Histogram of 1st-peak timestamp( ~ 8/31)",xlab=expression(t[peak_1]),xlim=c(0,250))
hist((df_result$control_index[df_result$control_index<=30]),main="Histogram of Control Index( ~ 8/31)",xlab=expression(C[t_max]),breaks = 20,xlim=c(0,30))

# Europe
europe <- data.frame(country=c(continent$country[continent$continent=="Europe"]))
df_result_eu<-left_join(europe,df_result)
hist((df_result_eu$current_state),main="Histogram of Current State( ~ 8/31, Europe)",xlab=expression(S[t_max]),breaks = 20,xlim=c(0,50),freq=FALSE,ylim=c(0,0.2))
barplot(table(df_result_eu$peak_number),main="Histogram of number of Peaks( ~ 8/14, Europe)",xlab=" # of peak")
hist(df_result_eu$peak1, main="Histogram of 1st-peak timestamp( ~ 8/14, Europe)",xlab=expression(t[peak_1]),xlim=c(0,250))
hist((df_result_eu$control_index[df_result_eu$control_index<=30]),main="Histogram of Control Index( ~ 8/14, Europe)",xlab=expression(C[t_max]),breaks = 20,xlim=c(0,10))





