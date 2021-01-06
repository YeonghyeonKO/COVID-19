derivative_analysis <- function(Country,criteria=0.1, 
                                start_date=1,end_date=max_date, 
                                exclusion_criteria=50,gkf_bandwidth=14,
                                save_image = FALSE,
                                save_excel= FALSE,
                                first_day_criteria = 50){
  
  # objects for storing results 
  current_state_box <- c()
  peak_box <- data.frame()
  peak_number_box <- c()
  max_spread_box <- c()
  max_suppression_box <- c()
  control_index_box <- c()
  df_result <- data.frame()
  country_list <- c()
  
  # criteria constant 
  c <- criteria
  par(mfrow=c(3,1))
  
  # Categorization Algorithm
  for(i in Country){
    
    ### Initialization
    temp <- which(colnames(df_sum)==i)
    
    # since 2020/01/01
    t=start_date
    x = t:end_date - (t-1)
    y <-df_sum[t:end_date,temp]
    
    # gaussian kernel estimation & 1st, 2nd derivative
    gkf <- ksmooth(x,y,kernel="normal",bandwidth = gkf_bandwidth)$y  # gaussian kernel 
    if(max(gkf)<exclusion_criteria) next
    deriv_gkf <- diff(gkf)                                # 1st derivative 
    deriv2_gkf <- diff(deriv_gkf)                         # 2nd derivative
    
    max_gkf <- max(gkf)
    max_deriv <- max(deriv_gkf)
    min_deriv <- -min(deriv_gkf)
    current_state <- sum(deriv_gkf)/max(deriv_gkf)
    
    # criteria constant 
    c0 <- 0.2
    c1 <- max(abs(deriv_gkf))*c
    c2 <- max(abs(deriv2_gkf))*c
    
    
    ### Peak Detection
    exp_point <- c()
    non_peak_point <- c()
    non_peak_point2 <- c()
    
    # condition_0 : local maximum
    for(j in 4:(length(x)-4)){
      if(deriv_gkf[j]*deriv_gkf[j+1]<=0 & deriv2_gkf[j+1]<= -c2) exp_point <-c(exp_point,j+1)
    }
    
    # condition_1 : exclude small peaks (peak value must be at least more than 0.2*max_gkf)
    if(length(exp_point)>0){
      for(j in 1:length(exp_point)){
        if(gkf[exp_point[j]] < c0*max_gkf){
          non_peak_point <- c(non_peak_point,j) 
        }
      }
      if(!is.null(non_peak_point))  exp_point <- exp_point[-(non_peak_point)]
    }
    
    # condition_2 : resolution criteria (able to distinguish between 2 peaks)  
    if(length(exp_point)>1){ 
      mark = 0
      k=2
      n=length(exp_point)
      for(j in 1:(n-1) ){
        
        z1 <- exp_point[k-1]
        z2 <- exp_point[k]
        
        small_peak <- min(gkf[z1],gkf[z2])
        off_peak <- min(gkf[z1:z2])
        
        if(off_peak>small_peak*0.8){
          if(gkf[z1]<=gkf[z2]){
            exp_point<-exp_point[-(k-1)]
          }else{
            exp_point<-exp_point[-k]
          }
          mark = mark+1
        }else{
          k = k+1
        }
      }
    }
    
    # condition_3 : exclude peaks which are not real peaks
    if(length(exp_point)>0){ 
      if(gkf[exp_point[length(exp_point)]]<max(gkf[(exp_point[length(exp_point)]+2):min((exp_point[length(exp_point)]+30),max_date)])){
        exp_point<-exp_point[-length(exp_point)]
      }
    }
    
    ### Save result
    country_list <- c(country_list,i)
    current_state_box <- c(current_state_box,current_state)
    peak_number_box <- c(peak_number_box,length(exp_point))
    max_spread_box <- c(max_spread_box,max_deriv)
    max_suppression_box <- c(max_suppression_box,min_deriv)
    peak_box <- rbind(peak_box,c(exp_point,rep(NA,5-length(exp_point))))
    
    
    # print/write result image
    if(save_image){
      ### Plotting
      png(filename=paste0("segmentation/image/",i,".png"))
      par(mfrow=c(3,1))
      
      plot(gkf,type="l",xlab="",ylim=c(0,max(y)))
      title(main=paste(i,ifelse(current_state<1," - Stabilized"," - Unstabilized"),
                       "\n 1) # of peaks : ",length(exp_point),"   2) current state : ",current_state,
                       "   3) control index : ",max_deriv/min_deriv))
      abline(h=0,lty=2,col=2)
      for(j in exp_point){
        abline(v=j,col=4)}
      lines(y,lty=2)
      
      plot(deriv_gkf,type="l",xlab="");  abline(h=0,col=2)
      abline(h=c1,lty=2,col=3)
      abline(h=-c1,lty=2,col=3)
      for(j in exp_point){
        abline(v=j,col=4)}
      
      plot(deriv2_gkf,type="l",xlab="Since 2020-01-01");  abline(h=0,col=2)
      abline(h=c2,lty=2,col=3)
      abline(h=-c2,lty=2,col=3)
      for(j in exp_point){
        abline(v=j,col=4)}
      
      dev.off()
    }else{
      
      plot(gkf,type="l",xlab="",ylim=c(0,max(y)))
      title(main=paste(i,ifelse(current_state<1," - Stabilized"," - Unstabilized"),
                       "\n 1) # of peaks : ",length(exp_point),"   2) current state : ",current_state,
                       "   3) control index : ",max_deriv/min_deriv))
      abline(h=0,lty=2,col=2)
      for(j in exp_point){
        abline(v=j,col=4)}
      lines(y,lty=2)
      
      plot(deriv_gkf,type="l",xlab="");  abline(h=0,col=2)
      abline(h=c1,lty=2,col=3)
      abline(h=-c1,lty=2,col=3)
      for(j in exp_point){
        abline(v=j,col=4)}
      
      plot(deriv2_gkf,type="l",xlab="Since 2020-01-01");  abline(h=0,col=2)
      abline(h=c2,lty=2,col=3)
      abline(h=-c2,lty=2,col=3)
      for(j in exp_point){
        abline(v=j,col=4)}
    }
  }    
  
  # dataframe for total results
  colnames(peak_box) <- c("peak1","peak2","peak3","peak4","peak5")
  peak_box 
  
  df_result <- data.frame(country=country_list,current_state=current_state_box,
                          max_spread=max_spread_box,
                          max_suppression=max_suppression_box,
                          control_index=max_spread_box/max_suppression_box,
                          peak_number=peak_number_box)
  
  df_result <- cbind(df_result,peak_box)
  
  if(save_excel){
    write.csv(df_result,"categorizaioin_result.csv")  
  }
  return(df_result)
}


df_result = derivative_analysis(Country = country,gkf_bandwidth = 14,first_day_criteria = 50,exclusion_criteria = 100,save_image = TRUE,save_excel = FALSE)













