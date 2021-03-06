#########################################
#### set_date : setting study period ####
#########################################

set_date <- function(x=0){
  if(x==0){
    max_date <- as.integer(as.Date(Sys.Date(),formate="%Y/%m/%d")-
                             as.Date("2019-12-31",formate="%Y/%m/%d"))
  }else{
    max_date <- as.integer(as.Date(x,formate="%Y/%m/%d")-
                             as.Date("2019-12-31",formate="%Y/%m/%d")) 
  }
  return(max_date)
}



#########################################
#### preprocessing data for analysis ####
#########################################

preprocessing_data <- function(){
  df_sum<-NULL
  for(i in country){
    df_sub<-select(as.tbl(df[df$countriesAndTerritories==i,]) ,
                   cases,dateRep) %>%
      mutate(Date=as.Date(as.character(dateRep),format='%d/%m/%Y')) %>%
      group_by(Date)%>%
      summarise(new_case=sum(cases))
    
    colnames(df_sub) <- c("Date",i)
    
    if(is.null(df_sum)){
      df_sum <- df_sub
    }else{
      df_sum <- full_join(df_sum, df_sub,by='Date')
    }
  }
  df_sum<-df_sum %>% arrange(Date)
  df_sum = as.data.frame(df_sum)
  
  # negative, NA value imputation
  for(i in 1:length(country)){
    df_sum[c(which(is.na(df_sum[,i+1])), which(df_sum[,i+1]<0)), i+1] <- 0
  }  
  return(df_sum)
}



###########################################################
#### categorizing/properties analysis using derivative ####
###########################################################

derivative_analysis <- function(Country,criteria=0.1, 
                        start_date=1,end_date=max_date, 
                        exclusion_criteria=50,gkf_bandwidth=14,
                        save_image = FALSE,
                        save_excel= FALSE){
  
  # objects for storing results 
  current_state_box <- c()
  peak_box <- data.frame()
  break_box <- data.frame()
  peak_number_box <- c()
  break_number_box <- c()
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
    if(max(y)<exclusion_criteria) next
    
    # gaussian kernel estimation & 1st, 2nd derivative
    gkf <- ksmooth(x,y,kernel="normal",bandwidth = gkf_bandwidth)$y  # gaussian kernel 
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
    
    ### Break points
    break_point <- c(min(which(y>=1)))
    
    if(length(exp_point)>1){
      for(j in 1:(length(exp_point)-1)){
        break_point <- c(break_point,exp_point[j]+which.min(gkf[exp_point[j]:exp_point[j+1]]))
      }  
    }
    
    ### Save result
    country_list <- c(country_list,i)
    current_state_box <- c(current_state_box,current_state)
    peak_number_box <- c(peak_number_box,length(exp_point))
    break_number_box <- c(break_number_box,length(break_point))
    max_spread_box <- c(max_spread_box,max_deriv)
    max_suppression_box <- c(max_suppression_box,min_deriv)
    peak_box <- rbind(peak_box,c(exp_point,rep(NA,5-length(exp_point))))
    break_box <- rbind(break_box,c(break_point,rep(NA,5-length(break_point))))
    
    # print/write result image
    if(save_image){
      ### Plotting
      png(filename=paste0("image/",i,".png"))
      par(mfrow=c(3,1))
      
      plot(gkf,type="l",xlab="",ylim=c(0,max(y)))
      title(main=paste(i,ifelse(current_state<1," - Stabilized"," - Unstabilized"),
                       "\n 1) # of peaks : ",length(exp_point),"   2) current state : ",current_state,
                       "   3) control index : ",max_deriv/min_deriv))
      abline(h=0,lty=2,col=2)
      for(j in exp_point){
        abline(v=j,col=4)}
      for(j in break_point){
        abline(v=j,col=5,lty=2)}
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
      for(j in break_point){
        abline(v=j,col=5,lty=2)}
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
  colnames(break_box) <- c("break1","break2","break3","break4","break5")
  
  df_result <- data.frame(country=country_list,current_state=current_state_box,
                          max_spread=max_spread_box,
                          max_suppression=max_suppression_box,
                          control_index=max_spread_box/max_suppression_box,
                          peak_number=peak_number_box,
                          break_number=break_number_box)
  
  df_result <- cbind(df_result,peak_box,break_box)
  
  if(save_excel){
    write.csv(df_result,"categorizaioin_result.csv")  
  }
  return(df_result)
}



#######################################################################
#### segPoisson : segmented poisson regression using design matrix ####
#######################################################################

segPoisson <- function(Country,break_point = c(),
                       start_date = 1, end_date = max_date){
  
  temp = which(colnames(df_sum) == Country)
  t = ifelse(start_date == 1, 1, which(daily > 0)[1])
  x = t : end_date - (t-1)
  y = df_sum[t : end_date, temp]
  n = length(x)
  
  # design matrix
  if(length(break_point) > 0){
    m = length(break_point)
    b = break_point
    X_design = cbind(x, log(x))
    for(j in 1:m){
      X_design = cbind(X_design, c(rep(0,b[j]-1), x[b[j]:n] - (b[j])),
                       c(rep(0, b[j]-1), log(x[b[j]:n] - (b[j]-1))))
    }
  } else {
    X_design = cbind(x, log(x))
  }
  
  X_design = data.frame(X_design)
  fit <- try(glm(y ~  . , family = poisson(), data = X_design))
  
  plot(x, y, 
       xlab = paste("Days since", as.Date(start_date,origin = "2019-12-31")),
       ylab = "Daily cases",
       main = Country)
  
  cat(paste0("-----------", Country, "-----------\n"))
  lines(fitted(fit), col=2, lty=2)
  return(fit)
}



################################################
#### cumulative analysis : cumulative graph ####
################################################

cumulative_analysis <- function(Country, criteria=0.1, end_date=max_date, 
                                exclusion_criteria=1000, gkf_bandwidth=14,
                                save_image = FALSE){
  
  # Categorization Algorithm
  for(i in Country){
    
    ### Initialization
    temp <- which(colnames(df_sum)==i)
    
    # since first case of covid_19
    t=min(which(df_sum[,i]>=1))
    x = t:end_date - (t-1)
    y <-cumsum(df_sum[t:end_date,temp])
    if(max(y)<exclusion_criteria) next
    
    # gaussian kernel estimation & 1st, 2nd derivative
    gkf <- ksmooth(x,y,kernel="normal",bandwidth = gkf_bandwidth)$y  # gaussian kernel 
    
    # print/write result image
    if(save_image){
      ### Plotting
      png(filename=paste0("image/",i,".png"))
      
      dev.off()
    }else{
      par(mfrow=c(1,1))
      plot(y,type="l",col=1)
      title(main=paste(i))
    }
    
    ## Logistric
    if(df_result$break_number[df_result$country==i]==1){
      
      st = data.frame(a = c(max(y), max(y)*10), b = c(-15, 15), c = c(0, 0.2))
      fit11 = nls2(y ~ a/(1+exp(b-c*x)), start = st,
                   algorithm = "plinear-random",control = nls.control(maxiter = 3000))
      initial1 = coef(fit11)
      print(i)
      ## fitting
      fit12 = try(nls(y ~ a/(1+exp(b-c*x)), 
                      start = list(a = initial1[1],b = initial1[2],c = initial1[3])))
      
      if(is.atomic(fit12)) next
      try(lines(x,(fitted(fit12)),col=2))
    }else{
      for(j in 1:df_result$break_number[df_result$country==i]){
        x_start = df_result[df_result$country==i,12+j]
        x_end = ifelse(is.na(df_result[df_result$country==i,13+j]),max(x),df_result[df_result$country==i,13+j])
        x_=x[x_start:x_end]
        y_=y[x_] - y[x_start-1]
        x_=x_-x_start
        st = data.frame(a = c(0.5*max(y_), max(y_)*30), b = c(0, 25), c = c(0, 0.5))
        
        fit11 = try(nls2(y_ ~ a/(1+exp(b-c*x_)), start = st,
                   algorithm = "plinear-random",control = nls.control(maxiter = 5000)))
        if(is.atomic(fit11)) next
        initial1 = coef(fit11)
        print(i)
        print(x_)
        print(y_)
        ## fitting
        fit12 = try(nls(y_ ~ a/(1+exp(b-c*x_)), 
                      start = list(a = initial1[1],b = initial1[2],c = initial1[3])))
      
        if(is.atomic(fit12)) next
        try(lines(x_+x_start,(fitted(fit12)+y[x_start]),col=2))
        try(lines((x_+x_start)[-1],10*diff(fitted(fit12)+y[x_start]),col=3))
        print(coef(fit12))
        
        lines(x[-1],10*diff(y),col=4)
      }
    }
  }
}


