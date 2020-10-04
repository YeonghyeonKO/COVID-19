#### set_date : setting study period ####
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

#### preprocessing_data : preprocessing data for analysis ####
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
    df_sum[ c(which(is.na(df_sum[,i+1])), which(df_sum[,i+1]<0)),i+1] <- 0
  }  
  return(df_sum)
}


#### derivative_analysis : categorizing/properties analysis using derivative ####
derivative_analysis <- function(Country,criteria=0.1, 
                        start_date="2020-01-01",end_date=max_date, 
                        exclusion_criteria=50,gkf_bandwidth=14,
                        save_image = FALSE,
                        save_excel= FALSE,
                        first_break = 50){
  
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
    t = as.integer(as.Date(start_date,format="%Y-%m-%d") - as.Date("2019/12/31",formate="%Y/%m/%d"))
    x = t:end_date - (t-1)
    y = df_sum[t:end_date,temp]
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
    n = length(exp_point)
    break_point <- c(min(which(cumsum(y)>=first_break)))
    
    if(n>1){
      for(j in 1:(n-1)){
        break_point <- c(break_point,exp_point[j]+which.min(gkf[exp_point[j]:exp_point[j+1]]))
      }  
    }
    
    if(n>=1){
      if(length(which(deriv_gkf[(exp_point[n]+1):length(gkf)] >= c1)) >= 5){
        last_break = exp_point[n] + which.min(gkf[exp_point[n]:length(gkf)])
        break_point = c(break_point,last_break)
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
    
    ### Plotting
    
    if(save_image){png(filename=paste0("image/",i,".png"))}
    par(mfrow=c(3,1))
    
    plot(gkf,type="l",xlab="",ylab="Daily cases",ylim=c(0,max(y)))
    title(main=paste(i,ifelse(current_state<1," - Stabilized"," - Unstabilized"),
                     "\n 1) # of peaks : ",length(exp_point),"   2) current state : ",current_state,
                     "   3) control index : ",max_deriv/min_deriv))
    abline(h=0,lty=2,col=2)
    for(j in exp_point){
      abline(v=j,col=4)}
    for(j in break_point){
      abline(v=j,col=5,lty=2)}
    lines(y,lty=2)
    
    plot(deriv_gkf,type="l",xlab="",ylab=expression(nabla*(Daily_cases)));  abline(h=0,col=2)
    abline(h=c1,lty=2,col=3)
    abline(h=-c1,lty=2,col=3)
    for(j in exp_point){
      abline(v=j,col=4)}
    
    plot(deriv2_gkf,type="l",xlab=paste("Since",start_date),ylab=expression(nabla^2*(Daily_cases)));  abline(h=0,col=2)
    abline(h=c2,lty=2,col=3)
    abline(h=-c2,lty=2,col=3)
    for(j in exp_point){
      abline(v=j,col=4)}
    
    if(save_image){dev.off()}
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

#### segPoisson : segmented poisson regression using design matrix ####
segPoisson <- function(Country,break_point=c(),
                       start_date=1,end_date=max_date,
                       save_image = FALSE){
  
  temp = which(colnames(df_sum)==Country)
  t = ifelse(start_date==1,1,which(daily>0)[1])
  x = t:end_date - (t-1)
  y = df_sum[t:end_date,temp]
  n = length(x)
  
  # design matrix
  if(length(break_point)>0){
    m = length(break_point)
    b = break_point
    X_design = cbind(x,log(x))
    for(j in 1:m){
      X_design = cbind(X_design,c(rep(0,b[j]-1),x[b[j]:n]-(b[j])),
                       c(rep(0,b[j]-1),log(x[b[j]:n]-(b[j]-1))))
    }
  }else{
    X_design = cbind(x,log(x))
  }
  
  X_design = data.frame(X_design)
 
  X_design = data.frame(X_design)
  fit <- try(glm(y ~  . , family = poisson(),data=X_design))
  
  fitted_y = c(rep(0,b[1]-1),fitted(fit)[b[1]:n])
  
  if(save_image){png(filename=paste0("poi_image/",i,".png"))}
  par(mfrow=c(1,1))
  plot(x,y,xlab=paste("Days since",as.Date(start_date,origin = "2019-12-31")),
       ylab="Daily cases",main=Country)
  cat(paste0("-----------",Country,"-----------\n"))
  lines(fitted_y,col=2,lty=2)
  if(save_image){dev.off()}
  
  return(summary(fit))
}


#### segmented Logistic : cumulative graph ####
segLogistic = function(Country,start_date=1,
                        end_date=max_date,
                        break_point=c(),
                        exclusion_criteria=1000,gkf_bandwidth=14,
                       max_iter = 500,
                        save_image = FALSE,
                       prediction_plot = FALSE,
                       daily_plot = FALSE){
  if(save_image){png(filename=paste0("cumulative/",i,"_Logistic.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] #min(which(df_sum[,temp]>0))
  x = t:end_date - (t-1)
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases",main=paste0(Country, " / Logistic"),xlim=c(0,n+30),ylim=c(0,2*max(y))) 
  }else{
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases",main=paste0(Country, " / Logistic")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  b = b - b[1]+1
  tryCatch(
    expr = {
      ## Logistic regression ##
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-1)]
        X_design_ =  X_design[1:(b[2]-1),]
      }
      
      fit11 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[3]-1]
          X_design_ = X_design[1:b[3]-1,] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
        } 
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[4]-1]
          X_design_ = X_design[1:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) + a3/(1+exp(b3-c3*x3)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) + a3/(1+exp(b3-c3*x3)), 
                         data=X_design_,start = initial_list) )
      
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
      }
    },
    error = function(e){
    cat("not fitted! \n")
   },
    finally = {
      # plotting
    if(length(coef_result[,1])==1){
      pred_y = predict(fit12,newdata=X_design_pred)
      lines(x_, pred_y, col=2)
    }
    if(length(coef_result[,1])==2){
      pred_y = predict(fit22,newdata=X_design_pred)
      lines(x_, pred_y, col=2)
    }
    if(length(coef_result[,1])==3){
      pred_y = predict(fit32,newdata=X_design_pred)
      lines(x_, pred_y, col=2)
    }
    abline(v=b[-1],lty=2)
    legend("topright",c("Real value","Estimated value"),pch=c(1,-1),lty=c(0,1),col=c(1,2))
    if(save_image){dev.off()}
    
    # return coefficient
    try(
      expr={
        rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
        print(coef_result)
        return(coef_result)
        }
      )
    }
  )
}



#### segmented Bertalanffy : cumulative graph ####
segBertalanffy = function(Country,start_date=1,
                       end_date=max_date,
                       break_point=c(),
                       exclusion_criteria=1000,gkf_bandwidth=14,
                       max_iter = 500,
                       save_image = FALSE,
                       prediction_plot = FALSE,
                       daily_plot = FALSE){
  if(save_image){png(filename=paste0("cumulative/",i,"_Bertalanffy.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] #min(which(df_sum[,temp]>0))
  x = t:end_date - (t-1)
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Bertalanffy"), xlim=c(0,n+30),ylim=c(0,2*max(y))) 
  }else{
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Bertalanffy")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  b = b - b[1]+1
  tryCatch(
    expr = {
      ## Bertalanffy regression ##
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-1)]
        X_design_ =  X_design[1:(b[2]-1),]
      }
      
      fit11 = try( nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random", control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1], b1=coef(fit11)[2], c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[3]-1]
          X_design_ = X_design[1:b[3]-1,] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
      } 
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[4]-1]
          X_design_ = X_design[1:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) + a3*((1-exp(-b3*x3))^c3) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) + a3*((1-exp(-b3*x3))^c3) , 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      # plotting
      if(length(coef_result[,1])==1){
        pred_y = predict(fit12,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==2){
        pred_y = predict(fit22,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==3){
        pred_y = predict(fit32,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","Estimated value"),pch=c(1,-1),lty=c(0,1),col=c(1,2))
      if(save_image){dev.off()}
      
      # return coefficient
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}


#### segmented Gompertz : cumulative graph ####
segGompertz = function(Country,start_date=1,
                          end_date=max_date,
                          break_point=c(),
                          exclusion_criteria=1000,gkf_bandwidth=14,
                          max_iter = 500,
                          save_image = FALSE,
                          prediction_plot = FALSE,
                          daily_plot = FALSE){
  if(save_image){png(filename=paste0("cumulative/",i,"_Gompertz.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] #min(which(df_sum[,temp]>0))
  x = t:end_date - (t-1)
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Gompertz"), xlim=c(0,n+30),ylim=c(0,2*max(y))) 
  }else{
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Gompertz")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  b = b - b[1]+1
  tryCatch(
    expr = {
      ## Gompertz regression ##
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-1)]
        X_design_ =  X_design[1:(b[2]-1),]
      }
      
      fit11 = try( nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) ,
                        data=X_design_, start = st,
                        algorithm = "plinear-random", control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1], b1=coef(fit11)[2], c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[3]-1]
          X_design_ = X_design[1:b[3]-1,] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
      } 
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[4]-1]
          X_design_ = X_design[1:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) + a3*exp(-b3*exp(-c3*x3)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) + a3*exp(-b3*exp(-c3*x3)) , 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      # plotting
      if(length(coef_result[,1])==1){
        pred_y = predict(fit12,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==2){
        pred_y = predict(fit22,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==3){
        pred_y = predict(fit32,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","Estimated value"),pch=c(1,-1),lty=c(0,1),col=c(1,2))
      if(save_image){dev.off()}
      
      # return coefficient
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}



#### segmented Logistic for daily graph ####
segLogistic_daily = function(Country,start_date=1,
                       end_date=max_date,
                       break_point=c(),
                       exclusion_criteria=1000,gkf_bandwidth=14,
                       max_iter = 500,
                       save_image = FALSE,
                       prediction_plot = FALSE){
  if(save_image){png(filename=paste0("daily/",i,"_Logistic.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] 
  x = t:end_date - (t-1)
  y_day = df_sum[t:end_date,temp]
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  b = b - b[1]+1
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y_day,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Logistic"), xlim=c(0,n+30),ylim=c(0,2*max(y_day))) 
  }else{
    plot(x,y_day,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Logistic")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  tryCatch(
    expr = {
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-b[1])]
        X_design_ =  X_design[1:(b[2]-b[1]),]
      }
      
      fit11 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[3]-1]
          X_design_ = X_design[1:b[3]-1,] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
        
      }
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[4]-1]
          X_design_ = X_design[1:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) + a3/(1+exp(b3-c3*x3)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) + a2/(1+exp(b2-c2*x2)) + a3/(1+exp(b3-c3*x3)), 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
        
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      if(length(coef_result[,1])==1){
        pred_y = predict(fit12,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==2){
        pred_y = predict(fit22,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==3){
        pred_y = predict(fit32,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","Estimated value"),pch=c(1,-1),lty=c(0,1),col=c(1,2))
      if(save_image){dev.off()}
      
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}


#### segmented Bertalanffy for daily graph ####
segBertalanffy_daily = function(Country,start_date=1,
                             end_date=max_date,
                             break_point=c(),
                             exclusion_criteria=1000,gkf_bandwidth=14,
                             max_iter = 500,
                             save_image = FALSE,
                             prediction_plot = FALSE){
  if(save_image){png(filename=paste0("daily/",i,"_Bertalanffy.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] 
  x = t:end_date - (t-1)
  y_day = df_sum[t:end_date,temp]
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  b = b - b[1]+1
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y_day,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Bertalanffy"), xlim=c(0,n+30),ylim=c(0,2*max(y_day))) 
  }else{
    plot(x,y_day,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Bertalanffy")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  tryCatch(
    expr = {
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-b[1])]
        X_design_ =  X_design[1:(b[2]-b[1]),]
      }

      fit11 = try(nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[3]-1]
          X_design_ = X_design[1:b[3]-1,] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
        
      }
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[4]-1]
          X_design_ = X_design[1:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) + a3*((1-exp(-b3*x3))^c3) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) + a2*((1-exp(-b2*x2))^c2) + a3*((1-exp(-b3*x3))^c3) , 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
        
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      if(length(coef_result[,1])==1){
        pred_y = predict(fit12,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==2){
        pred_y = predict(fit22,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==3){
        pred_y = predict(fit32,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","Estimated value"),pch=c(1,-1),lty=c(0,1),col=c(1,2))
      if(save_image){dev.off()}
      
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}


#### segmented Gompertz for daily graph ####
segGompertz_daily = function(Country,start_date=1,
                                end_date=max_date,
                                break_point=c(),
                                exclusion_criteria=1000,gkf_bandwidth=14,
                                max_iter = 500,
                                save_image = FALSE,
                                prediction_plot = FALSE){
  if(save_image){png(filename=paste0("daily/",i,"_Gompertz.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] 
  x = t:end_date - (t-1)
  y_day = df_sum[t:end_date,temp]
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  b = b - b[1]+1
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y_day,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Gompertz"), xlim=c(0,n+30),ylim=c(0,2*max(y_day))) 
  }else{
    plot(x,y_day,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Gompertz")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  tryCatch(
    expr = {
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-b[1])]
        X_design_ =  X_design[1:(b[2]-b[1]),]
      }
      
      fit11 = try(nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[3]-1]
          X_design_ = X_design[1:b[3]-1,] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
        
      }
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y
          X_design_ = X_design
        }else{
          y_ = y[1:b[4]-1]
          X_design_ = X_design[1:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) + a3*exp(-b3*exp(-c3*x3)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) + a2*exp(-b2*exp(-c2*x2)) + a3*exp(-b3*exp(-c3*x3)) , 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
        
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      if(length(coef_result[,1])==1){
        pred_y = predict(fit12,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==2){
        pred_y = predict(fit22,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])==3){
        pred_y = predict(fit32,newdata=X_design_pred)
        pred_y = c(diff(pred_y)[1],diff(pred_y))
        lines(x_, pred_y, col=2)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","Estimated value"),pch=c(1,-1),lty=c(0,1),col=c(1,2))
      if(save_image){dev.off()}
      
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}



#### separately Logistic : do not preserve continuity ####
segLogistic_separate = function(Country,start_date=1,
                       end_date=max_date,
                       break_point=c(),
                       exclusion_criteria=1000,gkf_bandwidth=14,
                       max_iter = 500,
                       save_image = FALSE,
                       prediction_plot = FALSE,
                       daily_plot = FALSE){
  if(save_image){png(filename=paste0("separate/",i,"_Logistic.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] #min(which(df_sum[,temp]>0))
  x = t:end_date - (t-1)
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Logistic"), xlim=c(0,n+30),ylim=c(0,2*max(y))) 
  }else{
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Logistic")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  b = b - b[1]+1
  tryCatch(
    expr = {
      ## Logistic regression ##
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-1)]
        X_design_ =  X_design[1:(b[2]-1),]
      }
      
      fit11 = try(nls2(y_ ~ a1/(1+exp(b1-c1*x1)) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1/(1+exp(b1-c1*x1)) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y[b[2]:n]
          X_design_ = X_design[b[2]:n,]
        }else{
          y_ = y[b[2]:(b[3]-1)]
          X_design_ = X_design[b[2]:(b[3]-1),] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ y_[1] + a2/(1+exp(b2-c2*x2)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ y_[1] + a2/(1+exp(b2-c2*x2)) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
      } 
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y[b[3]:n]
          X_design_ = X_design[b[3]:n]
        }else{
          y_ = y[b[3]:b[4]-1]
          X_design_ = X_design[b[3]:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ y_[1] + a3/(1+exp(b3-c3*x3)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ y_[1] + a3/(1+exp(b3-c3*x3)), 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      if(length(coef_result[,1])>=1){
        pred_y = predict(fit12,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])>=2){
        pred_y = predict(fit22,newdata=X_design_pred)
        lines(x_, pred_y, col=3)
      }
      if(length(coef_result[,1])>=3){
        pred_y = predict(fit32,newdata=X_design_pred)
        lines(x_, pred_y, col=4)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","fit1","fit2","fit3"),pch=c(1,-1,-1,-1),lty=c(0,1,1,1),col=c(1,2,3,4))
      if(save_image){dev.off()}
      
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}



#### separately Bertalanffy : do not preserve continuity ####
segBertalanffy_separate = function(Country,start_date=1,
                                end_date=max_date,
                                break_point=c(),
                                exclusion_criteria=1000,gkf_bandwidth=14,
                                max_iter = 500,
                                save_image = FALSE,
                                prediction_plot = FALSE,
                                daily_plot = FALSE){
  if(save_image){png(filename=paste0("separate/",i,"_Bertalanffy.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] #min(which(df_sum[,temp]>0))
  x = t:end_date - (t-1)
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Bertalanffy"), xlim=c(0,n+30),ylim=c(0,2*max(y))) 
  }else{
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Bertalanffy")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  b = b - b[1]+1
  tryCatch(
    expr = {
      ## Bertalanffy regression ##
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-1)]
        X_design_ =  X_design[1:(b[2]-1),]
      }
      
      fit11 = try(nls2(y_ ~ a1*((1-exp(-b1*x1))^c1) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1*((1-exp(-b1*x1))^c1) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y[b[2]:n]
          X_design_ = X_design[b[2]:n,]
        }else{
          y_ = y[b[2]:(b[3]-1)]
          X_design_ = X_design[b[2]:(b[3]-1),] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ y_[1] + a2*((1-exp(-b2*x2))^c2) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ y_[1] + a2*((1-exp(-b2*x2))^c2) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
      } 
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y[b[3]:n]
          X_design_ = X_design[b[3]:n]
        }else{
          y_ = y[b[3]:b[4]-1]
          X_design_ = X_design[b[3]:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ y_[1] + a3*((1-exp(-b3*x3))^c3) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ y_[1] + a3*((1-exp(-b3*x3))^c3), 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      if(length(coef_result[,1])>=1){
        pred_y = predict(fit12,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])>=2){
        pred_y = predict(fit22,newdata=X_design_pred)
        lines(x_, pred_y, col=3)
      }
      if(length(coef_result[,1])>=3){
        pred_y = predict(fit32,newdata=X_design_pred)
        lines(x_, pred_y, col=4)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","fit1","fit2","fit3"),pch=c(1,-1,-1,-1),lty=c(0,1,1,1),col=c(1,2,3,4))
      if(save_image){dev.off()}
      
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}



#### separately Gompertz : do not preserve continuity ####
segGompertz_separate = function(Country,start_date=1,
                                   end_date=max_date,
                                   break_point=c(),
                                   exclusion_criteria=1000,gkf_bandwidth=14,
                                   max_iter = 500,
                                   save_image = FALSE,
                                   prediction_plot = FALSE,
                                   daily_plot = FALSE){
  if(save_image){png(filename=paste0("separate/",i,"_Gompertz.png"))}
  
  temp = which(colnames(df_sum)==Country)
  t = break_point[1] #min(which(df_sum[,temp]>0))
  x = t:end_date - (t-1)
  y = cumsum(df_sum[t:end_date,temp])
  if(max(y)<exclusion_criteria) next
  n = length(x)
  
  # design matrix
  m = length(break_point)
  b = break_point
  coef_result = NULL
  
  if(m>1){
    X_design = cbind(x)
    for(j in 1:(m-1)){
      X_design = cbind( X_design,x)
    }
  }else{
    X_design = cbind(x)
  }
  X_design = data.frame(X_design)
  colnames(X_design) = c("x1","x2","x3","x4","x5")[1:m]
  
  x_ = 1:(n+30)
  if(m>1){
    X_design_pred = cbind(x_)
    for(j in 1:(m-1)){
      X_design_pred = cbind( X_design_pred,x_)
    }
  }else{
    X_design_pred = cbind(x_)
  }
  X_design_pred = data.frame(X_design_pred)
  colnames(X_design_pred) = c("x1","x2","x3","x4","x5")[1:m]
  
  if(prediction_plot){
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Gompertz"), xlim=c(0,n+30),ylim=c(0,2*max(y))) 
  }else{
    plot(x,y,xlab=paste("Days since",as.Date(t,origin = "2019-12-31")),
         ylab="Daily cases", main=paste0(Country, " / Gompertz")) 
  }
  
  cat(paste("-----  ", Country,"  -----\n"))
  b = b - b[1]+1
  tryCatch(
    expr = {
      ## Bertalanffy regression ##
      # 1st
      st = data.frame(a1 = c(0.1*max(y), max(y)*2), b1 = c(0, 25), c1 = c(0, 0.5))
      
      if(m==1){
        y_ = y
        X_design_ =  X_design
      }else{
        y_ = y[1:(b[2]-1)]
        X_design_ =  X_design[1:(b[2]-1),]
      }
      
      fit11 = try(nls2(y_ ~ a1*exp(-b1*exp(-c1*x1)) ,
                       data=X_design_, start = st,
                       algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
      initial_list = data.frame(a1=coef(fit11)[1],b1=coef(fit11)[2],c1=coef(fit11)[3])
      as.list(initial_list)
      fit12 = try( nls(y_ ~ a1*exp(-b1*exp(-c1*x1)) , 
                       data=X_design_,start = initial_list) )
      coef_result = data.frame(a=coef(fit12)[1],b=coef(fit12)[2],c=coef(fit12)[3])
      
      a1=coef(fit12)[1]; b1=coef(fit12)[2]; c1=coef(fit12)[3]
      
      #2nd
      if(m>=2){
        if(m==2){
          y_ = y[b[2]:n]
          X_design_ = X_design[b[2]:n,]
        }else{
          y_ = y[b[2]:(b[3]-1)]
          X_design_ = X_design[b[2]:(b[3]-1),] 
        }
        
        st = data.frame(a2 = c(0.1*max(y), max(y)*2), b2 = c(0, 25), c2 = c(0, 0.5)) 
        fit21 = try(nls2(y_ ~ y_[1] + a2*exp(-b2*exp(-c2*x2)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a2=coef(fit21)[1], b2=coef(fit21)[2], c2=coef(fit21)[3])
        as.list(initial_list)
        fit22 = try( nls(y_ ~ y_[1] + a2*exp(-b2*exp(-c2*x2)) , 
                         data=X_design_,start = initial_list) )
        
        a2=coef(fit22)[1]; b2=coef(fit22)[2]; c2=coef(fit22)[3]
        coef_result = rbind(coef_result,coef(fit22))
      } 
      
      #3rd
      if(m>=3){
        if(m==3){
          y_ = y[b[3]:n]
          X_design_ = X_design[b[3]:n]
        }else{
          y_ = y[b[3]:b[4]-1]
          X_design_ = X_design[b[3]:b[4]-1,] 
        }
        
        st = data.frame(a3 = c(0.1*max(y), max(y)*2), b3 = c(0, 25), c3 = c(0, 0.5)) 
        fit31 = try(nls2(y_ ~ y_[1] + a3*exp(-b3*exp(-c3*x3)) ,
                         data=X_design_, start = st,
                         algorithm = "plinear-random",control = nls.control(maxiter = max_iter)))
        initial_list = data.frame(a3=coef(fit31)[1], b3=coef(fit31)[2], c3=coef(fit31)[3])
        as.list(initial_list)
        fit32 = try( nls(y_ ~ y_[1] + a3*exp(-b3*exp(-c3*x3)) , 
                         data=X_design_,start = initial_list) )
        
        a3=coef(fit32)[1]; b3=coef(fit32)[2]; c3=coef(fit32)[3]
        coef_result = rbind(coef_result,coef(fit32))
      }
    },
    error = function(e){
      cat("not fitted! \n")
    },
    finally = {
      if(length(coef_result[,1])>=1){
        pred_y = predict(fit12,newdata=X_design_pred)
        lines(x_, pred_y, col=2)
      }
      if(length(coef_result[,1])>=2){
        pred_y = predict(fit22,newdata=X_design_pred)
        lines(x_, pred_y, col=3)
      }
      if(length(coef_result[,1])>=3){
        pred_y = predict(fit32,newdata=X_design_pred)
        lines(x_, pred_y, col=4)
      }
      abline(v=b[-1],lty=2)
      legend("topright",c("Real value","fit1","fit2","fit3"),pch=c(1,-1,-1,-1),lty=c(0,1,1,1),col=c(1,2,3,4))
      if(save_image){dev.off()}
      
      try(
        expr={
          rownames(coef_result) = c(1,2,3,4,5)[1:length(coef_result[,1])]
          print(coef_result)
          return(coef_result)
        }
      )
    }
  )
}

