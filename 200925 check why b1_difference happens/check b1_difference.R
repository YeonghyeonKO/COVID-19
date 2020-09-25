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

# Analysis Period( ~ 9/22)
max_date <- set_date("2020/9/22") 

#### Data import & Data Preprocessing & Peak Detection Analysis ####
# data import
df <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")

# generate country, population dataframe
country <- as.character(unique(df$countriesAndTerritories))
country <- country[-which(country=="Cases_on_an_international_conveyance_Japan")]

# preprocess data
df_sum = preprocessing_data()



start_data <- as.Date("2019-12-31")
first_case_Costa_Rica <- as.Date("2020-03-07")
first_case_Ukraine <- as.Date("2020-03-05")
first_case_Argentina <- as.Date("2020-03-04")
first_case_Paraguay <- as.Date("2020-03-09")

zero_Costa_Rica <- which(df_sum$Costa_Rica[(first_case_Costa_Rica-start_data+1):length(df_sum$Costa_Rica)] == 0)
zero_Ukraine <- which(df_sum$Ukraine[(first_case_Ukraine-start_data+1):length(df_sum$Ukraine)] == 0)
zero_Argentina <- which(df_sum$Argentina[(first_case_Argentina-start_data+1):length(df_sum$Argentina)] == 0)
zero_Paraguay <- which(df_sum$Paraguay[(first_case_Paraguay-start_data+1):length(df_sum$Paraguay)] == 0)


par(mfrow=c(2,2))
plot(1: (length(df_sum$Costa_Rica) - (first_case_Costa_Rica-start_data)), df_sum$Costa_Rica[(first_case_Costa_Rica-start_data+1):length(df_sum$Costa_Rica)],
     xlab = "Costa_Rica", ylab = "daily cases")
points(zero_Costa_Rica, rep(0, length(zero_Costa_Rica)), col=2)

plot(1: (length(df_sum$Ukraine) - (first_case_Ukraine-start_data)), df_sum$Ukraine[(first_case_Ukraine-start_data+1):length(df_sum$Ukraine)],
     xlab = "Ukraine", ylab = "daily cases")
points(zero_Ukraine, rep(0, length(zero_Ukraine)), col=2)

plot(1: (length(df_sum$Argentina) - (first_case_Argentina-start_data)), df_sum$Argentina[(first_case_Argentina-start_data+1):length(df_sum$Argentina)],
     xlab = "Argentina", ylab = "daily cases")
points(zero_Argentina, rep(0, length(zero_Argentina)), col=2)

plot(1: (length(df_sum$Paraguay) - (first_case_Paraguay-start_data)), df_sum$Paraguay[(first_case_Paraguay-start_data+1):length(df_sum$Paraguay)],
     xlab = "Paraguay", ylab = "daily cases")
points(zero_Paraguay, rep(0, length(zero_Paraguay)), col=2)

title("Details about b1_difference", outer=TRUE, line = -2)
