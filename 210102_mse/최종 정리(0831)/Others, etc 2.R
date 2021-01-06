





colnames(ourworld)
ourworld["stringency_index"]
for(i in country){
  
  
}


# number of test data
ourworld = read.csv("C:/Users/김학용/Desktop/코로나 연구 인턴/참고 문헌 및 기타 자료/our_world_in_data.csv")

US_test = select(ourworld,location,date,new_tests)%>%
  filter(location == "United States")

Korea_test = select(ourworld,location,date,new_tests)%>%
  filter(location == "South Korea")

France_test = select(ourworld,location,date,new_tests)%>%
  filter(location == "France")

par(mfrow=c(1,1))
plot(US_test[,c(2,3)],type='l',main = "Daily New Tests",ylab="",xlab="US")
plot(Korea_test[,c(2,3)],type='l',ylab="",xlab="South Korea")
plot(France_test[,c(2,3)],type='l',ylab="",xlab="France")

acf(France_test$new_tests[140:210])


# incorrectly fitted country
# scatter plot : logistic vs gompertz
par(mfrow=c(1,1))
plot(log(coefficient_result$a1_Logi),log(coefficient_result$a1_Gom))
text(log(coefficient_result$a1_Logi),log(coefficient_result$a1_Gom),
     labels=rownames(coefficient_result),cex=0.5,col='red',pos=2)
plot(log(coefficient_result$a2_Logi),log(coefficient_result$a2_Gom))

par(mfrow=c(2,1))

incorrect_list = c('Albania', 'Zambia', 'Rwanda', "Venezuela")
i = 'Albania'
for(i in "India"){
  temp = which(rownames(break_sum)==i)
  break_point = break_sum$break1[temp]
  for(j in 2:5){
    if(!is.na(break_sum[temp,j])){
      break_point = c(break_point,break_sum[temp,j])
    }
  }
  segLogistic(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = TRUE,
              max_iter = 1000)
  segGompertz_separate(Country=i,break_point = break_point,save_image = FALSE,prediction_plot = TRUE,
                       max_iter = 1000)
}