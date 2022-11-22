###Predict water temperature for each sampling site from air temperature data##
##Air temperature data are from Met Office MIDAS Database (https://catalogue.ceda.ac.uk/uuid/dbd451271eb04662beade68da43546e1)##
##for each sampling site we indicate here the name and code of Met Office sation used

##load packages
library(mgcv)

##data preparation of the recorded temperature in Loch Morie and Glenquey Reservoir 2016
##"temp_16" dataset
temp_data_16<-read.csv("C:/Users/Matteo/OneDrive/Documents/paper_CTmax/csv_scripts/temp_data_16.csv",
                       colClasses = c("numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric"))
tmp_16<-aggregate(cbind(WUS,WDS,AUS,ADS)~date+hour+day+site,data = temp_data_16,FUN=mean)
tmp_16$date.time<-strptime(with(tmp_16, paste(date,day,hour)), "%Y %j %H")
tmp_16$date.time<-as.POSIXct(tmp_16$date.time,format="%Y-%m-%d %T")
tmp_16$hour<-ifelse(tmp_16$hour==24,0,tmp_16$hour)

###LOCH MORIE sampling sites temperatures
tmp_16_m<-tmp_16[which(tmp_16$site=="morie"),]
tmp_16_m$dayf<-as.factor(tmp_16_m$day)

###GLENQUEY RESERVOIR sampling site temperatures
tmp_16_g<-tmp_16[which(tmp_16$site=="glenquey"),]
tmp_16_g$dayf<-as.factor(tmp_16_g$day)


##Tulloch 16idge station (code: 00105) air temperaure data for Loch Morie 2016
##data preparation
tl_bridge<-read.csv("C:/Users/Matteo/OneDrive/Documents/paper_CTmax/csv_scripts/tulloc_bridge_2016.csv",
                    skip=280,header = T)[(1:8784),c(1,36)]
tl_bridge$ob_time<-as.POSIXct(tl_bridge$ob_time,format="%Y-%m-%d %T")
tl_bridge$day<-rep(1:366,each=24)
tl_bridge$hour<-rep(seq(0,23),366)
colnames(tl_bridge)[c(1,2)]<-c("date.time","Air")

#extract the Air temperature data correspondent to sampling location available measurement 
tl_bridge_pr<-tl_bridge[c(which(tl_bridge$date.time==tmp_16_m$date.time[1]):which(tl_bridge$date.time==tmp_16_m$date.time[length(tmp_16_m$date)])),]
tl_bridge_pr<-merge(tl_bridge_pr,tmp_16_m,all=T,by="date.time")
tl_bridge_pr<-tl_bridge_pr[-which(is.na(tl_bridge_pr$date)|is.na(tl_bridge_pr$Air)),]
tl_bridge_pr$hourf<-as.factor(tl_bridge_pr$hour.x)
tl_bridge_pr$seq<-as.numeric(tl_bridge_pr$date.time)

##step1: predict WUS and WDS from Tulloch 16idge sation air temperatures

##create training (75%) and validation(25%) dataframes
##LASTS-FIRSTS method: the latest 75% of the data are used as training, the earliest 25% of the data as validation
tb_train<-tl_bridge_pr[c((round(length(tl_bridge_pr$date.time)*0.25):length(tl_bridge_pr$date.time))),]
tb_val<-tl_bridge_pr[c(1:round(length(tl_bridge_pr$date.time)*0.25)-1),]


###implementation of 10 GAM models for Upstream location
##running time ~10sec on ACER Aspire7 mod.A715-71G-743K,Intel Core i7-7700HQ 2.8 GHz, 16GB RAM 

##models to be implemented
a<-c('s(Air,k=30)',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)',
     'te(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,5,5))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,bs="re")+day.y',
     'te(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(day.y,bs="re")+day.y',
     'ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(Air,k=30)+s(hour.y,k=24)+s(day.y,bs="re")+day.y')


loc<-c("WUS","WDS")
WUS<-data.frame(model=a,MAE=rep(NA,length(a)),AEsd=rep(NA,length(a)),r=rep(NA,length(a)))
WDS<-data.frame(model=a,MAE=rep(NA,length(a)),AEsd=rep(NA,length(a)),r=rep(NA,length(a)))

system.time(
  for (j in 1:length(loc)) {
    show(loc[j])
    b<-paste0(paste(loc[j],"~"),a)
    for (i in 1:length(a)) {
      
      show(a[i])
      af<-as.formula(b[i]) 
      m1<-gam(af,data=tb_train)
      
      tb_val$p1<-c(predict.gam(m1,tb_val,type="response"))
      tb_val$diff1<-tb_val$p1-tb_val[,which(colnames(tb_val)==loc[j])]
      MAE<-mean(tb_val$diff1);AEsd<-sd(tb_val$diff1);r<-cor(tb_val[,which(colnames(tb_val)==loc[j])],tb_val$p1)
      if(loc[j]=="WUS"){WUS[i,c(2:4)]<-c(MAE,AEsd,r)}else{WDS[i,c(2:4)]<-c(MAE,AEsd,r)}
    }  
  })


##best model for upstream and downstream
##the selection here is based on the sum of the ranking position of each accuracy metric (MAE, AEsd, r)
##the lowest sum is the best model
##IT IS POSSIBLE TO EXPLORE WUS AND WDS DATA FRAMES TO SEE THE RESULT OF EACH MODEL

mWUS<-paste("WUS~",WUS$model[which.min(rank(round(abs(WUS$MAE),1))+rank(round(WUS$AEsd,1))+rank(ifelse(WUS$r>0,round(-WUS$r,2),round(WUS$r,2))))])
mWDS<-paste("WDS~",WDS$model[which.min(rank(round(abs(WDS$MAE),1))+rank(round(WDS$AEsd,1))+rank(ifelse(WDS$r>0,round(-WDS$r,2),round(WDS$r,2))))])
show(mWUS)
show(mWDS)

##retrain the two best models with the whole available data


mWUS2<-gam(as.formula(mWUS),data = tl_bridge_pr)
mWDS2<-gam(as.formula(mWDS),data = tl_bridge_pr)

#model checking (change model name for checking upstream/downstream models)
summary (mWUS2)
plot(mWUS2)
gam.check(mWUS2)

##predict missing data from 01/06/2016 to 28/07/2016
colnames(tl_bridge)[c(3,4)]<-c("day.y","hour.y")
morie<-tl_bridge[c(which(tl_bridge$date.time=="2016-06-01 06:00:00"):which(tl_bridge$date.time=="2016-08-03 06:00:00")),]

morie$WUSp<-predict.gam(mWUS2,morie,type = "response")
morie$WDSp<-predict.gam(mWDS2,morie,type = "response")

morie<-merge(morie,tmp_16_m,by="date.time",all=T)
morie$WUS2<-ifelse(is.na(morie$WUS),morie$WUSp,morie$WUS)
morie$WDS2<-ifelse(is.na(morie$WDS),morie$WDSp,morie$WDS)
morie$TFus<-ifelse(is.na(morie$WUS),"predicted","measured")
morie$TFds<-ifelse(is.na(morie$WDS),"predicted","measured")

##########Glenquey reserevoir 2017##########
##Strathallan airfield station (code 00212) air temperaure data for Glenquey reservoir 2016 
str_air<-read.csv("C:/Users/Matteo/OneDrive/Documents/paper_CTmax/csv_scripts/strath_air_2016.csv",
                  skip=280,header = T)[(1:8784),c(1,36)]
str_air$ob_time<-as.POSIXct(str_air$ob_time,format="%Y-%m-%d %T")
str_air$day<-rep(1:366,each=24)
str_air$hour<-rep(seq(0,23),366)
colnames(str_air)[c(1,2)]<-c("date.time","Air")

#extract the Air temperature data correspondent to sampling location available measurement
str_air_pr<-str_air[c(which(str_air$date.time==tmp_16_g$date.time[1]):which(str_air$date.time==tmp_16_g$date.time[length(tmp_16_g$date)])),]
str_air_pr<-merge(str_air_pr,tmp_16_g,all=T,by="date.time")
str_air_pr<-str_air_pr[-which(is.na(str_air_pr$date)|is.na(str_air_pr$Air)),]
str_air_pr$hourf<-as.factor(str_air_pr$hour.x)
str_air_pr$seq<-as.numeric(str_air_pr$date.time)

##create training (75%) and validation(25%) data frames
##LASTS-FIRSTS method: the latest 75% of the data are used as training, the earliest 25% of the data as validation
sa_train<-str_air_pr[c((round(length(str_air_pr$date.time)*0.25):length(str_air_pr$date.time))),]
sa_val<-str_air_pr[c(1:round(length(str_air_pr$date.time)*0.25)-1),]


###implementation of 10 GAM models for Upstream location
##running time ~5sec on ACER Aspire7 mod.A715-71G-743K,Intel Core i7-7700HQ 2.8 GHz, 16GB RAM 

##models to be implemented
a<-c('s(Air,k=30)',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)',
     'te(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,5,5))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,bs="re")+day.y',
     'te(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(day.y,bs="re")+day.y',
     'ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(Air,k=20)+s(hour.y,k=24)+s(day.y,bs="re")+day.y')


loc<-c("WUS","WDS")
WUS<-data.frame(model=a,MAE=rep(NA,length(a)),AEsd=rep(NA,length(a)),r=rep(NA,length(a)))
WDS<-data.frame(model=a,MAE=rep(NA,length(a)),AEsd=rep(NA,length(a)),r=rep(NA,length(a)))

system.time(
  for (j in 1:length(loc)) {
    show(loc[j])
    b<-paste0(paste(loc[j],"~"),a)
    for (i in 1:length(a)) {
      
      show(a[i])
      af<-as.formula(b[i]) 
      m1<-gam(af,data=sa_train)
      
      sa_val$p1<-c(predict.gam(m1,sa_val,type="response"))
      sa_val$diff1<-sa_val$p1-sa_val[,which(colnames(sa_val)==loc[j])]
      MAE<-mean(sa_val$diff1);AEsd<-sd(sa_val$diff1);r<-cor(sa_val[,which(colnames(sa_val)==loc[j])],sa_val$p1)
      if(loc[j]=="WUS"){WUS[i,c(2:4)]<-c(MAE,AEsd,r)}else{WDS[i,c(2:4)]<-c(MAE,AEsd,r)}
    }  
  })


##best model for upstream and downstream
##the selection here is based on the sum of the ranking position of each accuracy metric (MAE, AEsd, r)
##the lowest sum is the best model
##IT IS POSSIBLE TO EXPLORE WUS AND WDS DATA FRAME TO SEE THE RESULT OF EACH MODEL
mWUS<-paste("WUS~",WUS$model[which.min(rank(round(abs(WUS$MAE),1))+rank(round(WUS$AEsd,2))+rank(ifelse(WUS$r>0,round(-WUS$r,2),round(WUS$r,2))))])
mWDS<-paste("WDS~",WDS$model[which.min(rank(round(abs(WDS$MAE),1))+rank(round(WDS$AEsd,2))+rank(ifelse(WDS$r>0,round(-WDS$r,2),round(WDS$r,2))))])
show(mWUS)
show(mWDS)


##retrain the two best models with the whole available data
mWUS2<-gam(as.formula(mWUS),data = str_air_pr)
mWDS2<-gam(as.formula(mWDS),data = str_air_pr)

##predict missing data from 01/06/2016 to 28/07/2016
colnames(str_air)[c(3,4)]<-c("day.y","hour.y")
glenquey16<-str_air[c(which(str_air$date.time=="2016-06-01 06:00:00"):which(str_air$date.time=="2016-08-03 06:00:00")),]

glenquey16$WUSp<-predict.gam(mWUS2,glenquey16,type = "response")
glenquey16$WDSp<-predict.gam(mWDS2,glenquey16,type = "response")

glenquey16<-merge(glenquey16,tmp_16_g,by="date.time",all=T)
glenquey16$WUS2<-ifelse(is.na(glenquey16$WUS),glenquey16$WUSp,glenquey16$WUS)
glenquey16$WDS2<-ifelse(is.na(glenquey16$WDS),glenquey16$WDSp,glenquey16$WDS)
glenquey16$TFus<-ifelse(is.na(glenquey16$WUS),"predicted","measured")
glenquey16$TFds<-ifelse(is.na(glenquey16$WDS),"predicted","measured")



##########merge the two data sets#################
glenquey16$site<-rep("Glenquey",length(glenquey16$date.time))
glenquey16$type<-rep("dam",length(glenquey16$date.time))
glenquey16$sp<-rep("16",length(glenquey16$date.time))

morie$site<-rep("Morie",length(morie$date.time))
morie$type<-rep("natural",length(morie$date.time))
morie$sp<-rep("16",length(morie$date.time))

temp_2016<-rbind(morie,glenquey16)
temp_2016<-temp_2016[,c(4,3,10,12,11,17,16,20,1,18,19,21)]
colnames(temp_2016)[c(1,2,9)]<-c("hour","day","date.time")


