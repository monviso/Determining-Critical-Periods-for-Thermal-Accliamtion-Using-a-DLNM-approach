###Predict water temperature for each sampling site from air temperature data##
##Air temperature data are from Met Office MIDAS Database (https://catalogue.ceda.ac.uk/uuid/dbd451271eb04662beade68da43546e1)##
##for each sampling site we indicate here the name and code of Met Office sation used

##load packages
library(mgcv)

##data preparation of the recorded temperature in Loch Rescobie, Loch of Lintrathen and Glenquey Reservoir 2017
temp_data_17<-read.csv("temp_data_17.csv",
                    colClasses = c("NULL","numeric","numeric","factor","numeric","numeric","factor","POSIXct"))

g<-temp_data_17[which(temp_data_17$site=="Glenquey"),]
l<-temp_data_17[which(temp_data_17$site=="Lintrathen"),]
re<-temp_data_17[which(temp_data_17$site=="Rescobie"),]



############GELNQUEY RESERVOIR 2017#####################
##Strathallan airfield station (code 00212) air temperaure data for Glenquey reservoir 2017
str_air<-read.csv("strath_air_2017.csv",
             skip=280,header = T)[c(1:8754),c(1,36)]
str_air$ob_time<-as.POSIXct(str_air$ob_time,format="%Y-%m-%d %T")
str_air_m<-str_air[which(str_air$ob_time>=min(g$progr) & str_air$ob_time<=max(g$progr)),]
colnames(str_air_m)<-c("progr","air")
gm<-merge(str_air_m,g,by="progr",all=T)             
gm$seq<-as.numeric(gm$progr)
colnames(gm)<-c("date.time","Air","hour.y","day.y","site","WDS","WUS","type","seq")

##create training (75%) and validation(25%) data frames
##LASTS-FIRSTS method: the latest 75% of the data are used as training, the earliest 25% of the data as validation
gm_train<-gm[c((round(length(gm$date.time)*0.25):length(gm$date.time))),]
gm_val<-gm[c(1:round(length(gm$date.time)*0.25)-1),]

###implementation of 10 GAM models for Upstream location
##running time ~16sec on ACER Aspire7 mod.A715-71G-743K,Intel Core i7-7700HQ 2.8 GHz, 16GB RAM 

##models to be implemented
a<-c('s(Air,k=30)',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)',
     #'s(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,k=5)',
     'te(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,5,5))',
     #'s(Air,k=20)+s(hour.y,bs="cc",k=15)+s(day.y,k=5)+ti(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,7,4))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,bs="re")+day.y',
     'te(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(day.y,bs="re")+day.y',
     'ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(Air,k=30)+s(hour.y,k=24)+s(day.y,bs="re")+day.y')
     #'s(hour.y,bs="cc",k=24)+s(day.y,k=10)+ti(day.y,hour.y,bs=c("tp","cc"),k=c(5,24))')


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
      m1<-gam(af,data=gm_train)
      
      gm_val$p1<-c(predict.gam(m1,gm_val,type="response"))
      gm_val$diff1<-gm_val$p1-gm_val[,which(colnames(gm_val)==loc[j])]
      MAE<-mean(gm_val$diff1);AEsd<-sd(gm_val$diff1);r<-cor(gm_val[,which(colnames(gm_val)==loc[j])],gm_val$p1)
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
mWUS2<-gam(as.formula(mWUS),data = gm)
mWDS2<-gam(as.formula(mWDS),data = gm)

##predict missing data from 01/06/2016 to 21/07/2016
p<-data.frame(progr=seq(as.POSIXct("2017-06-01 01:00:00"),as.POSIXct("2017-07-21 00:00:00"),by="hour"),
              hour=rep(seq(1,24),50),day=rep(152:201,each=24))
str_air_short<-str_air[which(str_air$ob_time>=min(p$progr) & str_air$ob_time<=max(p$progr)),]
colnames(str_air_short)<-c("progr","air")

pm<-merge(p,str_air_short,by="progr",all=T)
colnames(pm)<-c("date.time","hour.y","day.y","Air")

pm$p_us<-predict.gam(mWUS2,pm,type = "response")
pm$p_ds<-predict.gam(mWDS2,pm,type = "response")

d<-merge(gm,pm,by="date.time",all=T)
glenquey17<-data.frame(hour=d$hour.y.y,day=d$day.y.y,site=rep("Glenquey",length(d$date.time)),WDS=d$WDS,WUS=d$WUS,
               WDS2=ifelse(is.na(d$WDS),d$p_ds,d$WDS),WUS2=ifelse(is.na(d$WUS),d$p_us,d$WUS),
               type=rep("dam",length(d$date.time)),date.time=d$date.time)

glenquey17$TFus<-ifelse(is.na(glenquey17$WUS),"predicted","measured")
glenquey17$TFds<-ifelse(is.na(glenquey17$WDS),"predicted","measured")


################ LOCH RESCOBIE##########################
##Cairnwell station  (code 00145) air temperature data for Loch Rescobie  2017
cairn<-read.csv("cairnwell_2017.csv",
         skip=280,header = T)[c(1:8016),c(1,36)]
cairn$ob_time<-as.POSIXct(cairn$ob_time,format="%Y-%m-%d %T")

cairn_m<-cairn[which(cairn$ob_time>=min(re$progr) & cairn$ob_time<=max(re$progr)),]
colnames(cairn_m)<-c("progr","air")
rm<-merge(cairn_m,re,by="progr",all=T)             
rm$seq<-as.numeric(rm$progr)
colnames(rm)<-c("date.time","Air","hour.y","day.y","site","WDS","WUS","type","seq")

##create training (75%) and validation(25%) data frames
##LASTS-FIRSTS method: the latest 75% of the data are used as training, the earliest 25% of the data as validation
rm_train<-rm[c((round(length(rm$date.time)*0.25):length(rm$date.time))),]
rm_val<-rm[c(1:round(length(rm$date.time)*0.25)-1),]

###implementation of 10 GAM models for Upstream location
##running time ~16sec on ACER Aspire7 mod.A715-71G-743K,Intel Core i7-7700HQ 2.8 GHz, 16GB RAM 

##models to be implemented
a<-c('s(Air,k=30)',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)',
     #'s(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,k=5)',
     'te(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,5,5))',
     #'s(Air,k=20)+s(hour.y,bs="cc",k=15)+s(day.y,k=5)+ti(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,7,4))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,bs="re")+day.y',
     'te(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(day.y,bs="re")+day.y',
     'ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(Air,k=30)+s(hour.y,k=24)+s(day.y,bs="re")+day.y')
    #'s(hour.y,bs="cc",k=24)+s(day.y,k=10)+ti(day.y,hour.y,bs=c("tp","cc"),k=c(5,24))')



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
      m1<-gam(af,data=rm_train)
      
      rm_val$p1<-c(predict.gam(m1,rm_val,type="response"))
      rm_val$diff1<-rm_val$p1-rm_val[,which(colnames(rm_val)==loc[j])]
      MAE<-mean(rm_val$diff1);AEsd<-sd(rm_val$diff1);r<-cor(rm_val[,which(colnames(rm_val)==loc[j])],rm_val$p1)
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
mWUS2<-gam(as.formula(mWUS),data = rm)
mWDS2<-gam(as.formula(mWDS),data = rm)

##predict missing data from 01/06/2016 to 20/07/2016
p<-data.frame(progr=seq(as.POSIXct("2017-06-01 01:00:00"),as.POSIXct("2017-07-20 00:00:00"),by="hour"),
              hour=rep(seq(1,24),49),day=rep(152:200,each=24))
cairn_short<-cairn[which(cairn$ob_time>=min(p$progr) & cairn$ob_time<=max(p$progr)),]
colnames(cairn_short)<-c("progr","air")

pm<-merge(p,cairn_short,by="progr",all=T)
colnames(pm)<-c("date.time","hour.y","day.y","Air")

pm$p_us<-predict.gam(mWUS2,pm,type = "response")
pm$p_ds<-predict.gam(mWDS2,pm,type = "response")

d<-merge(rm,pm,by="date.time",all=T)
rescobie<-data.frame(hour=d$hour.y.y,day=d$day.y.y,site=rep("Rescobie",length(d$date.time)),WDS=d$WDS,WUS=d$WUS,
                       WDS2=ifelse(is.na(d$WDS),d$p_ds,d$WDS),WUS2=ifelse(is.na(d$WUS),d$p_us,d$WUS),
                       type=rep("dam",length(d$date.time)),date.time=d$date.time)

rescobie$TFus<-ifelse(is.na(rescobie$WUS),"predicted","measured")
rescobie$TFds<-ifelse(is.na(rescobie$WDS),"predicted","measured")



################ LOCH OF LINTRATHEN ##########################

##Cairnwell station  (code 00145) air temperature data for Loch of lintrethen  2017
cairn<-read.csv("cairnwell_2017.csv",
                skip=280,header = T)[c(1:8016),c(1,36)]
cairn$ob_time<-as.POSIXct(cairn$ob_time,format="%Y-%m-%d %T")

cairn_m<-cairn[which(cairn$ob_time>=min(l$progr) & cairn$ob_time<=max(l$progr)),]
colnames(cairn_m)<-c("progr","air")
lm<-merge(cairn_m,l,by="progr",all=T)             
lm$seq<-as.numeric(lm$progr)
colnames(lm)<-c("date.time","Air","hour.y","day.y","site","WDS","WUS","type","seq")
lm<-lm[complete.cases(lm),]
##create training (75%) and validation(25%) data frames
##LASTS-FIRSTS method: the latest 75% of the data are used as training, the earliest 25% of the data as validation
lm_train<-lm[c((round(length(lm$date.time)*0.25):length(lm$date.time))),]
lm_val<-lm[c(1:round(length(lm$date.time)*0.25)-1),]

###implementation of 10 GAM models for Upstream location
##running time ~15sec on ACER Aspire7 mod.A715-71G-743K,Intel Core i7-7700HQ 2.8 GHz, 16GB RAM 

##models to be implemented
a<-c('s(Air,k=30)',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)',
     #'s(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,k=5)',
     'te(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,5,5))',
     #'s(Air,k=20)+s(hour.y,bs="cc",k=15)+s(day.y,k=5)+ti(Air,hour.y,day.y,bs=c("tp","cc","tp"),k=c(8,7,4))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))',
     's(Air,k=30)+s(hour.y,bs="cc",k=24)+s(day.y,bs="re")+day.y',
     'te(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(day.y,bs="re")+day.y',
     'ti(Air,hour.y,bs=c("tp","cc"),k=c(10,8))+s(Air,k=30)+s(hour.y,k=24)+s(day.y,bs="re")+day.y')
     #'s(hour.y,bs="cc",k=24)+s(day.y,k=10)+ti(day.y,hour.y,bs=c("tp","cc"),k=c(5,24))')


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
      m1<-gam(af,data=lm_train)
      
      lm_val$p1<-c(predict.gam(m1,lm_val,type="response"))
      lm_val$diff1<-lm_val$p1-lm_val[,which(colnames(lm_val)==loc[j])]
      MAE<-mean(lm_val$diff1);AEsd<-sd(lm_val$diff1);r<-cor(lm_val[,which(colnames(lm_val)==loc[j])],lm_val$p1)
      if(loc[j]=="WUS"){WUS[i,c(2:4)]<-c(MAE,AEsd,r)}else{WDS[i,c(2:4)]<-c(MAE,AEsd,r)}
    }  
  })

##best model for upstream and downstream
##the selection here is based on the sum of the ranking position of each accuracy metric (MAE, AEsd, r)
##the lowest sum is the best model
##IT IS POSSIBLE TO EXPLORE WUS AND WDS DATA FRAME TO SEE THE RESULT OF EACH MODEL
mWUS<-paste("WUS~",WUS$model[which.min(rank(abs(WUS$MAE))+rank(WUS$AEsd)+rank(-(WUS$r)))])
mWDS<-paste("WDS~",WDS$model[which.min(rank(abs(WDS$MAE))+rank(WDS$AEsd)+rank(-(WDS$r)))])
show(mWUS)
show(mWDS)


##retrain the two best models with the whole available data
mWUS2<-gam(as.formula(mWUS),data = lm)
mWDS2<-gam(as.formula(mWDS),data = lm)

##predict missing data from 01/06/2016 to 19/07/2016
p<-data.frame(progr=seq(as.POSIXct("2017-06-01 01:00:00"),as.POSIXct("2017-07-19 00:00:00"),by="hour"),
              hour=rep(seq(1,24),48),day=rep(152:199,each=24))
cairn_short<-cairn[which(cairn$ob_time>=min(p$progr) & cairn$ob_time<=max(p$progr)),]
colnames(cairn_short)<-c("progr","air")

pm<-merge(p,cairn_short,by="progr",all=T)
colnames(pm)<-c("date.time","hour.y","day.y","Air")

pm$p_us<-predict.gam(mWUS2,pm,type = "response")
pm$p_ds<-predict.gam(mWDS2,pm,type = "response")

d<-merge(lm,pm,by="date.time",all=T)
lintrathen<-data.frame(hour=d$hour.y.y,day=d$day.y.y,site=rep("Lintrathen",length(d$date.time)),WDS=d$WDS,WUS=d$WUS,
                     WDS2=ifelse(is.na(d$WDS),d$p_ds,d$WDS),WUS2=ifelse(is.na(d$WUS),d$p_us,d$WUS),
                     type=rep("dam",length(d$date.time)),date.time=d$date.time)

lintrathen$TFus<-ifelse(is.na(lintrathen$WUS),"predicted","measured")
lintrathen$TFds<-ifelse(is.na(lintrathen$WDS),"predicted","measured")


###########merge all datasets##############

temp_2017<-rbind(glenquey17,lintrathen,rescobie)
temp_2017$sp<-rep("si",length(temp_2017$hour))
