####Prepare dataframes and produce exposure matrix

#load libraries
library(tidyr)
library(plyr)
library(mgcv)

##run predtk and prettk_2017

#unique dataframe
temp_data<-rbind(temp_2017,temp_2016)
colnames(temp_data)<-c("hour","day","site","orig_DS","orig_US","Downstream","Upstream","type","date.time","TFus","TFds","sp")
head(temp_data)
##long format 
library(tidyr)
temp_data_long<-gather(temp_data,location, temperature,Downstream:Upstream,factor_key = T)
temp_data_long$loc_site<-paste(temp_data_long$location,temp_data_long$site,temp_data_long$sp)
temp_data_long$loc_site<-as.factor(temp_data_long$loc_site)

for (i in 1:length(temp_data_long$hour)) { ##fill evantual NA values 
  temp_data_long$temperature[i]<-ifelse(is.na(temp_data_long$temperature[i]),
                                        temp_data_long$temperature[i-1],
                                        temp_data_long$temperature[i])
}

###############Standardize Temperatures##################

###Z-score ()
Zscore<-by(temp_data_long$temperature,temp_data_long$loc_site,
           function(x)(scale(x)))
temp_data_long$Zscore<-c(Zscore$`Downstream Glenquey si`,Zscore$`Downstream Lintrathen si`,Zscore$`Downstream Rescobie si`,
                         Zscore$`Downstream Morie br`,Zscore$`Downstream Glenquey br`,
                         Zscore$`Upstream Glenquey si`,Zscore$`Upstream Lintrathen si`,Zscore$`Upstream Rescobie si`,
                         Zscore$`Upstream Morie br`,Zscore$`Upstream Glenquey br`)




###########Distance from Trand (TrDi)########

#compute the smooth trend for each sampling location temperatures exposure########
sample<-as.character(unique(temp_data_long$loc_site))
sm.eff_r<-data.frame(NA,NA,NA)
colnames(sm.eff_r)<-c("trend_T","upper_r","lower_r")

for (i in sample) {
  gamname<-paste0("gam",i)
  data.st<-subset.data.frame(temp_data_long,loc_site==i)
  data.st$lin<-seq(1,length(data.st$hour))
  g<-gam(temperature~s(lin,bs="cr",k=length(unique(data.st$day))+10),data=data.st, method="REML") ##k=n° of days +10
  want <- seq(1, nrow(data.st), length.out = length(data.st$hour))
  pdat<-with(data.st,data.frame(lin = lin[want]))
  df.res <- df.residual(g)
  p_r<-predict(g, newdata = pdat, type = "response", se.fit = TRUE)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  pdat_r <- transform(pdat, trend_T = p_r$fit, trend_se = p_r$se.fit)
  pdat_r <- transform(pdat_r,upper_r = trend_T + (crit.t * trend_se),lower_r = trend_T - (crit.t * trend_se))##compute the Confidence interval of the smooth (not used)
  sm_r<-pdat_r[,c(2,4,5)]
  sm.eff_r<-rbind(sm.eff_r,sm_r)

}
temp_data_long<-cbind(temp_data_long,sm.eff_r[-1,])

########Trend Distance (TrDi)#########
temp_data_long$TrDi<-temp_data_long$temperature - temp_data_long$trend_T


###Absoulte temperature rate of change (Atr)#########
Trate<-by(temp_data_long$temperature,temp_data_long$loc_site,
          function(x)(c(0,diff(x))))

temp_data_long$Trate<-c(Trate$`Downstream Glenquey si`,Trate$`Downstream Lintrathen si`,Trate$`Downstream Rescobie si`,
                        Trate$`Downstream Morie br`,Trate$`Downstream Glenquey br`,
                        Trate$`Upstream Glenquey si`,Trate$`Upstream Lintrathen si`,Trate$`Upstream Rescobie si`,
                        Trate$`Upstream Morie br`,Trate$`Upstream Glenquey br`)

############BUILD EXPOSURE MATRIX##########

###load Critical Thermal maxima data###

CT_max<-read.csv("C:/Users/Matteo/OneDrive/Documents/paper_CTmax/csv_scripts/allCT.csv",
                 colClasses = c("character","numeric","factor","numeric","factor",rep("numeric",2)))[c(1:830),]
CT_max$date<-as.POSIXct(CT_max$date, format="%d/%m/%Y")

#assign type (dam/nat) to each sampling location
CT_max$type<-ifelse(CT_max$site%in%c("Glenquey","Lintrathen")&CT_max$location=="Downstream","dam","natural")

##Computing Mean CTmax for each sample
CTmax_ag<-aggregate(CTmax~nsampl*site*location*sp,data=CT_max,FUN = mean)
date<-aggregate(date~nsampl*site*location*sp,data=CT_max,FUN = mean)
CTmax_ag$date<-date[,5]
CT<-CTmax_ag
CT$time<-rep("06:00:00",length(CT$nsampl))##assuming all sample being at 06:00 am
CT$date.time<-as.POSIXct(paste(CT$date,CT$time),format="%Y-%m-%d %T")
rownames(CT)<-seq(1:length(CT$nsampl))

##############number of days in the past over which we explore the delayed effect

d<-28

locations<-paste (CT$location,CT$site,CT$sp)


############exposure matrix of z-score (Zs)###########
Zs_wide<-temp_data_long[,c(7,13,14)]
Zs_wide<-spread(Zs_wide,loc_site,Zscore)
Zs<-matrix(rep(NA,(nrow(CT)*d*24)+25),nrow = 25)
for (i in 1:nrow(CT)) {
  exposure<-rev(Zs_wide[which(Zs_wide$date.time<=CT$date.time[i] & Zs_wide$date.time>=(CT$date.time[i]-(d*24*60*60))),which(colnames(Zs_wide)==locations[i])]) 
  Zs[i,]<-exposure
}

#############exposure matrix of Distance from Trend (TrDi)##########
TrDi_wide<-temp_data_long[,c(7,13,18)]
TrDi_wide<-spread(TrDi_wide,loc_site,TrDi)
TrDi<-matrix(rep(NA,(nrow(CT)*d*24)+25),nrow = 25)
locations<-paste (CT$location,CT$site,CT$sp)
for (i in 1:nrow(CT)) {
  exposure<-rev(TrDi_wide[which(TrDi_wide$date.time<=CT$date.time[i] & TrDi_wide$date.time>=(CT$date.time[i]-(d*24*60*60))),which(colnames(TrDi_wide)==locations[i])]) 
  TrDi[i,]<-exposure
}


###########exposure matrix of Absolute Temperature Rate of Change (Atr)#######
Atr_wide<-temp_data_long[,c(7,13,19)]
Atr_wide<-spread(Atr_wide,loc_site,Trate)
Atr<-matrix(rep(NA,(nrow(CT)*d*24)+25),nrow = 25)
for (i in 1:nrow(CT)) {
  exposure<-rev(Atr_wide[which(Atr_wide$date.time<=CT$date.time[i] & Atr_wide$date.time>=(CT$date.time[i]-(d*24*60*60))),which(colnames(Atr_wide)==locations[i])]) 
  Atr[i,]<-exposure
}




