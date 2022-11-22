##########Fit DLNM Models#######

########fitting DLNM
library(mgcv)
library(lme4)

CT$loc_site<-as.factor(paste0(CT$location,CT$site))
CT$sp<-as.factor(CT$sp)


### DLNM based on Zscore standardized temperatures

#build crossbasis object 
cb.Zs<-crossbasis(Zs,argvar =list(fun="cr",df=3), arglag=list(fun="cr",df=3,intercept=T))
summary(cb.Zs)
#default penalty matrix
Pen_Zs.d <- cbPen(cb.Zs)

#fitting GAM
gam.Zs<-gam(CTmax~cb.Zs+s(loc_site,bs="re")+ sp
            ,data = CT,method = "GCV.Cp",paraPen=list(cb.Zs=Pen_Zs.d),)
summary(gam.Zs)
AIC(gam.Zs)

#estimation of exposures effects
pred_gam_Zs<-crosspred(cb.Zs,gam.Zs,cumul=F,cen=0,ci.level = 0.95,
                       at=seq(-2,2,0.1))
#default plots
plot(pred_gam_Zs,theta=230,phi=20,ltheta=-80,cumul=F,ylab="lag (hour)",zlab="CTmax variation",xlab="Zs (°C)")
plot(pred_gam_Zs, lag=600,cumul=F, ylab="CTmax variation", xlab="Zs", main="lag= 10")

pred_Zs.at<-crosspred(cb.Zs,gam.Zs,cumul=F,cen=0,ci.level = 0.95, at=2)
plot(pred_Zs.at, var=2,cumul=F, ylab="CTmax variation", xlab="lag (hour)", main="Zs = 2")

###2D heat map plot with filter of all non-significative combiantions of lag-exposure (run "heatplot_function" script)
plot_heat(pred_gam_Zs,Zs)



##DLNM based on Distance from Trend standardized temperatures

#build crossbasis object 
cb.TrDi<-crossbasis(TrDi,argvar =list(fun="cr",df=3), arglag=list(fun="cr",df=3,intercept=T))
summary(cb.TrDi)
#default penalty matrix
Pen_TrDi.d <- cbPen(cb.TrDi) 

#fitting GAM
gam.TrDi<-gam(CTmax~cb.TrDi+s(loc_site,bs="re")+sp
            ,data = CT,method = "GCV.Cp",paraPen=list(cb.TrDi=Pen_TrDi.d))

summary(gam.TrDi)
AIC(gam.TrDi)

#estimation of exposures effects
pred_gam_TrDi<-crosspred(cb.TrDi,gam.TrDi,cumul=F,cen=0,ci.level = 0.95,
                        at=seq(min(TrDi,na.rm = TRUE),max(TrDi,na.rm = TRUE),0.1))

#default plots
plot(pred_gam_TrDi,theta=230,phi=20,ltheta=-80,cumul=F,ylab="lag (hour)",zlab="CTmax variation",xlab="Temperature (°C) -TrDi")
plot(pred_gam_TrDi, lag=300,cumul=F, ylab="CTmax variation", xlab="Temperature (°C) -TrDi", main="lag = 300",
     cex.main=1.75, cex.lab=1.5, cex.axis=1.25)

pred_TrDi.at<-crosspred(cb.TrDi,gam.TrDi,cumul=F,cen=0,ci.level = 0.95, at=-2)
plot(pred_TrDi.at, var=-2,cumul=F, ylab="CTmax variation", xlab="lag (hour)", main="Temperature - TrDi = -2 °C",
     cex.main=1.75, cex.lab=1.5, cex.axis=1.25)

###2D heat map plot with filter of all non-significative combiantions of lag-exposure (run "heatplot_function" script)
plot_heat(pred_gam_TrDi,TrDi)




##DLNM based on Temperatures rate of change

#build crossbasis object 
cb.Atr<-crossbasis(Atr,argvar =list(fun="cr",df=3), arglag=list(fun="cr",df=3,intercept=T))
summary(cb.Atr)
#default penalty matrix
Pen_Atr.d <- cbPen(cb.Atr)

#fitting GAM
gam.Atr<-gam(CTmax~cb.Atr+s(loc_site,bs="re")+sp
            ,data = CT,method = "GCV.Cp",paraPen=list(cb.Atr=Pen_Atr.d))

summary(gam.Atr)
AIC(gam.Atr)

#estimation of exposures effects
pred_gam_Atr<-crosspred(cb.Atr,gam.Atr,cumul=F,cen=0,ci.level = 0.95,
                       at=seq(min(Atr,na.rm = TRUE),max(Atr,na.rm = TRUE),0.1))

#default plots
plot(pred_gam_Atr,theta=230,phi=20,ltheta=-80,cumul=F,ylab="lag (hour)",zlab="CTmax variation",xlab="Temperature (°C)-Atr")
plot(pred_gam_Atr, lag=0,cumul=F, ylab="CTmax variation", xlab="Temperature (°C)-Atr", main="lag= 300",
     cex.main=1.75, cex.lab=1.5, cex.axis=1.25)

pred_Atr.at<-crosspred(cb.Atr,gam.Atr,cumul=F,cen=0,ci.level = 0.95, at=3)
plot(pred_Atr.at, var=3,cumul=F, ylab="CTmax variation", xlab="lag (hour)", main="Temperature (°C) - Atr = +3 °C/hour",
     cex.main=1.75, cex.lab=1.5, cex.axis=1.25)

###2D heat map plot with filter of all non-significant combiantions of lag-exposure (run "heatplot_function" script)
plot_heat(pred_gam_Atr,Atr)

