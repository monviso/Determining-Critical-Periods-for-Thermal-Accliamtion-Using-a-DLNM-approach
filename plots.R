##some plots for data visualization
library(ggplot2)

############## Water temperature 2016
#to explore other sampling location temperature possible values are: 
#"site"= Glenquey, Morie; 
#y and color in ggplot: WUS2+TFus, WDS2+TFds
site<-"Morie"


ggplot(temp_2016[which(temp_2016$site==site),])+
  geom_line(aes(x=date.time,y=WUS2,color=TFus,group=1),size=1)+
  labs(y="temperature (°C)",x="Date - year 2016",color="Data")+
  theme_bw()

############## Water temperature 2017
#to explore other sampling location temperature possible values are: 
#"site"= Glenquey, Lintrathen; Rescobie 
#y and color in ggplot: WUS2+TFus, WDS2+TFds
site<-"Glenquey"


ggplot(temp_2017[which(temp_2017$site==site),])+
  geom_line(aes(x=date.time,y=WDS2,color=TFds,group=1),size=1)+
  labs(y="temperature (°C)",x="Date - year 2017",color="Data")+
  theme_bw()
