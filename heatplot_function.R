##function to produce a ggplot with non-significant combinations of lag-exposure removed
#SEE JOURANL ARTICLE FOR THE DEFINITION HERE OF SIGNIFICAITIVE EFFECTS

#x: name of crosspred object (predictions from distributed lag linear (DLMs) and non-linear models (DLNMs))
#y: name of matrix object (exposure matrix)

plot_heat<-function(x,y){
  av <- require("ggplot2")
  if(av==F){plot("install and/or load GGplot")}else{
  m<-c(x[[7]])
mh<-c(x[[12]])
ml<-c(x[[11]])
pred2d<-data.frame(matfit=m,up=mh,low=ml)
pred2d$exp<-rep(seq(min(y,na.rm = TRUE),max(y,na.rm = TRUE),0.1),ncol(y))
pred2d$lag<-rep(1:ncol(y),each=nrow(x[[7]]))
pred2d$fit095<-ifelse(pred2d$up>0 & pred2d$low>0 |pred2d$up<0&pred2d$low<0 , pred2d$matfit,NA)

lg<- bquote(atop(CT[max]~ variation))
ggplot(data=pred2d)+
  geom_tile(aes(x=lag,y=exp,fill=fit095))+
  scale_fill_viridis_c(option="magma")+
  theme_bw()+
  labs(x="lag", y="Exposure", fill=lg)+
  theme(axis.text = element_text(size = 13),axis.title = element_text(size=15),legend.text =element_text(size=12),legend.title = element_text(size=15),
        legend.position = 'top', legend.spacing.x = unit(1.0, 'cm'))}
}


