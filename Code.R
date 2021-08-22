##########################################################
###########################Create world map###Fig.1#######
install.packages('tidyverse')
library(tidyverse)
install.packages('maps')
library(maps)
world <- map_data("world")
#Input data
Book1<-read.csv('BAAD.csv',header=TRUE)
Book2<-read.csv('ICP.csv',header=TRUE)
Book3<-read.csv('LAI.csv',header=TRUE)
Book3<-subset(Book3,LAI_observed!='NA')
ggplot() +
  geom_map(
    data = world, map = world,
    aes(x=Longitude, y=Latitude, map_id = region),
    color = "grey", fill = "darkgrey", size = 0.5
  ) +
  geom_point(
    data = Book1,
    aes(x=Longitude, y=Latitude),color='blue',size=0.5,
    alpha = 0.7
  )+  
  geom_point(
    data = Book2,
    aes(x=Longitude, y=Latitude),color='green', size=0.5,
    alpha = 0.7
  )+
  geom_point(
    data = Book3,
    aes(x=Longitude, y=Latitude),color='yellow', size=0.5,
    alpha = 0.7
  )+
  scale_colour_manual(values=c("blue","green",'yellow'))+
  theme_void() +####delete x and y axises
  theme(legend.position = c(0.15, 0.4),
        legend.box.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 10, colour = "black")
  )+
  guides(colour = guide_legend(override.aes = list(size = 3)))

############################################################################
###########################Analyses related to Hm###########################
#########Get data for all species
Input_Hm<-read.csv('Ac_D_H_Total_data.csv',header=TRUE)
Data_input1<-subset(Input,D!='NA'&H!='NA')
library(nlme)
datasetG <- groupedData(H~1|Site, Data_input1)
initVals <- getInitial(H~SSasympOrig(D, m, log(c)), data=datasetG)
X1<-initVals[1]
X2<-exp(initVals[2])
Modelzhou<-nlme(H~m*(1-exp(-c*D)),datasetG,fixed=list(m~1,c~1),random=pdDiag(list(m~1,c~1)),method='ML',start=c(m=X1,c=X2))
result <- data.frame(matrix(nrow =979, ncol = 3))
colnames(result) <- c("Number",'m','c')
Results_all<-cbind(result,coef(Modelzhou))
#########Get data for evergreen and deciduous/angiosperm and gymnosperm
Data_input2<-subset(Input_Hm,D!='NA'&H!='NA'&Type=='evergreen/deciduous/angiosperm/gymnosperm')###Substitute the evergreen to deciduous or angiosperm or gymnosperm
library(nlme)
datasetG <- groupedData(H~1|Site, Data_input2)
initVals <- getInitial(H~SSasympOrig(D, m, log(c)), data=datasetG)
X1<-initVals[1]
X2<-exp(initVals[2])
Modelzhou<-nlme(H~m*(1-exp(-c*D)),datasetG,fixed=list(m~1,c~1),random=pdDiag(list(m~1,c~1)),method='ML',start=c(m=X1,c=X2))
result <- data.frame(matrix(nrow =979, ncol = 3))
colnames(result) <- c("Number",'m','c')
Results_PFTs<-cbind(result,coef(Modelzhou))
###########Environmental effects on Hm or b#all species
Data_input_Hm<-subset(Results_all,Hm!='NA'&MI!='NA')
Data_input_b<-subset(Results_all,b!='NA'&MI!='NA')
##############OLS analyses
Model_b<-lm(b~log(MI),Data_input_b)
Model_Hm<-lm(log(Hm)~log(MI),Data_input_Hm)
summary(Model_b)
summary(Model_Hm)
###############OLS plot for Hm
visreg(fit=Model_Hm,xvar="MI",band=TRUE,rug=2,xtrans=log,jitter=TRUE,points=list(cex=1,col='#293352'),line=list(cex=0.5,col="red"),ylab=expression(bold('ln H '[m]*' (m)')),xlab=expression(bold('ln MI')))
text(x=-0.8,y=4.1,labels='(a)',font=2,cex=1.2)
###############Robust analyses
install.packages('mblm')
library(mblm)
Data_input_Hm$log_Hm<-log(Data_input_Hm$Hm)
Data_input_b$log_MI<-log(Data_input_b$MI)
fit_Hm<-mblm(log_Hm~log_MI,data=Data_input_Hm,repeated=FALSE)
summary(fit_Hm)
fit_Hm<-mblm(b~log_MI,data=Data_input_b,repeated=FALSE)
summary(fit_b)
################robust plot for Hm and b respectively
effects_fit1 <- effects::effect(term="log_MI", mod=fit_Hm)
summary(effects_fit1)
x_fit1 <- as.data.frame(effects_fit1)
ggplot(data=Data_input_Hm, aes(x=log(MI), y=log(Hm)))+geom_point(pch=16,size=2,col='black')+
  geom_smooth(method='lm',formula=y~x,se=TRUE,color='black',linetype=1)+
  geom_line(data=x_fit1, aes(x=log_MI, y=fit), color="blue") +
  geom_line(data=x_fit1, aes(x=log_MI, y=lower), color="blue",linetype='dashed')+
  geom_line(data=x_fit1, aes(x=log_MI, y=upper), color="blue",linetype='dashed')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(bold('ln MI')),y=expression(bold('ln H'[m]*'(m)')))+
  annotate ('text',x=-0.8, y=4.2,cex=5,label=expression(bold("(a)")))
ggplot(data=Data_input_b, aes(x=log(MI), y=log(b)))+geom_point(pch=16,size=2,col='black')+
  geom_smooth(method='lm',formula=y~x,se=TRUE,color='black',linetype=1)+
  geom_line(data=x_fit1, aes(x=log_MI, y=fit), color="blue") +
  geom_line(data=x_fit1, aes(x=log_MI, y=lower), color="blue",linetype='dashed')+
  geom_line(data=x_fit1, aes(x=log_MI, y=upper), color="blue",linetype='dashed')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(bold('ln MI')),y=expression(bold('ln b (m'^-1*')')))+
  annotate ('text',x=-0.78, y=1.59,cex=5,label=expression(bold("(b)")))
##################Histograms for all PFTs of fitted b values
Model<-read.csv('Book2.csv',header=TRUE)
Model$PFTs <- factor(Model$PFTs, levels = c('All species','Evergreen','Deciduous','Angiosperm','Gymnosperm'))
ggplot(data=Model, aes(x=PFTs, y=value,fill=PFTs)) +
  geom_bar(stat="identity",width=0.6)+
  geom_errorbar(data=Model,aes(ymin=value-se, ymax=value+se), width=.2,
                position=position_dodge(.9))+
  labs(x='',y=expression(bold('Fitted b')))+
  scale_fill_manual(values=c('steelblue3','steelblue3','steelblue3','steelblue3','steelblue3'))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate ('text',x=0.7, y=4.6,cex=5,label=expression(bold("(c)")))
####################H distribution for BAAD and ICP
Data<-read.csv('Book3.csv',header=TRUE)
Data_input<-subset(Data,H!='NA')
ggplot(Data_input, aes(H, fill = Database)) + geom_histogram(alpha = 0.5,bins=30,position = 'identity')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y='Count',x='Height')
##################################################################################
#####################Ac/DH related analyses########################################
Input_Ac<-read.csv('Ac_D_H_Total_data.csv',header=TRUE)
Data_input_Ac<-subset(Input_Ac,R!='NA'&MI!='NA')
Model<-lmer(log(Ac)~log(MI)+(1|species),Data_input_Ac)
summary(Model)
######################Separate environmental effects on Ac,D and H
ggplot()+
  geom_smooth(data=Data_input_Ac,aes(x=log(MI), y=log(Ac)),method='lm',formula=y~x,se=TRUE,color='darkgreen',linetype=1)+
  geom_smooth(data=Data_input_Ac,aes(x=log(MI), y=log(H)),method='lm',formula=y~x,se=TRUE,color='red',linetype=1)+
  geom_smooth(data=Data_input_Ac,aes(x=log(MI), y=log(D)),method='lm',formula=y~x,se=TRUE,color='blue',linetype=1)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y='',x='ln MI')+
  annotate ('text',x=-0.7, y=3.3,cex=5,label=expression(bold("(a)")))
########################Variation between PFTs##Ac/DH
summary(aov(log(R)~as.factor(Type2),Data_input_Ac))
summary(Model)
install.packages('gplots')
library(gplots)
plotmeans(log(R)~Type2,xlab='',ylab='ln Ac/DH',lwd=1.5,n.label=FALSE,barwidth=2,Data_input)
r.squaredGLMM(Model)
##################Environmental changes on Ac/DH
Model_Ac<-lm(log(R)~log(MI),Data_input_Ac)
visreg(fit=Model_Ac,xvar="MI",band=TRUE,rug=2,xtrans=log,jitter=TRUE,points=list(cex=1,col='#E7B800'),line=list(cex=0.5,col="#293352"),ylab=expression(bold('ln Ac/DH')),xlab=expression(bold('ln MI')))###angiosperm
text(x=-0.15,y=5.2,labels='(a)',font=2,cex=1.2)
visreg(fit=Model_Ac,xvar="MI",band=TRUE,rug=2,xtrans=log,jitter=TRUE,points=list(cex=1,col='darkblue'),line=list(cex=0.5,col="#293352"),ylab=expression(bold('ln Ac/DH')),xlab=expression(bold('ln MI')))###gymnosperm
text(x=-0.3,y=5.5,labels='(b)',font=2,cex=1.2)
###################robust analyses for Ac/DH
Input_Ac<-read.csv('Ac_D_H_Total_data.csv',header=TRUE)
Data_input_Ac<-subset(Input_Ac,R!='NA'&MI!='NA')
Data_input_Ac$log_MI<-log(Data_input_Ac$MI)
Data_input_Ac$log_R<-log(Data_input_Ac$R)
fit_Ac<-mblm(log_R~log_MI,data=Data_input_Ac,repeated=FALSE)
summary(fit_Ac)
confint.mblm(fit,level=0.95)
######################robust plot for Ac/DH
effects_fit1 <- effects::effect(term="log_MI", mod=fit_Ac)
summary(effects_fit1)
x_fit1 <- as.data.frame(effects_fit1)
ggplot(data=Data_input_Ac, aes(x=log(MI), y=log(R)))+geom_point(pch=16,size=2)+
  geom_smooth(method='lm',formula=y~x,se=TRUE,color='black',linetype=1)+
  geom_line(data=x_fit1, aes(x=log_MI, y=fit), color="blue") +
  geom_line(data=x_fit1, aes(x=log_MI, y=lower), color="blue",linetype='dashed')+
  geom_line(data=x_fit1, aes(x=log_MI, y=upper), color="blue",linetype='dashed')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(bold('ln MI')),y=expression(bold('ln Ac/DH')))+
  annotate ('text',x=-0.75, y=6,cex=5,label=expression(bold("(c)")))

##############################################################
########Huber###############################################
Input_Huber<-read.csv('BAAD.csv',header=TRUE)
Data_input_Huber<-subset(Input_Huber,Huber_crown!='NA')
####################OLS analyses for Huber value
Model_Huber<-lm(log(Huber_crown)~Gst+log(Gppfd),Data_input_Huber)
summary(Model)
####################OLS plot
visreg(fit=Model_Huber,xvar="Gst",band=TRUE,rug=2,jitter=TRUE,points=list(cex=1,col='#293352'),line=list(cex=0.5,col="red"),ylab=expression(bold('ln vH')),xlab=expression(bold('Gst (°C)')))###All species
text(x=11.5,y=-6.3,labels='(d)',font=2,cex=1.2)
visreg(fit=Model_Huber,xvar="Gppfd",band=TRUE,rug=2,xtrans=log,jitter=TRUE,points=list(cex=1,col='#293352'),line=list(cex=0.5,col="red"),ylab=expression(bold('ln vH')),xlab=expression(bold('ln Gppfd (mol m'^-2*'day'^-1*')')))###All species
text(x=2.97,y=-6.4,labels='(e)',font=2,cex=1.2)
#####################robust analyses for Huber value
Data_input_Huber$log_Gppfd<-log(Data_input_Huber$Gppfd)
Data_input_Huber$log_vH<-log(Data_input_Huber$Huber_crown)
fit_Huber<-rlm(log_vH~Gst+log_Gppfd,data=Data_input_Huber,method='MM')
library(MASS)
summary(fit)
install.packages('sfsmisc')
library(sfsmisc)###Calaculate p values for rlm
f.robftest(fit, var = "log_Gppfd")
f.robftest(fit, var = "Gst")
########################Robust plot for Huber value
visreg(fit=fit_Huber,xvar="Gst",band=TRUE,points=list(cex=1,col='#293352'),line=list(cex=0.5,col="#293352"),ylab=expression(bold('ln vH (m'^2*'/m'^2*')')),xlab=expression(bold('Gst (°C)')))###All species
abline(a=-8.95,b=0.058,col='blue')
text(x=12,y=-6.4,labels='(e)',font=2,cex=1.2)
visreg(fit=fit_Huber,xvar="log_Gppfd",band=TRUE,points=list(cex=1,col='#293352'),line=list(cex=0.5,col="#293352"),ylab=expression(bold('ln vH (m'^2*'/m'^2*')')),xlab=expression(bold('ln Gppfd (mol m'^-2*'day'^-1*')')))###All species
abline(a=-1.12,b=-2.13,col='blue')
text(x=2.95,y=-6.4,labels='(f)',font=2,cex=1.2)

######################################################################
###########################The ratio of fine roots to leaf area#######################
Input_roots<-read.csv('BAAD.csv',header=TRUE)
Data_input_roots<-subset(Input_roots,Root_removed!='NA')
##########################OLS analyses for roots ratio
Model_roots<-lm(log(Root_removed)~Gst,Data_input_roots)
summary(Model)
#########################OLS Plot
visreg(fit=Model_roots,xvar="Gst",band=TRUE,rug=2,jitter=TRUE,points=list(cex=1,col='#293352'),line=list(cex=0.5,col="red"),ylab='ln ?? ',xlab=expression(bold('Gst (°C)')))
text(x=11.3,y=-2.5,labels='(c)',font=2,cex=1.2)
#########################robust analyses for roots ratio
Data_input_roots$log_roots<-log(Data_input_roots$Root_removed)
fit_roots<-mblm(log_roots~Gst,data=Data_input_roots,repeated=FALSE)
summary(fit_roots)
confint.mblm(fit_roots,level=0.95)
########################robust plot
effects_fit1 <- effects::effect(term="Gst", mod=fit_roots)
summary(effects_fit1)
x_fit1 <- as.data.frame(effects_fit1)
ggplot(data=Data_input_roots, aes(x=Gst, y=log(Root_removed)))+geom_point(pch=16,size=2)+
  geom_smooth(method='lm',formula=y~x,se=TRUE,color='black',linetype=1)+
  geom_line(data=x_fit1, aes(x=Gst, y=fit), color="blue") +
  geom_line(data=x_fit1, aes(x=Gst, y=lower), color="blue",linetype='dashed')+
  geom_line(data=x_fit1, aes(x=Gst, y=upper), color="blue",linetype='dashed')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(bold('Gst (°C)')),y='ln ?? (g/m2)')+
  annotate ('text',x=11, y=-2.2,cex=5,label=expression(bold("(d)")))

#################################################################
######LAI model evaluation######################################
library(rpmodel)
pmodel<-read.csv('LAI_pmodel_input.csv',header=TRUE)
result <- data.frame(matrix(nrow = 2648, ncol = 8))
colnames(result) <- c("Site",'A0','chi','ca','vpd','Temperature','ppfd','Elevation')
for(i in 1:2648){
  R<-pmodel[(12*i-11):(12*i),5:9]
  R1<-pmodel[1:12,11:16]
  rpmodel_outputs<-rpmodel(tc = as.numeric(R1$MAT), 
                           vpd = as.numeric(R1$vpd), 
                           co2 = as.numeric(R1$CO2),
                           fapar =1.0, 
                           ppfd = as.numeric(R1$ppfd),
                           elv=R1$Elevation,
                           kphio=0.08718,   # quantum yield efficiency,
                           beta= 146)           # unit cost ratio a/b,
  RR<-data.frame(sapply(rpmodel_outputs,c))#####convert from list to matrix then to dataframe,'c'here means a generic function which combines its arguments
  result[i,1]<-i
  result[i,2]<-sum(RR$gpp)/12#per year
  result[i,4]<-sum(RR$ca)/144#per month
  result[i,3]<-sum(RR$chi)/144#per month
  result[i,5]<-mean(R1$vpd)
  result[i,6]<-mean(R$MAT)
  result[i,7]<-mean(R$ppfd)
  result[i,8]<-mean(R$Elevation)
}
Data1<-read.csv('LAI.csv',header=TRUE)
Result<-cbind(result,Data1$MP_predicted,Data1$Gst_predicted,Data1$MI_predicted,Data1$PFTs,Data1$LAI_observed,Data1$fapar_predicted)
colnames(Result) <- c("Site",'A0','chi','ca','D','Temperature','ppfd','Elevation','MP_predicted','Gst_predicted','MI_predicted','PFTs','LAI_observed','fapar_predicted')  
##########################calculate fapar_observed
Result$fapar_observed<-1-exp(-0.5*Result$LAI_observed)
fapar1<-nls(fapar_observed~1-z/I(0.5*A0),Result,start=c(z=267))
fapar2<-nls(fapar_observed~ca*(1-chi)*f0*MP_predicted/I(1.6*D*A0),Result,start=c(f0=74))
summary(fapar1)
summary(fapar2)
############################carbon-limited fapar predicted
Result$fapar_carbon<-1-267/I(0.5*Result$A0)
############################water-limited fapar predicted
Result$fapar_water<-Result$ca*(1-Result$chi)*89*Result$MP_predicted/I(1.6*Result$D*Result$A0)
library(matrixStats)
for (i in 1:length(Data$fapar_observed)){
  Result$fapar_final1[i]<--logSumExp(c(-1000*Result$fapar_carbon[i],-1000*Result$fapar_water[i]),na.rm=TRUE)
  Result$fapar_final2[i]<-Result$fapar_final1[i]/1000
}
#############################calculate predicted LAI
Result$LAI_predicted<--2*log(1-Result$fapar_final2)
############################The climate distribution of fapar differences (observed-predicted)
Data<-read.csv("LAI_combined_data.csv",header=TRUE)
Model<-lm(fapar_fraction~MI_predicted,Data)
summary(Model)
Plot<-subset(Data,fapar_difference!='NA'&Gst_predicted!='NA'&MI_predicted!='NA')
ggplot(data=Plot,aes(x=Gst_predicted,y=MI_predicted,col=fapar_difference))+
  geom_point(size=2)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_gradient(low = "gray1", high = "green")+  
  theme(legend.position = c(0.9, 0.85),
        legend.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 10, colour = "black"))+
  labs(colour='vH')+
  labs(x=expression(bold('Gst (°C)')),y=expression(bold(' MI')))+
  annotate ('text',x=2.2, y=9,cex=4,label=expression(bold('(a)')))
########################Differences between carbon and water in climate distribution
Plot1<-subset(Data,fapar_carbon!='NA')
Plot2<-subset(Data,fapar_water!='NA')
ggplot()+
  geom_point(data=Plot1,aes(x=Gst_predicted, y=MI_predicted),pch=16,size=2,col='brown1',alpha=0.3)+
  geom_point(data=Plot2,aes(x=Gst_predicted, y=MI_predicted),pch=16,size=2,col='slateblue',alpha=0.3)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression(bold('Gst (°C)')),y=expression(bold('MI')))+
  annotate ('text',x=2.2, y=9,cex=4,label=expression(bold('(b)')))
#######################Predicted and observed fAPAR comparison
ggplot()+
  geom_point(data=Data,aes(x=fapar_carbon, y=fapar_observed), color="brown1",alpha=0.3)+
  geom_point(data=Data,aes(x=fapar_water, y=fapar_observed), color="slateblue",alpha=0.3)+
  geom_smooth(data=Data,aes(x=fapar_predicted, y=fapar_observed),method='lm',formula=y~0+x,se=TRUE,linetype=1,color='black')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y='fAPAR(observed)',x='fAPAR(predicted)')+
  geom_abline(intercept=0,slope=1,col='black',linetype='longdash',size=1.05)+
  annotate ('text',x=0.78, y=0.16,cex=3,label=expression(bold("R"^2*'=0.847')))+
  annotate ('text',x=0.8, y=0.06,cex=3,label=expression(bold('RMSE=0.304')))+
  annotate ('text',x=0.04, y=1.16,cex=4,label=expression(bold('(a)')))
##############################Plot A0 and LAI
Data<-read.csv("LAI_combined_data.csv",header=TRUE)
ggplot()+
  geom_point(data=Data,aes(x=A0, y=LAI_used_carbon), color="brown1",alpha=0.3)+
  geom_point(data=Data,aes(x=A0, y=LAI_used_water), color="slateblue",alpha=0.3)+
  labs(y='Observed LAI',x=expression('A'[0]*'(g C'^-1*'m'^-2*'year'^-1*')'))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate ('text',x=700, y=12,cex=4,label=expression(bold('(b)')))

##########################################################
###Calculation number##Table S2############################
Data<-read.csv('BAAD.csv',header=TRUE)
Data<-read.csv("LAI_combined_data.csv",header=TRUE)
Data<-subset(Data[1:2648,])
Data_input<-subset(Data,LAI_used!='NA'&Type=='deciduous'&Type2=='angiosperm')###Change those PFTs names respectively to calculate the number
#Data number
length(Data_input$LAI_used)
#Site number
length(unique(Data_input$Site))
#Species number
length(unique(Data_input$Species))

##########################################################################
##Create the global distribution of each trait############################
install.packages('sf') # Core vector GIS data package
install.packages('sp') # Another core vector GIS package
install.packages('rgeos') # Extends vector data functionality
install.packages('rgdal') # Interface to the Geospatial Data Abstraction Library
install.packages('lwgeom') # Extends vector data functionality
library(rgeos)
library(rgdal)
library(lwgeom)
library(sf)
library(sp)
#############################LAI
Data_LAI<-read.csv("LAI.csv",header=TRUE)
Plot<-subset(Data_LAI,LAI!='NA')
library(viridis)
ggplot() +
  geom_map(
    data = world, map = world,
    aes(x=long, y=lat, map_id = region),
    color = "grey", fill = "darkgrey", size = 0.1
  ) +
  geom_point(
    data = Plot,
    aes(x=Longitude, y=Latitude,col=LAI),size=1,
    alpha = 0.7
  )+ scale_color_gradient(low = "blue", high = "red")+
  theme_void()+
  labs(colour='LAI')

##############################Ac/DH
Data_Ac<-read.csv("Ac_D_H_Total_data.csv",header=TRUE)
Plot<-subset(Data_Ac,R!='NA')
ggplot() +
  geom_map(
    data = world, map = world,
    aes(x=long, y=lat, map_id = region),
    color = "grey", fill = "darkgrey", size = 0.1
  ) +
  geom_point(
    data = Plot,
    aes(x=lon, y=lat,col=R),size=1,alpha = 0.7)+  
  scale_color_gradient(low = "blue", high = "red")+
  theme_void()+
  labs(colour='Ac/DH') 

###########################vH
Data_Huber<-read.csv("BAAD.csv",header=TRUE)
Plot<-subset(Data_Huber,Huber_crown!='NA')
ggplot() +
  geom_map(
    data = world, map = world,
    aes(x=long, y=lat, map_id = region),
    color = "grey", fill = "darkgrey", size = 0.1
  ) +
  geom_point(
    data = Plot,
    aes(x=lon, y=lat,col=Huber_crown),size=4,alpha = 0.7)+  
  scale_color_gradient(low = "blue", high = "red")+
  theme_void()+
  labs(colour='vH') 

###############################Hm
Plot<-subset(Data,b!='NA')
ggplot() +
  geom_map(
    data = world, map = world,
    aes(x=long, y=lat, map_id = region),
    color = "grey", fill = "darkgrey", size = 0.1
  ) +
  geom_point(
    data = Plot,
    aes(x=lon, y=lat,col=b),size=1,alpha = 0.7)+  
  scale_color_gradient(low = "blue", high = "red")+
  theme_void()+
  labs(colour=expression('b (m'^-1*')')) 

#################################################################
#########Dot plots###############################################
#########################LAI
Data<-read.csv("LAI.csv",header=TRUE)
Plot<-subset(Data,LAI_used!='NA'&PFTs!='NA')
ggplot(data=Plot,aes(x=Gst_predicted,y=MI_predicted,col=PFTs))+
  geom_point(size=2)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 10, colour = "black"))+
  labs(colour='vH')+
  labs(x=expression(bold('Gst (°C)')),y=expression(bold('ln MI')))+
  annotate ('text',x=3, y=8.5,cex=5,label=expression(bold("(a)")))
###########################Ac/DH
Data<-read.csv("Ac_D_H_Total_data.csv",header=TRUE)
Plot<-subset(Data,R!='NA'&Gst!='NA'&MI!='NA'&PFTs!='NA')
ggplot(data=Plot,aes(x=Gst,y=MI,col=PFTs))+
  geom_point(size=2)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 10, colour = "black"))+
  labs(colour='R')+
  labs(x=expression(bold('Gst (°C)')),y=expression(bold('ln MI')))+
  annotate ('text',x=3, y=15,cex=5,label=expression(bold("(b)"))) 
################################Hm/b
Data<-read.csv("Results_combined_evergreen.csv",header=TRUE)
Plot<-subset(Data,Hm!='NA'&Gst!='NA'&MI!='NA')
ggplot(data=Plot,aes(x=Gst,y=MI,col=Type))+
  geom_point(size=2)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 10, colour = "black"))+
  labs(colour='vH')+
  labs(x=expression(bold('Gst (°C)')),y=expression(bold('ln MI')))+
  annotate ('text',x=3, y=37,cex=5,label=expression(bold("(d)")))  

#######################################################################
###Bubble plots########################################################
###########################roots
Data<-read.csv('BAAD.csv',header=TRUE)
Data_input<-subset(Data,Root_removed!='NA')
ggplot(data=Data_input,aes(x=MT,y=MP,size=Root_removed,col=PFT))+
  geom_point()+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.text = element_text(size = 10, colour = "black"))+
  labs(y=expression('MP (mm)'),size='?? (kg/m2)')+
  labs(x=expression(bold('MT (°C)')),y=expression(bold('MP (mm)')))+
  annotate ('text',x=6, y=270,cex=5,label=expression(bold("(a)")))  
############################huber value
Data<-read.csv('BAAD.csv',header=TRUE)
Data_input<-subset(Data,Huber_crown!='NA')
ggplot(data=Data_input,aes(x=Gst,y=Huber_crown,col=PFT,size=Elevation))+
  geom_point()+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.background = element_blank(),
        legend.key.size = unit(4,"mm"),
        legend.text = element_text(size = 10, colour = "black"))+
  labs(x=expression('Gppfd (mol m'^-2*'day'^-1*')'),y='vH',size='Elevation (m)')+
  labs(x=expression(bold('Gst (°C)')),y=expression(bold('ln MI')))+
  annotate ('text',x=11.5, y=0.00185,cex=5,label=expression(bold("(b)")))  



