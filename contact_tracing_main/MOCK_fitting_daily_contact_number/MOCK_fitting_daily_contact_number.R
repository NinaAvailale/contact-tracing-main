rm(list = ls())
Sys.setlocale("LC_TIME","english")
setwd("C:\\contact_tracing_main\\MOCK_fitting_daily_contact_numbers")
getwd()


# library("writexl")

###Load samples of daily contact numbers in different locations.
ctnum_part <- read.csv("mock_data_ct_nocontacts_for_fit.csv")

# install.packages("fitdistrplus")
library(fitdistrplus)
unique(ctnum_part$patt)
par(mfrow=c(1,4))


gamma_pa<-as.data.frame(matrix(data=NA,ncol=3,nrow = 0))
colnames(gamma_pa)<-c("patt","shape","rate")
gamma_pa[1:5,1]<-c("Dwelling","Work","Cospace-time","Community","Other")


com<-subset(ctnum_part,ctnum_part$patt == "Dwelling")
model1<-fitdist(com$no.contacts,"gamma")
qqcomp(model1, legendtext = "Gamma")
model2<-fitdist(com$no.contacts,"weibull")
qqcomp(model2,legendtext = "Weibull")
model3<-fitdist(com$no.contacts,"lnorm")
qqcomp(model3,legendtext = "Lognormal") 
model4<-fitdist(com$no.contacts,"exp")
qqcomp(model4,legendtext = "Exponential")
gamma_pa$shape[1]<-model1$estimate['shape']
gamma_pa$rate[1]<-model1$estimate['rate']
par(mfrow=c(1,4))
com<-subset(ctnum_part,ctnum_part$patt == "Work")
model1<-fitdist(com$no.contacts,"gamma")
qqcomp(model1, legendtext = "Gamma")
model2<-fitdist(com$no.contacts,"weibull")
qqcomp(model2,legendtext = "Weibull")
model3<-fitdist(com$no.contacts,"lnorm")
qqcomp(model3,legendtext = "Lognormal") 
model4<-fitdist(com$no.contacts,"exp")
qqcomp(model4,legendtext = "Exponential")
gamma_pa$shape[2]<-model1$estimate['shape']
gamma_pa$rate[2]<-model1$estimate['rate']
par(mfrow=c(1,4))
com<-subset(ctnum_part,ctnum_part$patt == "Cospace-time")
model1<-fitdist(com$no.contacts,"gamma")
qqcomp(model1, legendtext = "Gamma")
model2<-fitdist(com$no.contacts,"weibull")
qqcomp(model2,legendtext = "Weibull")
model3<-fitdist(com$no.contacts,"lnorm")
qqcomp(model3,legendtext = "Lognormal") 
model4<-fitdist(com$no.contacts,"exp")
qqcomp(model4,legendtext = "Exponential")
gamma_pa$shape[3]<-model1$estimate['shape']
gamma_pa$rate[3]<-model1$estimate['rate']
par(mfrow=c(1,4))
com<-subset(ctnum_part,ctnum_part$patt == "Community")
model1<-fitdist(com$no.contacts,"gamma")
qqcomp(model1, legendtext = "Gamma")
model2<-fitdist(com$no.contacts,"weibull")
qqcomp(model2,legendtext = "Weibull")
model3<-fitdist(com$no.contacts,"lnorm")
qqcomp(model3,legendtext = "Lognormal") 
model4<-fitdist(com$no.contacts,"exp")
qqcomp(model4,legendtext = "Exponential")
gamma_pa$shape[4]<-model1$estimate['shape']
gamma_pa$rate[4]<-model1$estimate['rate']
par(mfrow=c(1,4))
com<-subset(ctnum_part,ctnum_part$patt == "Other")
model1<-fitdist(com$no.contacts,"gamma")
qqcomp(model1, legendtext = "Gamma")
model2<-fitdist(com$no.contacts,"weibull")
qqcomp(model2,legendtext = "Weibull")
model3<-fitdist(com$no.contacts,"lnorm")
qqcomp(model3,legendtext = "Lognormal") 
model4<-fitdist(com$no.contacts,"exp")
qqcomp(model4,legendtext = "Exponential")
gamma_pa$shape[5]<-model1$estimate['shape']
gamma_pa$rate[5]<-model1$estimate['rate']

# write.table(gamma_pa,"gamma_fitting_para_nocontacts.csv",row.names = F,
            # col.names = T,sep = ",")


shape_para<-gamma_pa[,2]
rate_para<-gamma_pa[,3]

##caculate the mean and variance of daily contact number distribution of 
##fitted gamma distribution under each contact pattern
gamma_sta<-gamma_pa
gamma_sta[,2:3]<-NA
gam_fit_sta<-gamma_sta

colnames(gamma_sta)<-c("patt","mean","variance")
colnames(gam_fit_sta)<-c("patt","mean","variance")

for (i in 1:5)
{
  
  gam_fit_sta$mean[i]<-shape_para[i]/rate_para[i]
  gam_fit_sta$variance[i]<-shape_para[i]/(rate_para[i])^2
}

# write.table(gam_fit_sta,"gam_fit_statistics_nocontacts.csv",row.names = FALSE,col.names = TRUE,sep = ",")






