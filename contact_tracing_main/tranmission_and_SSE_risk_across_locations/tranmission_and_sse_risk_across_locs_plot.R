rm(list=ls())
Sys.setlocale("LC_TIME","english")
setwd("C:\\contact_tracing_main")
getwd()
library(ggplot2)
library(RColorBrewer)


mytheme2<- theme(axis.text.x = element_text(size = 6, color="black",vjust = 0.5, angle = 90,hjust = 1),
                 plot.margin=unit(rep(0.5,4),'lines'),
                 axis.text.y = element_text(size = 6,color = "black",vjust = 0.5,angle = 0),
                 axis.title.y = element_text(size = 7,color="black",vjust=2),
                 axis.title.x = element_text(size = 7, color="black",vjust = 0),
                 panel.border = element_rect(colour = NA, fill=NA, size=0)
)

###load solved p(tau) under dwelling location during model solving process
df.p.tau<-read.csv("model/ptau_by_time.csv")
df.p.tau[which.max(df.p.tau$beta),1]
p_constant<-sum(df.p.tau$beta)
ci_perstep_constant<-1/0.05


###load average transmission rate constant vs home for each contact patterns (derived from epidemiological data and contact tracing collected)
con_trans_vs_home<-read.csv("data_processed/contacts_by_patt_transmission_vs_dwelling.csv")



 scientific_10 <- function(x) {
   gsub("e", " x 10^", scientific_format()(x))
 }
 
###load gamma fitted parameters for distributions of daily contact numbers for each contact pattern 
###(derived from contact tracing data collected)
 
 gamma_pa<-read.csv("data_processed/gamma_fitting_para_nocontacts.csv")
 gamma_pa<-as.data.frame(gamma_pa)
 shape_para<-gamma_pa[,2]
 rate_para<-gamma_pa[,3]
 gamma_sta<-read.csv("data_processed/gam_fit_statistics_nocontacts.csv")
 gamma_cases_sta<-gamma_sta
con_trans_vs_home<-con_trans_vs_home[c(3,5,2,1,4),]
for (i in 1:5)
{

  gamma_cases_sta$mean[i]<-shape_para[i]*(p_constant*con_trans_vs_home[i,2])/(rate_para[i]*ci_perstep_constant)
  gamma_cases_sta$variance[i]<-shape_para[i]*((p_constant*con_trans_vs_home[i,2])/(rate_para[i]*ci_perstep_constant))^2
}

write.table(gamma_cases_sta,"tranmission_and_SSE_risk_across_locations/gammafit_secondary_cases_statistics_by_location.csv",row.names = FALSE,
            col.names = TRUE,sep = ",")


abtenpro<-as.data.frame(matrix(data = NA,nrow = 5,ncol = 2))
colnames(abtenpro)<-c("patt","probability")

abtenpro[,1]<-c('Dwelling', 'Work', 'Cospace-time','Community',"Other")
x<-seq(0,1000,.01)
con2<-con_trans_vs_home
con2[,3]<-shape_para
con2[,4]<-rate_para*ci_perstep_constant
####find the 99 quantile of possion distribution for each patt
quan_r<-abtenpro
colnames(quan_r)[2]<-"R0"

frac<-read.csv("model/frac_of_R_across_all_locs.csv")
frac$patt<-factor(frac$patt,c('Dwelling', 'Work', 'Cospace-time','Community',"Other"))
R0<-10
quan_r[,2] <- round(R0*frac$contribution,3)
quan_r[,3]<-qpois(0.99,quan_r[,2])
x<-seq(0,1000,.01)

for (i in 1:5)
{abtenpro[i,2]<-1-pgamma(quan_r[i,3],con2[i,3],con2[i,4]/(p_constant *con2[i,2]))

}

abtenpro$patt<-factor(abtenpro$patt,c('Dwelling', 'Work', 'Cospace-time','Community',"Other"))
a<-ggplot(abtenpro, aes(patt, probability, fill= patt)) +
  geom_col(position = "dodge",width = 0.5)+
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(2,8,6,1,4)])+
  ylab(paste0("Probability of SSE"))+
  xlab("")+
  theme_test(base_size = 13)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.key.size = unit(0.35, 'cm'))+
  theme( legend.text=element_text(size=13))+
  theme_classic()+
  theme(legend.position = "none")

ggsave(a, file=paste("Figures/Figure3B",".pdf",sep=""),width=4.5,height=3,dpi = 500) 


b<-ggplot(frac, aes(patt, contribution, fill= patt)) +
  geom_col(position = "dodge",width = 0.5)+
  scale_fill_manual(values = brewer.pal(8, "Set2")[c(2,8,6,1,4)])+
  ylab(paste0("Contribution to transmission"))+
  xlab("")+
  theme_test(base_size = 13)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.key.size = unit(0.35, 'cm'))+
  theme( legend.text=element_text(size=13))+
  theme_classic()+
  theme(legend.position = "none")
ggsave(b, file=paste("Figures/Figure3A",".pdf",sep=""),width=4.5,height=3,dpi = 500) 


