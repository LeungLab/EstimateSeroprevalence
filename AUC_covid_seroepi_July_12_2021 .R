library(tidyverse)
library(party)
library(permimp)
library(furrr)
library(cvAUC)
library(pROC)

setwd("/Users/asbhuiya/Documents/R/Covid_19_data/Data_last_time_points/last time/Data")
getwd()

# Make training and testing data 
make_data=function(dat){
  nrows=dim(dat)[1]
  train_idx=sample(1:nrows,round(nrows*.8)) # this is 5 fold cross-validation 
  train=dat[train_idx,]
  test=dat[-train_idx,]
  return(list(train,test))
}

#Function for cvAUC from classificaiton
Cross_val=function(train,test,nvars,ntrees=1000){ #added ntrees and mtry as settings here.
  out=cforest(case_control~.,data=train,control=cforest_unbiased(ntree=ntrees,mtry=ceiling((ncol(train)-1)/3)))#mtry=ceiling(sqrt(ncol(train)-1))))
  vimp=permimp(out, conditional = TRUE, asParty = TRUE)
  df=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% arrange(desc(Cond_Importance))
  mod=cforest(as.formula(paste('case_control~',paste(df$Biomarker[1:nvars],collapse="+"),sep="")),data=train,control=cforest_unbiased(ntree=ntrees,mtry=ceiling(nvars/3)))#ceiling(sqrt(nvars))))#
  return(predict(mod,newdata=test))
}

#Function for calculation of MAE
Cross_val1=function(train,test,nvars=3,ntrees=1000,mtry=1){
  out=cforest(stday2~.,data=train,control=cforest_unbiased(ntree=ntrees,mtry=mtry))
  vimp=permimp(out, conditional = TRUE, asParty = TRUE)
  df=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% arrange(desc(Cond_Importance))
  mod=cforest(as.formula(paste('stday2~',paste(df$Biomarker[1:nvars],collapse="+"),sep="")),
              data=train,control=cforest_unbiased(ntree=ntrees,mtry=mtry))#ceiling(sqrt(nvars))))
  return(mean(abs(test$stday2-predict(mod,newdata=test))))
}

#Function to calculate AUC
auc=function(dat,niter=niter,nvars=nvars){
  train_test=1:niter %>% purrr::map(~make_data(dat))
  res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=nvars,ntrees=1000))
  test_df=bind_rows(train_test %>% purrr::map(~.[[2]]),.id="Iter")
  test_predictions=unlist(res)
  test_df$prediction=test_predictions
  test_df=test_df %>% mutate(Iter=as.numeric(Iter))
  AUC=with(test_df,ci.cvAUC(prediction, case_control, folds = Iter))
  return(AUC)
}

#Load data
dat_Dan <- read.csv("Dan_et_al_last_time_d30_pos.csv")
dat_Peluso <- read.csv("Peluso_et_al_last_time_d30_pos_neg.csv")
dat_Whitcombe <- read.csv("Whitcombe_et_al_last_time_d30_pos_neg.csv")
dat_Markmann<-read.csv("Markmann_et_al_last_time_d30_pos_neg.csv")
dat_Isho<- read.csv("Isho_et_al_last_time_d30_pos_neg.csv")

### Peluso data
dat_Peluso=dat_Peluso %>% dplyr::mutate_at(vars(c( "age","sex", "hosp_status", "case_control")), ~as.factor(.)) %>%
  dplyr::mutate_at(vars(c("stday2","N_Abbott", "S_DiaSorin", "N_LIPS","S_LIPS", "N.full._Lum","N.frag._Lum","RBD_Lum", "S_Lum", "N_Roche", 
                          "N_Split_Luc","RBD_Split_Luc","S_Ortho_IgG", "S_Ortho_Ig")),~as.numeric(.))

str(dat_Peluso)
dat_Peluso %>% summarize_all(~sum(is.na(.)))

# cvAUC for
dat=dat_Peluso %>% dplyr::select(case_control=case_control,N.full._Lum,N.frag._Lum,RBD_Lum, S_Lum) %>% 
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
auc1_P=auc(dat,niter=100,nvars=1)#AUC for top one biomarker
auc2_P=auc(dat,niter=100,nvars=2)#AUC for top two biomarker
auc3_P=auc(dat,niter=100,nvars=3)#AUC for top three biomarker

#cvAUC for single biomarker
dat=dat_Peluso %>% dplyr::select(case_control=case_control,N.full._Lum) %>% 
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
auc_N.full_lum=auc(dat,niter=100,nvars=1)
#Nucleocapsid IgG
dat=dat_Peluso %>% dplyr::select(case_control=case_control,"N.frag._Lum") %>% 
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
auc_N.frag_lum=auc(dat,niter=100,nvars=1)
#RBD IgG
dat=dat_Peluso %>% dplyr::select(case_control=case_control,"RBD_Lum") %>% 
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
auc_RBD_lum=auc(dat,niter=100,nvars=1)
#Spike IgG
dat=dat_Peluso %>% dplyr::select(case_control=case_control, "S_Lum") %>% 
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
auc_S_lum=auc(dat,niter=100,nvars=1)

#AUC for saturated model
dat=dat_Peluso  %>% dplyr::select(case_control=case_control, N.full._Lum, N.frag._Lum, RBD_Lum, S_Lum)%>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
auc_Full_Peluso=auc(dat,niter = 100,nvars=4)

### Whitcombe data
dat_Whitcombe=dat_Whitcombe %>% dplyr::mutate_at(vars(c("age", "sex")), ~as.factor(.)) %>% 
  dplyr::mutate_at(vars(c("stday2", ends_with('IgA'),ends_with('IgG'),ends_with('IgM'))),~as.numeric(.))
dat_Whitcombe %>% summarize_all(~sum(is.na(.)))
str(dat_Whitcombe)

#cvAUC for Whitcombe

niter=100
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,ends_with ('IgA'), ends_with('IgG'), ends_with('IgM')) %>% 
  mutate(case_control=as.numeric(case_control))
auc1_W=auc(dat,niter = 100,nvars=1)
auc2_W=auc(dat,niter = 100,nvars=2)
auc3_W=auc(dat,niter = 100,nvars=3)

#single biomarker AUC
#Spike IgG
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,S_IgG) %>% 
  mutate(case_control=as.numeric(case_control))
auc_S_IgG=auc(dat,niter = 100,nvars=1)
#RBD IgG
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,RBD_IgG) %>% 
  mutate(case_control=as.numeric(case_control))
auc_R_IgG=auc(dat,niter = 100,nvars=1)
#Nucleocapsid IgG
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,NP_IgG) %>% 
  mutate(case_control=as.numeric(case_control))
auc_N_IgG=auc(dat,niter = 100,nvars=1)

#best two from importance list
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,NP_IgA, RBD_IgA) %>% 
  mutate(case_control=as.numeric(case_control))
auc_N_R_aa=auc(dat,niter = 100,nvars=2)
# check best two NP
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,NP_IgA, NP_IgG) %>% 
  mutate(case_control=as.numeric(case_control))
auc_N_ag=auc(dat,niter = 100,nvars=2)
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,NP_IgM, NP_IgG) %>% 
  mutate(case_control=as.numeric(case_control))
auc_N_gm=auc(dat,niter = 100,nvars=2)
#saturated model
dat=dat_Whitcombe %>% dplyr::select(case_control=case_control,ends_with ('IgA'), ends_with('IgG'), ends_with('IgM')) %>% 
  mutate(case_control=as.numeric(case_control))
auc_Full_Whit=auc(dat,niter=100,nvars=9)

#Isho dataset
dat_Isho=dat_Isho %>% dplyr::mutate_at(vars(c( "age","sex", "case_control")), ~as.factor(.)) %>%
  dplyr::mutate_at(vars(c("stday2",ends_with ('iga'), ends_with('igg'), ends_with('igm'))),~as.numeric(.))

#Calculate the AUC for classification
dat=dat_Isho %>% dplyr::select(case_control=case_control,spike_rob_igg,spike_rob_igm,spike_rob_iga,rbd_rob_igg,rbd_rob_igm,rbd_rob_iga,np_rob_igg,np_rob_igm,np_rob_iga) %>% 
mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a")
auc1_I=auc(dat,niter=100, nvars = 1)#AUC for top one biomarker
auc2_I=auc(dat, niter=100, nvars = 2)#AUC for top two biomarker
auc3_I=auc(dat, niter=100, nvars = 3)#AUC for top three biomarker

#specific single marker AUC
#Nucleocapsid IgG
dat=dat_Isho %>% dplyr::select(case_control=case_control,np_rob_igg) %>% mutate(case_control=as.numeric(case_control)) 
auc_np_igg=auc(dat,niter=100, nvars = 1) 
#RBD IgG
dat=dat_Isho %>% dplyr::select(case_control=case_control,rbd_rob_igg) %>% mutate(case_control=as.numeric(case_control)) 
auc_rbd_igg=auc(dat,niter=100, nvars = 1)
#Spike IgG
dat=dat_Isho %>% dplyr::select(case_control=case_control,spike_rob_igg) %>% mutate(case_control=as.numeric(case_control)) 
auc_spike_igg=auc(dat,niter=100, nvars = 1)
#Check for best two np marker
dat=dat_Isho %>% dplyr::select(case_control=case_control,np_rob_igg, np_rob_igm) %>% mutate(case_control=as.numeric(case_control)) 
auc_np_gm=auc(dat,niter=100, nvars = 2) 
dat=dat_Isho %>% dplyr::select(case_control=case_control,np_rob_igg, np_rob_iga) %>% mutate(case_control=as.numeric(case_control)) 
auc_np_ga=auc(dat,niter=100, nvars = 2) 
#Saturated model for AUC
dat=dat_Isho %>% dplyr::select(case_control=case_control,spike_rob_igg,spike_rob_igm,spike_rob_iga,rbd_rob_igg,rbd_rob_igm,rbd_rob_iga,np_rob_igg,np_rob_igm,np_rob_iga) %>% mutate(case_control=as.numeric(case_control)) %>%  filter_all(~.!='n.a.')
aucFull=auc(dat,niter=100, nvars=9)


#Markmann data
#get auc (we have only two markers for pre-pandemic controls)
niter=100
dat=dat_Markmann %>% dplyr::select(case_control,RBD_IgG, NP_IgG) %>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!='n.a.')
np=auc(dat,niter = 100,nvars=1)
rp=auc(dat,niter = 100,nvars=2)

# RBD IgG
dat=dat_Markmann %>% dplyr::select(case_control,RBD_IgG) %>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!='n.a.')
auc_R_igg_Mark=auc(dat,niter = 100,nvars=1)
#Nucleocapsid IgG
dat=dat_Markmann %>% dplyr::select(case_control, NP_IgG) %>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!='n.a.')
auc_N_igg_Mark=auc(dat,niter = 100,nvars=1)

### all summary statistics
#Load data again
dat_Dan <- read.csv("Dan_et_al_last_time_d30_pos.csv")
dat_Peluso <- read.csv("Peluso_et_al_last_time_d30_pos_neg.csv")
dat_Whitcombe <- read.csv("Whitcombe_et_al_last_time_d30_pos_neg.csv")
dat_Markmann<-read.csv("Markmann_et_al_last_time_d30_pos_neg.csv")
dat_Isho<- read.csv("Isho_et_al_last_time_d30_pos_neg.csv")
#for median time since infection
dat_Dan  %>% summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))
dat_Peluso %>% filter(case_control=="1") %>% summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))
dat_Whitcombe %>% filter(case_control=="1") %>%  summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))
dat_Markmann %>% filter(case_control=="1") %>% summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))
dat_Isho %>% filter(pn=="1") %>% summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))

#for age (for Peluso data age is binned)
dat_Dan  %>% summarize(quantile(age,.25), median(age), quantile(age,.75))
dat_Whitcombe %>% filter(case_control=="1") %>% mutate_at(vars(c("age")), ~as.numeric(.)) %>% summarize(quantile(age,.25), median(age), quantile(age,.75))
dat_Markmann %>% filter(case_control=="1") %>% mutate_at(vars(c("age")), ~as.numeric(.)) %>% filter_all(~.!="n.a.") %>% summarize(quantile(age,.25), median(age), quantile(age,.75))
dat_Isho %>% filter(pn=="1") %>% mutate_at(vars(c("age")), ~as.numeric(.)) %>% summarize(quantile(age,.25), median(age), quantile(age,.75))

#for sex
dat_Dan %>% group_by(sex) %>% summarise(n=n()) %>% summarise(pc=paste0((n[sex=="M"]/sum(n)*100 )))
dat_Peluso %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) %>% summarise(pc=paste0((n[sex=="Male"]/sum(n)*100 )))
dat_Whitcombe %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) %>% summarise(pc=paste0((n[sex=="M"]/sum(n)*100 )))
dat_Markmann %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) %>% summarise(pc=paste0((n[sex=="Male"]/sum(n)*100 )))
dat_Isho %>% filter(pn=="1") %>% group_by(sex) %>% summarise(n=n()) %>% summarise(pc=paste0((n[sex=="1 Male"]/sum(n)*100 )))

# number of patients or subject
D=dat_Dan  %>% nrow()
P=dat_Peluso %>% filter(case_control=="1") %>% nrow() 
W=dat_Whitcombe %>% filter(case_control=="1") %>% nrow()
M=dat_Markmann %>% filter(case_control=="1") %>% nrow()
I=dat_Isho %>% filter(pn=="1") %>% nrow()
D+P+W+M+I
# over all age quantile
Dan_age=dat_Dan$age
Whit_age=(dat_Whitcombe %>% filter(case_control=="pos"))$age
Mark_age=(dat_Markmann %>% filter(case_control=="1"))$age
Isho_age=(dat_Isho %>% filter(pn=="1"))$age
write.csv(Dan_age,"Dan_age.csv", row.names = F)
write.csv(Whit_age,"Whit_age.csv", row.names = F)
write.csv(Mark_age,"Mark_age.csv", row.names = F)
write.csv(Isho_age,"Isho_age.csv", row.names = F)
Dan_age=read.csv("Dan_age.csv") %>% mutate(age=x)
Whit_age=read.csv("Whit_age.csv") %>% mutate(age=x)
Mark_age=read.csv("Mark_age.csv") %>% mutate(age=x)
Isho_age=read.csv("Isho_age.csv") %>% mutate(age=x)
age_quantile=as.data.frame(rbind(Dan_age,Whit_age, Mark_age,Isho_age)) %>% filter_all(~.!='n.a.') %>% summarize(quantile(age,.25), median(age), quantile(age,.75))
age_quantile
#over all stday2 quantile
Dan_stday2=dat_Dan$stday2
Pelu_stday2=(dat_Peluso %>% filter(case_control=="1"))$stday2
Whit_stday2=(dat_Whitcombe %>% filter(case_control=="pos"))$stday2
Mark_stday2=(dat_Markmann %>% filter(case_control=="1"))$stday2
Isho_stday2=(dat_Isho %>% filter(pn=="1"))$stday2
write.csv(Dan_stday2,"Dan_stday2.csv", row.names = F)
write.csv(Pelu_stday2, "Pelu_stday2.csv",row.names = F)
write.csv(Whit_stday2,"Whit_stday2.csv", row.names = F)
write.csv(Mark_stday2,"Mark_stday2.csv", row.names = F)
write.csv(Isho_stday2,"Isho_stday2.csv", row.names = F)
Dan_stday2=read.csv("Dan_stday2.csv") %>% mutate(stday2=x)
Pelu_stday2=read.csv("Pelu_stday2.csv") %>% mutate(stday2=x)
Whit_stday2=read.csv("Whit_stday2.csv") %>% mutate(stday2=x)
Mark_stday2=read.csv("Mark_stday2.csv") %>% mutate(stday2=x)
Isho_stday2=read.csv("Isho_stday2.csv") %>% mutate(stday2=x)
stday2_quantile=as.data.frame(rbind(Dan_stday2,Pelu_stday2,Whit_stday2, Mark_stday2,Isho_stday2)) %>% filter_all(~.!='n.a.') %>% summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))
stday2_quantile
Dan_sex=dat_Dan %>% group_by(sex) %>% summarise(n=n()) 
Pelu_sex=dat_Peluso %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) 
Whit_sex=dat_Whitcombe %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) 
Mark_sex=dat_Markmann %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) 
Isho_sex=dat_Isho %>% filter(case_control=="1") %>% group_by(sex) %>% summarise(n=n()) 
sex_all=as.data.frame(rbind(Dan_sex,Pelu_sex,Whit_sex,Mark_sex,Isho_sex))
sex_all$sex1=0
sex_all$sex1[sex_all$sex %in% "M"]<- 1
sex_all$sex1[sex_all$sex %in% "F"]<- 2
sex_all$sex1[sex_all$sex %in% "1 Male"]<- 1
sex_all$sex1[sex_all$sex %in% "1 Female"]<- 2
sex_all$sex1[sex_all$sex %in% "Male"]<- 1
sex_all$sex1[sex_all$sex %in% "Female"]<- 2

sex_quantile=sex_all %>% group_by(sex1) %>% summarise(n=n()) %>% summarise(pc=paste0((n[sex1=="1"]/sum(n)*100 )))
sex_quantile

#Conditional importance plot 

dfDan=read.csv("Dan_cond_imp1.csv")
dfIsho=read.csv("Isho_cond_imp1.csv")
dfWhitcombe=read.csv("Whitcombe_cond_imp1.csv")
dfMarkmann=read.csv("Markmann_cond_imp1.csv")
dfPeluso=read.csv("Peluso_cond_imp1.csv")
g1=ggplot(dfDan, aes(x=reorder(Biomarker,Cond_Importance), y=Cond_Importance,fill=Cond_Importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("")+
  xlab("")+
  ggtitle("Dan et. al.")+
  guides(fill=F)+
  scale_fill_gradient(low="black", high="black")
g2=ggplot(dfPeluso, aes(x=reorder(Biomarker,Cond_Importance), y=Cond_Importance,fill=Cond_Importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("")+
  xlab("")+
  ggtitle("Peluso et. al.")+
  guides(fill=F)+
  scale_fill_gradient(low="black", high="black")
g3=ggplot(dfWhitcombe, aes(x=reorder(Biomarker,Cond_Importance), y=Cond_Importance,fill=Cond_Importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("")+
  xlab("")+
  ggtitle("Whitcombe et. al.")+
  guides(fill=F)+
  scale_fill_gradient(low="black", high="black")
g4=ggplot(dfIsho, aes(x=reorder(Biomarker,Cond_Importance), y=Cond_Importance,fill=Cond_Importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("")+
  xlab("")+
  ggtitle("Isho et. al.")+
  guides(fill=F)+
  scale_fill_gradient(low="black", high="black")
g5=ggplot(dfMarkmann, aes(x=reorder(Biomarker,Cond_Importance), y=Cond_Importance,fill=Cond_Importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("")+
  xlab("")+
  ggtitle("Markmann et. al.")+
  guides(fill=F)+
  scale_fill_gradient(low="black", high="black")
library(ggpubr)
figure=ggarrange(g1, g2,g3,g4,g5, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)

gg=annotate_figure(figure, bottom=text_grob("Variable Importance"), 
                right= "Importance for Time-since-infection")
gg
ggsave("conditional_imp_Plot.pdf", width=20, height=20, unit="cm")
ggsave("conditional_imp_Plot.tiff", width=20, height=20, unit="cm")
ggsave("conditional_imp_Plot", device="pdf",width=20, height=20, unit="cm")
