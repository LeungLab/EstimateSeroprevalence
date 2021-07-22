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

#Function for calculation of MAE
Cross_val1=function(train,test,nvars=3,ntrees=1000,mtry=1){
  out=cforest(stday2~.,data=train,control=cforest_unbiased(ntree=ntrees,mtry=mtry))
  vimp=permimp(out, conditional = TRUE, asParty = TRUE)
  df=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% arrange(desc(Cond_Importance))
  mod=cforest(as.formula(paste('stday2~',paste(df$Biomarker[1:nvars],collapse="+"),sep="")),
              data=train,control=cforest_unbiased(ntree=ntrees,mtry=mtry))#ceiling(sqrt(nvars))))
  return(mean(abs(test$stday2-predict(mod,newdata=test))))
}

#Load data
dat_Dan <- read.csv("Dan_et_al_last_time_d30_pos.csv")
dat_Peluso <- read.csv("Peluso_et_al_last_time_d30_pos_neg.csv")
dat_Whitcombe <- read.csv("Whitcombe_et_al_last_time_d30_pos_neg.csv")
dat_Markmann<-read.csv("Markmann_et_al_last_time_d30_pos_neg.csv")
dat_Isho<- read.csv("Isho_et_al_last_time_d30_pos_neg.csv")

# Dan et al importance
dat_Dan=dat_Dan %>% mutate_at(vars(c("sex","ethnicity","severity")),~as.factor(.)) %>%
  mutate_at(vars("stday2","age",ends_with('igg'), ends_with('iga')),~as.numeric(.))
str(dat_Dan)
dat_Dan %>% summarize(quantile(stday2,.25), median(stday2), quantile(stday2,.75))# for median (IQR)

dat_Dan=dat_Dan %>% filter_all(~.!="n.a.")
dat=dat_Dan %>% 
  dplyr::select(age,sex,ethnicity,severity,stday2,ends_with('igg'), ends_with('iga'))
out=cforest(stday2~.,data=dat,control=cforest_unbiased(ntree=1000,mtry=ceiling((ncol(dat)-1)/3))) 
vimp=permimp(out, conditional = TRUE, asParty = TRUE)
df=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% 
  arrange(desc(Cond_Importance))
df
# Dan et al MAE
niter=100
dat=dat_Dan %>% dplyr::select(stday2,ends_with('igg'),ends_with('iga'))
train_test_Dan=1:niter %>% purrr::map(~make_data(dat))
#MAE for single biomarker
#plan(multisession)
res_Dan1=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000,mtry=1),.progress=T)
#MAE for two biomarker
res_Dan2_mtry1=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000,mtry=1),.progress=T)
res_Dan2_mtry2=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000,mtry=2),.progress=T)
#MAE for three biomarker
res_Dan3_mtry1=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000,mtry=1),.progress=T)
res_Dan3_mtry2=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000,mtry=2),.progress=T)
res_Dan3_mtry3=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000,mtry=3),.progress=T)

dfDan=rbind(data.frame(Nvars=1,mtry=1,Dan=unlist(res_Dan1)),
            data.frame(Nvars=2,mtry=1,Dan=unlist(res_Dan2_mtry1)),
            data.frame(Nvars=2,mtry=2,Dan=unlist(res_Dan2_mtry2)),
            data.frame(Nvars=3,mtry=1,Dan=unlist(res_Dan3_mtry1)),
            data.frame(Nvars=3,mtry=2,Dan=unlist(res_Dan3_mtry2)),
            data.frame(Nvars=3,mtry=3,Dan=unlist(res_Dan3_mtry3)))

dfDan_out=dfDan %>% dplyr::group_by(Nvars, mtry) %>% dplyr::summarize(mean(Dan),sd(Dan))

dfDan_out
write.csv(dfDan_out,"Generated_data/MAE_dfDan_out.csv",row.names=F)
#single specific biomarker MAE
dat=dat_Dan %>% dplyr::select(stday2,rbd_igg) %>% mutate_all(as.numeric) %>% filter_all(~.!="n.a.")
train_test_Dan=1:niter %>% purrr::map(~make_data(dat)) 
res_Dan1=train_test_Dan %>% future_map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_RBD_igg_Dan=data.frame(Nvars=1, Dan=unlist(res_Dan1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Dan),sd(Dan)) 

dat=dat_Dan %>% dplyr::select(stday2,rbd_iga) %>% mutate_all(as.numeric) %>% filter_all(~.!="n.a.")
train_test_Dan=1:niter %>% purrr::map(~make_data(dat)) 
res_Dan1=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_RBD_iga_Dan=data.frame(Nvars=1, Dan=unlist(res_Dan1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Dan),sd(Dan)) 
#Spike IgG
dat=dat_Dan %>% dplyr::select(stday2,spike_igg) %>% mutate_all(as.numeric) %>% filter_all(~.!="n.a.")
train_test_Dan=1:niter %>% purrr::map(~make_data(dat)) 
res_Dan1=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_S_igg_Dan=data.frame(Nvars=1, Dan=unlist(res_Dan1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Dan),sd(Dan)) 
#Nucleocapsid IgG 
dat=dat_Dan %>% dplyr::select(stday2, np_igg) %>% mutate_all(as.numeric) %>% filter_all(~.!="n.a.")
train_test_Dan=1:niter %>% purrr::map(~make_data(dat)) 
res_Dan1=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_N_igg_Dan=data.frame(Nvars=6, Dan=unlist(res_Dan1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Dan),sd(Dan)) 
#saturated model for MAE
dat=dat_Dan %>% 
  dplyr::select(age,sex,ethnicity,severity,stday2,ends_with('igg'), ends_with('iga'))
train_test_Dan=1:niter %>% purrr::map(~make_data(dat)) 
res_Dan9=train_test_Dan %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=9,ntrees=1000, mtry=9),.progress=T)
MAE_Full_Dan=data.frame(Nvars=9, Dan=unlist(res_Dan9)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Dan),sd(Dan)) 
MAE_Full_Dan
### Peluso data
dat_Peluso=dat_Peluso %>% dplyr::mutate_at(vars(c( "age","sex", "hosp_status", "case_control")), ~as.factor(.)) %>%
  dplyr::mutate_at(vars(c("stday2","N_Abbott", "S_DiaSorin", "N_LIPS","S_LIPS", "N.full._Lum","N.frag._Lum","RBD_Lum", "S_Lum", "N_Roche", 
                          "N_Split_Luc","RBD_Split_Luc","S_Ortho_IgG", "S_Ortho_Ig")),~as.numeric(.))

str(dat_Peluso)
dat_Peluso %>% summarize_all(~sum(is.na(.)))
#Calculate variable importance
dat=dat_Peluso %>% 
  dplyr::select(age, sex, hosp_status,stday2,N_Abbott, S_DiaSorin,ends_with('LIPS'),
                ends_with('Lum'),ends_with('Luc'),S_Ortho_IgG, S_Ortho_Ig)
out=cforest(stday2~.,data=dat,control=cforest_unbiased(ntree=1000,mtry=ceiling((ncol(dat)-1)/3))) 
vimp=permimp(out, conditional = TRUE, asParty = TRUE)
dfPeluso=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% 
  arrange(desc(Cond_Importance))
dfPeluso

# Peluso et al MAE
niter=100
dat_Peluso=dat_Peluso %>% dplyr::mutate_at(vars(c( "age","sex", "hosp_status", "case_control")), ~as.factor(.)) %>%
  dplyr::mutate_at(vars(c("stday2","N_Abbott", "S_DiaSorin", "N_LIPS","S_LIPS", "N.full._Lum","N.frag._Lum","RBD_Lum", "S_Lum", "N_Roche", 
                          "N_Split_Luc","RBD_Split_Luc","S_Ortho_IgG", "S_Ortho_Ig")),~as.numeric(.))

dat=dat_Peluso %>% filter(case_control=="1") %>%  dplyr::select(stday2,N_Abbott, S_DiaSorin, N_LIPS,S_LIPS,N.full._Lum,N.frag._Lum,RBD_Lum, S_Lum, N_Roche, 
                                                                N_Split_Luc,RBD_Split_Luc,S_Ortho_IgG, S_Ortho_Ig)%>% filter_all(~.!="n.a.")
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
res_Peluso1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000, mtry=1),.progress=T)
#MAE for two biomarkers
res_Peluso2_mtry1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=1),.progress=T)
res_Peluso2_mtry2=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=2),.progress=T)
#MAE for three biomarkers
res_Peluso3_mtry1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=1),.progress=T)
res_Peluso3_mtry2=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=2),.progress=T)
res_Peluso3_mtry3=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=3),.progress=T)


dfPeluso=rbind(data.frame(Nvars=1,mtry=1,Peluso=unlist(res_Peluso1)),
               data.frame(Nvars=2,mtry=1,Peluso=unlist(res_Peluso2_mtry1)),
               data.frame(Nvars=2,mtry=2,Peluso=unlist(res_Peluso2_mtry2)),
               data.frame(Nvars=3,mtry=1,Peluso=unlist(res_Peluso3_mtry1)),
               data.frame(Nvars=3,mtry=2,Peluso=unlist(res_Peluso3_mtry2)),
               data.frame(Nvars=3,mtry=3,Peluso=unlist(res_Peluso3_mtry3)))

dfPeluso_out=dfPeluso %>% dplyr::group_by(Nvars, mtry) %>% dplyr::summarize(mean(Peluso),sd(Peluso))

dfPeluso_out
write.csv(dfPeluso_out,"Generated_data/MAE_dfPeluso_out.csv",row.names=F)
#single biomarker 
dat=dat_Peluso %>% filter(case_control=="1") %>%  dplyr::select(stday2,N_Abbott, S_DiaSorin, N_LIPS,S_LIPS,N.full._Lum,N.frag._Lum,RBD_Lum, S_Lum, N_Roche, 
                                                                N_Split_Luc,RBD_Split_Luc,S_Ortho_IgG, S_Ortho_Ig) %>% filter_all(~.!='n.a.')
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
res_Peluso1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_N_Abbot_Peluso=data.frame(Nvars=1,Peluso=unlist(res_Peluso1))%>% dplyr::summarize(mean(Peluso), sd(Peluso))
MAE_N_Abbot_Peluso
#S-ortho-total IgG
dat=dat_Peluso %>% filter(case_control=="1") %>%  dplyr::select(stday2, S_Ortho_Ig) %>% filter_all(~.!='n.a.')
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
res_Peluso1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_S_Ortho_tot_Peluso=data.frame(Nvars=1,Peluso=unlist(res_Peluso1))%>% dplyr::group_by(Nvars) %>% dplyr::summarize(mean(Peluso), sd(Peluso))
MAE_S_Ortho_tot_Peluso
#RBD split Luc (IgG)
dat=dat_Peluso %>% filter(case_control=="1") %>%  dplyr::select(stday2, RBD_Split_Luc) %>% filter_all(~.!='n.a.')
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
res_Peluso1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_R_split_luc_Peluso=data.frame(Nvars=1,Peluso=unlist(res_Peluso1))%>% dplyr::group_by(Nvars) %>% dplyr::summarize(mean(Peluso), sd(Peluso))
MAE_R_split_luc_Peluso
#two nucleocapsid 
dat=dat_Peluso %>% filter(case_control=="1") %>%  dplyr::select(stday2,N_Abbott,  N_LIPS) %>% filter_all(~.!='n.a.')
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
res_Peluso1=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000),.progress=T)
MAE_N_2_Peluso=data.frame(Nvars=2,Peluso=unlist(res_Peluso1))%>% dplyr::group_by(Nvars) %>% dplyr::summarize(mean(Peluso), sd(Peluso))
MAE_N_2_Peluso
#saturated model
dat=dat_Peluso %>% filter(case_control=="1") %>%  dplyr::select(stday2,age, sex,hosp_status, S_DiaSorin, N_LIPS,S_LIPS,N.full._Lum,N.frag._Lum,RBD_Lum, S_Lum, N_Roche, 
                                                                N_Split_Luc,RBD_Split_Luc,S_Ortho_IgG, S_Ortho_Ig) %>% filter_all(~.!='n.a.')
train_test_Peluso=1:niter %>% purrr::map(~make_data(dat)) 
res_Peluso15=train_test_Peluso %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=15,ntrees=1000, mtry=15),.progress=T)
MAE_Full_Peluso=data.frame(Nvars=15,Peluso=unlist(res_Peluso15))%>% dplyr::group_by(Nvars) %>% dplyr::summarize(mean(Peluso), sd(Peluso))
MAE_Full_Peluso
### Whitcombe data
dat_Whitcombe=dat_Whitcombe %>% dplyr::mutate_at(vars(c("age", "sex")), ~as.factor(.)) %>% 
  dplyr::mutate_at(vars(c("stday2", ends_with('IgA'),ends_with('IgG'),ends_with('IgM'))),~as.numeric(.))
dat_Whitcombe %>% summarize_all(~sum(is.na(.)))
str(dat_Whitcombe)
#conditional importance for Whitcombe (time-since-infection)
dat=dat_Whitcombe %>% filter(case_control=='1') %>% dplyr::select(stday2, age, sex, ends_with('IgA'),ends_with('IgG'),ends_with('IgM'))
out=cforest(stday2~.,data=dat,control=cforest_unbiased(ntree=1000,mtry=ceiling((ncol(dat)-1)/3))) 
vimp=permimp(out, conditional = TRUE, asParty = TRUE)
dfWhitcombe=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% 
  arrange(desc(Cond_Importance))
dfWhitcombe

# Whitecombe et al MAE
niter=100
dat=dat_Whitcombe %>% filter(case_control=="1") %>% dplyr::select(stday2,ends_with('IgA'),ends_with('IgG'),ends_with('IgM')) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
# MAE for single biomarker
res_Whitcombe1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000, mtry = 1),.progress=T)
#MAE for two biomarker
res_Whitcombe2_mtry1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=1),.progress=T)
res_Whitcombe2_mtry2=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=2),.progress=T)
#MAE for three biomarker
res_Whitcombe3_mtry1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=1),.progress=T)
res_Whitcombe3_mtry2=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=2),.progress=T)
res_Whitcombe3_mtry3=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=3),.progress=T)


dfWhitcombe=rbind(data.frame(Nvars=1,mtry=1,Whitcombe=unlist(res_Whitcombe1)),
                  data.frame(Nvars=2,mtry=1,Whitcombe=unlist(res_Whitcombe2_mtry1)),
                  data.frame(Nvars=2,mtry=2,Whitcombe=unlist(res_Whitcombe2_mtry2)),
                  data.frame(Nvars=3,mtry=1,Whitcombe=unlist(res_Whitcombe3_mtry1)),
                  data.frame(Nvars=3,mtry=2,Whitcombe=unlist(res_Whitcombe3_mtry2)),
                  data.frame(Nvars=3,mtry=3,Whitcombe=unlist(res_Whitcombe3_mtry3)))

dfWhitcombe_out=dfWhitcombe %>% dplyr::group_by(Nvars, mtry) %>% dplyr::summarize(mean(Whitcombe),sd(Whitcombe))

dfWhitcombe_out
write.csv(dfWhitcombe_out,"Generated_data/MAE_dfWhitcombe_out.csv",row.names=F)
#single biomarker 
#Nucleocapsid IgG
dat=dat_Whitcombe %>% dplyr::select(stday2,NP_IgG) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
res_Whitcombe1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_Whit_N_IgG_out=data.frame(Nvars=1,Whitcombe=unlist(res_Whitcombe1))
MAE_Whit_N_IgG_out=MAE_Whit_N_IgG_out %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Whitcombe),sd(Whitcombe)) 
#RBD IgG
dat=dat_Whitcombe %>% dplyr::select(stday2,RBD_IgG) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
res_Whitcombe1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_Whit_R_IgG_out=data.frame(Nvars=1,Whitcombe=unlist(res_Whitcombe1))
MAE_Whit_R_IgG_out=MAE_Whit_R_IgG_out %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Whitcombe),sd(Whitcombe)) 
#Spike IgG
dat=dat_Whitcombe %>% dplyr::select(stday2,S_IgG) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
res_Whitcombe1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_Whit_S_IgG_out=data.frame(Nvars=1,Whitcombe=unlist(res_Whitcombe1))
MAE_Whit_S_IgG_out=MAE_Whit_S_IgG_out %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Whitcombe),sd(Whitcombe)) 

#Check for best two NP
dat=dat_Whitcombe %>% dplyr::select(stday2,NP_IgA, NP_IgG) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
res_Whitcombe1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_Whit_N_ag_out=data.frame(Nvars=1,Whitcombe=unlist(res_Whitcombe1))
MAE_Whit_N_ag_out=MAE_Whit_N_ag_out %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Whitcombe),sd(Whitcombe)) 

dat=dat_Whitcombe %>% dplyr::select(stday2,NP_IgM, NP_IgG) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
res_Whitcombe1=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_Whit_N_mg_out=data.frame(Nvars=1,Whitcombe=unlist(res_Whitcombe1))
MAE_Whit_N_mg_out=MAE_Whit_N_mg_out %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Whitcombe),sd(Whitcombe)) 
#saturated model
dat=dat_Whitcombe %>% filter(case_control=="1") %>% dplyr::select(stday2,age, sex,ends_with('IgA'),ends_with('IgG'),ends_with('IgM')) %>% filter_all(~.!="n.a.")
train_test_Whitcombe=1:niter %>% purrr::map(~make_data(dat)) 
res_Whitcombe12=train_test_Whitcombe %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=11,ntrees=1000, mtry=12),.progress=T)
MAE_Whit_Full_out=data.frame(Nvars=12,Whitcombe=unlist(res_Whitcombe12))
MAE_Whit_Full_out=MAE_Whit_Full_out %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Whitcombe),sd(Whitcombe)) 
MAE_Whit_Full_out
#Isho dataset
dat_Isho=dat_Isho %>% dplyr::mutate_at(vars(c( "age","sex", "case_control")), ~as.factor(.)) %>%
  dplyr::mutate_at(vars(c("stday2",ends_with ('iga'), ends_with('igg'), ends_with('igm'))),~as.numeric(.))

##Variable importance 
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,age,sex,spike_rob_igg,spike_rob_igm,spike_rob_iga,rbd_rob_igg,rbd_rob_igm,rbd_rob_iga,np_rob_igg,np_rob_igm,np_rob_iga) %>% filter_all(~.!='n.a.')
out=cforest(stday2~.,data=dat,control=cforest_unbiased(ntree=1000,mtry=ceiling((ncol(dat)-1)/3)))
vimp=permimp(out, conditional = TRUE, asParty = TRUE)
dfIsho=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% arrange(desc(Cond_Importance))
dfIsho

###Calculate the MAE for time-since-infection
niter=100
dat=dat_Isho %>% filter (case_control=='1') %>% select(stday2, spike_rob_igg,spike_rob_igm,spike_rob_iga,rbd_rob_igg,rbd_rob_igm,rbd_rob_iga,np_rob_igg,np_rob_igm,np_rob_iga)%>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
# MAE for single biomarker
res_Isho1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000, mtry = 1),.progress=T)
#MAE for two biomarker
res_Isho2_mtry1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=1),.progress=T)
res_Isho2_mtry2=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=2),.progress=T)
#MAE for three biomarker
res_Isho3_mtry1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=1),.progress=T)
res_Isho3_mtry2=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=2),.progress=T)
res_Isho3_mtry3=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=3),.progress=T)


dfIsho=rbind(data.frame(Nvars=1,mtry=1,Isho=unlist(res_Isho1)),
             data.frame(Nvars=2,mtry=1,Isho=unlist(res_Isho2_mtry1)),
             data.frame(Nvars=2,mtry=2,Isho=unlist(res_Isho2_mtry2)),
             data.frame(Nvars=3,mtry=1,Isho=unlist(res_Isho3_mtry1)),
             data.frame(Nvars=3,mtry=2,Isho=unlist(res_Isho3_mtry2)),
             data.frame(Nvars=3,mtry=3,Isho=unlist(res_Isho3_mtry3)))

dfIsho_out=dfIsho %>% dplyr::group_by(Nvars, mtry) %>% dplyr::summarize(mean(Isho),sd(Isho))

dfIsho_out
write.csv(dfIsho_out,"Generated_data/MAE_dfIsho_out.csv",row.names=F)
#MAE for time-since-infection
#Specific single marker
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,np_rob_igg) %>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
res_Isho1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_np_igg=data.frame(Nvars=2, Isho=unlist(res_Isho1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Isho),sd(Isho)) 
# RBD IgG
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,rbd_rob_igg) %>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
res_Isho1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_rbd_igg=data.frame(Nvars=2, Isho=unlist(res_Isho1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Isho),sd(Isho)) 
#Spike IgG
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,spike_rob_igg) %>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
res_Isho1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_spike_igg=data.frame(Nvars=2, Isho=unlist(res_Isho1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Isho),sd(Isho)) 
#Check for best two np marker
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,np_rob_igg, np_rob_igm) %>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
res_Isho1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000),.progress=T)
MAE_np_gm=data.frame(Nvars=2, Isho=unlist(res_Isho1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Isho),sd(Isho)) 
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,np_rob_igg, np_rob_iga) %>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
res_Isho1=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000),.progress=T)
MAE_np_ga=data.frame(Nvars=2, Isho=unlist(res_Isho1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Isho),sd(Isho)) 
#saturated model 
dat=dat_Isho %>% filter(case_control=='1') %>% dplyr::select(stday2,age, sex,spike_rob_igg,spike_rob_igm,spike_rob_iga,rbd_rob_igg,rbd_rob_igm,rbd_rob_iga,np_rob_igg,np_rob_igm,np_rob_iga) %>% filter_all(~.!='n.a.')
train_test_Isho=1:niter %>% purrr::map(~make_data(dat)) 
res_Isho11=train_test_Isho %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=11,ntrees=1000, mtry=11),.progress=T)
MAE_Full_Isho=data.frame(Nvars=11, Isho=unlist(res_Isho11)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Isho),sd(Isho)) 
MAE_Full_Isho

#Markmann data
##Variable importance 
dat_Markmann=dat_Markmann %>% dplyr::mutate_at(vars(c("age", "sex", "Ethnicity")), ~as.factor(.)) %>% mutate_at(vars(c("stday2","RBD_.IgM","RBD_IgA", "RBD_tot_Ig", "RBD_IgG","NP_IgG","S1..N.terminal.domain.Ig.P.N")), ~as.numeric(.)) %>% filter_all(~.!='n.a.')
dat=dat_Markmann %>% dplyr::select(stday2,age,sex,Ethnicity,RBD_.IgM,RBD_IgA, RBD_tot_Ig, RBD_IgG,NP_IgG,S1..N.terminal.domain.Ig.P.N)  %>% filter_all(~.!='n.a.')
out=cforest(stday2~.,data=dat,control=cforest_unbiased(ntree=1000,mtry=ceiling((ncol(dat)-1)/3)))
vimp=permimp(out, conditional = TRUE, asParty = TRUE)
dfMarkmann=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% arrange(desc(Cond_Importance))
dfMarkmann

# MAE for time-since-infection
dat=dat_Markmann %>% filter(case_control=='1') %>% dplyr::select(stday2,RBD_.IgM,RBD_IgA, RBD_tot_Ig, RBD_IgG,NP_IgG,S1..N.terminal.domain.Ig.P.N)  %>% filter_all(~.!='n.a.')
niter=100
train_test_Markmann=1:niter %>% purrr::map(~make_data(dat)) 
# MAE for single biomarker
res_Markmann1=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000, mtry = 1),.progress=T)
#MAE for two biomarker
res_Markmann2_mtry1=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=1),.progress=T)
res_Markmann2_mtry2=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=2,ntrees=1000, mtry=2),.progress=T)
#MAE for three biomarker
res_Markmann3_mtry1=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=1),.progress=T)
res_Markmann3_mtry2=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=2),.progress=T)
res_Markmann3_mtry3=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=3,ntrees=1000, mtry=3),.progress=T)


dfMarkmann=rbind(data.frame(Nvars=1,mtry=1,Markmann=unlist(res_Markmann1)),
                 data.frame(Nvars=2,mtry=1,Markmann=unlist(res_Markmann2_mtry1)),
                 data.frame(Nvars=2,mtry=2,Markmann=unlist(res_Markmann2_mtry2)),
                 data.frame(Nvars=3,mtry=1,Markmann=unlist(res_Markmann3_mtry1)),
                 data.frame(Nvars=3,mtry=2,Markmann=unlist(res_Markmann3_mtry2)),
                 data.frame(Nvars=3,mtry=3,Markmann=unlist(res_Markmann3_mtry3)))

dfMarkmann_out=dfMarkmann %>% dplyr::group_by(Nvars, mtry) %>% dplyr::summarize(mean(Markmann),sd(Markmann))

dfMarkmann_out
write.csv(dfMarkmann_out,"Generated_data/MAE_dfMarkmann_out.csv",row.names=F)
#Single biomarker MAE
dat=dat_Markmann %>% filter(case_control=='1') %>% dplyr::select(stday2,RBD_.IgM) %>% filter_all(~.!='n.a.')
train_test_Markmann=1:niter %>% purrr::map(~make_data(dat)) 
res_Markmann1=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_R_igm_Mark=data.frame(Nvars=1, Markmann=unlist(res_Markmann1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Markmann),sd(Markmann)) 

dat=dat_Markmann %>% filter(case_control=='1') %>% dplyr::select(stday2,NP_IgG) %>% filter_all(~.!='n.a.')
train_test_Markmann=1:niter %>% purrr::map(~make_data(dat)) 
res_Markmann1=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=1,ntrees=1000),.progress=T)
MAE_N_igg_Mark=data.frame(Nvars=1, Markmann=unlist(res_Markmann1)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Markmann),sd(Markmann)) 
#Saturated model 
dat=dat_Markmann %>% filter(case_control=='1') %>% dplyr::select(stday2,age, sex, Ethnicity,RBD_.IgM,RBD_IgA, RBD_tot_Ig, RBD_IgG,NP_IgG,S1..N.terminal.domain.Ig.P.N)  %>% filter_all(~.!='n.a.')
train_test_Markmann=1:niter %>% purrr::map(~make_data(dat)) 
res_Markmann6=train_test_Markmann %>% purrr::map(~Cross_val1(.[[1]],.[[2]],nvars=9,ntrees=1000, mtry=6),.progress=T)
MAE_Full_Mark=data.frame(Nvars=6, Markmann=unlist(res_Markmann6)) %>% dplyr::group_by(Nvars) %>%  dplyr::summarize(mean(Markmann),sd(Markmann)) 
MAE_Full_Mark
