library(tidyverse)
library(party)
library(permimp)
library(furrr)
library(cvAUC)
library(pROC)
plan(multiprocess)
setwd("/Users/asbhuiya/Documents/R/Covid_19_data/Data_last_time_points/last time")
getwd()

#To make training and testing dataset
make_data=function(dat){
  nrows=dim(dat)[1]
  train_idx=sample(1:nrows,round(nrows*.8)) # this is 5 fold cross-validation
  train=dat[train_idx,]
  test=dat[-train_idx,]
  return(list(train,test))
}

#Function for cvAUC from classificaiton
Cross_val=function(train,test,nvars,ntrees=1000){
  out=cforest(case_control~.,data=train,control=cforest_unbiased(ntree=ntrees,mtry=ceiling((ncol(train)-1)/3)))#mtry=ceiling(sqrt(ncol(train)-1))))
  vimp=permimp(out, conditional = TRUE, asParty = TRUE)
  df=data.frame(Biomarker=names(vimp$values),Cond_Importance=as.numeric(vimp$values)) %>% arrange(desc(Cond_Importance))
  mod=cforest(as.formula(paste('case_control~',paste(df$Biomarker[1:nvars],collapse="+"),sep="")),data=train,control=cforest_unbiased(ntree=ntrees,mtry=ceiling(nvars/3)))#ceiling(sqrt(nvars))))#
  return(predict(mod,newdata=test))
}

#Load all data for ROC comparison
# Note: rows are removed to get same length for biomarkers
dat_Peluso <- read.csv("Data/Peluso_et_al_last_time_d30_pos_neg.csv")
dat_Whitcombe <- read.csv("Data/Whitcombe_et_al_last_time_d30_pos_neg.csv")
dat_Markmann<-read.csv("Data/Markmann_et_al_last_time_d30_pos_neg_roc_test.csv")
dat_Isho<- read.csv("Data/Isho_et_al_last_time_d30_pos_neg.csv")

#Isho data set
#for Single biomarker
# RBD IgG
niter=100
dat=dat_Isho %>% dplyr::select(case_control=case_control,rbd_rob_igg) %>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test=1:niter %>% purrr::map(~make_data(dat))
Isho_R_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#Spike IgG
dat<-dat_Isho %>% dplyr::select(case_control=case_control,spike_rob_igg) %>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test<-1:niter %>% purrr::map(~make_data(dat))
Isho_S_res<-train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#NP IgG
dat<-dat_Isho %>% dplyr::select(case_control=case_control,np_rob_igg) %>%
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test<-1:niter %>% purrr::map(~make_data(dat))
Isho_N_res<-train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))

#NP vs Spike
# p-value is summarized (mean, sd)
N_vs_S<-pmap(list(Isho_N_res,Isho_S_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                              pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(N_vs_S))
sd(unlist(N_vs_S))
#for two best biomarker - RBD IgG & Spike IgG
dat<-dat_Isho %>% dplyr::select(case_control=case_control,rbd_rob_igg, spike_rob_igg) %>%  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test<-1:niter %>% purrr::map(~make_data(dat))
Isho_best_2<-train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=2,ntrees=1000))
#Spike vs best two biomarkers
S_vs_2<-pmap(list(Isho_S_res,Isho_best_2,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                               pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(S_vs_2))
sd(unlist(S_vs_2))
# Best three biomarkers
dat<-dat_Isho %>% dplyr::select(case_control=case_control,spike_rob_igm, spike_rob_igg, rbd_rob_igg) %>%  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test<-1:niter %>% purrr::map(~make_data(dat))
Isho_best_3<-train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=3,ntrees=1000))
#Best two vs best three
Isho_2_vs_3<-pmap(list(Isho_best_2,Isho_best_3,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                                     pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(Isho_2_vs_3))
sd(unlist(Isho_2_vs_3))

### Peluso data set
dat<-dat_Peluso %>% dplyr::select(case_control=case_control,N.full._Lum,N.frag._Lum,RBD_Lum, S_Lum) %>%
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")

# RBD IgG
dat=dat_Peluso %>% dplyr::select(case_control=case_control,RBD_Lum) %>%
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test=1:niter %>% purrr::map(~make_data(dat))
Pelu_R_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#Spike IgG
dat=dat_Peluso %>% dplyr::select(case_control=case_control,S_Lum) %>%
  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test=1:niter %>% purrr::map(~make_data(dat))
Pelu_S_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#NP IgG
dat=dat_Peluso %>% dplyr::select(case_control=case_control,N.full._Lum) %>% mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test=1:niter %>% purrr::map(~make_data(dat))
Pelu_N_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#Spike vs RBD
P_R_vs_S<-pmap(list(Pelu_R_res,Pelu_S_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                                pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(P_R_vs_S))
sd(unlist(P_R_vs_S))
#Spike vs NP
P_N_vs_S<-pmap(list(Pelu_N_res,Pelu_S_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                                pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(P_N_vs_S))
sd(unlist(P_N_vs_S))
#RBD vs NP
P_R_vs_N<-pmap(list(Pelu_R_res,Pelu_N_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                                pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(P_R_vs_N))
sd(unlist(P_R_vs_N))
#for two - RBD & Spike
dat=dat_Peluso %>% dplyr::select(case_control=case_control,RBD_Lum, S_Lum) %>%  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test=1:niter %>% purrr::map(~make_data(dat))
Pelu_best_1=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))

Pelu_best_2=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=2,ntrees=1000))
#Spike/RBD vs best two
P_S_vs_2=pmap(list(Pelu_best_1,Pelu_best_2,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(P_S_vs_2))
sd(unlist(P_S_vs_2))

#RBD vs best two
P_R_vs_2=pmap(list(Pelu_R_res,Pelu_best_2,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.vlaue})
mean(unlist(P_R_vs_2))
sd(unlist(P_R_vs_2))
#NP vs best two
P_N_vs_2=pmap(list(Pelu_N_res,Pelu_best_2,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(P_N_vs_2))
sd(unlist(P_N_vs_2))
#Best 3
dat=dat_Peluso %>% dplyr::select(case_control=case_control,RBD_Lum, S_Lum,N.full._Lum) %>%  mutate(case_control=as.numeric(case_control)) %>% filter_all(~.!="n.a.")
train_test=1:niter %>% purrr::map(~make_data(dat))
Pelu_best_3=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=3,ntrees=1000))
#best 2 vs 3
P_2_vs_3<-pmap(list(Pelu_best_2,Pelu_best_3,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                                  pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(P_2_vs_3))
sd(unlist(P_2_vs_3))
####
dat=dat_Whitcombe %>%dplyr::select(case_control=case_control,ends_with("IgA"), ends_with("IgG"), ends_with("IgM")) %>% mutate(case_control=as.numeric(case_control))
#Signle biomarker
#Spike IgG
dat=dat_Whitcombe %>%dplyr::select(case_control=case_control,S_IgG) %>% mutate(case_control=as.numeric(case_control))
train_test=1:niter %>% purrr::map(~make_data(dat))
W_S_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#RBD IgG
dat=dat_Whitcombe %>%dplyr::select(case_control=case_control,RBD_IgG) %>% mutate(case_control=as.numeric(case_control))
train_test=1:niter %>% purrr::map(~make_data(dat))
W_R_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#NP IgG
dat=dat_Whitcombe %>%dplyr::select(case_control=case_control,NP_IgG) %>% mutate(case_control=as.numeric(case_control))
train_test=1:niter %>% purrr::map(~make_data(dat))
W_N_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))

#Spike vs NP
W_S_vs_N<-pmap(list(W_R_res,W_N_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                          pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(W_S_vs_N))
sd(unlist(W_S_vs_N))
#RBD vs NP
W_R_vs_N<-pmap(list(W_R_res,W_N_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                          pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(W_R_vs_N))
sd(unlist(W_R_vs_N))
#best two
dat=dat_Whitcombe %>%dplyr::select(case_control=case_control,NP_IgA, RBD_IgA) %>% mutate(case_control=as.numeric(case_control))
train_test=1:niter %>% purrr::map(~make_data(dat))
W_best_2=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=2,ntrees=1000))
#Best three
dat=dat_Whitcombe %>%dplyr::select(case_control=case_control,NP_IgA, RBD_IgA, S_IgM) %>% mutate(case_control=as.numeric(case_control))
train_test=1:niter %>% purrr::map(~make_data(dat))
W_best_3=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=3,ntrees=1000))

#Spike vs best two
W_S_vs_2<-pmap(list(W_best_2,W_S_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                           pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(W_S_vs_2))
sd(unlist(W_S_vs_2))
#best 2 vs 3
W_2_vs_3<-pmap(list(W_best_2,W_best_3,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                            pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(W_2_vs_3))
sd(unlist(W_2_vs_3))
##Markman
dat=dat_Markmann %>% dplyr::select(case_control=case_control, RBD_.IgM,RBD_IgA, RBD_tot_Ig, RBD_IgG,NP_IgG,S1..N.terminal.domain.Ig.P.N)  %>% filter_all(~.!='n.a.')
#Single biomrker
#RBD IgG
dat=dat_Markmann %>% dplyr::select(case_control=case_control, RBD_IgG)  %>% filter_all(~.!='n.a.')
train_test=1:niter %>% purrr::map(~make_data(dat))
M_R_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#NP
dat=dat_Markmann %>% dplyr::select(case_control=case_control,NP_IgG)  %>% filter_all(~.!='n.a.')
train_test=1:niter %>% purrr::map(~make_data(dat))
M_N_res=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=1,ntrees=1000))
#two biomarker
dat=dat_Markmann %>% dplyr::select(case_control=case_control,RBD_IgG,NP_IgG)  %>% filter_all(~.!='n.a.')
train_test=1:niter %>% purrr::map(~make_data(dat))
M_best_2=train_test %>% purrr::map(~Cross_val(.[[1]],.[[2]],nvars=2,ntrees=1000))
#RBD vs NP
M_R_vs_N<-pmap(list(M_R_res,M_N_res,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                          pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(M_R_vs_N))
sd(unlist(M_R_vs_N))
#RBD vs two best marker
M_S_vs_2<-pmap(list(M_R_res,M_best_2,train_test),function(x,y,z) {roc.test(pROC::roc(z[[2]]$case_control,as.numeric(x)),
                                                                           pROC::roc(z[[2]]$case_control,as.numeric(y)))$p.value})
mean(unlist(M_S_vs_2))
sd(unlist(M_S_vs_2))


