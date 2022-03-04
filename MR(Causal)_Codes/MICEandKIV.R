setwd("C:/Users/main/Desktop/FINAL/brain")
#Import libraries
library(dplyr)
library(data.table)
library(nlmr)
library(mice)

#Imputation of missing values (metabolomics)
a<-fread('meta_blood_bio.txt')
a<-as.data.frame(a)
blood<-a[c(1:14,184:245)]
meta<-a[c(1:14)]  
b<-fread('249.txt',fill=T)  
b<-as.data.frame(b)
meta<-left_join(meta,b,by=c('IID'='eid'))
na_count=rep(0,nrow(meta))
for(i in 1:nrow(meta)){
  na_count[i]=sum(is.na(meta[i,]))
}
meta_full<-meta[which(na_count!=249),]
colnames(meta_full)
mat<-meta_full[15:263]
mat.imp<-mice(mat)
mat.comp<-complete(mat.imp,3)
mat_fin<-cbind(meta_full[14],mat.comp)
write.table(meta_full,file='meta_full.txt',row.names = F,quote=F)
write.table(mat_fin,file='meta_fin.txt',row.names = F,quote=F)
a=fread('meta_fin.txt')

#Imputation of missing values (blood)
colnames(blood2)
blood2<-blood[15:ncol(blood)]
blood2.imp<-mice(blood2)
blood2.comp<-complete(blood2.imp,3)
head(blood)
blood0<-blood[c(1:14)]
blood3<-cbind(blood0,blood2.comp)
meta2<-meta[complete.cases(meta),]
head(blood)
write.table(blood3,file='FINAL_blood_filled.txt',row.names = F,quote=F)
write.table(meta2,file='FINAL_meta.txt',row.names = F,quote=F)

a<-fread('FINAL_blood_filled.txt')
a<-as.data.frame(a)
a2=a[2:14]
mod=lm(delta_no_INT~.,data=a2)
a$resid=mod$residuals
a=a[-c(2:13)]
write.table(a,file='blood_for_MR.txt',row.names=F,quote=F)


#Kernel IV with Lymphocyte percentage
kiv<-fread('lymp_per_KIV.txt')
b<-fread('FINAL_blood_filled.txt')
b<-b[1:nrow(b),c(1,63)]
kiv<-left_join(kiv,b,by=c('name'='IID'))
kiv<-kiv[1:nrow(kiv),c(1,2,3,5)]
colnames(kiv)<-c('z','id','y','x')
sum(is.na(kiv))
kiv<-kiv[complete.cases(kiv),]

ex<-kiv[sample(1:34129,2000),]
K_z=matrix(0,2000,2000)
K_x=matrix(0,2000,2000)

sig_x<-var(ex$x)
sig_z<-var(ex$z)
for(i in 1:2000){
  for(j in 1:2000){
    K_x[i,j]<-exp(-((ex$x[i]-ex$x[j])^2)/(2))
    K_z[i,j]<-exp(-((ex$z[i]-ex$z[j])^2)/(2))
    #K_Zz[i,j]<-exp(-((ex2$z[i]-ex$z[j])^2)/(2))
  }
}

W=K_x%*%solve(K_z+2000*0.01*diag(2000))%*%K_z
alpha=solve(W%*%t(W)+2000*0.01*K_x+0.00001*diag(2000))%*%W%*%matrix(ex$y)
dim(alpha)

K_Xx=matrix(0,2000,1000)
test<-seq(from=1,to=100.9,by=0.1)
length(test)
for(i in 1:2000){
  for(j in 1:1000){
    K_Xx[i,j]=exp(-((ex$x[i]-test[j])^2)/(2))
  }
}
res<-t(alpha)%*%K_Xx
plot(test,res,type='l')

#blood 따로 meta 따로 

a<-fread('test_23545exp.txt')
ex<-a[sample(1:8486,2000),]
K_z=matrix(0,2000,2000)
K_x=matrix(0,2000,2000)
sig_x<-var(ex$X)
sig_z<-var(ex$Z)
for(i in 1:2000){
  for(j in 1:2000){
    K_x[i,j]<-exp(-((ex$X[i]-ex$X[j])^2)/(2))
    K_z[i,j]<-exp(-((ex$Z[i]-ex$Z[j])^2)/(2))
    #K_Zz[i,j]<-exp(-((ex2$z[i]-ex$z[j])^2)/(2))
  }
}

W=K_x%*%solve(K_z+2000*0.01*diag(2000))%*%K_z
alpha=solve(W%*%t(W)+2000*0.01*K_x+0.00001*diag(2000))%*%W%*%matrix(ex$Y)

K_Xx=matrix(0,2000,2000)
test<-seq(from=min(ex$X),to=max(ex$X),length.out=2000)
for(i in 1:2000){
  for(j in 1:2000){
    K_Xx[i,j]=exp(-((ex$X[i]-test[j])^2)/(2))
  }
}
res<-t(alpha)%*%K_Xx
plot(test,res,type='l',main='Total Lipids in Small LDL', xlab='exposure', ylab='delta age')


#Nonlinear MR results
b<-read.csv('NONLINEAR_MR.csv')

list<-c(30660,30790,23610, 23631, 23514, 23466, 23449, 23614, 23436, 23545, 23593, 23467)
tab<-as.data.frame(list)
colnames(tab)[1]<-'code'
tab$p_het<-0
tab$p_het_trend<-0

tab$p_nonlinear_quad<-0
tab$p_nonlinear_Q<-0

for(i in c(1)){
  q_num=10
  db<-fread(paste0('test_',as.character(list[i]),'exp.txt'))
  mod=piecewise_mr(y=db$Y,x=db$X,g=db$Z,fig=T,q=q_num)
  if(all(db$X>1)){
    mod2<-fracpoly_mr(y=db$Y,x=db$X,g=db$Z,fig=T,q=q_num,d='both')
  }else{
    db$X<-db$X+1.00001-min(db$X)
    mod2<-fracpoly_mr(y=db$Y,x=db$X,g=db$Z,fig=T,q=q_num,d='both')
  }
  tab$p_nonlinear_quad[i]<-mod$p_tests[1]
  tab$p_nonlinear_Q[i]<-mod$p_tests[2]
  tab$p_het[i]<-mod2$p_heterogeneity[1]
  tab$p_het_trend[i]<-mod2$p_heterogeneity[2]
  
  if(mod2$p_heterogeneity[1]>0.05 & mod2$p_heterogeneity[2]>0.05){
    print('het')
    print(list[i])
  }
  if(any(mod$p_tests<0.05)){
    print('nonlinear')
    print(list[i])
  }
}
trace(piecewise_mr,edit=T)
dim(db)
View(tab)
write.csv(tab,file='NONLINEAR_MR_RESULT.csv',row.names = F)
mod$p_tests
mod2$p_heterogeneity
summary(mod)


for(i in 1:1){
  q_num=10
  db<-fread(paste0('test_',as.character(list[i]),'exp.txt'))
  mod=piecewise_mr(y=db$Y,x=db$X,g=db$Z,fig=T,q=q_num)
  if(all(db$X>1)){
    mod2<-fracpoly_mr(y=db$Y,x=db$X,g=db$Z,fig=T,q=q_num,d='both')
  }else{
    db$X<-db$X+1.00001-min(db$X)
    mod2<-fracpoly_mr(y=db$Y,x=db$X,g=db$Z,fig=T,q=q_num,d='both')
  }

  if(mod2$p_heterogeneity[1]>0.05 & mod2$p_heterogeneity[2]>0.05){
    print('het')
    print(list[i])
  }
  if(any(mod2$p_tests[1:2]<0.05)){
    print('nonlinear')
    print(list[i])
  }
}
summary(mod)
mod2$p_tests



