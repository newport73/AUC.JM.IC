
# AUC.JM.IC

<!-- badges: start -->
<!-- badges: end -->

The goal of AUC.JM.IC is to ...

## Installation

You can install the development version of AUC.JM.IC like so:

``` r
library(devtools)

install.packages("BiocManager")
BiocManager::install("Icens", force = TRUE)

install_github("newport73/AUC.JM.IC") # GitHub에서 패키지 설치
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(AUC.JM.IC)
## basic example code
library(AUC.JM.IC)
## basic example code
nt=nrow(data2)
xxx=cbind(data2$age_init,data2$CEP, data2$male) # covariate vector for Longitudinal submodel
xx=cbind(data2id$age_init, data2id$CEP, data2id$male) ## covariate vector for hazard submodel
yyy=data2$IST ## response varible at longitudinal data
ttt=data2$age ## observation time at longitudinal data
######################################################################################
fit1<-ic_joint(data2,data2id,xx,ttt,yyy,xxx)
gamma=fit1$gamma; betas<-fit1$beta; survt=fit1$time; ee=fit1$ee; clamb=fit1$clamb
## regression coeffcients (gamma, beta) at joint model
nn=nrow(ee)
xp=matrix(0,nn,2)
survp=matrix(0,nn,length(survt))
## joint AUC(s,s+th_v) , BS(s,s+th_v) : three s points(75,80,85), th_v=(3,5,7,9)
s1=75; s2=80; s3=85
ss=c(s1,s2,s3)
th_v=c(3,5,7,9)
aucp_joint=brp_joint<-matrix(0,length(ss), length(th_v))

for(i in 1:nn){
for(j in 1:length(survt)) {
if(data2id$TR[i]<=survt[1]) { xp0=fit1$zmat[i,1]}
else { xp0=fit1$zmat[i,j]}
survp[i,j]=exp(-clamb[j]*exp(c(xx[i,],xp0)%*%gamma))
} }


for(kkk in 1:3){
  s=ss[kkk]
  aucp_joint[kkk,]<-AUC_ic3(s,th_v,data2.id,survp,survt)$auc1[1:4]
  brp_joint[kkk,]<-brier3(data2.id,survt,survp,s,th_v)$bs[1:4] 
  }

auc_jj<-data.frame(round(aucp_joint,4))
rownames(auc_jj)<-c("s=75", "s=80","s=85")

colnames(auc_jj)<-c("AUC(s,s+3)","AUC(s,s+5)","AUC(s,s+7)","AUC(s,s+9)")
auc_jj
bs_jj<-data.frame(round(brp_joint,4))
rownames(bs_jj)<-c("s=75", "s=80","s=85")
colnames(bs_jj)<-c("BS(s,s+3)","BS(s,s+5)","BS(s,s+7)","BS(s,s+9)")
bs_jj
```

