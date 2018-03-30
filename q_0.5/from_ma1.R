rm(list=ls())
setwd('/Users/shangshuhan/Downloads/程序与论文/code/Part-I-median/q_0.5')
source("sigma_prepare.R")
#####  main code ######
library(mvtnorm)
library(MASS)
library(SparseM)
library(xtable)
library(HyperbolicDist)
library(drc)

n<-50
T<-10
R<-50
rho<-0.8 #rho的三种可能情况
q<-0.5 #中位数

##真实值
betai1<- 10
betai2<- 20
betai3<- 5
betai4<- -0.5
bb<-c(betai1,betai2,betai3,betai4)

# 变量个数
p_num=length(bb)


#生成两种不同类型的时间数据：x1表示非随机时间点，x2表示随机时间点
#x1=matrix(rep(10*c(1:T/T),n),nc=T,nr=n,byrow=T)
# x2=matrix(10*rbeta(n*T,2,2),nc=T,nr=n)
X=array(NA,c(n,T,R))
Y=array(NA,c(n,T,R))
Epsilon = array(NA,c(n,T,R))

for (h in 1:R){
  ni = rep(T,n)
  ep=NULL
  for(i in 1:n){
    k=ni[i]
    out=diag(k)/2
    out[row(out)==col(out)+1]=rho
    sigma=out+t(out)
    ep=rbind(ep,rmvnorm(1, rep(0,k),sigma))
  }
  Epsilon[,,h]  = matrix(ep,ncol=T,byrow = TRUE)
  X[,,h]=matrix(10*rbeta(n*T,2,2),nc=T,nr=n)
  Y[,,h]=f(X[,,h],bb)+Epsilon[,,h]
}
####################################################
#############beta1 Laplace分布，x1非随机情况下###
beta1<-array(0,dim=c(7,p_num,R))
beta2<-array(0,dim=c(4,p_num,R))

timestart = Sys.time()
for (r in 1:R){  
  
  print(r)
  
  x=X[,,r]
  y=Y[,,r]
  
  
  xx = c(t(x))
  yy = c(t(y))
  
  x1 = exp(xx)
  
  
  fit0<-drm(yy~x1,fct=LL.4())
  beta_new<-fit0$fit$par
  
  beta_n=c(beta_new[2],beta_new[3],log(beta_new[4]),beta_new[1])
  
  
  fit1<-optim(beta_n,sump1)
  beta0.5<-fit1$par
  
  e<-y-f(x,beta0.5) 
  
  St<-cov(e)
  
  Si<-cov(t(e))
  
  ut<-c(1:T)
  for (t in 1:T)
    ut[t]<- mean(abs(e[,t]-mean(e[,t])))
  
  ui<- c(1:n)
  for (i in 1:n)
    ui[i]<- mean(abs(e[i,]-mean(e[i,])))
  ## 7个权重
  
  w1<- matrix(1,n,T)
  
  w2<- matrix(1,n,T)
  for(t in 1:T)
    w2[,t]<- rep(1/ut[t],n)
  
  w3<- matrix(1,n,T)
  for(i in 1:n)
    w3[i,]<- rep(1/ui[i],T)
  
  w4<- matrix(1,n,T)
  for(t in 1:T)
    w4[,t]<- rep(1/sqrt(St[t,t]),n)
  
  w5<- matrix(1,n,T)
  for(i in 1:n)
    w5[i,]<- rep(1/sqrt(Si[i,i]),T)
  
  w6<-1/abs(St)
  
  w7<-abs(ginv(St))
  
  ### first method  
  fit1<-optim(beta_n,sumw1)
  beta1[1,,r]<-fit1$par
  fit2<-optim(beta_n,sumw2)
  beta1[2,,r]<-fit2$par
  fit3<-optim(beta_n,sumw3)
  beta1[3,,r]<-fit3$par
  fit4<-optim(beta_n,sumw4)
  beta1[4,,r]<-fit4$par
  fit5<-optim(beta_n,sumw5)
  beta1[5,,r]<-fit5$par
  fit6<-optim(beta_n,sumw6)
  beta1[6,,r]<-fit6$par
  fit7<-optim(beta_n,sumw7)
  beta1[7,,r]<-fit7$par
  
  
  ### second method ####
  
  id = rep(1:n,rep(T,n)) 
  dx=abs(outer(id,id,'-'))
  
  ####
  beta2[1,,r]=optim(beta1[1,,r],Ut4_cs)$par
  beta2[2,,r]=optim(beta1[1,,r],Ut4_ar1)$par
  beta2[3,,r]=optim(beta1[1,,r],Ut4_ma1)$par
  beta2[4,,r]=optim(beta1[1,,r],Ut4_toep)$par
  
}

timeend = Sys.time()
runningtime = timeend - timestart

## MSE和BIAS
MA1=MAPE1(beta1,7)
rownames(MA1)<-c('w1','w2','w3','w4','w5','w6','w7')
colnames(MA1)<-c('b1','b2','b3','b4')
BA1=BIAS1(beta1,7)
rownames(BA1)<-c('w1','w2','w3','w4','w5','w6','w7')
colnames(BA1)<-c('b1','b2','b3','b4')


MA2=MAPE1(beta2,4)
rownames(MA2)<-c('cs','ar1','ma1','toep')
colnames(MA2)<-c('b1','b2','b3','b4')
BA2=BIAS1(beta2,4)
rownames(BA2)<-c('cs','ar1','ma1','toep')
colnames(BA2)<-c('b1','b2','b3','b4')


print(rbind(BA1,BA2))
print(rbind(MA1,MA2))

save.image('from_ma1.RData')