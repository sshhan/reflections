#四参数的logistic函数
f<- function(x,b){
  result<-matrix(NA,nrow=n,ncol=T)
  for(t in 1:T) result[,t]<-b[1]+(b[2]-b[1])/(1+exp(b[4]*(x[,t]-b[3])))
  return(result)
}
#构造损失函数
I<-function(x){
  ifelse(x<0,1,0)
}
cost<-function(b){
  (y-f(x,b))*(q-I(y-f(x,b)))
}
#######误差和函数#######
sump1<-function(b){
  sum(cost(b))
}
sumw1<-function(b){
  sum(w1*cost(b))
}
sumw2<-function(b){
  sum(w2*cost(b))
}

sumw3<-function(b){
  sum(w3*cost(b))
}

sumw4<-function(b){
  sum(w4*cost(b))
}

sumw5<-function(b){
  sum(w5*cost(b))
}

sumw6<-function(b){
  sum(rep(1,T)%*%w6%*%t(cost(b)))
}
sumw7<-function(b){
  sum(rep(1,T)%*%w7%*%t(cost(b)))
}

################### four parameter Lt function ########################
Lt4<-function(b){
  epsilon=yy-(b[1]+(b[2]-b[1])/(1+exp(b[4]*(xx-b[3]))))
  sum(epsilon*(q-I(epsilon)))
}
################### Ut_cs function - cs 结构 ########################
Ut4_cs<-function(b){
  epsilon=yy-(b[1]+(b[2]-b[1])/(1+exp(b[4]*(xx-b[3]))))
  e=I(epsilon)
  S=q-e
  delta=outer(c(e),c(e),"*")
  a=delta*(dx==0)
  M2=sum (ni*(ni-1))             
  det=(sum(a)-sum(diag(a)))/M2 #\sum ni(ni-1) 等相关的相关系数估计
  r=(det-q^2)/(q-q^2) 
  r=max(0, min(r,0.95)) # 相关系数非负
  
  ut<-array(0,dim=c(4,1))
  for(i in 1:n){
    J=matrix(1,nr=ni[i],nc=ni[i])
    VI=(1/(ni[i]*(1+(ni[i]-1)*r))*J+(diag(ni[i])-1/ni[i]*J)/(1-r))/(q-q^2)
    
    xi=xx[id==i]
    Di<-array(0,dim=c(ni[i],4))
    for (j in 1:ni[i]){
      Di[j,1]= -exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,2]= -1/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,3]=-(b[2]-b[1])*b[4]*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
      Di[j,4]=(b[2]-b[1])*(xi[j]-b[3])*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
    }
    Si=as.matrix(S[id==i],nr=ni[i],nc=1)
    DI=t(Di)%*%VI
    ut=ut+DI%*%Si
  }
  t(ut)%*%ut
}
#### Ut_ar1 - ar1结构 ###########################################
ar1=function(rho,m)
{
  e=1:m
  rho^abs(outer(e,e,"-"))
  
}
Ut4_ar1<-function(b){
  epsilon=yy-(b[1]+(b[2]-b[1])/(1+exp(b[4]*(xx-b[3]))))
  e=I(epsilon)
  
  S=q-e
  delta=outer(c(e),c(e),"*")
  no = rep(c(1:ni[1]),n)
  a2=(outer(no,no,"-"))
  a3=(dx==0)*(a2==1)
  rho=sum(delta*a3)/(sum(diag(delta))+0.5*sum(((e^2*(no==1))+(e^2*(no==T)))))
  
  ut<-array(0,dim=c(4,1))
  for(i in 1:n){
    xi=xx[id==i]
    
    #J=matrix(1,nr=ni[i],nc=ni[i])
    #V = ar1(rho,ni[1])*(q-q^2)
    #VI = solve(V)
    VI = 1/(1-rho^2)*(ma1(-rho,ni[i])+diag(c(0,rep(rho^2,ni[i]-2),0)))/(q-q^2)
    
    Di<-array(0,dim=c(ni[i],4))
    for (j in 1:ni[i]){
      Di[j,1]= -exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,2]= -1/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,3]=-(b[2]-b[1])*b[4]*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
      Di[j,4]=(b[2]-b[1])*(xi[j]-b[3])*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
    }
    Si=as.matrix(S[id==i],nr=ni[i],nc=1)
    DI=t(Di)%*%VI
    ut=ut+DI%*%Si
  }
  t(ut)%*%ut
}
#### Ut4_ma1 - MA(1)结构##############################
ma1=function(rho,m){
  out=diag(m)/2
  #Ctau= pmvnorm(lower=c(-Inf, -Inf),upper=c(qnorm(q),qnorm(q)),corr=matrix(c(1,rho,rho,1),2,2))
  
  #rk=(Ctau-q^2)/(q-q^2)
  
  out[row(out)==col(out)+1]=rho
  out+t(out)
}
Ut4_ma1<-function(b){
  epsilon=yy-(b[1]+(b[2]-b[1])/(1+exp(b[4]*(xx-b[3]))))
  e=I(epsilon)
  
  S=q-e
  delta=outer(c(e),c(e),"*")
  no = rep(c(1:ni[1]),n)
  a2=(outer(no,no,"-"))
  a3=(dx==0)*(a2==1)
  
  rho=sum(delta*a3)/(sum(diag(delta))+0.5*sum(((e^2*(no==1))+(e^2*(no==T)))))
  
  ut<-array(0,dim=c(4,1))
  for(i in 1:n){
    xi=xx[id==i]
    
    
    V = ma1(rho,ni[1])*(q-q^2)
    VI = solve(V)
    
    Di<-array(0,dim=c(ni[i],4))
    for (j in 1:ni[i]){
      Di[j,1]= -exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,2]= -1/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,3]=-(b[2]-b[1])*b[4]*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
      Di[j,4]=(b[2]-b[1])*(xi[j]-b[3])*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
    }
    Si=as.matrix(S[id==i],nr=ni[i],nc=1)
    DI=t(Di)%*%VI
    ut=ut+DI%*%Si
  }
  t(ut)%*%ut
}

#### Ut4_toep - toeplitz结构#############################################
top=function(h_rho){
  
  a=c(1,h_rho)
  toeplitz(a)
  
}
Ut4_toep<-function(b){
  epsilon=yy-(b[1]+(b[2]-b[1])/(1+exp(b[4]*(xx-b[3]))))
  e=I(epsilon)
  S=q-e
  
  S_tilde = S/(q*(1-q))^0.5
  denominator = sum(S_tilde^2)/(n*ni[1])
  
  h_rho = rep(0, ni[1]-1)
  for(k in 1:(ni[1]-1)){
    for(i in 1:n){
      Si_tilde=as.matrix(S_tilde[id==i],nr=ni[i],nc=1)
      for(j in 1:(ni[i]-k)){
        h_rho[k]= h_rho[k]+ Si_tilde[j]*Si_tilde[j+k]
      }
    }
    h_rho[k] = h_rho[k]/(n*(ni[1]-k))
  }
  h_rho = h_rho/denominator
  
  ut<-array(0,dim=c(4,1))
  for(i in 1:n){
    xi=xx[id==i]
    
    V = top(h_rho)*(q-q^2)
    VI = solve(V)
    
    Di<-array(0,dim=c(ni[i],4))
    for (j in 1:ni[i]){
      Di[j,1]= -exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,2]= -1/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,3]=-(b[2]-b[1])*b[4]*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
      Di[j,4]=(b[2]-b[1])*(xi[j]-b[3])*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
    }
    Si=as.matrix(S[id==i],nr=ni[i],nc=1)
    DI=t(Di)%*%VI
    ut=ut+DI%*%Si
  }
  t(ut)%*%ut
}
################### Ut_zi function 估计方差的函数 ########################
Ut4_zi<-function(b){
  
  zi = rexp(n,1)
  epsilon=yy-(b[1]+(b[2]-b[1])/(1+exp(b[4]*(xx-b[3]))))
  e=I(epsilon)
  S=q-e
  delta=outer(e,e,"*")
  a=delta*(dx==0)
  M2=sum (ni*(ni-1))             
  det=(sum(a)-sum(diag(a)))/M2 #\sum ni(ni-1) 等相关的相关系数估计
  r=(det-q^2)/(q-q^2) 
  r=max(0, min(r,0.95)) # 相关系数非负
  
  ut<-array(0,dim=c(4,1))
  for(i in 1:n){
    J=matrix(1,nr=ni[i],nc=ni[i])
    VI=(1/(ni[i]*(1+(ni[i]-1)*r))*J+(diag(ni[i])-1/ni[i]*J)/(1-r))/(q-q^2)
    
    xi=xx[id==i]
    Di<-array(0,dim=c(ni[i],4))
    for (j in 1:ni[i]){
      Di[j,1]= -exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,2]= -1/(1+exp(b[4]*(xi[j]-b[3])))
      Di[j,3]=-(b[2]-b[1])*b[4]*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
      Di[j,4]=(b[2]-b[1])*(xi[j]-b[3])*exp(b[4]*(xi[j]-b[3]))/(1+exp(b[4]*(xi[j]-b[3])))^2
    }
    Si=as.matrix(S[id==i],nr=ni[i],nc=1)
    DI=t(Di)%*%VI*zi[i]
    ut=ut+DI%*%Si
  }
  #ut1 = ut*n^(-0.5)
  t(-ut)%*%(-ut)
}
######## MAPE function - 绝对最小偏差 ###############
MAPE<-function(beta,num_b){
  bb=c(betai1,betai2,betai3,betai4)
  MAPEbeta<-matrix(0,num_b,p_num)
  for (p in 1:num_b)
    for (h in 1:p_num)
    {  MAPEbeta[p,h]<-mean(abs(beta[p,h,]-bb[h])/abs(bb[h]))*100
    }
  return(MAPEbeta)
}

######## BIAS 偏的函数 ############
BIAS<-function(beta,num_b){
  bb=c(betai1,betai2,betai3,betai4)
  biasbeta<-matrix(0,num_b,p_num)
  for (p in 1:num_b)
    for (h in 1:p_num)
    {  biasbeta[p,h]<-(mean(beta[p,h,])-bb[h])/abs(bb[h])*100
    }
  return(biasbeta)
}
######## MAPE1 function - 均方误差 ###############
MAPE1<-function(beta,num_b){
  bb=c(betai1,betai2,betai3,betai4)
  MAPEbeta<-matrix(0,num_b,p_num)
  for (p in 1:num_b)
    for (h in 1:p_num)
    {  MAPEbeta[p,h]<-mean((beta[p,h,]-bb[h])^2)
    }
  return(MAPEbeta)
}

######## BIAS1 偏的函数 ############
BIAS1<-function(beta,num_b){
  bb=c(betai1,betai2,betai3,betai4)
  biasbeta<-matrix(0,num_b,p_num)
  for (p in 1:num_b)
    for (h in 1:p_num)
    {  biasbeta[p,h]<-(mean(beta[p,h,])-bb[h])
    }
  return(biasbeta)
}