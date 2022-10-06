setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("resp1new_SS600.csv",row.names = 1)
rep=11
N1=N2=N3=200 
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
N=N1+N2+N3
resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #
Unif=T
setwd('/Users/ruoyizhu/Desktop')
resp=data.matrix(read.csv("/Users/ruoyizhu/Desktop/resp.csv"))
Group=read.csv("/Users/ruoyizhu/Desktop/Group.csv")[,1]
indic=data.matrix(read.csv("/Users/ruoyizhu/Desktop/indic.csv"))



write.csv(resp,file = "resp.csv", row.names = F)
write.csv(Group,file = "Group.csv", row.names = F)
write.csv(indic,file = "indic.csv", row.names = F)
COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2) #
write.csv(COV,file = "COV.csv")

result0=LRT_function(resp=resp2,Group=Group,indic=indic,Unif=Unif)
init=DIF_init(resp=resp2,Group=Group,indic=indic,Unif=Unif)
result01=LRT_function(resp=resp2,Group=Group,indic=indic,Unif=F)
result0=LRT_function(resp,Group,indic,Unif=F)

result0=vemm_dif(u=resp,G=G,y=y,domain=domain,new_a0=gra00,new_b0=grd00,new_beta0=grbeta00,new_gamma0=grgamma00,Sigma0=Sigma0,eta=eta,indic=indic,lbd=lbd)




result1=reg_DIF_alllbd(resp=data.matrix(resp),indic=indic,Group=Group,Method='LRT',Unif=T)
data.matrix(resp)

domain=result0$domain
y=result0$y
m<-cbind(result0$Amat,result0$Dmat)
for (r in 1:domain){
  m<-cbind(m,t(result0$Gamma[r,,]))
}
m<-cbind(m,result0$Beta)
gp=NULL
for(yy in 1:(y-1)){
  gp=c(gp,rep(yy,domain))
}
colnames(m)<-c(paste("a",1:(result0$domain),sep=""),"d",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))


domain=result0$domain
y=result0$y
m<-cbind(result0$Amat,result0$Dmat)
if (y==2){
  for (r in 1:domain){
    m<-cbind(m,result0$Gamma[,r,])
  }
} else {
  for (r in 1:domain){
    m<-cbind(m,t(result0$Gamma[,r,]))
  }
}
m<-cbind(m,result0$Beta)

m1<-matrix(result0$Mu)

m2<-result0$Sig
m=cbind(m1,m2)
colnames(m)<-c(paste("Mean"),paste("Var dimension",1:(result0$domain),sep=""))
rownames(m)<-paste("group",rep(1:(result0$y),each=result0$domain)," dimension",rep(1:(result0$domain),result0$y),sep="")



plot(seq(10,20,5),result0$ICs,type = "b",xlab="Tuning parameter",ylab="Information Criteria")



setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("resp1new_SS600.csv",row.names = 1)
rep=1
N1=N2=N3=200 
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
N=N1+N2+N3
resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #
setwd('/Users/ruoyizhu/Desktop/DIF_sims')
write.csv(resp,file = "resp.csv", row.names = F)
write.csv(Group,file = "Group.csv", row.names = F)
write.csv(indic,file = "indic.csv", row.names = F)
COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2) #
write.csv(COV,file = "COV.csv")



# 2-group
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("resp1new_SS600.csv",row.names = 1)
rep=19
N1=N2=N3=200 
Group=c(rep('G1', N1), rep('G2', N2))
N=N1+N2
resp=responses[c(((rep-1)*N+1):((rep-1)*N+N1),((rep-1)*N+N1+N2+1):((rep-1)*N+N1+N2+N3)),]
indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #
setwd('/Users/ruoyizhu/Desktop/DIF_sims')
write.csv(resp,file = "resp2grp3.csv", row.names = F)
write.csv(Group,file = "Group2grp3.csv", row.names = F)
write.csv(indic,file = "indic2grp3.csv", row.names = F)

resp=data.matrix(read.csv("resp2grp3.csv"))
Group=read.csv("Group2grp3.csv")[,1]
indic=data.matrix(read.csv("indic2grp3.csv"))


resp=data.matrix(read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF/DIF_sims/resp2grp.csv"))
Group=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF/DIF_sims/Group2grp.csv")[,1]
indic=data.matrix(read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF/DIF_sims/indic.csv"))


resp=data.matrix(read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF/DIF_sims/NonUnif2grp/NonUnif_resp2grp.csv"))
Group=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF/DIF_sims/NonUnif2grp/Group2grp500.csv")[,1]
indic=data.matrix(read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF/DIF_sims/NonUnif2grp/indic.csv"))


result0=reg_DIF_alllbd(resp=resp,indic=indic,Group=Group,Method='LRT',Unif=T)

result0=LRT_function(resp=resp,indic=indic,Group=Group,Unif=T)

domain=result0$domain
y=result0$y
m<-cbind(result0$Amat,result0$Dmat)
if (y==2){
  for (r in 1:domain){
    m<-cbind(m,result0$Gamma[,r,])
  }
} else {
  for (r in 1:domain){
    m<-cbind(m,t(result0$Gamma[,r,]))
  }
}
m<-cbind(m,result0$Beta)
gp=NULL
#for(yy in 1:(y-1)){
#  gp=c(gp,rep(yy,domain))
#}
#colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))
#rownames(m)<-c(paste("Item",1:ncol(indic),sep=""))
for(r in 1:domain){
  gp=c(gp,rep(r,y-1))
}
colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:(y-1),sep=""),result0()$domain),gp,sep=""),paste("beta",1:(y-1),sep=""))



# 4-group
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
A11=read.csv("Para1.csv",row.names = 1)[,1]
A12=A11;A13=A11;A14=A11
A12[c(4:5)]=A11[c(4:5)]-0.8
A13[c(4:5)]=A11[c(4:5)]-1.1
A21=read.csv("Para1.csv",row.names = 1)[,2]
A22=A21;A23=A21;A24=A21
A22[c(12:13)]=A21[c(12:13)]-0.8
A23[c(12:13)]=A21[c(12:13)]-1.1
D11=read.csv("Para1.csv",row.names = 1)[,7]
D12=D11;D13=D11;D14=D11
D12[c(4:5,12:13)]=D12[c(4:5,12:13)]+0.75
D13[c(4:5,12:13)]=D13[c(4:5,12:13)]+1
D14[c(4:5,12:13)]=D14[c(4:5,12:13)]+0.75

Amat1=cbind(A11,A21);Amat2=cbind(A12,A22);Amat3=cbind(A13,A23);Amat4=cbind(A14,A24)
Dmat1=cbind(D11)
Dmat2=cbind(D12);Dmat3=cbind(D13);Dmat4=cbind(D14)

#write.csv(cbind(Amat1,Amat2,Amat3,Dmat1,Dmat2,Dmat3), file = "Para4new.csv")
N1=N2=N3=N4=1000
N=N1+N2+N3+N4
set.seed(10)
#Theta=mvrnorm(N*50,c(0,0),matrix(c(1,0.25,0.25,1),2,2))
Theta=mvrnorm(N,c(0,0),matrix(c(1,0.25,0.25,1),2,2))
resp=matrix(0,N,J)
for (rep in 1:1){
  Theta1=Theta[((rep-1)*N+1):((rep-1)*N+N1),]
  Theta2=Theta[((rep-1)*N+N1+1):((rep-1)*N+N1+N2),]
  Theta3=Theta[((rep-1)*N+N1+N2+1):((rep-1)*N+N1+N2+N3),]
  #Theta4=Theta[((rep-1)*N+N1+N2+N3+1):(rep*N),]
  datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
  datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
  datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
  #datasetF3=simdata(Amat4,Dmat4,itemtype = "dich",Theta=Theta4)
  #resp[((rep-1)*N+1):(rep*N),]=rbind(datasetR,datasetF1,datasetF2,datasetF3)
  resp[((rep-1)*N+1):(rep*N),]=rbind(datasetR,datasetF1,datasetF2)
}
resp=as.data.frame(resp)
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3), rep('G4', N4))
indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #
setwd('/Users/ruoyizhu/Desktop')
write.csv(resp,file = "NonUnif_resp3.csv", row.names = F)
write.csv(Group,file = "Group4grp.csv", row.names = F)
write.csv(indic,file = "indic4grp.csv", row.names = F)



result0=reg_DIF_alllbd(resp=resp,indic=indic,Group=Group,Method='EM',Unif=T)

m<-cbind(result0$Amat,result0$Dmat)
for (r in 1:domain){
  m<-cbind(m,result0$Gamma[,r,,1])
}
m<-cbind(m,result0$Beta)


#username: uwpmetrics  (or uwpmetrics@gmail.com)
# password: COE_psychometrics21

#remote desktop
# pmetrics@uw.edu
# COE_psychometrics19
#  pin:123456

#Miller321!!








