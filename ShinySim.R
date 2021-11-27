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
write.csv(resp,file = "resp.csv")
write.csv(Group,file = "Group.csv")
write.csv(indic,file = "indic.csv")
COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2) #
write.csv(COV,file = "COV.csv")


result0=reg_DIF_alllbd(resp=resp,indic=indic,Group=Group,Method='EMM',Unif=T)


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


m1<-matrix(result0$Mu)

m2<-result0$Sig
m=cbind(m1,m2)
colnames(m)<-c(paste("Mean"),paste("Var dimension",1:(result0$domain),sep=""))
rownames(m)<-paste("group",rep(1:(result0$y),each=result0$domain)," dimension",rep(1:(result0$domain),result0$y),sep="")






plot(seq(10,20,5),result0$ICs,type = "b",xlab="Tuning parameter",ylab="Information Criteria")




