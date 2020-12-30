dataset=read.csv('/Users/ruoyizhu/Documents/MGRM Proj/PROMIS.csv')
Sex=dataset$Sex
Race=dataset$RACE2
Edu=dataset$EDUCAT2
Age=dataset$AGE_CAT
Lan=dataset$Survey_language
Resp1 <- dataset[,8:17]
Resp2 <- dataset[,18:28]
Resp=cbind(Resp1,Resp2)

#delete all missing (no need for lordif, required by mirt)
dataset1=na.omit(dataset)
summary(dataset1)
Sex=dataset1$Sex
Race=dataset1$RACE2
Edu=dataset1$EDUCAT2
Age=dataset1$AGE_CAT
table(Age)
Lan=dataset1$Survey_language
table(Lan)
Resp1 <- dataset1[,c(6,8:17)]
Resp2 <- dataset1[,18:28]
Resp=cbind(Resp1,Resp2)
class(Resp)

table(Resp$var1_dep)
Resp[which(Resp$var1_dep==2),"var1_dep"]=Resp[which(Resp$var1_dep==2),"var1_dep"]-1
Resp[which(Resp$var1_dep==3),"var1_dep"]=Resp[which(Resp$var1_dep==3),"var1_dep"]-2
Resp[which(Resp$var1_dep==4),"var1_dep"]=Resp[which(Resp$var1_dep==4),"var1_dep"]-3
table(Resp$var1_dep)

table(Resp$var2_dep)
Resp[which(Resp$var2_dep==2),"var2_dep"]=Resp[which(Resp$var2_dep==2),"var2_dep"]-1
Resp[which(Resp$var2_dep==3),"var2_dep"]=Resp[which(Resp$var2_dep==3),"var2_dep"]-2
Resp[which(Resp$var2_dep==4),"var2_dep"]=Resp[which(Resp$var2_dep==4),"var2_dep"]-3
table(Resp$var2_dep)

table(Resp$var3_dep)
Resp[which(Resp$var3_dep==2),"var3_dep"]=Resp[which(Resp$var3_dep==2),"var3_dep"]-1
Resp[which(Resp$var3_dep==3),"var3_dep"]=Resp[which(Resp$var3_dep==3),"var3_dep"]-2
Resp[which(Resp$var3_dep==4),"var3_dep"]=Resp[which(Resp$var3_dep==4),"var3_dep"]-3
table(Resp$var3_dep)

table(Resp$var4_dep)
Resp[which(Resp$var4_dep==2),"var4_dep"]=Resp[which(Resp$var4_dep==2),"var4_dep"]-1
Resp[which(Resp$var4_dep==3),"var4_dep"]=Resp[which(Resp$var4_dep==3),"var4_dep"]-2
Resp[which(Resp$var4_dep==4),"var4_dep"]=Resp[which(Resp$var4_dep==4),"var4_dep"]-3
table(Resp$var4_dep)

table(Resp$var5_dep)
Resp[which(Resp$var5_dep==2),"var5_dep"]=Resp[which(Resp$var5_dep==2),"var5_dep"]-1
Resp[which(Resp$var5_dep==3),"var5_dep"]=Resp[which(Resp$var5_dep==3),"var5_dep"]-2
Resp[which(Resp$var5_dep==4),"var5_dep"]=Resp[which(Resp$var5_dep==4),"var5_dep"]-3
table(Resp$var5_dep)

table(Resp$var6_dep)
Resp[which(Resp$var6_dep==2),"var6_dep"]=Resp[which(Resp$var6_dep==2),"var6_dep"]-1
Resp[which(Resp$var6_dep==3),"var6_dep"]=Resp[which(Resp$var6_dep==3),"var6_dep"]-2
Resp[which(Resp$var6_dep==4),"var6_dep"]=Resp[which(Resp$var6_dep==4),"var6_dep"]-3
table(Resp$var6_dep)

table(Resp$var7_dep)
Resp[which(Resp$var7_dep==2),"var7_dep"]=Resp[which(Resp$var7_dep==2),"var7_dep"]-1
Resp[which(Resp$var7_dep==3),"var7_dep"]=Resp[which(Resp$var7_dep==3),"var7_dep"]-2
Resp[which(Resp$var7_dep==4),"var7_dep"]=Resp[which(Resp$var7_dep==4),"var7_dep"]-3
table(Resp$var7_dep)

table(Resp$var8_dep)
Resp[which(Resp$var8_dep==2),"var8_dep"]=Resp[which(Resp$var8_dep==2),"var8_dep"]-1
Resp[which(Resp$var8_dep==3),"var8_dep"]=Resp[which(Resp$var8_dep==3),"var8_dep"]-2
Resp[which(Resp$var8_dep==4),"var8_dep"]=Resp[which(Resp$var8_dep==4),"var8_dep"]-3
table(Resp$var8_dep)

table(Resp$var9_dep)
Resp[which(Resp$var9_dep==2),"var9_dep"]=Resp[which(Resp$var9_dep==2),"var9_dep"]-1
Resp[which(Resp$var9_dep==3),"var9_dep"]=Resp[which(Resp$var9_dep==3),"var9_dep"]-2
Resp[which(Resp$var9_dep==4),"var9_dep"]=Resp[which(Resp$var9_dep==4),"var9_dep"]-3
table(Resp$var9_dep)

table(Resp$var10_dep)
Resp[which(Resp$var10_dep==2),"var10_dep"]=Resp[which(Resp$var10_dep==2),"var10_dep"]-1
Resp[which(Resp$var10_dep==3),"var10_dep"]=Resp[which(Resp$var10_dep==3),"var10_dep"]-2
Resp[which(Resp$var10_dep==4),"var10_dep"]=Resp[which(Resp$var10_dep==4),"var10_dep"]-3
table(Resp$var10_dep)


table(Resp$var1_anx)
Resp[which(Resp$var1_anx==2),"var1_anx"]=Resp[which(Resp$var1_anx==2),"var1_anx"]-1
Resp[which(Resp$var1_anx==3),"var1_anx"]=Resp[which(Resp$var1_anx==3),"var1_anx"]-2
Resp[which(Resp$var1_anx==4),"var1_anx"]=Resp[which(Resp$var1_anx==4),"var1_anx"]-3
table(Resp$var1_anx)

table(Resp$var2_anx)
Resp[which(Resp$var2_anx==2),"var2_anx"]=Resp[which(Resp$var2_anx==2),"var2_anx"]-1
Resp[which(Resp$var2_anx==3),"var2_anx"]=Resp[which(Resp$var2_anx==3),"var2_anx"]-2
Resp[which(Resp$var2_anx==4),"var2_anx"]=Resp[which(Resp$var2_anx==4),"var2_anx"]-3
table(Resp$var2_anx)

table(Resp$var3_anx)
Resp[which(Resp$var3_anx==2),"var3_anx"]=Resp[which(Resp$var3_anx==2),"var3_anx"]-1
Resp[which(Resp$var3_anx==3),"var3_anx"]=Resp[which(Resp$var3_anx==3),"var3_anx"]-2
Resp[which(Resp$var3_anx==4),"var3_anx"]=Resp[which(Resp$var3_anx==4),"var3_anx"]-3
table(Resp$var3_anx)

table(Resp$var4_anx)
Resp[which(Resp$var4_anx==2),"var4_anx"]=Resp[which(Resp$var4_anx==2),"var4_anx"]-1
Resp[which(Resp$var4_anx==3),"var4_anx"]=Resp[which(Resp$var4_anx==3),"var4_anx"]-2
Resp[which(Resp$var4_anx==4),"var4_anx"]=Resp[which(Resp$var4_anx==4),"var4_anx"]-3
table(Resp$var4_anx)

table(Resp$var5_anx)
Resp[which(Resp$var5_anx==2),"var5_anx"]=Resp[which(Resp$var5_anx==2),"var5_anx"]-1
Resp[which(Resp$var5_anx==3),"var5_anx"]=Resp[which(Resp$var5_anx==3),"var5_anx"]-2
Resp[which(Resp$var5_anx==4),"var5_anx"]=Resp[which(Resp$var5_anx==4),"var5_anx"]-3
table(Resp$var5_anx)

table(Resp$var6_anx)
Resp[which(Resp$var6_anx==2),"var6_anx"]=Resp[which(Resp$var6_anx==2),"var6_anx"]-1
Resp[which(Resp$var6_anx==3),"var6_anx"]=Resp[which(Resp$var6_anx==3),"var6_anx"]-2
Resp[which(Resp$var6_anx==4),"var6_anx"]=Resp[which(Resp$var6_anx==4),"var6_anx"]-3
table(Resp$var6_anx)

table(Resp$var7_anx)
Resp[which(Resp$var7_anx==2),"var7_anx"]=Resp[which(Resp$var7_anx==2),"var7_anx"]-1
Resp[which(Resp$var7_anx==3),"var7_anx"]=Resp[which(Resp$var7_anx==3),"var7_anx"]-2
Resp[which(Resp$var7_anx==4),"var7_anx"]=Resp[which(Resp$var7_anx==4),"var7_anx"]-3
table(Resp$var7_anx)

table(Resp$var8_anx)
Resp[which(Resp$var8_anx==2),"var8_anx"]=Resp[which(Resp$var8_anx==2),"var8_anx"]-1
Resp[which(Resp$var8_anx==3),"var8_anx"]=Resp[which(Resp$var8_anx==3),"var8_anx"]-2
Resp[which(Resp$var8_anx==4),"var8_anx"]=Resp[which(Resp$var8_anx==4),"var8_anx"]-3
table(Resp$var8_anx)

table(Resp$var9_anx)
Resp[which(Resp$var9_anx==2),"var9_anx"]=Resp[which(Resp$var9_anx==2),"var9_anx"]-1
Resp[which(Resp$var9_anx==3),"var9_anx"]=Resp[which(Resp$var9_anx==3),"var9_anx"]-2
Resp[which(Resp$var9_anx==4),"var9_anx"]=Resp[which(Resp$var9_anx==4),"var9_anx"]-3
table(Resp$var9_anx)

table(Resp$var10_anx)
Resp[which(Resp$var10_anx==2),"var10_anx"]=Resp[which(Resp$var10_anx==2),"var10_anx"]-1
Resp[which(Resp$var10_anx==3),"var10_anx"]=Resp[which(Resp$var10_anx==3),"var10_anx"]-2
Resp[which(Resp$var10_anx==4),"var10_anx"]=Resp[which(Resp$var10_anx==4),"var10_anx"]-3
table(Resp$var10_anx)

table(Resp$var11_anx)
Resp[which(Resp$var11_anx==2),"var11_anx"]=Resp[which(Resp$var11_anx==2),"var11_anx"]-1
Resp[which(Resp$var11_anx==3),"var11_anx"]=Resp[which(Resp$var11_anx==3),"var11_anx"]-2
Resp[which(Resp$var11_anx==4),"var11_anx"]=Resp[which(Resp$var11_anx==4),"var11_anx"]-3
table(Resp$var11_anx)

Resp
Resp_ordered <- Resp[order(Resp$AGE_CAT),] 
Resp_ordered2<- Resp_ordered[,-1]
table(Resp$AGE_CAT)


J=21

N1=1143
N2=1935
N3=2141
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
N=N1+N2+N3

m=2
r=2
y1=c(0,0)
y2=c(1,0)
y3=c(0,1)

#starting values estimated by multiplegroup
s <- 'D1 = 1-10
          D2 = 11-21
          COV = D1*D2'
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
#omnibus
md.noncons0 <- multipleGroup(Resp_ordered2, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes'))
coef(md.noncons0,simplify=T)$G1$items[,3]
(coef(md.noncons0,simplify=T)$G2$items-coef(md.noncons0,simplify=T)$G1$items)[,3]
(coef(md.noncons0,simplify=T)$G3$items-coef(md.noncons0,simplify=T)$G1$items)[,3]

A11=c(3.74, 4.78, 4.37, 5.02, 4.70, 5.08, 5.22, 6.63, 4.22, 3.51,rep(0,11))
A12=A11
A21=c(rep(0,10),3.17, 4.08, 4.32, 5.12, 4.82, 5.45, 4.75, 4.49, 4.36, 3.91, 4.04)
A22=A21
Amat1=Amat2=cbind(A11,A21)

# MLE starting values
gra00=Amat1
grd00=matrix(c(-0.6373672, -0.9818470, -0.1000673,  3.0270722, -0.8048397,  1.4439140,  2.3179328, -0.5755666,  0.9272971,  0.1702086,  0.9673069,  1.8925296,  3.5242466, -0.1254136 ,1.3748429,  1.5009523,  1.8913342,  0.5744377, -0.6449433,  1.5288644, -0.4072117 ),J,1)
grbeta00=matrix(c(-0.96290193, -1.08642646, -1.48517328, -1.83900230, -1.55070890, -1.53020888, -1.47102188, -2.41962502, -1.36065259, -0.98405835, -0.57950287, -0.54938674, -0.52748067, -0.32662476, -0.20029678, -0.04582263, -0.27723828, -0.60073243, -0.26727144, -0.28619881, -0.02209986,-0.158152324,  0.007947862, -0.640897426, -1.384699620, -0.820888148, -1.076224065, -0.936424413, -1.320135561, -0.655502818, -0.545332623, -1.166835401, -1.009112450,-1.275761280, -0.794785532, -0.617675917, -0.513221353, -0.707333325, -1.205034151, -0.820128767, -0.530055644, -0.293250002  ),J,2)
mu100=c(0,0)
mu200=c(-0.14208419, -0.06937317)
mu300=c(-0.34790747, -0.29761618)
Sig100=matrix(c(1,0.9,0.9,1),2,2)
Sig200=matrix(c(1.109854,0.985275,0.985275,1.125161),2,2)
Sig300=matrix(c(1.055942,0.939480,0.939480,1.084270),2,2)
  
r=2
m=2
y=3
eta.vec=seq(1,20,1)
bics=rep(0,length(eta.vec))
ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
#Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
biass=matrix(0,length(eta.vec),3)
RMSEs=matrix(0,length(eta.vec),3)
theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
for (k in 1:length(eta.vec))
{
  eta=eta.vec[k]
  ptm <- proc.time()
  sim=Reg_DIF(resp=Resp_ordered2,m=2,r=2,y=3,N.vec=c(N1,N2,N3),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
  print(proc.time() - ptm)
  bics[k]=sim$bic
  #Gammas[,,,k]=sim$Gamma
  ADmat[,,k]=sim$est
  Betas[,,k]=sim$Beta
  theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
}

kk=which.min(bics[7:14])

eta.2[rep]=eta.vec[kk]
#Gammas.13[,,,i]=Gammas[,,,kk]
ADmat.2[,,rep]=ADmat[,,kk]
Betas.2[,,rep]=Betas[,,kk]
biass.2[rep,]=biass[kk,]
RMSEs.2[rep,]=RMSEs[kk,]
print(ADmat.2[,,rep])
print(eta.2[rep])
print(Betas.2[,,rep])
print(biass.2[rep,])
print(RMSEs.2[rep,])
write.csv(eta.2[rep],file = paste("eta2EM_",rep))
write.csv(ADmat.2[,,rep],file = paste("ADmat2EM_",rep))
write.csv(Betas.2[,,rep],file = paste("Beta2EM_",rep))
write.csv(theta.dist[,,kk],file = paste("theta2EM_",rep))


for (k in 1:length(eta.vec))
{
  eta=eta.vec[k]
  ptm <- proc.time()
  sim=Reg_EMM_DIF(resp=Resp_ordered2,m=2,r=2,y=3,N.vec=c(N1,N2,N3),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
  print(proc.time() - ptm)
  bics[k]=sim$bic
  #Gammas[,,,k]=sim$Gamma
  ADmat[,,k]=sim$est
  Betas[,,k]=sim$Beta
  theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
}


for (k in 1:length(eta.vec))
{
  eta=eta.vec[k]
  ptm <- proc.time()
  sim=Reg_Adaptive_DIF(resp=Resp_ordered2,m=2,r=2,y=3,N.vec=c(N1,N2,N3),eta=eta,lam=1,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
  print(proc.time() - ptm)
  bics[k]=sim$bic
  #Gammas[,,,k]=sim$Gamma
  ADmat[,,k]=sim$est
  Betas[,,k]=sim$Beta
  theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
}

