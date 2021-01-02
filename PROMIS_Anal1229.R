dataset=read.csv('/Users/ruoyizhu/Documents/MGRM Proj/PROMIS.csv')
dataset=read.csv('/Users/zhux0445/Documents/GitHub/RegDIF/PROMIS copy.csv')
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
coef(md.noncons0,simplify=T)$G1$items[,1]
coef(md.noncons0,simplify=T)$G1$items[,2]
coef(md.noncons0,simplify=T)$G1$items[,3]
(coef(md.noncons0,simplify=T)$G2$items-coef(md.noncons0,simplify=T)$G1$items)[,3]
(coef(md.noncons0,simplify=T)$G3$items-coef(md.noncons0,simplify=T)$G1$items)[,3]

A11=c(3.749509,  4.906052,  4.403274,  5.322981,  4.817137,  5.151824,  5.302439,  6.867881,  4.115629,  3.436797,rep(0,11))
A12=A11
A21=c(rep(0,10),2.985399,  3.880630,  4.295031,  4.929496,  4.552499 , 5.338969 , 4.617473 , 4.214355 , 4.150518 , 3.789304 , 3.849753 )
A22=A21
Amat1=Amat2=cbind(A11,A21)
gra00=Amat1
grd00=matrix(c(-0.6373672, -0.9818470, -0.1000673,  3.0270722, -0.8048397,  1.4439140,  2.3179328, -0.5755666 , 0.9272971 , 0.1702086 , 0.9673069 , 1.8925296 , 3.5242466 ,-0.1254136 , 1.3748429 , 1.5009523,  1.8913342,  0.5744377, -0.6449433,  1.5288644,-0.4072117 ),J,1)
grbeta00=matrix(c(-0.96290193, -1.08642646, -1.48517328, -1.83900230, -1.55070890,-1.53020888, -1.47102188, -2.41962502, -1.36065259, -0.98405835, -0.57950287, -0.54938674, -0.52748067, -0.32662476, -0.20029678, -0.04582263, -0.27723828, -0.60073243, -0.26727144, -0.28619881, -0.02209986, 
                  -0.158152324,  0.007947862, -0.640897426, -1.384699620, -0.820888148, -1.076224065, -0.936424413, -1.320135561, -0.655502818, -0.545332623, -1.166835401, -1.009112450, -1.275761280, -0.794785532, -0.617675917, -0.513221353, -0.707333325, -1.205034151, -0.820128767, -0.530055644, -0.293250002),J,2)
mu100=c(0,0)
mu200=c(0.134, -0.089 )
mu300=c(-0.243, -0.275 )
Sig100=matrix(c(1,0.911,0.911,1),2,2)
Sig200=matrix(c(1.222,1.220,1.220,1.435),2,2)
Sig300=matrix(c(1.069,1.094,1.094,1.314),2,2)

# MLE starting values

A11=c(4.727485, 6.182298, 5.622785, 6.915568, 6.097589, 6.582838, 6.879913, 8.736050, 5.296018, 4.384396,rep(0,11))
A12=A11
A21=c(rep(0,10),3.791942, 4.965282, 5.532676, 6.246125, 5.801095, 6.731871, 5.869306, 5.347072, 5.210306, 4.812828, 4.835288)
A22=A21
Amat1=Amat2=cbind(A11,A21)
gra00=Amat1
grd00=matrix(c(1.428995, 1.760881, 2.220391,  6.021333, 1.773182,  4.383001,  5.384768, 3.093435,  3.224947,  2.086015,  2.527317,  3.943190,  5.925238, 2.534635 ,3.817292,  4.512750,  4.440518,  2.805363, 1.567375,  3.577428, 1.678282 ),J,1)
grbeta00=matrix(c(0.256558933,  0.533222654,  0.063720755, -0.128242789, -0.003074121, -0.106987036,  0.147042629, -0.058750516,  0.037792447,  0.088173378, -0.311043154, -0.056780899, -0.275643371,  0.213744399,  0.375508140 , 0.290335854,  0.196523029, -0.229861046, 0.144789323,  0.095624015,  0.369718651,0.55112795,  0.87749027,  0.37315662, -0.42836366,  0.23498341, -0.20916744, -0.05180618,  0.06780513,  0.02165227,  0.05863596, -0.53335920, -0.17821646,-0.37374232,  0.19790551,  0.36318861,  0.44515429,  0.12494803, -0.26507690,  0.04064621,  0.28452507,  0.50321574),J,2)
mu100=c(0,0)
mu200=c(-0.5722514, -0.5662688)
mu300=c(-0.7567240, -0.7919464)
Sig100=matrix(c(1,0.9385956,0.9385956,1),2,2)
Sig200=matrix(c(0.7030104,0.7033863,0.7033863,0.8323150),2,2)
Sig300=matrix(c(0.6054572,0.6203933,0.6203933,0.7527996),2,2)

r=2
m=2
y=3
eta.vec=seq(1,25,1)
bics=rep(0,length(eta.vec))
ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
#Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
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

kk=which.min(bics)

eta.vec[kk]

ADmat[,,kk]
Betas[,,kk]
theta.dist[,,kk]



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


###########################
####    Non-Uniform    ####
###########################

#starting values estimated by multiplegroup
s <- 'D1 = 1-10
          D2 = 11-21
          COV = D1*D2'
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
#omnibus
md.noncons0 <- multipleGroup(Resp_ordered2, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts'))
coef(md.noncons0,simplify=T)$G1$items[,1]
coef(md.noncons0,simplify=T)$G1$items[,2]
coef(md.noncons0,simplify=T)$G1$items[,3]
(coef(md.noncons0,simplify=T)$G2$items-coef(md.noncons0,simplify=T)$G1$items)[,3]
(coef(md.noncons0,simplify=T)$G3$items-coef(md.noncons0,simplify=T)$G1$items)[,3]

A11=c(3.710687,  4.786595,  4.644714,  4.948962,  5.037059,  4.990929,  5.281461,  6.785962,  4.200693,  3.640178,rep(0,11))
A12=A11
A21=c(rep(0,10),2.608567,  3.756078,  4.166061,  5.791455 , 5.100036,  6.259032,  4.916154,  3.956793 ,4.507497,  4.036874,  4.136797 )
A22=A21
Amat1=Amat2=cbind(A11,A21)

gra00=Amat1
grd00=matrix(c(-0.279390198, -0.360158604,  0.007719413,  2.902898338, -0.707744688,  1.526532465,  2.519992471, -0.563860414,  1.029708300, 0.304492148,  0.773647067,  1.909845753,  3.509072163,  0.213329105,  1.763505494,  2.098557268,  2.232402942,  0.546085858, -0.419384693,  1.798717819,  0.059779088 ),J,1)
grgamma00=array(0,dim=c(r,r,J))
grgamma00[1,1,1:10]=(coef(md.noncons0,simplify=T)$G2$items[,c("a1","a2")]-coef(md.noncons0,simplify=T)$G1$items[,c("a1","a2")])[1:10,1]
grgamma00[2,1,1:10]=(coef(md.noncons0,simplify=T)$G3$items[,c("a1","a2")]-coef(md.noncons0,simplify=T)$G1$items[,c("a1","a2")])[1:10,1]
grgamma00[1,2,11:21]=(coef(md.noncons0,simplify=T)$G2$items[,c("a1","a2")]-coef(md.noncons0,simplify=T)$G1$items[,c("a1","a2")])[11:21,2]
grgamma00[2,2,11:21]=(coef(md.noncons0,simplify=T)$G3$items[,c("a1","a2")]-coef(md.noncons0,simplify=T)$G1$items[,c("a1","a2")])[11:21,2]

mu100=c(0,0)
mu200=c(-0.285, -0.284)
mu300=c(-0.579, -0.641 )
Sig100=matrix(c(1,0.914,0.914,1),2,2)
Sig200=matrix(c(2.194,2.057,2.057, 2.269),2,2)
Sig300=matrix(c(2.073,1.899,1.899, 2.044),2,2)

# MLE starting values
A11=c(3.887732, 5.476158, 4.598185, 5.070688, 4.839723, 5.038220, 5.371222, 6.180719, 4.198201, 3.616412,rep(0,11))
A12=A11
A21=c(rep(0,10),2.647928, 3.841790, 4.287103, 5.772638, 5.093658, 6.259422, 4.975799, 4.004270, 4.485499, 4.091389, 4.401138)
A22=A21
Amat1=Amat2=cbind(A11,A21)

gra00=Amat1
grd00=matrix(c(0.1815997,  0.2354828,  0.5585105,  3.5826869, -0.1031907,  2.1672437,  3.1935803,  0.2964340,  1.5510633,  0.7371429,  1.1500392,  2.4024150,  4.0687241,  0.8382791,  2.3306916,  2.7713551,2.8106932,  1.0783690,  0.1018597,  2.2761272,  0.5398323),J,1)
grgamma00[1,1,1:10]=c(0.3680776, -0.1459668,  0.4013647,  0.7255812,  0.7677556,  0.3720315,  0.3290386,  1.7985989,  0.5161178,  0.1114731)
grgamma00[2,1,1:10]=c( 0.03458868, -0.48083793,  0.10340428,  1.29530066,  0.30145681,  1.14761525,  0.86473662,  1.76403916,  0.37589465,  0.26740820)
grgamma00[1,2,11:21]=c( 0.85237225,  0.95222932,  0.83949648,  0.08927024,  0.25678010, -0.12443061,  0.50140125,  1.27521006,  0.57021110,  0.52155279,  0.45802679)
grgamma00[2,2,11:21]=c( 1.50987868,  0.99831730,  1.25167376,  0.01358535,  0.30033031,  0.02360375,  0.61014055,  1.29257372,  0.37820736,  0.36470687, -0.20463719)
mu100=c(0,0)
mu200=c(-0.2974597, -0.2869540)
mu300=c(-0.4826673, -0.5171873)
Sig100=matrix(c(1,0.9158126,0.9158126,1),2,2)
Sig200=matrix(c(0.9353043,0.8660313,0.8660313, 0.9510791),2,2)
Sig300=matrix(c(0.8409400,0.7709742,0.7709742,0.8362425),2,2)

r=2
m=2
y=3
eta.vec=seq(1,25,1)
bics=rep(0,length(eta.vec))
ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
#Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
for (k in 1:length(eta.vec))
{
  eta=eta.vec[k]
  ptm <- proc.time()
  sim=Reg_DIF(resp=Resp_ordered2,m=2,r=2,y=3,N.vec=c(N1,N2,N3),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=matrix(0,J,y-1),grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=T)
  print(proc.time() - ptm)
  bics[k]=sim$bic
  #Gammas[,,,k]=sim$Gamma
  ADmat[,,k]=sim$est
  Gammas[,,,k]=sim$Gamma
  theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  print(Gammas)
}

kk=which.min(bics)

eta.vec[kk]
Gammas[,,,kk]
ADmat[,,kk]
theta.dist[,,kk]


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

kk=which.min(bics)

eta.vec[kk]
Gammas[,,,kk]
ADmat[,,kk]
theta.dist[,,kk]

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


kk=which.min(bics)

eta.vec[kk]
Gammas[,,,kk]
ADmat[,,kk]
theta.dist[,,kk]