library(mirt) 
mirtCluster(4)
setwd('/Users/hyzhu27/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para4.csv",row.names = 1)
responses=read.csv("RESP6.csv",row.names = 1)
Amat1=params[,1:2]
Amat2=params[,3:4]
Amat3=params[,5:6]
Dmat1=matrix(params[,7],J,(m-1))
Dmat2=matrix(params[,8],J,(m-1))
Dmat3=matrix(params[,9],J,(m-1))
Amat1=as.matrix(Amat1)
Amat2=as.matrix(Amat2)
Amat3=as.matrix(Amat3)
rownames(Amat1) <- c()
rownames(Amat2) <- c()
rownames(Amat3) <- c()
J=20
N1=N2=N3=500 
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
N=N1+N2+N3
m=2
r=2

mirt.p.mat1=matrix(0,50,18) # 16 is the number of replications
bias.mirt1=matrix(0,50,3)
rmse.mirt1=matrix(0,50,3)
difrec.mirt.fn1=matrix(0,50,3) # DIF magnitude recovery with false negative included (only need the first column)
difrec.mirt.fn1.r=matrix(0,50,3) # DIF magnitude recovery for regularization method omnibus (only need the first column)
mirt.p.mat2=matrix(0,50,18) # 16 is the number of replications
bias.mirt2=matrix(0,50,3)
rmse.mirt2=matrix(0,50,3)
difrec.mirt.fn2=matrix(0,50,3) #(only need the second column)
mirt.p.mat3=matrix(0,50,18) # 16 is the number of replications
bias.mirt3=matrix(0,50,3)
rmse.mirt3=matrix(0,50,3)
difrec.mirt.fn3=matrix(0,50,3) #(only need the third column)
ADmat=matrix(0,20*50,21)
Traitdistmat=matrix(0,21*50,2)
for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:r]))
  mirt.p.mat1[(rep),]=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-(which(mirt.p.mat1[rep,]<0.05)+2)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-c(4,5,12,13)]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[rep,2:3]= c(mean(abs((coef(md.refit0,simplify=T)$G2$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[(which(mirt.p.mat1[rep,]<0.05)+2)]-0.5)),mean(abs((coef(md.refit0,simplify=T)$G3$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[(which(mirt.p.mat1[rep,]<0.05)+2)]-1)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3])) #omnibus DIF recovery in report (only need this)
  #difrec.mirt.fn1.r[rep,2:3]= c(mean(abs((coef(md.refit.r,simplify=T)$G2$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit.r,simplify=T)$G3$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  #difrec.mirt.fn1.r[rep,1]=mean(c(difrec.mirt.fn1.r[rep,2],difrec.mirt.fn1.r[rep,3])) #omnibus DIF recovery in report (only need this)
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:r]))
  mirt.p.mat2[(rep),]=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  if ((sum(mirt.p.mat2[rep,]<0.05))==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[which(mirt.p.mat2[rep,]<0.05)+2]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:r]))
  mirt.p.mat3[(rep),]=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[-(which(mirt.p.mat3[rep,]<0.05)+2)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[which(mirt.p.mat3[rep,]<0.05)+2]-1))
  ADmat[((rep-1)*20+1):(rep*20),]=cbind(coef(md.refit0,simplify=T)$G1$items[,1:3],coef(md.refit0,simplify=T)$G2$items[,1:3],coef(md.refit0,simplify=T)$G3$items[,1:3],coef(md.refit01,simplify=T)$G1$items[,1:3],coef(md.refit01,simplify=T)$G2$items[,1:3],coef(md.refit02,simplify=T)$G1$items[,1:3],coef(md.refit02,simplify=T)$G3$items[,1:3])
  Traitdistmat[((rep-1)*21+1):(rep*21),]=rbind(coef(md.refit0,simplify=T)$G1$means,coef(md.refit0,simplify=T)$G2$means,coef(md.refit0,simplify=T)$G3$means,coef(md.refit0,simplify=T)$G1$cov,coef(md.refit0,simplify=T)$G2$cov,coef(md.refit0,simplify=T)$G3$cov,coef(md.refit01,simplify=T)$G1$means,coef(md.refit01,simplify=T)$G2$means,coef(md.refit01,simplify=T)$G1$cov,coef(md.refit01,simplify=T)$G2$cov,coef(md.refit02,simplify=T)$G1$means,coef(md.refit02,simplify=T)$G3$means,coef(md.refit02,simplify=T)$G1$cov,coef(md.refit02,simplify=T)$G3$cov)
  }
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(na.omit(difrec.mirt.fn2))
colMeans(difrec.mirt.fn3)
write(ADmat,file = "Sim1lowcor_LRT_ADmat.csv")
write(Traitdistmat,file = "Sim1lowcor_LRT_Traitdistmat.csv")

for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #md.noncons0 <- multipleGroup(resp, cmodel, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  #omnibus
  # From Phil, when you test DIF on slope you should also test the intercept at the same time
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  #md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:r]))
  mirt.p.mat1[(rep),1:9]=DIF(md.noncons0, which.par = c('a1'), p.adjust = 'none',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat1[(rep),10:18]=DIF(md.noncons0, which.par = c('a2'), p.adjust = 'none',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-(which(mirt.p.mat1[rep,]<0.05)+2)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-c()]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[rep,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[(which(mirt.p.mat1[rep,(4:5)]<0.05)+3)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[(which(mirt.p.mat1[rep,(12:13)]<0.05)+11)]-0.5)))
  difrec.mirt.fn1[rep,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[(which(mirt.p.mat1[rep,(4:5)]<0.05)+3)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[(which(mirt.p.mat1[rep,(12:13)]<0.05)+11)]-1)))
  #difrec.mirt.fn1[rep,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[c(4,5,6,7,8,9)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[c(12,13,14,15,16,17)]-0.5)))
  #difrec.mirt.fn1[rep,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[c(4,5,6,7,8,9)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[c(12,13,14,15,16,17)]-1)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3]))
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[1:r]))
  mirt.p.mat2[(rep),1:9]=DIF(md.noncons01, which.par = c('a1'), p.adjust = 'none',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat2[(rep),10:18]=DIF(md.noncons01, which.par = c('a2'), p.adjust = 'none',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  #md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  if ((sum(mirt.p.mat2[rep,]<0.05))==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(c(abs((coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[(which(mirt.p.mat1[rep,(4:5)]<0.05)+3)]-0.5),abs((coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[(which(mirt.p.mat1[rep,(12:13)]<0.05)+11)]-0.5)))
  #difrec.mirt.fn2[rep,2]= mean(c(abs((coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[c(4,5,6,7,8,9)]-0.5),abs((coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[c(12,13,14,15,16,17)]-0.5)))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp02)[1:r]))
  mirt.p.mat3[(rep),1:9]=DIF(md.noncons02, which.par = c('a1'), p.adjust = 'none',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat3[(rep),10:18]=DIF(md.noncons02, which.par = c('a2'), p.adjust = 'none',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp02)[-(which(mirt.p.mat3[rep,]<0.05)+2)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]= mean(c(abs((coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[(which(mirt.p.mat1[rep,(4:5)]<0.05)+3)]-1),abs((coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[(which(mirt.p.mat1[rep,(12:13)]<0.05)+11)]-1)))
  #difrec.mirt.fn3[rep,3]= mean(c(abs((coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[c(4,5,6,7,8,9)]-1),abs((coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[c(12,13,14,15,16,17)]-1)))
  ADmat[((rep-1)*20+1):(rep*20),]=cbind(coef(md.refit0,simplify=T)$G1$items[,1:3],coef(md.refit0,simplify=T)$G2$items[,1:3],coef(md.refit0,simplify=T)$G3$items[,1:3],coef(md.refit01,simplify=T)$G1$items[,1:3],coef(md.refit01,simplify=T)$G2$items[,1:3],coef(md.refit02,simplify=T)$G1$items[,1:3],coef(md.refit02,simplify=T)$G3$items[,1:3])
  Traitdistmat[((rep-1)*21+1):(rep*21),]=rbind(coef(md.refit0,simplify=T)$G1$means,coef(md.refit0,simplify=T)$G2$means,coef(md.refit0,simplify=T)$G3$means,coef(md.refit0,simplify=T)$G1$cov,coef(md.refit0,simplify=T)$G2$cov,coef(md.refit0,simplify=T)$G3$cov,coef(md.refit01,simplify=T)$G1$means,coef(md.refit01,simplify=T)$G2$means,coef(md.refit01,simplify=T)$G1$cov,coef(md.refit01,simplify=T)$G2$cov,coef(md.refit02,simplify=T)$G1$means,coef(md.refit02,simplify=T)$G3$means,coef(md.refit02,simplify=T)$G1$cov,coef(md.refit02,simplify=T)$G3$cov)
  print(rep)
}
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt1)
sqrt(colMeans(rmse.mirt1^2))
colMeans(difrec.mirt.fn1)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(18*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(18*14)
colMeans(bias.mirt2)
sqrt(colMeans(rmse.mirt2^2))
colMeans(difrec.mirt.fn2)

sum(mirt.p.mat3[1:24,c(2,3,10,11)]<0.05)/(24*4)
sum(mirt.p.mat3[1:24,-c(2,3,10,11)]<0.05)/(24*14)
colMeans(bias.mirt3)
sqrt(colMeans(rmse.mirt3^2))
colMeans(difrec.mirt.fn3)


for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #md.noncons0 <- multipleGroup(resp, cmodel, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  #omnibus
  # From Phil, when you test DIF on slope you should also test the intercept at the same time
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  #md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:r]))
  mirt.p.mat1[(rep),1:9]=DIF(md.noncons0, which.par = c('a1'), p.adjust = 'none',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat1[(rep),10:18]=DIF(md.noncons0, which.par = c('a2'), p.adjust = 'none',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-(which(mirt.p.mat1[rep,]<0.05)+2)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-c()]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[rep,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[(which(mirt.p.mat1[rep,(4:9)]<0.05)+3)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[(which(mirt.p.mat1[rep,(12:17)]<0.05)+11)]-0.5)))
  difrec.mirt.fn1[rep,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[(which(mirt.p.mat1[rep,(4:9)]<0.05)+3)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[(which(mirt.p.mat1[rep,(12:17)]<0.05)+11)]-1)))
  #difrec.mirt.fn1[rep,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[c(4,5,6,7,8,9)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[c(12,13,14,15,16,17)]-0.5)))
  #difrec.mirt.fn1[rep,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[c(4,5,6,7,8,9)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[c(12,13,14,15,16,17)]-1)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3]))
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[1:r]))
  mirt.p.mat2[(rep),1:9]=DIF(md.noncons01, which.par = c('a1'), p.adjust = 'none',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat2[(rep),10:18]=DIF(md.noncons01, which.par = c('a2'), p.adjust = 'none',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  #md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  if ((sum(mirt.p.mat2[rep,]<0.05))==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(c(abs((coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[(which(mirt.p.mat1[rep,(4:9)]<0.05)+3)]-0.5),abs((coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[(which(mirt.p.mat1[rep,(12:17)]<0.05)+11)]-0.5)))
  #difrec.mirt.fn2[rep,2]= mean(c(abs((coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[c(4,5,6,7,8,9)]-0.5),abs((coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[c(12,13,14,15,16,17)]-0.5)))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp02)[1:r]))
  mirt.p.mat3[(rep),1:9]=DIF(md.noncons02, which.par = c('a1'), p.adjust = 'none',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat3[(rep),10:18]=DIF(md.noncons02, which.par = c('a2'), p.adjust = 'none',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp02)[-(which(mirt.p.mat3[rep,]<0.05)+2)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]= mean(c(abs((coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[(which(mirt.p.mat1[rep,(4:9)]<0.05)+3)]-1),abs((coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[(which(mirt.p.mat1[rep,(12:17)]<0.05)+11)]-1)))
  #difrec.mirt.fn3[rep,3]= mean(c(abs((coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[c(4,5,6,7,8,9)]-1),abs((coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[c(12,13,14,15,16,17)]-1)))
  ADmat[((rep-1)*20+1):(rep*20),]=cbind(coef(md.refit0,simplify=T)$G1$items[,1:3],coef(md.refit0,simplify=T)$G2$items[,1:3],coef(md.refit0,simplify=T)$G3$items[,1:3],coef(md.refit01,simplify=T)$G1$items[,1:3],coef(md.refit01,simplify=T)$G2$items[,1:3],coef(md.refit02,simplify=T)$G1$items[,1:3],coef(md.refit02,simplify=T)$G3$items[,1:3])
  Traitdistmat[((rep-1)*21+1):(rep*21),]=rbind(coef(md.refit0,simplify=T)$G1$means,coef(md.refit0,simplify=T)$G2$means,coef(md.refit0,simplify=T)$G3$means,coef(md.refit0,simplify=T)$G1$cov,coef(md.refit0,simplify=T)$G2$cov,coef(md.refit0,simplify=T)$G3$cov,coef(md.refit01,simplify=T)$G1$means,coef(md.refit01,simplify=T)$G2$means,coef(md.refit01,simplify=T)$G1$cov,coef(md.refit01,simplify=T)$G2$cov,coef(md.refit02,simplify=T)$G1$means,coef(md.refit02,simplify=T)$G3$means,coef(md.refit02,simplify=T)$G1$cov,coef(md.refit02,simplify=T)$G3$cov)
  print(rep)
}
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt1)
sqrt(colMeans(rmse.mirt1^2))
colMeans(difrec.mirt.fn1)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(18*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(18*14)
colMeans(bias.mirt2)
sqrt(colMeans(rmse.mirt2^2))
colMeans(difrec.mirt.fn2)

sum(mirt.p.mat3[1:24,c(2,3,10,11)]<0.05)/(24*4)
sum(mirt.p.mat3[1:24,-c(2,3,10,11)]<0.05)/(24*14)
colMeans(bias.mirt3)
sqrt(colMeans(rmse.mirt3^2))
colMeans(difrec.mirt.fn3)

