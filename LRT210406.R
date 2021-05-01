library(MASS)
library(mirt) 
mirtCluster(8)
library(cacIRT)
library(mvtnorm)
library(graphics)

setwd('/Users/zhux0445/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("RESP1.csv",row.names = 1)

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


#####  No anchor

rep=10

resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
resp01=resp[1:(N1+N2),]
resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
#select anchor
md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop')
ratio1=d$X2/d$df
(anchor1=which(ratio1==sort((d$X2/d$df)[c(1,3:11)])[1]))
(anchor2=which(ratio1==sort((d$X2/d$df)[c(2,12:20)])[1]))

#omnibus dif
md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor1,anchor2)]))
dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor1,anchor2)])
dif1
#dif1.t=dif1[which(dif1$adj_pvals<0.05),]
#power1[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif1.t))
#tpI1[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif1.t))


#Use two different anchors 


anchor11=1
anchor21=2
#omnibus dif
md.noncons01 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor11,anchor21)]))
dif11=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor11,anchor21)])
dif11

#Use other two different anchors 

anchor12=3
anchor22=20
#omnibus dif
md.noncons02 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor12,anchor22)]))
dif12=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor12,anchor22)])
dif12

#Use other two different anchors 


anchor13=10
anchor23=20
#omnibus dif
md.noncons03 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor13,anchor23)]))
dif13=DIF(md.noncons03, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor13,anchor23)])
dif13

# Another replication
rep=20

resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
resp01=resp[1:(N1+N2),]
resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
#select anchor
md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop')
ratio1=d$X2/d$df
(anchor1=which(ratio1==sort((d$X2/d$df)[c(1,3:11)])[1]))
(anchor2=which(ratio1==sort((d$X2/d$df)[c(2,12:20)])[1]))

#omnibus dif
md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor1,anchor2)]))
dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor1,anchor2)])
dif1
#dif1.t=dif1[which(dif1$adj_pvals<0.05),]
#power1[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif1.t))
#tpI1[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif1.t))


#Use two different anchors 


anchor11=1
anchor21=2
#omnibus dif
md.noncons01 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor11,anchor21)]))
dif11=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor11,anchor21)])
dif11

#Use other two different anchors 

anchor12=3
anchor22=20
#omnibus dif
md.noncons02 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor12,anchor22)]))
dif12=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor12,anchor22)])
dif12

#Use other two different anchors 


anchor13=10
anchor23=20
#omnibus dif
md.noncons03 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor13,anchor23)]))
dif13=DIF(md.noncons03, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor13,anchor23)])
dif13

