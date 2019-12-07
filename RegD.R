# DIF identification
# Regularized MGRM

library(MASS)
library(mirt)
# num of items is fixed to 20
J=20

# Dataset #1 (2 Uniform DIF items per scale)
# item (4,5,12,13) have DIF and item (1,11) are anchor (DIF free)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A13=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A23=A21
D11=rnorm(J,0,1)
D12=D13=D11
# focal group has more extreme intercepts
D12[4]=D11[4]+ 0.5
D12[5]=D11[5]+ 0.5
D13[4]=D11[4]+ 1
D13[5]=D11[5]+ 1
D12[12]=D11[12]+ 0.5
D12[13]=D11[13]+ 0.5
D13[12]=D11[12]+ 1
D13[13]=D11[13]+ 1
Amat01=cbind(A11,A21)
Amat1=Amat01
Amat1[2,]=Amat01[11,]
Amat1[11,]=Amat01[2,]
Amat02=cbind(A12,A22)
Amat2=Amat02
Amat2[2,]=Amat02[11,]
Amat2[11,]=Amat02[2,]
Amat03=cbind(A13,A23)
Amat3=Amat03
Amat3[2,]=Amat03[11,]
Amat3[11,]=Amat03[2,]
Dmat1=cbind(D11)
Dmat2=cbind(D12)
Dmat3=cbind(D13)

# Dataset #2 (4 Uniform DIF items per scale)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A21
D11=rnorm(J,0,1)
D12=D11
# focal group has more extreme intercepts
D12[3]=D11[3]+ 0.25
D12[4]=D11[4]+ 0.25 
D12[5]=D11[5]+ 0.5
D12[6]=D11[6]+ 0.5
D12[11]=D11[11]+ 0.25
D12[12]=D11[12]+ 0.25
D12[13]=D11[13]+ 0.5
D12[14]=D11[14]+ 0.5
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Dmat1=cbind(D11)
Dmat2=cbind(D12)

# Dataset #3 (2 Non-uniform DIF items per scale)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A21
D11=rnorm(J,0,1)
D12=D11
# focal group has smaller slopes and more extreme intercepts
A12[4]=A11[4]- 0.3
A12[5]=A11[5]- 0.6
A22[12]=A21[12]- 0.3
A22[13]=A21[13]- 0.6
D12[4]=D11[4]+ 0.25
D12[5]=D11[5]+ 0.5
D12[12]=D11[12]+ 0.25
D12[13]=D11[13]+ 0.5
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Dmat1=cbind(D11)
Dmat2=cbind(D12)

# Dataset #4 (4 Non-uniform DIF items per scale)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A21
D11=rnorm(J,0,1)
D12=D11
# focal group has smaller slopes and more extreme intercepts
A12[3]=A11[3]-0.3
A12[4]=A11[4]-0.3
A12[5]=A11[5]-0.6
A12[6]=A11[6]-0.6
A22[11]=A21[11]- 0.3
A22[12]=A21[12]- 0.3
A22[13]=A21[13]- 0.6
A22[14]=A21[14]- 0.6
D12[3]=D11[3]+ 0.25
D12[4]=D11[4]+ 0.25
D12[5]=D11[5]+ 0.5
D12[6]=D11[6]+ 0.5
D12[11]=D11[11]+ 0.25
D12[12]=D11[12]+ 0.25
D12[13]=D11[13]+ 0.5
D12[14]=D11[14]+ 0.5
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Dmat1=cbind(D11)
Dmat2=cbind(D12)



Group=c(rep('R', N1), rep('F1', N2),rep('F2', N3))
N1=N2=N3=500
y1=c(0,0)
y2=matrix(c(1,0),1,2)
y3=matrix(c(0,1),1,2)
y=matrix(0,N1+N2+N3,2)
y[1:N1,]=y1
for (i in (N1+1):(N1+N2)){
  y[i,]=y2
}
for (i in (N1+N2+1):(N1+N2+N3)){
  y[i,]=y3
}

set.seed(40)
Theta=mvrnorm(n=N1+N2+N3,mu=c(0,0),Sigma = matrix(c(1,0.85,0.85,1),2,2))
Theta1=Theta[1:N1,]
Theta2=Theta[(N1+1):(N1+N2),]
Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
resp=rbind(datasetR,datasetF1,datasetF2)
r=2
m=2

eta=0

colSums(gra-Amat1)/10
sqrt(colSums((gra-Amat1)^2)/10)
colMeans(grd-Dmat1)
sqrt(colMeans((grd-Dmat1)^2))
Sig

colMeans(grd+grbeta[,1]-Dmat2)
sqrt(colMeans((grd+grbeta[,1]-Dmat2)^2))

colMeans(grd+grbeta[,2]-Dmat3)
sqrt(colMeans((grd+grbeta[,2]-Dmat3)^2))

grbeta=grbeta0
gra=gra0
grd=grd0



