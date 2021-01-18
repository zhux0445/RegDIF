reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}
colSums(power)/(12*reps)
colSums(tpI)/(6*reps)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para4.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*50)
absbias2/(12*50)

###############
# Refit model #
###############
#grgamma=array(0,dim=c(r,r,J))
#grgamma[c(1,2),1,3:11]=matrix(t(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[1:9,]),2,9)
#grgamma[c(1,2),2,12:20]=matrix(t(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[10:18,]),2,9)
grgamma=matrix(0,J,y-1)
grgamma[3:20,]=as.matrix(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1))

sparsity=grgamma
for (j in 1:J){
  for (rr in 1:2){
    sparsity[j,rr]=ifelse(grgamma[j,rr]==0,0,1)
  }
} 

# use mirt for re-estimate
anchor1=which(sparsity[,1]==0)
anchor2=which(sparsity[,2]==0)
model3 <- '
    F1 = 1,3-11
    F2 = 2,12-20
    COV = F1*F2'
Group=c(rep('G1', N1), rep('G2', N2),rep('G3', N3))
#values <- multipleGroup(resp, model3, group = Group, pars = 'values')
#values
N1=N2=N3=500 
N=N1+N2+N3
mirtCluster(8)


rep=1#rep in 1:50
resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]),pars = 'values')
constrainA<- list(c(1,106,211), c(7,112,217),c(11,116,221),c(16,121,226),c(21,126,231),c(26,131,236),c(31,136,241),c(36,141,246),c(41,146,251),c(46,151,256),c(51,156,261),c(57,162,267),c(62,167,272),c(67,172,277),c(72,177,282),c(77,182,287),c(82,187,292),c(87,192,297),c(92,197,302),c(97,202,307))
Dgp1=seq(3,98,5)
Dgp2=seq(108,203,5)
Dgp3=seq(213,308,5)
Dgp23=cbind(Dgp2,Dgp3)
invsparsity=1-sparsity
nonzerobeta=which(rowSums(invsparsity)!=0)
constrainD23<- cbind( Dgp1,invsparsity* Dgp23)
constrainD<- vector("list", length(nonzerobeta))
for(nn in 1:length(nonzerobeta)){
  constrainD[[nn]]=c(fun.zero.omit((constrainD23[nonzerobeta,])[nn,]))
}
equalslopes <- multipleGroup(resp, model3, group = Group, constrain = c(constrainA,constrainD),invariance=c('free_means', 'free_var'))
dif.mag.est=cbind((coef(equalslopes,simplify=T)$G2$items-coef(equalslopes,simplify=T)$G1$items)[,3],(coef(equalslopes,simplify=T)$G3$items-coef(equalslopes,simplify=T)$G1$items)[,3])
