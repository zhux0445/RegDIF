# sim5

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(2,3,10,11),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(2,3,10,11),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3new.csv",row.names = 1)
Amat1=params[,1:2]
Amat2=params[,3:4]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)
for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAADmat5_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAADmat5_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAADmat5_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAADmat5_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(2,3,10,11),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(2,3,10,11),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim5 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMM_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMM_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5EMM_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5EMM_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim5 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adapt_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adapt_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5Adapt_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5Adapt_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim6

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6_",rep),row.names = 1)[c(2:7,10:15),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6_",rep),row.names = 1)[c(2:7,10:15),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6_",rep),row.names = 1)[-c(2:7,10:15),1])!=0)
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6_",rep),row.names = 1)[-c(2:7,10:15),2])!=0)
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4



# sim6 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMM_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMM_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6EMM_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6EMM_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim6 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adapt_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adapt_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6Adapt_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6Adapt_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:5){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[c(2,3,10,11),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[c(2,3,10,11),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

for (rep in 6:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 6:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim7 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMM_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMM_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7EMM_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7EMM_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim7 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adapt_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adapt_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7Adapt_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7Adapt_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim8

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7


absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim8 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7


absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMM_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMM_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8EMM_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8EMM_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim8 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adapt_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adapt_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8Adapt_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8Adapt_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4



#############
# 
#   Low Cor 
#
#############

# sim5

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5lowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5lowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5lowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5lowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5lowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5lowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5lowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5lowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5lowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5lowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5lowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5lowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim5 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim5 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  #(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  print(c(rep,a1[which(a1!=0)],a2[which(a2!=0)]))
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta5Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim6

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6lowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6lowcor_",rep),row.names = 1)[c(2:7,10:15),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6lowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6lowcor_",rep),row.names = 1)[c(2:7,10:15),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6lowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6lowcor_",rep),row.names = 1)[-c(2:7,10:15),1])!=0)
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6lowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=c(FALSE, FALSE,(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6lowcor_",rep),row.names = 1)[-c(2:7,10:15),2])!=0)
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6lowcor_",rep),row.names = 1)[c(2:7,10:15),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6lowcor_",rep),row.names = 1)[c(2:7,10:15),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6lowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6lowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim6 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim6 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta6Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}
colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta6Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4



# sim7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)

for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7lowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7lowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7lowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7lowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7lowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7lowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7lowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7lowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7lowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7lowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  #print(c(rep,a1[which(a1!=0)],a2[which(a2!=0)]))
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7lowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7lowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim7 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)

for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMMlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMMlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  #print(c(rep,a1[which(a1!=0)],a2[which(a2!=0)]))
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7EMMlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4


# sim7 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)

for (rep in 1:50){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adaptlowcor_",rep),row.names = 1)[c(4,5,12,13),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta7Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adaptlowcor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*4)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma7Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  #print(c(rep,a1[which(a1!=0)],a2[which(a2!=0)]))
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta7Adaptlowcor_",rep),row.names = 1)[c(4:5,12:13),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim8

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8lowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8lowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8lowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8lowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8lowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8lowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8lowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8lowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8lowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8lowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8lowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8lowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim8 EMM

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:47)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMMlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(47*12)
colSums(tpI)/(8*47)
sd((power[1:47,]/12)[,1])/7#sqrt(50-1)=7
sd((power[1:47,]/12)[,2])/7
sd((power[1:47,]/12)[,3])/7
sd((tpI[1:47,]/8)[,1])/7#sqrt(50-1)=7
sd((tpI[1:47,]/8)[,2])/7
sd((tpI[1:47,]/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  print(c(rep,a1[which(a1!=0)],a2[which(a2!=0)]))
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8EMMlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4

# sim8 Adapt

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  a1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1])!=0
  a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2])!=0
  a5=(a1|a2)!=0
  a6=(a3|a4)!=0
  a7=(a1|a2|a3|a4)!=0
  power[rep,c(2,3)]=c(sum(a5),sum(a6))
  power[rep,1]=sum(a7)
  b1=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0
  b3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta8Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adaptlowcor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0
  b5=(b1|b2)!=0
  b6=(b3|b4)!=0
  b7=(b1|b2|b3|b4)!=0
  tpI[rep,c(2,3)]=c(sum(b5),sum(b6))
  tpI[rep,1]=sum(b7) 
}

colSums(power)/(50*12)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

absbias1=0
absbias2=0
length1=0
length2=0
for (rep in 1:reps){
  a1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  #print(c(rep,a1[which(a1!=0)],a2[which(a2!=0)]))
  absbias1=absbias1+sum(abs(a1[which(a1!=0)]+0.4))
  absbias2=absbias2+sum(abs(a2[which(a2!=0)]+0.6))
  length1=length1+length(a1[which(a1!=0)])
  length2=length2+length(a2[which(a2!=0)])
}

absbias1/length1
absbias2/length2


absbias3=0
absbias4=0
length3=0
length4=0
for (rep in 1:reps){
  a3=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),1]
  a4=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NABeta8Adaptlowcor_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias3=absbias3+sum(abs(a3[which(a3!=0)]-0.25))
  absbias4=absbias4+sum(abs(a4[which(a4!=0)]-0.6))
  length3=length3+length(a3[which(a3!=0)])
  length4=length4+length(a4[which(a4!=0)])
}

absbias3/length3
absbias4/length4
