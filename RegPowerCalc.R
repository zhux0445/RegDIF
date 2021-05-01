for(rep in 1:50){
  print(read.csv(paste("/Users/ruoyizhu/Desktop/DIF replication/Beta2_",rep),row.names = 1))
}
read.csv(paste("/Users/ruoyizhu/Desktop/DIF replication/Beta2_",rep),row.names = 1)

# sim1

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/sim1results/Beta1_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(14*50)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:50){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:50){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:50){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/200
absbias2/200

# sim3
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
  }

colSums(power)/200
colSums(tpI)/(14*50)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:50){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)
  }
biassum/50

for (rep in 1:50){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:50){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),2]-1))
  }

absbias1/200
absbias2/200

# 60% sim2

reps=20
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:20)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*20)
colSums(tpI)/(6*20)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:20,22:29,31:35,37:50)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/ADmat2_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/ADmat2_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/47

for (rep in c(1:20,22:29,31:35,37:50)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/ADmat2_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/ADmat2_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*47))
sqrt(rmsesum[3]/(20*47))

absbias1=0
absbias2=0
for (rep in c(1:20,22:29,31:35,37:50)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/new DIF results/Beta2_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2]-1))
}

absbias1/(12*47)
absbias2/(12*47)

# 60% sim4

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(6*50)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:50){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/ADmat4_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/ADmat4_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:50){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/ADmat4_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/ADmat4_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:50){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta4_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2]-1))
}

absbias1/(12*50)
absbias2/(12*50)


####### Non-uniform

reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
  }

colSums(power)/200
colSums(tpI)/(14*50)

reps=50
Gamma.6=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*50)
colSums(tpI)/(6*50)

reps=50
Gamma.7=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma7_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(14*50)

reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*50)
colSums(tpI)/(6*50)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:50){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:50){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:50){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/200
absbias2/200



reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/non-uniform DIF_1500/Gamma5_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}



# 60% uniform

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:26){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(14*50)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:50){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:50){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:50){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/200
absbias2/200


# group-lasso sim9

reps=50
Betas.9=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:43){
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta9mm_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta9_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta9mm_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta9_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*43)
colSums(tpI)/(14*43)

# group-lasso sim11

reps=50
Betas.11=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:50){
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta11_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta11_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta11_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta11_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*50)
colSums(tpI)/(14*50)


# group-lasso sim10

reps=50
Betas.10=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:50){
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta10_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta10_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta10_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta10_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(6*50)

# group-lasso sim12

reps=50
Betas.12=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:21,26:50)){
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta12_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta12_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta12_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Beta12_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*46)
colSums(tpI)/(6*46)







sim4lrt1=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Sim2_LRTfn1.csv",row.names = 1)
sim4lrt2=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Sim2_LRTfn2.csv",row.names = 1)
sim4lrt3=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Sim2_LRTfn3.csv",row.names = 1)
colMeans(sim4lrt1)
colMeans(sim4lrt2)
colMeans(sim4lrt3)

sim9lrt1=read.csv("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Sim9_LRTpvs.csv",row.names = 1)
sum(sim9lrt1[,c(2,3,10,11)]<0.05)/(50*4)
sum(sim9lrt1[,-c(2,3,10,11)]<0.05)/(50*14)


reps=48
Gamma.6=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)

reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)


reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=49
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)


reps=40
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)


reps=39
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)

reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:23,25:37,39:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*48)
colSums(tpI)/(14*48)

reps=29
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)


reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*50)
colSums(tpI)/(14*50)



reps=29
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)


reps=33
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)


reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=39
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=40
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=47
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=30
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 31:49){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*19)
colSums(tpI)/(14*19)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:30){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMBIC_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/30

for (rep in 1:30){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMBIC_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*30))
sqrt(rmsesum[3]/(20*30))

absbias1=0
absbias2=0
for (rep in 1:30){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*30)
absbias2/(4*30)

# 60% sim2 EMM lowcor 
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 6:49){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*44)
colSums(tpI)/(6*44)

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(6*50)

##########
##########
# Dec 23
##########
#sim1 adapt

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:11,13:29,31:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*48)
colSums(tpI)/(14*48)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:11,13:29,31:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/48

for (rep in c(1:11,13:29,31:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*48))
sqrt(rmsesum[3]/(20*48))

absbias1=0
absbias2=0
for (rep in c(1:11,13:29,31:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1adapt_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*48)
absbias2/(4*48)

#sim1 EMM

reps=31
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:4,6:31)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*30)
colSums(tpI)/(14*30)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:4,6:31)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/30

for (rep in c(1:4,6:31)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*30))
sqrt(rmsesum[3]/(20*30))

absbias1=0
absbias2=0
for (rep in c(1:4,6:31)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*30)
absbias2/(4*30)


# sim2 EM
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 6:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*45)
colSums(tpI)/(6*45)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 6:50){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMHighCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMHighCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/45

for (rep in 6:50){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMHighCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMHighCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*45))
sqrt(rmsesum[3]/(20*45))

absbias1=0
absbias2=0
for (rep in 6:50){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMHighCor_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*45)
absbias2/(12*45)



# sim2 Adapt
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:11,13:29,31:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*48)
colSums(tpI)/(6*48)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:11,13:29,31:50)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/48

for (rep in c(1:11,13:29,31:50)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*48))
sqrt(rmsesum[3]/(20*48))

absbias1=0
absbias2=0
be1.length=0
be2.length=0
for (rep in c(1:11,13:29,31:50)){
  be1=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[c(4:9,12:17),1]
  be2=read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2adapt_",rep),row.names = 1)[c(4:9,12:17),2]
  absbias1=absbias1+sum(abs(be1[which(be1!=0)]-0.5))
  absbias2=absbias2+sum(abs(be2[which(be2!=0)]-1))
  be1.length=be1.length+length(be1[which(be1!=0)])
  be2.length=be2.length+length(be2[which(be2!=0)])
}

absbias1/be1.length
absbias2/be2.length

# sim 3 adapt 
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:13,16:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*48)
colSums(tpI)/(14*48)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:13,16:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/48

for (rep in c(1:13,16:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*48))
sqrt(rmsesum[3]/(20*48))

absbias1=0
absbias2=0
for (rep in c(1:13,16:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3adapt_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*48)
absbias2/(4*48)


# sim4 EM
reps=12
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:5,8:12)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*10)
colSums(tpI)/(6*10)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:5,8:12)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EM_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EM_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/10

for (rep in c(1:5,8:12)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EM_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EM_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*10))
sqrt(rmsesum[3]/(20*10))

absbias1=0
absbias2=0
for (rep in c(1:5,8:12)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EM_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*10)
absbias2/(12*10)

# sim4 Adapt
reps=45
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:20,26:45)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*40)
colSums(tpI)/(6*40)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:20,26:45)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/40

for (rep in c(1:20,26:45)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*40))
sqrt(rmsesum[3]/(20*40))

absbias1=0
absbias2=0
for (rep in c(1:20,26:45)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adapt_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*40)
absbias2/(12*40)

#sim1 EMM lowcor

reps=31
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:4,6:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:4,6:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/30

for (rep in c(1:4,6:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1EMMBIC_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*30))
sqrt(rmsesum[3]/(20*30))

absbias1=0
absbias2=0
for (rep in c(1:4,6:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1EMMBIC_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*30)
absbias2/(4*30)

#sim1 adapt lowcor

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1AdaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1AdaptLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in c(1:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1AdaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat1AdaptLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in c(1:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*50)
absbias2/(4*50)

# sim 2 EM lowcor 

reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 6:49){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*44)
colSums(tpI)/(6*44)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 6:49){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EM_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EM_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/44

for (rep in 6:49){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EM_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EM_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*44))
sqrt(rmsesum[3]/(20*44))

absbias1=0
absbias2=0
for (rep in 6:49){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EM_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*44)
absbias2/(12*44)

# sim2 EMM LowCor
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(6*50)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:50)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMMLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMMLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in c(1:50)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMMLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2EMMLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in c(1:50)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*50)
absbias2/(12*50)

# sim2 Adapt LowCor
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:6,10:25,31:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*42)
colSums(tpI)/(6*42)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:6,10:25,31:50)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2AdaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2AdaptLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/42

for (rep in c(1:6,10:25,31:50)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2AdaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat2AdaptLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*42))
sqrt(rmsesum[3]/(20*42))

absbias1=0
absbias2=0
for (rep in c(1:6,10:25,31:50)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta2AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*42)
absbias2/(12*42)

# sim3 EM lowcor
reps=50
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3LowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3LowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in c(1:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3LowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3LowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in c(1:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*50)
absbias2/(4*50)

# sim3 EMM lowcor

reps=40
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/40

for (rep in c(1:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3EMMLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*40))
sqrt(rmsesum[3]/(20*40))

absbias1=0
absbias2=0
for (rep in c(1:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*40)
absbias2/(4*40)

# sim3 adapt lowcor
reps=47
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para1.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3AdaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3AdaptLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/47

for (rep in c(1:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3AdaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat3AdaptLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*47))
sqrt(rmsesum[3]/(20*47))

absbias1=0
absbias2=0
for (rep in c(1:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2]-1))
}

absbias1/(4*47)
absbias2/(4*47)

# sim4 EM lowcor
reps=35
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4LowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4LowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/35

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4LowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4LowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*35))
sqrt(rmsesum[3]/(20*35))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4LowCor_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*35)
absbias2/(12*35)

# sim4 EMM lowcor
reps=45
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:9,13:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*42)
colSums(tpI)/(6*42)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:9,13:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EMMLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EMMLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/42

for (rep in c(1:9,13:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EMMLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4EMMLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*42))
sqrt(rmsesum[3]/(20*42))

absbias1=0
absbias2=0
for (rep in c(1:9,13:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*42)
absbias2/(12*42)

# sim4 Adapt lowcor
reps=38
Betas.3=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:22,26:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)) 
}

colSums(power)/(12*35)
colSums(tpI)/(6*35)


biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para2.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:22,26:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adaptLowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/35

for (rep in c(1:22,26:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adaptLowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat4adaptLowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*35))
sqrt(rmsesum[3]/(20*35))

absbias1=0
absbias2=0
for (rep in c(1:22,26:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[c(4:9,12:17),1]-0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1)[c(4:9,12:17),2]-1))
}

absbias1/(12*35)
absbias2/(12*35)

# sim5 EMM
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIC_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIC_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIC_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)


# sim5 adapt
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIC_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIC_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIC_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIC_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIC_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)

# sim 6 EM
reps=36
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6Dec_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6Dec_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/36

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6Dec_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6Dec_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*36))
sqrt(rmsesum[3]/(20*36))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6Dec_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*36)
absbias2/(12*36)


# sim 6 Adapt
reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*50)
absbias2/(12*50)

# sim7 EMM

reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIC_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIC_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIC_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIC_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)

# sim7 adapt
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7adapt_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)

# sim8 EM
reps=4
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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

for (rep in c(1:4,26:29)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/8

for (rep in c(1:4,26:29)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*8))
sqrt(rmsesum[3]/(20*8))

absbias1=0
absbias2=0
for (rep in c(1:4,26:29)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*8)
absbias2/(12*8)
# sim8 EMM
reps=16
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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

for (rep in c(1:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMM_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMM_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/16

for (rep in c(1:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMM_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMM_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*16))
sqrt(rmsesum[3]/(20*16))

absbias1=0
absbias2=0
for (rep in c(1:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMM_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*reps)
absbias2/(12*reps)

# sim8 adapt
reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8adapt_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8adapt_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/reps

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8adapt_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8adapt_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*reps))
sqrt(rmsesum[3]/(20*reps))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*reps)
absbias2/(12*reps)



# sim5 EM lowcor
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5BIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5BIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5BIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5BIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)


# sim5 EMM lowcor
reps=49
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/49

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*49))
sqrt(rmsesum[3]/(20*49))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*49)
absbias2/(4*49)


# sim5 adapt lowcor
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat5adaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma5adaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)


# sim 6 EM lowcor
reps=40
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:17,19:reps)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*39)
colSums(tpI)/(6*39)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para4.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:17,19:reps)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6LowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6LowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/39

for (rep in c(1:17,19:reps)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6LowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6LowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*39))
sqrt(rmsesum[3]/(20*39))

absbias1=0
absbias2=0
for (rep in c(1:17,19:reps)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*39)
absbias2/(12*39)

# sim 6 EMM lowcor
reps=48
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/48

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*48))
sqrt(rmsesum[3]/(20*48))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*48)
absbias2/(12*48)

# sim 6 Adapt lowcor
reps=49
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/49

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*49))
sqrt(rmsesum[3]/(20*49))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*49)
absbias2/(12*49)

# sim7 EM lowcor
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:23,25:37,39:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*48)
colSums(tpI)/(14*48)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:23,25:37,39:50)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7BIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7BIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/48

for (rep in c(1:23,25:37,39:50)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7BIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7BIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*48))
sqrt(rmsesum[3]/(20*48))

absbias1=0
absbias2=0
for (rep in c(1:23,25:37,39:50)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7BIClowcor_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*48)
absbias2/(4*48)

# sim7 EMM lowcor
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*reps)
colSums(tpI)/(14*reps)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7EMMBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)

# sim7 adapt lowcor
reps=50
Gamma.5=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,10,11),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(4*50)
colSums(tpI)/(14*50)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para3.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat7adaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma7AdaptBIClowcor_",rep),row.names = 1)[c(2,3,10,11),2]+1))
}

absbias1/(4*50)
absbias2/(4*50)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para4.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/49

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat6AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*49))
sqrt(rmsesum[3]/(20*49))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma6AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*49)
absbias2/(12*49)

# sim 8 EM lowcor
reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in c(1:22,31:37,40,41,46,47)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*33)
colSums(tpI)/(6*33)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/zhux0445/Documents/Github/RegDIF_SimData/Para4.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in c(1:22,31:37,40,41,46,47)){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8LowCor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8LowCor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/33

for (rep in c(1:22,31:37,40,41,46,47)){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8LowCor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8LowCor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*33))
sqrt(rmsesum[3]/(20*33))

absbias1=0
absbias2=0
for (rep in c(1:22,31:37,40,41,46,47)){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8LowCor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*33)
absbias2/(12*33)

# sim 8 EMM lowcor
reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMMBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8EMMBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8EMMBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*50)
absbias2/(12*50)

# sim 8 Adapt lowcor
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
