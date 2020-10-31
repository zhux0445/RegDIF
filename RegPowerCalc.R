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
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
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

reps=20
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:20){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8adapt_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}

colSums(power)/(12*reps)
colSums(tpI)/(6*reps)

