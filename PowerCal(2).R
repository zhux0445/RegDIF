read.csv(paste("/Users/ruoyizhu/Desktop/DIF replication/Beta2_",rep),row.names = 1)

# sim1

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta1LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7



# sim2

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7



reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)

sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:46){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta2LowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*46)
colSums(tpI)/(8*46)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7





# sim3

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:42,45:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*48)
colSums(tpI)/(16*48)
sd((power/4)[,1])/sqrt(47)#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3Adapt_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:38,40:50)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3AdaptLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*49)
colSums(tpI)/(16*49)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMM_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/200
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3EMMLowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*50)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[c(4,5,12,13),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta3LowCor_",rep),row.names = 1)[-c(4,5,12,13),2])!=0)) 
}

colSums(power)/(4*50)
colSums(tpI)/(16*50)
sd((power/4)[,1])/7#sqrt(50-1)=7
sd((power/4)[,2])/7
sd((power/4)[,3])/7
sd((tpI/16)[,1])/7#sqrt(50-1)=7
sd((tpI/16)[,2])/7
sd((tpI/16)[,3])/7



# sim4

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4Adapt_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in c(1:16,26:37)){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4AdaptLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*length(c(1:16,26:37)))
colSums(tpI)/(8*length(c(1:16,26:37)))
sd((power/12)[c(1:16,26:37),1])/7#sqrt(50-1)=7
sd((power/12)[c(1:16,26:37),2])/7
sd((power/12)[c(1:16,26:37),3])/7
sd((tpI/8)[c(1:16,26:37),1])/7#sqrt(50-1)=7
sd((tpI/8)[c(1:16,26:37),2])/7
sd((tpI/8)[c(1:16,26:37),3])/7


reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMM_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:50){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4EMMLowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*50)
colSums(tpI)/(8*50)
sd((power/12)[,1])/7#sqrt(50-1)=7
sd((power/12)[,2])/7
sd((power/12)[,3])/7
sd((tpI/8)[,1])/7#sqrt(50-1)=7
sd((tpI/8)[,2])/7
sd((tpI/8)[,3])/7

reps=50
Betas.1=array(double(20*2*reps),dim = c(20,2,reps))
power=matrix(0,50,3)
tpI=matrix(0,50,3)
for (rep in 1:47){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[c(4:9,12:17),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[-c(4:9,12:17),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta4LowCor_",rep),row.names = 1)[-c(4:9,12:17),2])!=0)) 
}

colSums(power)/(12*47)
colSums(tpI)/(8*47)
sd((power/12)[1:47,1])/7#sqrt(50-1)=7
sd((power/12)[1:47,2])/7
sd((power/12)[1:47,3])/7
sd((tpI/8)[1:47,1])/7#sqrt(50-1)=7
sd((tpI/8)[1:47,2])/7
sd((tpI/8)[1:47,3])/7



