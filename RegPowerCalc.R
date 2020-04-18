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





