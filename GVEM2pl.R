library(Matrix) #calculate rank
identify<-function(u){
  scale=rowSums(u) #num of item responsed correctly per examinee
  p=apply(u,2,mean) #the frequency per item
  q=1-p
  y=qnorm(p,0,1) #inverse of the standard normal CDF
  y=dnorm(y,0,1) #the density function
  s=sd(scale)
  #scale=do.call(cbind,rep(list(scale),ncol(u)))
  #x1=scale[a==1]
  r=NULL
  for (i in 1:dim(u)[2]) {
    x1=scale[u[,i]==1]
    x2=scale[u[,i]==0]
    r[i]=(mean(x1)-mean(x2))/s*p[i]*q[i]/y[i]
  }
  return(r)
}

#initialization
init<-function(u,domain,indic){
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  a0=t(rep(1,domain)%o%(r/sqrt(1-r^2)))*indic
  b0=-qnorm(colSums(u)/person,0,1)/r
  Sigma = diag(domain)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(b0)-theta%*%t(a0) #-a\theta+b=-(a\theta-b)
  eta0=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)    
  eps0=xi
  return(list(a0,b0,eta0,eps0,Sigma))
}


#GVEM-2PL CFA
vem_2PLCFA <- function(u, domain, indic,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 0
  initial=init(u,domain,indic)
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  new_a = initial[[1]] * indic
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  par_Sigma = Sigma
  while(converge==0 && rankMatrix(Sigma) == domain){
    par_Sigma=Sigma
    Spart=matrix(0,nrow=domain,ncol=domain)
    for (i in 1:person){
      sigma_part=matrix(0,nrow=domain,ncol=domain)
      mu_part= matrix(0,nrow=domain,ncol=1)
      for (j in 1:item) {
        sigma_part=sigma_part+eta[i,j]*new_a[j,]%*%t(new_a[j,])
        mu_part=mu_part+(2*eta[i,j]*new_b[j]+u[i,j]-0.5)*new_a[j,]
      }
      sigmahat=solve(solve(Sigma)+2*sigma_part)
      muhat=sigmahat%*%mu_part
      SIGMA[,,i]=sigmahat
      MU[,i]=muhat
      xi[i,]=sqrt(new_b^2-2*new_b*(new_a%*%muhat)+diag(new_a%*%(sigmahat+muhat%*%t(muhat))%*%t(new_a)))
      Spart=Spart+sigmahat+muhat%*%t(muhat)
    }
    
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(eta)< 0.01,0.125)
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    b_part = matrix(0,nrow=person,ncol=item)
    for(i in 1:person){
      b_part[i,]=new_a%*%MU[,i]
    }
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)
    par_a=new_a
    
    for(j in 1:item){
      a_nu=0
      a_de=0
      Ind=indic[j,]
      iind=which(Ind==1)
      for(i in 1:person){
        sigma=SIGMA[,,i]
        sigma=sigma[iind,iind]
        mu=MU[iind,i]
        a_de=a_de+eta[i,j]*sigma+eta[i,j]*(mu%*%t(mu))
        a_nu=a_nu+(u[i,j]-0.5+2*new_b[j]*eta[i,j])*mu
      }
      new_a[j,iind]=solve(a_de)%*%a_nu/2
    }
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2") <0.0001){
      converge=1
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
  }else{
    rsigma = Sigma
  }
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n))
}