#9/28/2020 mirt re-estimate step (using different modeling method, cannot be used)

###############
# Refit model #
###############
sparsity=grbeta
for (j in 1:J){
  for (rr in 1:2){
    sparsity[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
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
