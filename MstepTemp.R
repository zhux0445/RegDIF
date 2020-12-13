M_step=function(max.tol,m,eta)
  // [[Rcpp::export]]
arma::vec M_step(int j, double ng, arma::mat X, int y, int G, arma::mat yallgroup, double maxtol, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int eta)
{
  arma::vec d=grd.row(j);
  arma::vec a=gra.row(j);
  arma::vec gam=grgamma.slice(j);
  arma::vec bet=grbeta.row(j);
  arma::cube Pstar=zeros<cube>(G,(m-1),y);
  arma::cube Qstar=zeros<cube>(G,(m-1),y);
  arma::cube P=zeros<cube>(G,m,y); 
  int miter=0;
  for(miter=0; miter<500; miter++){
    for(int yy=0; yy<y; yy++){
      for(int g = 0; g < G; g++){
        (Pstar.slice(yy)).row(g)=1/(1+exp(-(d+a*X.row(g)+ yallgroup.row(yy)*bet+yallgroup.row(yy)*gam*X.row(g))));
        (P.slice(yy)).row(g)=-diff(join_horiz(1,(Pstar.slice(yy)).row(g),0),1,1);
      }
    }
    Qstar=ones<cube>(G,(m-1),y)-Pstar;
    arma::vec Dsco=zeros<vec>(m-1);
    for (int yy=0; yy<y; yy++){
      Dsco += sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy));
    }
    arma::cube PQdif=zeros<cube>(G,m,y); 
    for (int yy=0; yy<y; yy++){
      PQdif.slice(yy)=-diff(join_horiz(zeros<colvec>(G),Pstar.slice(yy)%Qstar.slice(yy),zeros<colvec>(G)),1,1).t();
    }
    
    int len=0;
    for (int kk=0; kk<r; kk++){
      if (a(kk)!=0){
        len++;
      }
    }
    arma::vec Asco=zeros<vec>(len);
    for (int yy=0; yy<y; yy++){
      Asco += sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1)*
    }
    for (yy in 1:y){
      Asco=Asco+(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy]*PQdif[,,yy],1,sum))%*%(X%*%ifelse(a==0,0,1))
    }    
    if(sum(add)<maxtol){
      break;
    }
  }
}