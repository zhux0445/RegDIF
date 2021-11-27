#include<RcppArmadillo.h>
#include <RcppEigen.h>

using namespace arma;



// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


//[[Rcpp::export]]
arma::mat scocal(int j, arma::rowvec ng, arma::mat rgk, arma::rowvec a, arma::rowvec d, arma::rowvec bet, arma::mat gam, double maxtol,arma::mat X, arma::mat yallgroup, int y, int G, int r, int m, int eta)
{
  arma::cube Pstar=zeros<cube>(G,(m-1),y);
  arma::cube Qstar=zeros<cube>(G,(m-1),y);
  arma::cube P=zeros<cube>(G,m,y); 
  

  for(int yy=0; yy<y; yy++){
    for(int g = 0; g < G; g++){
      (Pstar.slice(yy)).row(g)=1/(1+exp(-(d+a*(X.row(g)).t()+ yallgroup.row(yy)*bet.t()+yallgroup.row(yy)*gam*X.row(g).t())));
      (P.slice(yy)).row(g)=-diff(join_horiz(ones<colvec>(1),(Pstar.slice(yy)).row(g),zeros<colvec>(1)),1,1);
    }
  }
  Qstar=ones<cube>(G,(m-1),y)-Pstar;
  arma::rowvec Dsco=zeros<rowvec>(m-1);
  for (int yy=0; yy<y; yy++){
    Dsco += sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy));
  }
  arma::cube PQdif=zeros<cube>(G,m,y); 
  for (int yy=0; yy<y; yy++){
    PQdif.slice(yy)=-diff(join_horiz(zeros<colvec>(G),Pstar.slice(yy)%Qstar.slice(yy),zeros<colvec>(G)),1,1);
  }
  
  int len=0;
  for (int kk=0; kk<r; kk++){
    if (a(kk)!=0){
      len++;
    }
  }
  
  arma::rowvec Ascoall=zeros<rowvec>(a.n_elem);
  arma::rowvec a01=a;
  for (int kk=0; kk<r; kk++){
    if (a(kk)!=0){
      a01(kk)=1;
    }
  }
  for (int yy=0; yy<y; yy++){
    Ascoall += (sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1).t()*X);
  }
  arma::vec Asco=Ascoall(find(a!=0));
  int len2=0;
  for (int kk=0; kk<r; kk++){
    for (int mm=0; mm<2; mm++){
      if (gam(kk,mm)!=0){
        len2++;
      }
    }
  }
  arma::mat Gamscoall=gam;
  for (int yy=1; yy<y; yy++){
    Gamscoall.row(yy-1)=(sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1)).t()*X;
  }
  for (int kk=0; kk<r; kk++){
    for (int mm=0; mm<(y-1); mm++){
      if (gam(kk,mm)==0){
        Gamscoall(kk,mm)=0;
      }
    }
  }
  arma::vec Gamsco=Gamscoall.elem(find(gam!=0));
  int len3=0;
  for (int mm=0; mm<2; mm++){
    if (bet(mm)!=0){
      len3++;
    }
  }
  arma::rowvec Betscoall=bet;
  for (int yy=1; yy<y; yy++){
    Betscoall.subvec(yy-1,yy-1) = sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy),0);
  }
  for (int mm=0; mm<2; mm++){
    if (bet(mm)==0){
      Betscoall(mm)=0;
    }
  }
  arma::vec Betsco=Betscoall.elem(find(Betscoall!=0));
  arma::rowvec minusgrad = -join_horiz(Dsco,Asco.t(),Gamsco.t(),Betsco.t());
  
  
  
  
  return(minusgrad);
}
