#include<RcppArmadillo.h>
#include <RcppEigen.h>

using namespace arma;



// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


//[[Rcpp::export]]
arma::mat M_step(int j, arma::rowvec ng, arma::mat rgk, arma::mat X, int y, int G, arma::mat yallgroup, double maxtol, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int eta)
{
  arma::rowvec d=grd.row(j-1);
  arma::rowvec a=gra.row(j-1);
  arma::mat gam=grgamma.slice(j-1);
  arma::rowvec bet=grbeta.row(j-1);
  arma::cube Pstar=zeros<cube>(G,(m-1),y);
  arma::cube Qstar=zeros<cube>(G,(m-1),y);
  arma::cube P=zeros<cube>(G,m,y); 
  int miter=0;
  for(miter=0; miter<500; miter++){
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
    arma::rowvec Asco=zeros<rowvec>(len);
    arma::rowvec a01=a;
    for (int kk=0; kk<r; kk++){
      if (a(kk)!=0){
        a01(kk)=1;
      }
    }
    for (int yy=0; yy<y; yy++){
      Asco += sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1).t()*(X*a01.t());
    }
    int len2=0;
    for (int kk=0; kk<r; kk++){
      for (int mm=0; mm<2; mm++){
        if (gam(kk,mm)!=0){
          len2++;
        }
      }
    }
    arma::vec Gamsco=zeros<vec>(len2);
    arma::mat Gamscoall=gam;
    for (int yy=1; yy<y; yy++){
      Gamscoall.row(yy-1)=(sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1)).t()*X;
    }
    for (int kk=0; kk<r; kk++){
      for (int mm=0; mm<2; mm++){
        if (gam(kk,mm)==0){
          Gamscoall(kk,mm)=0;
        }
      }
    }
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
    arma::rowvec minusgrad = join_horiz(Dsco,Asco,Gamscoall.elem(find(Gamscoall!=0)),Betscoall(find(Betscoall!=0)));
    
    arma::mat FI =zeros<arma::mat>( minusgrad.n_elem, minusgrad.n_elem);
    
    for (int yy=0; yy<y; yy++){
      ( FI.row(0)).subvec(0,0) = ( FI.row(0)).subvec(0,0)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%square(Pstar.slice(yy)%Qstar.slice(yy))%(1/(P.slice(yy)).col(1)+1/(P.slice(yy)).col(2))));
    }
    
    for (int yy=0; yy<y; yy++){
      ( FI.row(1)).subvec(1,1) = ( FI.row(1)).subvec(1,1)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X*a01.t())%(X*a01.t())%sum(square(PQdif.slice(yy))/P.slice(yy),1)));
      ( FI.row(1)).subvec(0,0) = ( FI.row(1)).subvec(0,0)+ sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X*a01.t())%Pstar.slice(yy)%Qstar.slice(yy)%((PQdif.slice(yy)).col(1)/(P.slice(yy)).col(1)-(PQdif.slice(yy)).col(2)/(P.slice(yy)).col(2)));
      FI(0,1)=FI(1,0);
    }
    
    for (int kk=2; kk<(y+1); kk++){
      ( FI.row(kk)).subvec(kk,kk) = (-sum(ng.subvec((kk*G),((kk+1)*G-1)).t()%(X*a01.t())%(X*a01.t())%sum(square(PQdif.slice(kk-1))/P.slice(kk-1),1)));
      FI(kk,1)=FI(kk,kk);
      FI(1,kk)=FI(kk,kk);
      ( FI.row(kk)).subvec(0,0)=sum(ng.subvec((kk*G),((kk+1)*G-1)).t()%(X*a01.t())%Pstar.slice(kk-1)%Qstar.slice(kk-1)%((PQdif.slice(kk-1)).col(1)/(P.slice(kk-1)).col(1)-(PQdif.slice(kk-1)).col(2)/(P.slice(kk-1)).col(2)));
      FI(0,kk)=FI(kk,0);
    }
    
    for (int kk=(2+len2); kk<(4+len2); kk++){
      ( FI.row(kk)).subvec(kk,kk) = (-sum(ng.subvec(((kk-2)*G),((kk-1)*G-1)).t()%square(Pstar.slice(kk-3)%Qstar.slice(kk-3))%(1/(P.slice(kk-3)).col(1)+1/(P.slice(kk-3)).col(2))));
      FI(kk,0)=FI(kk,kk);
      FI(0,kk)=FI(kk,kk);
      ( FI.row(kk)).subvec(1,1)=sum(ng.subvec(((kk-2)*G),((kk-1)*G-1)).t()%(X*a01.t())%Pstar.slice(kk-3)%Qstar.slice(kk-3)%((PQdif.slice(kk-3)).col(1)/(P.slice(kk-3)).col(1)-(PQdif.slice(kk-3)).col(2)/(P.slice(kk-3)).col(2)));
      FI(1,kk)=FI(kk,1);
      
    }
    
    arma::rowvec scoall01=ones<arma::rowvec>(2+2*(y-1));
    arma::rowvec gam01=gam*a01.t();
    gam01.elem(find(gam01!=0)).ones();
    scoall01.subvec(2,y)= gam01;
    arma::rowvec bet01=bet;
    bet01.elem(find(bet01!=0)).ones();
    scoall01.subvec(y+1,2*y-1)= bet01;
    arma::mat FI2= scoall01*FI*scoall01.t();
    
    
    
    
    if(sum(add)<maxtol){
      break;
    }
  }
}