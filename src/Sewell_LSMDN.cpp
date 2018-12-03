
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
SEXP c_initialize1( SEXP Data, SEXP DIMS, SEXP Yy,
                    SEXP XSCALE, SEXP BETAIN,
                    SEXP BETAOUT, SEXP WW ) {
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double Xscale= Rcpp::as<double>(XSCALE);
  
  double ret =0,dx=0, eta=0;
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      for(int j = 0; j < dims(0); j++)
      {
        if(i != j)
        {
          dx = Xscale*arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
          eta = (BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)));
          ret += Y.slice(tt)(i,j)*eta-
            log(1+exp(eta));
        }
      }
    }
  }
  
  return wrap(ret);
  
  END_RCPP
}

// [[Rcpp::export]]
SEXP c_initialize1_grad( SEXP Data, SEXP DIMS, SEXP Yy,
                         SEXP XSCALE, SEXP BETAIN,
                         SEXP BETAOUT, SEXP WW ){
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double Xscale= Rcpp::as<double>(XSCALE);
  
  double dx=0, eta=0;
  Rcpp::NumericVector ret(3);
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      for(int j = 0; j < dims(0); j++)
      {
        if(i != j)
        {
          dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
          eta = BIN*(1-Xscale*dx/ww(j))+BOUT*(1-Xscale*dx/ww(i));
          ret(0) = ret(0)+ dx*(BIN/ww(j)+BOUT/ww(i))*(1/(1+exp(-eta))-Y(i,j,tt));
          ret(1) = ret(1)+ (1-Xscale*dx/ww(j))*(Y(i,j,tt)-1/(1+exp(-eta)));
          ret(2) = ret(2)+ (1-Xscale*dx/ww(i))*(Y(i,j,tt)-1/(1+exp(-eta)));
        }
      }
    }
  }
  
  return wrap(ret);
  
  END_RCPP
}

// [[Rcpp::export]]
SEXP c_update2( SEXP Xitm1, SEXP DIMS, SEXP TUNEX, SEXP Yy,
                SEXP BETAIN, SEXP BETAOUT, SEXP TUNEBIO,
                SEXP WW, SEXP t2X, SEXP s2X,
                SEXP xiBIN, SEXP xiBOUT, SEXP nuBIN,
                SEXP nuBOUT, SEXP CAUCHY,
                SEXP RNORMS, SEXP RNORMSBIO,
                SEXP ELOUT, SEXP ELIN, SEXP SUBSEQ, SEXP DEG ){
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT,dinmax,doutmax)
   Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  int Cauchy = Rcpp::as<int>(CAUCHY);
  
  Rcpp::NumericVector XOLD(Xitm1);
  arma::cube Xold(XOLD.begin(),dims[1],dims[2],dims[0]);
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);
  for(int i=0;i<dims(0);i++)
  {
    Xnew.slice(i)=Xold.slice(i);
  }
  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); 
  arma::colvec centerings = arma::zeros(dims(1),1);
  
  Rcpp::NumericVector rnormsVec(RNORMS);
  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  arma::colvec rnormsBIO = Rcpp::as<arma::colvec>(RNORMSBIO);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  Rcpp::NumericVector ELinVec(ELIN);
  arma::cube ELin(ELinVec.begin(),dims(0),dims(3),dims(2));
  Rcpp::NumericVector ELoutVec(ELOUT);
  arma::cube ELout(ELoutVec.begin(),dims(0),dims(4),dims(2));
  Rcpp::IntegerMatrix subseq(SUBSEQ);
  int n0 = subseq.ncol(); 
  Rcpp::NumericVector DegrVec(DEG);
  arma::cube degr(DegrVec.begin(),dims(0),2,dims(2));
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  double xiBin = Rcpp::as<double>(xiBIN);
  double xiBout = Rcpp::as<double>(xiBOUT);
  double nuBin = Rcpp::as<double>(nuBIN);
  double nuBout = Rcpp::as<double>(nuBOUT);
  double BinNew =0, BoutNew =0;
  double tunex = Rcpp::as<double>(TUNEX);
  double tuneBIO = Rcpp::as<double>(TUNEBIO);
  
  double t2 = Rcpp::as<double>(t2X);
  double s2 = Rcpp::as<double>(s2X);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double AccProb =0, contrOld=0, contrNew=0;
  double dz=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+3,1);
  
  
  
  //-------------------- Latent Positions-------------------
  
  for(int tt=0;tt < dims[2]; tt++)
  {
    for(int i=0;i<dims[0];i++)
    {
      AccProb=0;
      if(Cauchy<0.5)
      {
        //---Normal Random Walk---
        //  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*arma::randn(dims(1),1);
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*rnorms.slice(i).col(tt);
      }else{
        //---Cauchy Random Walk---
        for(int ell=0;ell<dims(1);ell++)
        {
          uu = arma::randu();
          Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + tunex*tan(PI*(uu-0.5));
        }
      }
      
      //IN edges
      if(degr(i,0,tt)>0)
      {
        for(int j=0;j<degr(i,0,tt);j++)
        {
          dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELin(i,j,tt)-1).col(tt),2);
          dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(ELin(i,j,tt)-1).col(tt),2);
          AccProb += (dx-dz)*(BIN/ww(i)+BOUT/ww(ELin(i,j,tt)-1));
        }
      }
      //OUT edges
      if(degr(i,1,tt)>0)
      {
        for(int j=0;j<degr(i,1,tt);j++)
        {
          dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
          dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
          AccProb += (dx-dz)*(BIN/ww(ELout(i,j,tt)-1)+BOUT/ww(i));
        }
      }
      //Control estimate
      
      contrOld=0;contrNew=0;
      for(int j=0;j<n0;j++)
      {
        dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
        dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
        contrNew += log(1+exp(BIN*(1-dz/ww(i))+BOUT*(1-dz/ww(subseq(i,j)-1))))+
          log(1+exp(BIN*(1-dz/ww(subseq(i,j)-1))+BOUT*(1-dz/ww(i))));
        contrOld += log(1+exp(BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(subseq(i,j)-1))))+
          log(1+exp(BIN*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
      }
      AccProb += (dims(0)-1)/n0*(contrOld-contrNew);
      
      
      if(tt==0)
      {
        insides = trans(Xnew.slice(i).col(tt))*(Xnew.slice(i).col(tt))/t2;
        AccProb += -0.5*insides(0,0);
        insides = trans(Xold.slice(i).col(tt))*(Xold.slice(i).col(tt))/t2;
        AccProb -= -0.5*insides(0,0);
      }
      if(tt>0)
      {
        muit = Xnew.slice(i).col(tt-1);
        insides = trans(Xnew.slice(i).col(tt)-muit)*(Xnew.slice(i).col(tt)-muit)/s2;
        AccProb += -0.5*insides(0,0);
        insides = trans(Xold.slice(i).col(tt)-muit)*(Xold.slice(i).col(tt)-muit)/s2;
        AccProb -= -0.5*insides(0,0);  
      }
      if(tt <dims[2]-1)
      {
        muit = Xnew.slice(i).col(tt);
        insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
        AccProb += -0.5*insides(0,0);
        muit = Xold.slice(i).col(tt);
        insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
        AccProb -= -0.5*insides(0,0);  
      }
      
      uu= arma::randu();
      if(uu<exp(AccProb))
      {
        AccRate(3+tt*dims(0)+i) = 1;
      }else
      {
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
      }
      
    }
  }
  
  /*---Centering---*/
  for(int i=0;i<dims(0);i++)
  {
    for(int tt=0;tt<dims(2);tt++)
    {
      Xconc.row(i*dims(2)+tt) = trans(Xnew.slice(i).col(tt));
    }
  }
  for(int ell=0;ell<dims(1);ell++)
  {
    centerings(ell) = sum(Xconc.col(ell))/(dims(0)*dims(2));
  }
  for(int i=0;i<dims(0);i++)
  {
    for(int tt=0;tt<dims(2);tt++)
    {
      for(int ell=0;ell<dims(1);ell++)
      {
        Xnew.slice(i)(ell,tt) = Xnew.slice(i)(ell,tt) - centerings(ell);
      }
    }
  }
  
  
  /*-------------------- BetaIn and BetaOut-------------------*/
  
  AccProb=0;
  if(Cauchy<0.5)
  {
    //  BinNew = BIN + tuneBIO*arma::randn();
    BinNew = BIN + tuneBIO*rnormsBIO(0);
  }else{
    uu = arma::randu();
    BinNew = BIN + tuneBIO*tan(PI*(uu-0.5));
  }
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      
      //OUT edges
      if(degr(i,1,tt)>0)
      {
        for(int j=0;j<degr(i,1,tt);j++)
        {
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
          AccProb += (BinNew-BIN)*(1-dx/ww(ELout(i,j,tt)-1));
        }
      }
      //Control estimate
      contrOld=0;contrNew=0;
      for(int j=0;j<n0;j++)
      {
        dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
        contrNew += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
        contrOld += log(1+exp(BIN*(1-dx/ww(subseq(i,j)-1))   +BOUT*(1-dx/ww(i))));
      }
      AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);
      
    }
  }
  
  AccProb += -0.5*(BinNew-nuBin)*(BinNew-nuBin)/xiBin;
  AccProb -= -0.5*(BIN-nuBin)*(BIN-nuBin)/xiBin;
  
  uu= arma::randu();
  if(uu<exp(AccProb))
  {
    AccRate(0) = 1;
  }else
  {
    BinNew = BIN;
  }
  
  
  AccProb=0;
  if(Cauchy<0.5)
  {
    //  BoutNew = BOUT + tuneBIO*arma::randn();
    BoutNew = BOUT + tuneBIO*rnormsBIO(1);
  }else{
    uu = arma::randu();
    BoutNew = BOUT + tuneBIO*tan(PI*(uu-0.5));
  }
  
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      
      //OUT edges
      if(degr(i,1,tt)>0)
      {
        for(int j=0;j<degr(i,1,tt);j++)
        {
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
          AccProb += (BoutNew-BOUT)*(1-dx/ww(i));
        }
      }
      //Control estimate
      contrOld=0;contrNew=0;
      for(int j=0;j<n0;j++)
      {
        dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
        contrNew += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BoutNew*(1-dx/ww(i))));
        contrOld += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
      }
      AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);
      
    }
  }
  
  AccProb += -0.5*(BoutNew-nuBout)*(BoutNew-nuBout)/xiBout;
  AccProb -= -0.5*(BOUT-nuBout)*(BOUT-nuBout)/xiBout;
  
  uu= arma::randu();
  if(uu<exp(AccProb))
  {
    AccRate(1) = 1;
  }else
  {
    BoutNew = BOUT;
  }
  
  
  return Rcpp::List::create(Xnew,BinNew,BoutNew,AccRate);
  
  END_RCPP
}

// [[Rcpp::export]]
SEXP c_update1( SEXP Xitm1, SEXP DIMS, SEXP TUNEX, SEXP Yy,
                SEXP BETAIN, SEXP BETAOUT, SEXP TUNEBIO,
                SEXP WW, SEXP t2X, SEXP s2X,
                SEXP xiBIN, SEXP xiBOUT, SEXP nuBIN,
                SEXP nuBOUT, SEXP CAUCHY,
                SEXP RNORMS, SEXP RNORMSBIO){
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT,K)
   Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  int Cauchy = Rcpp::as<int>(CAUCHY);
  
  Rcpp::NumericVector XOLD(Xitm1);
  arma::cube Xold(XOLD.begin(),dims[1],dims[2],dims[0]);
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);
  for(int i=0;i<dims(0);i++)
  {
    Xnew.slice(i)=Xold.slice(i);
  }
  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); 
  arma::colvec centerings = arma::zeros(dims(1),1);
  
  Rcpp::NumericVector rnormsVec(RNORMS);
  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  arma::colvec rnormsBIO = Rcpp::as<arma::colvec>(RNORMSBIO);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  double xiBin = Rcpp::as<double>(xiBIN);
  double xiBout = Rcpp::as<double>(xiBOUT);
  double nuBin = Rcpp::as<double>(nuBIN);
  double nuBout = Rcpp::as<double>(nuBOUT);
  double BinNew =0, BoutNew =0;
  double tunex = Rcpp::as<double>(TUNEX);
  double tuneBIO = Rcpp::as<double>(TUNEBIO);
  
  double t2 = Rcpp::as<double>(t2X);
  double s2 = Rcpp::as<double>(s2X);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double AccProb =0;
  double dz=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+3,1);
  
  
  /*-------------------- Latent Positions-------------------*/
  
  for(int tt=0;tt < dims[2]; tt++)
  {
    for(int i=0;i<dims[0];i++)
    {
      AccProb=0;
      if(Cauchy<0.5)
      {
        /*---Normal Random Walk---*/
        //  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*arma::randn(dims(1),1);
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*rnorms.slice(i).col(tt);
      }else{
        /*---Cauchy Random Walk---*/
        for(int ell=0;ell<dims(1);ell++)
        {
          uu = arma::randu();
          Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + tunex*tan(PI*(uu-0.5));
        }
      }
      
      
      for(int j=0;j<dims[0];j++)
      {
        if(j != i){
          dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
          dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
          AccProb += (dx-dz)*(Y.slice(tt)(j,i)*(BIN/ww(i)+BOUT/ww(j))+
            Y.slice(tt)(i,j)*(BIN/ww(j)+BOUT/ww(i)));
          AccProb += log(1+exp(BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(j))));
          AccProb += log(1+exp(BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i))));
          AccProb -= log(1+exp(BIN*(1-dz/ww(i))+BOUT*(1-dz/ww(j))));
          AccProb -= log(1+exp(BIN*(1-dz/ww(j))+BOUT*(1-dz/ww(i))));
        }
      }
      if(tt==0)
      {
        insides = trans(Xnew.slice(i).col(tt))*(Xnew.slice(i).col(tt))/t2;
        AccProb += -0.5*insides(0,0);
        insides = trans(Xold.slice(i).col(tt))*(Xold.slice(i).col(tt))/t2;
        AccProb -= -0.5*insides(0,0);
      }
      if(tt>0)
      {
        muit = Xnew.slice(i).col(tt-1);
        insides = trans(Xnew.slice(i).col(tt)-muit)*(Xnew.slice(i).col(tt)-muit)/s2;
        AccProb += -0.5*insides(0,0);
        insides = trans(Xold.slice(i).col(tt)-muit)*(Xold.slice(i).col(tt)-muit)/s2;
        AccProb -= -0.5*insides(0,0);  
      }
      if(tt <dims[2]-1)
      {
        muit = Xnew.slice(i).col(tt);
        insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
        AccProb += -0.5*insides(0,0);
        muit = Xold.slice(i).col(tt);
        insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
        AccProb -= -0.5*insides(0,0);  
      }
      
      uu= arma::randu();
      if(uu<exp(AccProb))
      {
        AccRate(3+tt*dims(0)+i) = 1;
      }else
      {
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
      }
      
    }
  }
  
  /*---Centering---*/
  for(int i=0;i<dims(0);i++)
  {
    for(int tt=0;tt<dims(2);tt++)
    {
      Xconc.row(i*dims(2)+tt) = trans(Xnew.slice(i).col(tt));
    }
  }
  for(int ell=0;ell<dims(1);ell++)
  {
    centerings(ell) = sum(Xconc.col(ell))/(dims(0)*dims(2));
  }
  for(int i=0;i<dims(0);i++)
  {
    for(int tt=0;tt<dims(2);tt++)
    {
      for(int ell=0;ell<dims(1);ell++)
      {
        Xnew.slice(i)(ell,tt) = Xnew.slice(i)(ell,tt) - centerings(ell);
      }
    }
  }
  
  
  /*-------------------- BetaIn and BetaOut-------------------*/
  AccProb=0;
  if(Cauchy<0.5)
  {
    //  BinNew = BIN + tuneBIO*arma::randn();
    BinNew = BIN + tuneBIO*rnormsBIO(0);
  }else{
    uu = arma::randu();
    BinNew = BIN + tuneBIO*tan(PI*(uu-0.5));
  }
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      for(int j = 0; j < dims(0); j++)
      {
        if(i != j)
        {
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
          AccProb += Y.slice(tt)(i,j)*(BinNew-BIN)*(1-dx/ww(j)) + 
            log(1+exp(BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)))) -
            log(1+exp(BinNew*(1-dx/ww(j))+BOUT*(1-dx/ww(i))));
        }
      }
    }
  }
  
  AccProb += -0.5*(BinNew-nuBin)*(BinNew-nuBin)/xiBin;
  AccProb -= -0.5*(BIN-nuBin)*(BIN-nuBin)/xiBin;
  
  uu= arma::randu();
  if(uu<exp(AccProb))
  {
    AccRate(0) = 1;
  }else
  {
    BinNew = BIN;
  }
  
  AccProb=0;
  if(Cauchy<0.5)
  {
    //  BoutNew = BOUT + tuneBIO*arma::randn();
    BoutNew = BOUT + tuneBIO*rnormsBIO(1);
  }else{
    uu = arma::randu();
    BoutNew = BOUT + tuneBIO*tan(PI*(uu-0.5));
  }
  
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      for(int j = 0; j < dims(0); j++)
      {
        if(i != j)
        {
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
          AccProb += Y.slice(tt)(i,j)*(BoutNew-BOUT)*(1-dx/ww(i)) + 
            log(1+exp(BinNew*(1-dx/ww(j))+BOUT*(1-dx/ww(i)))) -
            log(1+exp(BinNew*(1-dx/ww(j))+BoutNew*(1-dx/ww(i))));
        }
      }
    }
  }
  
  AccProb += -0.5*(BoutNew-nuBout)*(BoutNew-nuBout)/xiBout;
  AccProb -= -0.5*(BOUT-nuBout)*(BOUT-nuBout)/xiBout;
  
  uu= arma::randu();
  if(uu<exp(AccProb))
  {
    AccRate(1) = 1;
  }else
  {
    BoutNew = BOUT;
  }
  
  
  return Rcpp::List::create(Xnew,BinNew,BoutNew,AccRate);
  
  END_RCPP
}

// [[Rcpp::export]]
SEXP c_t2s2Parms( SEXP DATA, SEXP DIMS, SEXP THETAT,
                  SEXP THETAS, SEXP PHIT, SEXP PHIS ){
  BEGIN_RCPP
  
  // Dims is c(n,p,TT,K)
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(DATA);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  double thetaT=Rcpp::as<double>(THETAT);
  double phiT=Rcpp::as<double>(PHIT);
  double thetaS=Rcpp::as<double>(THETAS);
  double phiS=Rcpp::as<double>(PHIS);
  double shapeT=0, scaleT=0, shapeS=0, scaleS=0;
  
  arma::mat insides = arma::zeros(1,1);
  
  //---------------------------------------
  
  shapeT = thetaT + 0.5*dims(0)*dims(1);
  scaleT = phiT;
  shapeS = thetaS + 0.5*dims(0)*dims(1)*(dims(2)-1);
  scaleS = phiS;
  for(int i =0;i<dims(0);i++)
  {
    insides = 0.5*trans(X.slice(i).col(0))*(X.slice(i).col(0));
    scaleT += insides(0,0);
    for(int tt=1;tt<dims(2);tt++)
    {
      insides = 0.5*trans(X.slice(i).col(tt)-X.slice(i).col(tt-1))*
        (X.slice(i).col(tt)-X.slice(i).col(tt-1));
      scaleS += insides(0,0);
    }
  }
  
  
  return Rcpp::List::create(shapeT,scaleT,shapeS,scaleS);
  
  END_RCPP
  }

// [[Rcpp::export]]
SEXP c_WAccProb2( SEXP Data, SEXP DIMS, SEXP Yy,
                  SEXP BETAIN, SEXP BETAOUT, SEXP TUNEW,
                  SEXP WWOld, SEXP WWNew,
                  SEXP ELOUT, SEXP ELIN, SEXP SUBSEQ, SEXP DEG ){
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT,dinmax,doutmax)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec wwOld = Rcpp::as<arma::colvec>(WWOld);
  arma::colvec wwNew = Rcpp::as<arma::colvec>(WWNew);
  double tuneW = Rcpp::as<double>(TUNEW);
  
  double AccProb =0, contrOld=0, contrNew=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
  
  Rcpp::NumericVector ELinVec(ELIN);
  arma::cube ELin(ELinVec.begin(),dims(0),dims(3),dims(2));
  Rcpp::NumericVector ELoutVec(ELOUT);
  arma::cube ELout(ELoutVec.begin(),dims(0),dims(4),dims(2));
  Rcpp::IntegerMatrix subseq(SUBSEQ);
  int n0 = subseq.ncol(); 
  Rcpp::NumericVector DegrVec(DEG);
  arma::cube degr(DegrVec.begin(),dims(0),2,dims(2));
  
  
  /*---------------------------------------*/
  
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      
      //OUT edges
      if(degr(i,1,tt)>0)
      {
        for(int j=0;j<degr(i,1,tt);j++)
        {
          dx = arma::norm(X.slice(i).col(tt)-X.slice(ELout(i,j,tt)-1).col(tt),2);
          AccProb += dx*(BIN*(1/wwOld(ELout(i,j,tt)-1)-1/wwNew(ELout(i,j,tt)-1)) + 
            BOUT*(1/wwOld(i)-1/wwNew(i)));
        }
      }
      //Control estimate
      
      contrOld=0;contrNew=0;
      for(int j=0;j<n0;j++)
      {
        dx = arma::norm(X.slice(i).col(tt)-X.slice(subseq(i,j)-1).col(tt),2);
        contrNew += log(1+exp(BIN*(1-dx/wwNew(subseq(i,j)-1))+BOUT*(1-dx/wwNew(i))));
        contrOld += log(1+exp(BIN*(1-dx/wwOld(subseq(i,j)-1))+BOUT*(1-dx/wwOld(i))));
      }
      AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);
      
    }
  }
  
  
  for(int i =0;i<dims(0);i++)
  {
    AccProb += (tuneW*wwNew(i)-1)*log(wwOld(i))- (tuneW*wwOld(i)-1)*log(wwNew(i)) -
      (tuneW*wwNew(i)-0.5)*log(tuneW*wwNew(i))-tuneW*wwNew(i) +
      (tuneW*wwOld(i)-0.5)*log(tuneW*wwOld(i))+tuneW*wwOld(i);
  }
  
  uu= arma::randu();
  if(uu<exp(AccProb))
  {
    AccRate = 1;
  }else
  {
    wwNew = wwOld;
  }
  
  return Rcpp::List::create(wwNew,AccRate);
  
  END_RCPP
}

// [[Rcpp::export]]
SEXP c_WAccProb1( SEXP Data, SEXP DIMS, SEXP Yy,
                  SEXP BETAIN, SEXP BETAOUT, SEXP TUNEW,
                  SEXP WWOld, SEXP WWNew ){
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT,K)
   Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec wwOld = Rcpp::as<arma::colvec>(WWOld);
  arma::colvec wwNew = Rcpp::as<arma::colvec>(WWNew);
  double tuneW = Rcpp::as<double>(TUNEW);
  
  double AccProb =0,dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
  
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
  {
    for(int i = 0; i < dims(0); i++)
    {
      for(int j = 0; j < dims(0); j++)
      {
        if(i != j)
        {
          dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
          AccProb += Y.slice(tt)(i,j)*dx*(BIN*(1/wwOld(j)-1/wwNew(j)) + 
            BOUT*(1/wwOld(i)-1/wwNew(i))) +
            log(1+exp(BIN*(1-dx/wwOld(j))+BOUT*(1-dx/wwOld(i)))) -
            log(1+exp(BIN*(1-dx/wwNew(j))+BOUT*(1-dx/wwNew(i))));
        }
      }
    }
  }
  
  for(int i =0;i<dims(0);i++)
  {
    AccProb += (tuneW*wwNew(i)-1)*log(wwOld(i))- (tuneW*wwOld(i)-1)*log(wwNew(i)) -
      (tuneW*wwNew(i)-0.5)*log(tuneW*wwNew(i))-tuneW*wwNew(i) +
      (tuneW*wwOld(i)-0.5)*log(tuneW*wwOld(i))+tuneW*wwOld(i);
  }
  
  uu= arma::randu();
  if(uu<exp(AccProb))
  {
    AccRate = 1;
  }else
  {
    wwNew = wwOld;
  }
  
  return Rcpp::List::create(wwNew,AccRate);
  
  END_RCPP
  }

// [[Rcpp::export]]
SEXP c_missing( SEXP Data, SEXP DIMS, SEXP MMM, SEXP Yy, SEXP Ttt,
                SEXP BETAIN, SEXP BETAOUT, SEXP WW ){
  BEGIN_RCPP
  
  /*Dims is c(n,p,TT,K); 
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  int ttt = Rcpp::as<int>(Ttt)-1;
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::IntegerVector MM(MMM);
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double dx=0, uu=0, Prob=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  
  /*---------------------------------------*/
  
  for(int i = 0; i < dims(0); i++)
  {
    if(std::find(MM.begin(),MM.end(),i) !=MM.end())
    {
      for(int j = 0; j < dims(0); j++)
      {
        if(i != j)
        {
          dx = arma::norm(X.slice(i).col(ttt)-X.slice(j).col(ttt),2);
          Prob = BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i));
          Prob = 1/(1+exp(-Prob));
          uu= arma::randu();
          if(uu<Prob)
          {
            Y(i,j,ttt) = 1;
          }else{
            Y(i,j,ttt) = 0;
          }
          /*  if(std::find(MM.begin(),MM.end(),j) ==MM.end())
          {
          Prob = BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(j));
          Prob = 1/(1+exp(-Prob));
          uu= arma::randu();
          if(uu<Prob)
          {
          Y(j,i,ttt) = 1;
          }else{
          Y(j,i,ttt) = 0;
          }
          }*/
          
        }
    }
  }
}
  
  return Rcpp::wrap(Y);
  
  END_RCPP
    }

// [[Rcpp::export]]
SEXP c_postzeroprob( SEXP Xi1, SEXP Xi2, SEXP Xj1, SEXP Xj2, SEXP SS2, SEXP LAM, SEXP PP0 ){
  BEGIN_RCPP
  
  double lam = Rcpp::as<double>(LAM);
  double pp0 = Rcpp::as<double>(PP0);
  Rcpp::NumericMatrix xi1(Xi1);
  Rcpp::NumericMatrix xi2(Xi2);
  Rcpp::NumericMatrix xj1(Xj1);
  Rcpp::NumericMatrix xj2(Xj2);
  Rcpp::NumericVector ss2(SS2);
  Rcpp::NumericVector normCDF(1);
  Rcpp::NumericVector temp1(1);
  Rcpp::NumericVector ret(1);
  
  double angle = 0, temp = 0, x = 0, y = 0;
  int nrows = xi1.nrow();
  int ncols = xi1.ncol();
  
  for(int iter = 0; iter < ncols; iter++)
  {
    temp = 0;
    for(int time = 1; time < nrows; time++)
    {
      x = xj1(time,iter)-xi1(time-1,iter);
      y = xj2(time,iter)-xi2(time-1,iter);
      angle = atan2(y,x);
      temp += (xi1(time,iter)-xi1(time-1,iter))*cos(angle)+(xi2(time,iter)-xi2(time-1,iter))*sin(angle);
    }
    temp1(0) = (lam*temp-ss2(iter))/(lam*sqrt(ss2(iter)*(nrows-1)));
    normCDF = pnorm(temp1);
    ret = ret + normCDF*sqrt(2*PI*ss2(iter)/(nrows-1))*exp(pow(lam*temp-ss2(iter),2)/(2*lam*lam*ss2(iter)*(nrows-1)));
  }
  ret = ret*(1-pp0)/(ncols*pp0*lam);
  
  return Rcpp::wrap(ret);
  
  END_RCPP
}

// [[Rcpp::export]]
SEXP c_prediction( SEXP EX, SEXP SIG2, SEXP X1T, SEXP X2T, SEXP BIN, SEXP BOUT, SEXP WW ){
  BEGIN_RCPP
  
  
  Rcpp::NumericMatrix Ex(EX);
  Rcpp::NumericMatrix x1(X1T);
  Rcpp::NumericMatrix x2(X2T);
  Rcpp::NumericMatrix ww(WW);
  Rcpp::NumericVector sig2(SIG2);
  Rcpp::NumericVector Bin(BIN);
  Rcpp::NumericVector Bout(BOUT);
  
  int nrows = Ex.nrow();
  int nIter = x1.ncol();
  Rcpp::NumericMatrix yhat(nrows,nrows);
  double dx, tempyhat, sumX, tempX1, tempX2;
  
  
  for(int i = 0; i < nrows-1; i++)
  {
    for( int j = i+1; j < nrows; j++)
    {
      sumX=0;
      dx = sqrt( pow(Ex(i,0)-Ex(j,0),2) + pow(Ex(i,1)-Ex(j,1),2) );
      for(int l = 0; l < nIter; l++)
      {
        
        tempX1 = 1/(2*PI*sig2(l))*exp(-1/(2*sig2(l))*(pow(Ex(i,0)-x1(i,l),2)+pow(Ex(i,1)-x2(i,l),2)) );
        tempX2 = 1/(2*PI*sig2(l))*exp(-1/(2*sig2(l))*(pow(Ex(j,0)-x1(j,l),2)+pow(Ex(j,1)-x2(j,l),2)) );
        sumX = sumX + tempX1*tempX2/10000;
        
        tempyhat = tempX1*tempX2/10000/(1+exp(Bin(l)*(dx/ww(j,l)-1)+Bout(l)*(dx/ww(i,l)-1)));
        yhat(i,j) = yhat(i,j) + tempyhat;
        tempyhat = tempX1*tempX2/10000/(1+exp(Bin(l)*(dx/ww(i,l)-1)+Bout(l)*(dx/ww(j,l)-1)));
        yhat(j,i) = yhat(j,i) + tempyhat;
      }
      yhat(i,j) = yhat(i,j)/sumX;
      yhat(j,i) = yhat(j,i)/sumX;
    }
  }
  
  return yhat;
  
  END_RCPP
}

