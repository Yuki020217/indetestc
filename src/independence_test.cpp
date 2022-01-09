//independence_test.cpp

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//'@title Independence test of categorical variables
//'@description
//'@param a a matrix
//'@examples
//'chisq_test(matrix(c(98,38,289,67,41,262,13,8,57,18,12,30),4,3,b=T))
//'chisq_test(matrix(c(341,405,105,103,11,15),2,3,b=T))
//'@import Rcpp
//'@import RcppArmadillo
//'@export
// [[Rcpp::export]]
double fac(int n)
{
  if(n==1 || n==0){
    return 1;
  }
  else
  {
    return n*fac(n-1);
  }
}
// [[Rcpp::export]]
double temp(arma::mat a){
  int r = a.n_rows;
  int s = a.n_cols;

  arma::vec nj(s);
  arma::vec ni(r);
  int n = accu(a);
  double fisher_p_value = 0.0;

  int i,j;
  double temp1 = 1;
  double temp2 = 1;
  double temp3 = 1;
  double temp4 = 1;
  for(i=0;i<r;i++){
    ni(i) = accu(a.row(i));
    temp1 = temp1*fac(ni(i));
  }
  for(j=0;j<s;j++){
    nj(j) = accu(a.col(j));
    temp2 = temp2*fac(nj(j));
  }
  for(i=0;i<r;i++){
    for(j=0;j<s;j++){
      temp3 = temp3*fac(a(i,j));
    }
  }
  temp4 = fac(n);
  fisher_p_value = temp1*temp2/temp3/temp4;
  return fisher_p_value;
}

// [[Rcpp::export]]
double fisher(arma::mat a){
  int r = a.n_rows;
  int s = a.n_cols;
  arma::vec nj(s);
  arma::vec ni(r);
  double temp0 = 0;
  int min = 0;

  int i,j,k;
  for(i=0;i<r;i++){
    ni(i) = accu(a.row(i));
  }
  for(j=0;j<s;j++){
    nj(j) = accu(a.col(j));
  }
  min = std::min(ni(0),nj(0));
  for(k=a(0,0);k<min;k++){
    a(0,0) = k;
    a(0,1) = ni(0) - k;
    a(1,0) = nj(0) - k;
    a(1,1) = ni(1) - a(1,0);
    temp0 = temp0 + temp(a);
  }
  return temp0;
}

// [[Rcpp::export]]
List chisq_test(arma::mat a) {
  int r = a.n_rows;
  int s = a.n_cols;

  arma::vec nj(s);
  arma::vec ni(r);
  int n = accu(a);
  arma::mat mij(r,s);
  double X_squared = 0.0;
  double p_value = 0.0;
  double Likelihood_X_squared = 0.0;
  double Likelihood_p_value = 0.0;
  int df1;
  int df2;

  int i,j;
  double temp = 0;
  for(i=0;i<r;i++){
    ni(i) = accu(a.row(i));
  }
  for(j=0;j<s;j++){
    nj(j) = accu(a.col(j));
  }

  for(i=0;i<r;i++){
    for(j=0;j<s;j++){
      mij(i,j) = ni(i)*nj(j)/n;
      temp = (a(i,j)-mij(i,j))*(a(i,j)-mij(i,j))/mij(i,j);
      X_squared = X_squared + temp;
    }
    df1 = (r-1)*(s-1);
    p_value = R::pchisq(X_squared, df1, 0, 0);
  }
  for(i=0;i<r;i++){
    for(j=0;j<s;j++){
      temp = a(i,j)*log(a(i,j)/mij(i,j));
      Likelihood_X_squared = Likelihood_X_squared + temp;
    }
    df2 = (r-1)*(s-1);
    Likelihood_p_value = R::pchisq(Likelihood_X_squared, df2, 0, 0);
  }

  DataFrame result1= DataFrame::create(Named("X_squared") = X_squared, _["df"] = df1, _["p_value"] = p_value);
  DataFrame result2= DataFrame::create(Named("Likelihood_X_squared") = Likelihood_X_squared, _["df"] = df2, _["p_value"] = Likelihood_p_value);
  DataFrame result3= DataFrame::create(Named("p_value") = fisher(a));
  if(r == 2 & s == 2){
    return List::create(_["卡方检验"] = result1, _["似然比卡方检验"] = result2, _["Fisher精确性检验"] = result3);
  }
  else return List::create(_["卡方检验"] = result1, _["似然比卡方检验"] = result2);
}


