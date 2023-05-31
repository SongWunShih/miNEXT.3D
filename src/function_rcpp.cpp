#include <Rcpp.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
using namespace Rcpp;


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
// this is created by Song Wun, Shih

// FD
// [[Rcpp::export]]
double FD_h_fn(int m1, int m2, int n1, int n2, int k1, int k2, double vi, double ai1, double ai2){
  double result;
  if (ai1 >= k1 && ai1 <= n1-m1+k1 && ai2 >= k2 && ai2 <= n2-m2+k2) {
    result = vi * exp(Rf_lchoose(ai1, k1) + Rf_lchoose(n1-ai1, m1-k1) - Rf_lchoose(n1, m1) + Rf_lchoose(ai2, k2) + Rf_lchoose(n2-ai2, m2-k2) - Rf_lchoose(n2, m2));
  } else {
    result = 0;
  }
  return result;
}
// [[Rcpp::export]]
double FD_h_fn_all(double m1, double m2, double n1, double n2, double vi, double ai1, double ai2, double q) {
    double result = 0;
    if (q == 0) {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += FD_h_fn(m1, m2, n1, n2, k1, k2, vi, ai1, ai2);
                }
            }
        }
    } else if (q == 1) {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += (k1 + k2) / (m1 + m2) * log((k1 + k2) / (m1 + m2)) * FD_h_fn(m1, m2, n1, n2, k1, k2, vi, ai1, ai2);
                }
            }
        }
        result = -result;
    } else {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += pow((k1 + k2) / (m1 + m2), q) * FD_h_fn(m1, m2, n1, n2, k1, k2, vi, ai1, ai2);
                }
            }
        }
    }
    return result;
}

// [[Rcpp::export]]
double FD_h_hat_fn(int m1, int m2, int n1, int n2, NumericVector vi, NumericVector ai1, NumericVector ai2, double q, int S) {
  double output = 0;
  for(int i = 0; i < S; i++) {
    output += FD_h_fn_all(m1,m2,n1,n2,vi[i],ai1[i],ai2[i],q);
  }
  return output;
}

// PD
// [[Rcpp::export]]
double PD_g_fn(int m1, int m2, int n1, int n2, int k1, int k2, double Li, double ai1, double ai2){
  double result;
  if (ai1 >= k1 && ai1 <= n1-m1+k1 && ai2 >= k2 && ai2 <= n2-m2+k2) {
    result = Li * exp(Rf_lchoose(ai1, k1) + Rf_lchoose(n1-ai1, m1-k1) - Rf_lchoose(n1, m1) + Rf_lchoose(ai2, k2) + Rf_lchoose(n2-ai2, m2-k2) - Rf_lchoose(n2, m2));
  } else {
    result = 0;
  }
  return result;
}
// [[Rcpp::export]]
double PD_g_fn_all(double m1, double m2, double n1, double n2, double Li, double ai1, double ai2, double q, double Tbar) {
    double result = 0;
    if (q == 0) {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += PD_g_fn(m1, m2, n1, n2, k1, k2, Li, ai1, ai2);
                }
            }
        }
    } else if (q == 1) {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += (k1 + k2) / ((m1 + m2) * Tbar) * log((k1 + k2) / ((m1 + m2) * Tbar)) * PD_g_fn(m1, m2, n1, n2, k1, k2, Li, ai1, ai2);
                }
            }
        }
        result = -result;
    } else {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += pow((k1 + k2) / ((m1 + m2) * Tbar), q) * PD_g_fn(m1, m2, n1, n2, k1, k2, Li, ai1, ai2);
                }
            }
        }
    }
    return result;
}

// [[Rcpp::export]]
double PD_g_hat_fn(int m1, int m2, int n1, int n2, NumericVector Li, NumericVector ai1, NumericVector ai2, double q, int S, double Tbar) {
  double output = 0;
  for(int i = 0; i < S; i++) {
    output += PD_g_fn_all(m1,m2,n1,n2,Li[i],ai1[i],ai2[i],q,Tbar);
  }
  return output;
}

// TD
// [[Rcpp::export]]
double TD_f_fn(int m1, int m2, int n1, int n2, int k1, int k2, double ai1, double ai2){
  double result;
  if (ai1 >= k1 && ai1 <= n1-m1+k1 && ai2 >= k2 && ai2 <= n2-m2+k2) {
    result = exp(Rf_lchoose(ai1, k1) + Rf_lchoose(n1-ai1, m1-k1) - Rf_lchoose(n1, m1) + Rf_lchoose(ai2, k2) + Rf_lchoose(n2-ai2, m2-k2) - Rf_lchoose(n2, m2));
  } else {
    result = 0;
  }
  return result;
}
// [[Rcpp::export]]
double TD_f_fn_all(double m1, double m2, double n1, double n2, double ai1, double ai2, double q) {
    double result = 0;
    if (q == 0) {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += TD_f_fn(m1, m2, n1, n2, k1, k2, ai1, ai2);
                }
            }
        }
    } else if (q == 1) {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += (k1 + k2) / (m1 + m2) * log((k1 + k2) / (m1 + m2)) * TD_f_fn(m1, m2, n1, n2, k1, k2, ai1, ai2);
                }
            }
        }
        result = -result;
    } else {
        for (int k2 = 0; k2 <= m2; k2++) {
            for (int k1 = 0; k1 <= m1; k1++) {
                if (k1 == 0 && k2 == 0) {
                    result += 0;
                } else {
                    result += pow((k1 + k2) / (m1 + m2), q) * TD_f_fn(m1, m2, n1, n2, k1, k2, ai1, ai2);
                }
            }
        }
    }
    return result;
}

// [[Rcpp::export]]
double TD_f_hat_fn(int m1, int m2, int n1, int n2, NumericVector ai1, NumericVector ai2, double q, int S) {
  double output = 0;
  for(int i = 0; i < S; i++) {
    output += TD_f_fn_all(m1,m2,n1,n2,ai1[i],ai2[i],q);
  }
  return output;
}

// q = 0 ext
// [[Rcpp::export]]
double h0_cpp(double pi1, double pi2, int m1,  int m2s, int n2){
  double output = 0;
  output = pow((1-pi1),m1)*pow((1-pi2),n2)*(1-pow((1-pi2),m2s));
  return output;
}

// [[Rcpp::export]]
double h0_cpp_PD(double pi1, double pi2, double Li, int m1, int m2s, int n2){
  double output = 0;
  output = Li*pow((1-pi1),m1)*pow((1-pi2),n2)*(1-pow((1-pi2),m2s));
  return output;
}

// [[Rcpp::export]]
double h0_cpp_FD(double pi1, double pi2, double vi, int m1, int m2s, int n2){
  double output = 0;
  output = vi*pow((1-pi1),m1)*pow((1-pi2),n2)*(1-pow((1-pi2),m2s));
  return output;
}

// [[Rcpp::export]]
double h0_hat_cpp(NumericVector pi1, NumericVector pi2, int m1, int m2s, int n1, int n2){
  double output_all= 0; 
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      // sumsh = sumsh +  h0(pi1_tmp[i],pi2_tmp[i],m1,m2s,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      sumsh = sumsh +  h0_cpp(pi1_tmp[i],pi2_tmp[i],m1,m2s,n2)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h0_cpp(0,pi2_tmp[i],m1,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sumsh+sumx0;
    //    output_sh = sumsh;
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(pi2>0)];
    double sum2 = 0;
    for(int i=0; i < pi2_tmp.size(); i++){
      sum2 = sum2 + h0_cpp(0,pi2_tmp[i],0,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// [[Rcpp::export]]
double h0_hat_cpp_PD(NumericVector pi1, NumericVector pi2, NumericVector Li_v, int m1, int m2s, int n1, int n2){
  double output_all= 0; 
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)];
    NumericVector Li_tmp = Li_v[(pi1>0) & (pi2>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      // sumsh = sumsh +  h0(pi1_tmp[i],pi2_tmp[i],m1,m2s,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      sumsh = sumsh +  h0_cpp_PD(pi1_tmp[i],pi2_tmp[i],Li_tmp[i],m1,m2s,n2)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)];
    Li_tmp = Li_v[(pi1==0) & (pi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h0_cpp_PD(0,pi2_tmp[i],Li_tmp[i],m1,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sumsh+sumx0;
    //    output_sh = sumsh;
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(pi2>0)];
    NumericVector Li_tmp = Li_v[(pi2>0)];
    double sum2 = 0;
    for(int i=0; i < pi2_tmp.size(); i++){
      sum2 = sum2 + h0_cpp_PD(0,pi2_tmp[i],Li_tmp[i],0,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// [[Rcpp::export]]
double h0_hat_cpp_FD(NumericVector pi1, NumericVector pi2, NumericVector vi_v, int m1, int m2s, int n1, int n2){
  double output_all= 0; 
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)];
    NumericVector vi_tmp = vi_v[(pi1>0) & (pi2>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      // sumsh = sumsh +  h0(pi1_tmp[i],pi2_tmp[i],m1,m2s,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      sumsh = sumsh +  h0_cpp_FD(pi1_tmp[i],pi2_tmp[i],vi_tmp[i],m1,m2s,n2)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)];
    vi_tmp = vi_v[(pi1==0) & (pi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h0_cpp_FD(0,pi2_tmp[i],vi_tmp[i],m1,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sumsh+sumx0;
    //    output_sh = sumsh;
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(pi2>0)];
    NumericVector vi_tmp = vi_v[(pi2>0)];
    double sum2 = 0;
    for(int i=0; i < pi2_tmp.size(); i++){
      sum2 = sum2 + h0_cpp_FD(0,pi2_tmp[i],vi_tmp[i],0,m2s,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// q = 1 ext
// [[Rcpp::export]]
double TD_theo(double pi1, double pi2,int m1,int m2,int k1,int k2){
  double result;
  if((pi1 == 0) & (k1==0)){
    result = Rf_dbinom(k2, m2 , pi2, 0);
    
  }
  else if ((pi2 == 0) & (k2==0)){
    result = Rf_dbinom(k1, m1 , pi1, 0);
  }
  else{
    result = Rf_dbinom(k1, m1 , pi1, 0)*Rf_dbinom(k2, m2 , pi2, 0);
  }
  if((k1 == 0) & (k2==0)){
    result = 0;
  }
  return result;
}
// [[Rcpp::export]]
double PD_theo(double pi1, double pi2, double Li, int m1,int m2,int k1,int k2){
  double result;
  if((pi1 == 0) & (k1==0)){
    result = Li * Rf_dbinom(k2, m2 , pi2, 0);
    
  }
  else if ((pi2 == 0) & (k2==0)){
    result = Li * Rf_dbinom(k1, m1 , pi1, 0);
  }
  else{
    result =Li * Rf_dbinom(k1, m1 , pi1, 0)*Rf_dbinom(k2, m2 , pi2, 0);
  }
  if((k1 == 0) & (k2==0)){
    result = 0;
  }
  return result;
}
// [[Rcpp::export]]
double FD_theo(double pi1, double pi2, double vi, int m1,int m2,int k1,int k2){
  double result;
  if((pi1 == 0) & (k1==0)){
    result = vi * Rf_dbinom(k2, m2 , pi2, 0);
    
  }
  else if ((pi2 == 0) & (k2==0)){
    result = vi * Rf_dbinom(k1, m1 , pi1, 0);
  }
  else{
    result = vi * Rf_dbinom(k1, m1 , pi1, 0)*Rf_dbinom(k2, m2 , pi2, 0);
  }
  if((k1 == 0) & (k2==0)){
    result = 0;
  }
  return result;
}

//TD
// [[Rcpp::export]]
double h1(double pi1, double pi2, double m1, double m2, double n1, double n2){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=0; k2 <= m2; k2++){
    for(int k1=0; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp1 = 0; }
      else{ tmp1 = tmp1 + (k1+k2)/(m1+m2)*log((k1+k2)/(m1+m2))*TD_theo(pi1,pi2,m1,m2,k1,k2); }
    }
  }
  tmp1 = -tmp1;
  
  for(int k2=0; k2 <= n2; k2++){
    for(int k1=0; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp2 = 0; }
      else{ tmp2 = tmp2 + (k1+k2)/(m1+n2)*log((k1+k2)/(m1+n2))*TD_theo(pi1,pi2,m1,n2,k1,k2); }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}

// [[Rcpp::export]]
double h1_assem2( double pi2,double m2, double n2){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=1; k2 <= m2; k2++){
    
    tmp1 = tmp1 + (k2)/(m2)*log((k2)/(m2))*TD_theo(0,pi2,0,m2,0,k2);
  }
  tmp1 = -tmp1;
  
  for(int k2=1; k2 <= n2; k2++){
    
    tmp2 = tmp2 + (k2)/(n2)*log((k2)/(n2))*TD_theo(0,pi2,0,n2,0,k2); 
    
  }
  tmp2 = -tmp2;
  double result = tmp1-tmp2;
  //Rcout << "tmp1  is " << tmp1 << std::endl;
  // Rcout << "tmp2  is " << tmp2 << std::endl;
  //Rcout << "h1_assem2  is " << result << std::endl;
  return result;
}

// [[Rcpp::export]]
double h1_hat(NumericVector pi1, NumericVector pi2, NumericVector xi1, NumericVector xi2, int m1, int m2, int n1, int n2){
  double output_all = 0;
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(xi1>0) & (xi2>0)];
    NumericVector pi2_tmp = pi2[(xi1>0) & (xi2>0)];
    
    double sumsh = 0;
    //Rcout << "h1_hat size=" << pi1_tmp.size() << std::endl;
    //    double sumsh_p = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1(pi1_tmp[i],pi2_tmp[i],m1,m2,n1,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      //      sumsh_p = sumsh_p +  h1(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,1)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(xi1==0) & (xi2>0)];
    pi2_tmp = pi2[(xi1==0) & (xi2>0)];
    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h1(0,pi2_tmp[i],m1,m2,n1,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(xi1>0) & (xi2==0)];
    pi2_tmp = pi2[(xi1>0) & (xi2==0)];
    double sumy0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumy0 = sumy0 +  h1(pi1_tmp[i],0,m1,m2,n1,n2)/(1-pow(1-pi1_tmp[i], n1));
    }
    output_all = sumsh+sumx0+sumy0;
    //    output_sh = sumsh_p;
    
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(xi2>0)];
    //Rcout << "The value pi2_tmp is " << pi2_tmp << std::endl;
    double sum2 = 0;
    
    for(int i=0; i < pi2_tmp.size(); i++){
      
      sum2 = sum2 +  h1_assem2(pi2_tmp[i],m2,n2)/(1-pow(1-pi2_tmp[i], n2));
      //Rcout << "The value sum2 z is " << z << std::endl;
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// PD
// [[Rcpp::export]]
double h1_PD(double pi1, double pi2, double Li, double Tbar, double m1, double m2, double n1, double n2){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=0; k2 <= m2; k2++){
    for(int k1=0; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp1 = 0; }
      else{ tmp1 = tmp1 + (k1+k2)/(m1+m2)/Tbar*log((k1+k2)/(m1+m2)/Tbar)*PD_theo(pi1,pi2,Li,m1,m2,k1,k2); }
    }
  }
  tmp1 = -tmp1;
  
  for(int k2=0; k2 <= n2; k2++){
    for(int k1=0; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp2 = 0; }
      else{ tmp2 = tmp2 + (k1+k2)/(m1+n2)/Tbar*log((k1+k2)/(m1+n2)/Tbar)*PD_theo(pi1,pi2,Li,m1,n2,k1,k2); }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}

  
// [[Rcpp::export]]
double h1_assem2_PD( double pi2, double Li, double Tbar, double m2, double n2){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=1; k2 <= m2; k2++){
    
    tmp1 = tmp1 + (k2)/(m2*Tbar)*log((k2)/(m2*Tbar))*PD_theo(0,pi2,Li,0,m2,0,k2);
  }
  tmp1 = -tmp1;
  
  for(int k2=1; k2 <= n2; k2++){
    
    tmp2 = tmp2 + (k2)/(n2*Tbar)*log((k2)/(n2*Tbar))*PD_theo(0,pi2,Li,0,n2,0,k2); 
    
  }
  tmp2 = -tmp2;
  double result = tmp1-tmp2;
  //Rcout << "tmp1  is " << tmp1 << std::endl;
  // Rcout << "tmp2  is " << tmp2 << std::endl;
  //Rcout << "h1_assem2  is " << result << std::endl;
  return result;
}

// [[Rcpp::export]]
double h1_hat_PD(NumericVector pi1, NumericVector pi2, NumericVector Li_v, double Tbar, NumericVector xi1, NumericVector xi2, int m1, int m2, int n1, int n2){
  double output_all = 0;
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(xi1>0) & (xi2>0)];
    NumericVector pi2_tmp = pi2[(xi1>0) & (xi2>0)];
    NumericVector Li_tmp = Li_v[(xi1>0) & (xi2>0)];
    
    double sumsh = 0;
    //Rcout << "h1_hat size=" << pi1_tmp.size() << std::endl;
    //    double sumsh_p = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1_PD(pi1_tmp[i],pi2_tmp[i],Li_tmp[i],Tbar,m1,m2,n1,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      //      sumsh_p = sumsh_p +  h1(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,1)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(xi1==0) & (xi2>0)];
    pi2_tmp = pi2[(xi1==0) & (xi2>0)];
    Li_tmp = Li_v[(xi1==0) & (xi2>0)];

    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h1_PD(0,pi2_tmp[i],Li_tmp[i],Tbar,m1,m2,n1,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(xi1>0) & (xi2==0)];
    pi2_tmp = pi2[(xi1>0) & (xi2==0)];
    Li_tmp = Li_v[(xi1>0) & (xi2==0)];

    double sumy0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumy0 = sumy0 +  h1_PD(pi1_tmp[i],0,Li_tmp[i],Tbar,m1,m2,n1,n2)/(1-pow(1-pi1_tmp[i], n1));
    }
    output_all = sumsh+sumx0+sumy0;
    //    output_sh = sumsh_p;
    
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(xi2>0)];
    NumericVector Li_tmp = Li_v[(xi2>0)];
    //Rcout << "The value pi2_tmp is " << pi2_tmp << std::endl;
    double sum2 = 0;
    
    for(int i=0; i < pi2_tmp.size(); i++){
      
      sum2 = sum2 +  h1_assem2_PD(pi2_tmp[i],Li_tmp[i],Tbar,m2,n2)/(1-pow(1-pi2_tmp[i], n2));
      //Rcout << "The value sum2 z is " << z << std::endl;
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// FD
// [[Rcpp::export]]
double h1_FD(double pi1, double pi2, double vi, double m1, double m2, double n1, double n2){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=0; k2 <= m2; k2++){
    for(int k1=0; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp1 = 0; }
      else{ tmp1 = tmp1 + (k1+k2)/(m1+m2)*log((k1+k2)/(m1+m2))*FD_theo(pi1,pi2,vi,m1,m2,k1,k2); }
    }
  }
  tmp1 = -tmp1;
  
  for(int k2=0; k2 <= n2; k2++){
    for(int k1=0; k1 <= m1; k1++){
      if((k1 == 0) & (k2 == 0)){ tmp2 = 0; }
      else{ tmp2 = tmp2 + (k1+k2)/(m1+n2)*log((k1+k2)/(m1+n2))*FD_theo(pi1,pi2,vi,m1,n2,k1,k2); }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}

  
// [[Rcpp::export]]
double h1_assem2_FD(double pi2, double vi, double m2, double n2){
  
  double tmp1 = 0;
  double tmp2 = 0;
  
  for(int k2=1; k2 <= m2; k2++){
    
    tmp1 = tmp1 + (k2)/(m2)*log((k2)/(m2))*FD_theo(0,pi2,vi,0,m2,0,k2);
  }
  tmp1 = -tmp1;
  
  for(int k2=1; k2 <= n2; k2++){
    
    tmp2 = tmp2 + (k2)/(n2)*log((k2)/(n2))*FD_theo(0,pi2,vi,0,n2,0,k2); 
    
  }
  tmp2 = -tmp2;
  double result = tmp1-tmp2;
  //Rcout << "tmp1  is " << tmp1 << std::endl;
  // Rcout << "tmp2  is " << tmp2 << std::endl;
  //Rcout << "h1_assem2  is " << result << std::endl;
  return result;
}

// [[Rcpp::export]]
double h1_hat_FD(NumericVector pi1, NumericVector pi2, NumericVector vi_v, NumericVector xi1, NumericVector xi2, int m1, int m2, int n1, int n2){
  double output_all = 0;
  //  double output_sh = 0;
  if(m1>0){
    NumericVector pi1_tmp = pi1[(xi1>0) & (xi2>0)];
    NumericVector pi2_tmp = pi2[(xi1>0) & (xi2>0)];
    NumericVector vi_tmp = vi_v[(xi1>0) & (xi2>0)];
    
    double sumsh = 0;
    //Rcout << "h1_hat size=" << pi1_tmp.size() << std::endl;
    //    double sumsh_p = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1_FD(pi1_tmp[i],pi2_tmp[i],vi_tmp[i],m1,m2,n1,n2)/(1-pow(1-pi1_tmp[i], n1))/(1-pow(1-pi2_tmp[i], n2));
      //      sumsh_p = sumsh_p +  h1(pi1_tmp[i],pi2_tmp[i],t1,t2,UT1,UT2,T1,T2,1)/(1-pow(1-pi1_tmp[i], T1))/(1-pow(1-pi2_tmp[i], T2));
    }
    
    pi1_tmp = pi1[(xi1==0) & (xi2>0)];
    pi2_tmp = pi2[(xi1==0) & (xi2>0)];
    vi_tmp = vi_v[(xi1==0) & (xi2>0)];

    double sumx0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumx0 = sumx0 +  h1_FD(0,pi2_tmp[i],vi_tmp[i],m1,m2,n1,n2)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(xi1>0) & (xi2==0)];
    pi2_tmp = pi2[(xi1>0) & (xi2==0)];
    vi_tmp = vi_v[(xi1>0) & (xi2==0)];

    double sumy0 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumy0 = sumy0 +  h1_FD(pi1_tmp[i],0,vi_tmp[i],m1,m2,n1,n2)/(1-pow(1-pi1_tmp[i], n1));
    }
    output_all = sumsh+sumx0+sumy0;
    //    output_sh = sumsh_p;
    
  }
  else if(m1 == 0){
    NumericVector pi2_tmp = pi2[(xi2>0)];
    NumericVector vi_tmp = vi_v[(xi2>0)];
    //Rcout << "The value pi2_tmp is " << pi2_tmp << std::endl;
    double sum2 = 0;
    
    for(int i=0; i < pi2_tmp.size(); i++){
      
      sum2 = sum2 +  h1_assem2_FD(pi2_tmp[i],vi_tmp[i],m2,n2)/(1-pow(1-pi2_tmp[i], n2));
      //Rcout << "The value sum2 z is " << z << std::endl;
    }
    output_all = sum2;
    //    output_sh = sum2;
  }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// [[Rcpp::export]]
NumericVector un_abun(NumericVector xi,int n, int m){
  int s = xi.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = 1-exp(Rf_lchoose((n-xi[i]),m)-Rf_lchoose(n,m));
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector sh_abun(NumericVector xi1, NumericVector xi2, int n1, int m1, int n2, int m2){
  int s = xi1.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = (1-exp(Rf_lchoose((n1-xi1[i]),m1)-Rf_lchoose(n1,m1)) * exp(Rf_lchoose((n2-xi2[i]),m2)-Rf_lchoose(n2,m2)) );
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector un_abun_PD(NumericVector Li,NumericVector xi,int n, int m){
  int s = xi.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = Li[i]*(1-exp(Rf_lchoose((n-xi[i]),m)-Rf_lchoose(n,m)));
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector sh_abun_PD(NumericVector Li, NumericVector xi1, NumericVector xi2, int n1, int m1, int n2, int m2){
  int s = xi1.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = Li[i]*(1-exp(Rf_lchoose((n1-xi1[i]),m1)-Rf_lchoose(n1,m1)) * exp(Rf_lchoose((n2-xi2[i]),m2)-Rf_lchoose(n2,m2)));
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector un_abun_FD(NumericVector vi,NumericVector xi,int n, int m){
  int s = xi.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = vi[i]*(1-exp(Rf_lchoose((n-xi[i]),m)-Rf_lchoose(n,m)));
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector sh_abun_FD(NumericVector vi, NumericVector xi1, NumericVector xi2, int n1, int m1, int n2, int m2){
  int s = xi1.size();
  NumericVector out(s, 0.0);
  
  for(int i=0;i<s;i++){
    out[i] = vi[i]*(1-exp(Rf_lchoose((n1-xi1[i]),m1)-Rf_lchoose(n1,m1)) * exp(Rf_lchoose((n2-xi2[i]),m2)-Rf_lchoose(n2,m2)));
  }
  return(out);
}