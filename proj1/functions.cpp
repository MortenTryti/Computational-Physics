#include <iostream>
#include <math.h> 
#include <iomanip>
#include <armadillo>
#include <string>
#include <fstream>
#include "functions.hpp"

double f(double x)
{
  return 100*exp(-10*x);
}

double gfunc(double x, double h)
{
  double b = h*h*f(x);
  return b;
}

// declaring the forward function
int forwardback(arma::vec a, arma::vec b, arma::vec c, arma::vec g,arma::vec btilde,arma::vec gtilde, arma::vec &v){

    int N = a.size();

    btilde[0] = b[0];
    gtilde[0] = g[0];
    
    for (int i = 1; i<= N-1; i++){

        btilde[i] = b[i] - (a[i]/btilde[i-1])*c[i-1];
        gtilde[i] = g[i] - (a[i]/btilde[i-1])*gtilde[i-1];
        

    
    }


    //Back iterations, first declare the last point in v, then for looping us to the end 
    v[N-1] = gtilde[N-1]/btilde[N-1];
    
    
    for (int i = N-2; i>=0; i--){
      v[i] = (gtilde[i] - c[i]*v[i+1])/btilde[i];
      
    }

    return 0;
}

// declaring the special forward function from task 9
int forwardback_special(arma::vec g, arma::vec btilde, arma::vec gtilde, arma::vec &v){
    int N = g.size();

    btilde[0] = 2.0;
    gtilde[0] = g[0];
    for (int i = 1; i<= N-1; i++){
      
      btilde[i] = (i+2)/(i+1.0);
      gtilde[i] = g[i] + gtilde[i-1]/btilde[i-1];
      
     }


    //Back iterations, first declare the last point in v, then for looping us to the end 
    v[N-1] = gtilde[N-1]/btilde[N-1];
    
    
    for (int i = N-2; i>=0; i--){
      v[i] = (gtilde[i] + v[i+1])/btilde[i];
      
    }
  
    return 0;

}

