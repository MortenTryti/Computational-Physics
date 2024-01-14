#include <armadillo>
#include <random>
#include <algorithm>
#include "ising.hpp"

int L = 20;
double T1 = 1;
double T2 = 2.4;
int n = 10000000;
int n_burn = 0.1*n;
int main(){
  
  Ising ising1 = Ising(L,T1,n,55,false);
  Ising ising2 = Ising(L,T2,n,55,false);  
  
  ising1.MonteCarlo(n_burn);
  ising2.MonteCarlo(n_burn);  
  
  arma::vec eps1 = ising1.E_/(L*L);
  arma::vec eps2 = ising2.E_/(L*L);
  
  arma::mat epsilon = arma::zeros(n, 2);
  
  for (int i=0; i<n; i++){
    epsilon(i,0) = eps1(i);
    epsilon(i,1) = eps2(i);
  }

  epsilon.save("eps.bin");
  
  return 0;
}
