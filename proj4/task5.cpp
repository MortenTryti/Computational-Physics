#include <armadillo>
#include <random>
#include <algorithm>
#include "ising.hpp"

int L = 20;
double T1 = 1;
double T2 = 2.4;
int n = 10e4;
int n_burn = 0;
int main(){

std::vector<bool> boolvec = {true,false};
  
  
 //The following code finds the things of interest for an aligned spinsystem and a random spin system
  
for (int i=0; i<boolvec.size();i++){
  bool booli = boolvec.at(i); 
  //Prints out the current bool variable
  std::cout<<booli<<std::endl;
  //Defines the two systems for different T
  Ising ising1 = Ising(L,T1,n,booli);
  Ising ising2 = Ising(L,T2,n,booli);  
  
  //Does the MC algo with a set burn in
  ising1.MonteCarlo(n_burn);
  ising2.MonteCarlo(n_burn);  

  arma::vec m1_cumsum = arma::cumsum(arma::abs(ising1.M_));
  arma::vec eps1_cumsum = cumsum(ising1.E_);

  arma::vec m2_cumsum = arma::cumsum(arma::abs(ising2.M_));
  arma::vec eps2_cumsum = arma::cumsum(ising2.E_);
  
  arma::mat m_eps = arma::zeros(n, 5);
  //Saves the calculated variables
  for (int i=0; i<n; i++){
    m_eps(i,0) = i+1;
    m_eps(i,1) = eps1_cumsum(i)/(L*L*(i+1));
    m_eps(i,2) = m1_cumsum(i)/(L*L*(i+1));    
    m_eps(i,3) = eps2_cumsum(i)/(L*L*(i+1));    
    m_eps(i,4) = m2_cumsum(i)/(L*L*(i+1));
  }

//Saves to correct .bin file
if (booli == true){
  m_eps.save("task5True.bin");
  }
  
if (booli == false){
  m_eps.save("task5false.bin");
  }
  
  }
  
  return 0;
}
