#include <armadillo>
#include <random>
#include <algorithm>
#include "ising.hpp"

int L = 2;
double T = 2.4;
int n = 100000;
int main(){
  Ising isingT = Ising(L,T,n,-1,false);
  isingT.one_cycle();
  isingT.A_.print();
  std::cout<<isingT.magnetization()<<std::endl;
  std::cout<<isingT.M_A<<std::endl;    
  return 0;
}
