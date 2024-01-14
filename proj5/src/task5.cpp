#include <armadillo>
#include <iostream>
#include <complex>
#include <cmath>
#include "../include/header.hpp"

int main(int argc, char* argv[]){
  //Nr holes
  int nh = std::atoi(argv[1]);
  std::string nhS = argv[1];
  //Setting up the potential
  arma::mat M = pot_V(201,1.0,nh);
  //Saving the potential to then be plotted in plot5.py
  M.save("data/M"+nhS+".bin");
  return 0;
}
