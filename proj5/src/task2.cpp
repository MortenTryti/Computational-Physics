#include <armadillo>
#include <iostream>
#include <complex>
#include <cmath>
#include "../include/header.hpp"






int main() {
  //Random V matrise for å sjekke koden funke
  arma::mat V = arma::mat("1,2,3,4; 6,7,8,9 ; 11,12,13,14; 20,21,22,23");
  int M = 6;
  double h = 1;
  double dt = 1;
  arma::cx_mat U = arma::cx_mat(V,V);

  //Test if the A matrix is implimented correctly by looking at diagonals
  arma::cx_mat A_test = create_nonsparse_Amatrix(V,M,h,dt);
  arma::cx_mat B_test = create_nonsparse_Bmatrix(V,M,h,dt);
  arma::abs(A_test.diag()).print("A-diagonal");
  arma::abs(A_test.diag(1)).print("A-superdiagonal");
  arma::abs(A_test.diag(-1)).print("A-subdiagonal");
  A_test.diag(M-2).print("A-(M-2)-diagonal");
  A_test.diag(-(M-2)).print("A-(2-M)-diagonal");
  
  return 0;
}
