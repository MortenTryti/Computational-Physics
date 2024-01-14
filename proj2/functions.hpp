#include <armadillo>


// Create a tridiagonal matrix 
arma::mat tridiagmat(int n, int m, double sub, double main, double super);

double max_offdiag_symmetric(arma::mat A, int &i, int &j);

void Jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l);
