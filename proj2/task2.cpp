#include <armadillo>
#include <iostream>
#include "functions.hpp"
#include <cmath>


int main()
{
	//Dimension
	int N = 6;
	double d = 2.0;
	double a= 1.0;

	//Makes a tridiagonal matrix
	arma::mat A = tridiagmat(N,N,a,d,a);

	//Initiates eigenvector and eigenvalues vectors and matrices
	arma::vec eigval;
	arma::mat eigvec;
	
	//
	arma::eig_sym(eigval, eigvec, A);

	//Prints the eigenvalues
	eigval.print("Eigenvalues");
	arma::normalise(eigvec,1).print("The eigenvectors");

	//Need to compare to the analytic eigenvectors
	arma::vec analytic_eigval= arma::vec(N);
	arma::mat analytic_eigvec = arma::mat(N,N); 



	for (int i=1; i<=N; i++){
	analytic_eigval(i-1) = d + 2*a*std::cos(i*M_PI/(N+1));
	
		for (int j=1;j<=N;j++){
			analytic_eigvec(j-1,i-1) = std::sin(i*j*M_PI/(N+1));

		}

	}
	analytic_eigvec = arma::normalise(analytic_eigvec,1);	
	arma::sort( analytic_eigval).print("Analytical eigenvalues");
	analytic_eigvec.print("Analytical eigenvectors");
	return 0;
}