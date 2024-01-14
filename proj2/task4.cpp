#include <iostream>
#include <armadillo>
#include "functions.hpp"

//Tol
double eps = 1e-70;






int main(int argc, char* argv[]){

    // Datapoints N
    int N = std::atoi(argv[1]);
    
    //Steplength
    double h = 1.0/(N+1);

    //Super and subdiagonal
    double a = -1.0/(h*h);

    // Main diagonal
    double b = 2.0/(h*h);

    //Making the matrix
    arma::mat A = tridiagmat(N,N,a,b,a);


    //Setting up R
    arma::mat R = tridiagmat(N,N,0.0,1.0,0.0);

    int k;
    int l;

    while (max_offdiag_symmetric(A,k,l)>eps){
        double Loffdiag = max_offdiag_symmetric(A,k,l);
        Jacobi_rotate(A,R,k,l);
    
    
    }
    
    //Printing the result of the iterative method
    arma::vec diagA = A.diag(); 
    //arma::normalise(R,1).print("The eigenvectors are");
    A.print("Change of A after iterative method");
    arma::sort(diagA).print("Eigenvalues of the iterative method");

    //Need to compare to the analytic eigenvectors
	arma::vec analytic_eigval= arma::vec(N);
	arma::mat analytic_eigvec = arma::mat(N,N); 


    //This is the analytical expressions
	for (int i=1; i<=N; i++){
	analytic_eigval(i-1) = b + 2*a*std::cos(i*M_PI/(N+1));
	
		for (int j=1;j<=N;j++){
			analytic_eigvec(j-1,i-1) = std::sin(i*j*M_PI/(N+1));

		}

	}
	analytic_eigvec = arma::normalise(analytic_eigvec,1);	
	arma::sort( analytic_eigval).print("Analytical eigenvalues");

    (arma::sort(diagA)-arma::sort(analytic_eigval)).print("Diff between the sorted analytical and iterative eigenvalues");
	
    return 0;
}
