#include <iostream>
#include <math.h> 
#include <iomanip>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <armadillo>
#include "functions.hpp"

//Tol
double eps = 1e-70;






int main(int argc, char* argv[]){

    // Number of steps
    int n = std::atoi(argv[1]);

    // Datapoints N
    int N = n - 1;
    
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

    // sort the eigenvalues to find the lowest eigenvalue
    arma::vec diagA_sorted = arma::sort(diagA);

    // Find the three lowest eigenvalues
    double first_eigval = diagA_sorted(0);
    double second_eigval = diagA_sorted(1);
    double third_eigval = diagA_sorted(2);

    std::cout << "The three lowest eigenvalues from the Jacobi rotation algotithm are: " << std::endl;
    std::cout << first_eigval << std::endl;
    std::cout << second_eigval << std::endl;
    std::cout << third_eigval << std::endl;

    // Find the index of the lowest eigenvalues in the original vector of eigenvalues
    arma::uvec index_first_eigval = find(diagA == first_eigval);
    arma::uvec index_second_eigval = find(diagA == second_eigval);
    arma::uvec index_third_eigval = find(diagA == third_eigval);

    // The three eigenvectors from the Jacobi rotation algorithm corresponding to the three lowest eigenvalues    
    arma::vec j1 = R.col(index_first_eigval[0]);
    arma::vec j2 = R.col(index_second_eigval[0]);
    arma::vec j3 = R.col(index_third_eigval[0]);

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

    // sort the eigenvalues to find the lowest eigenvalue    
    arma::vec analytic_eigval_sorted = arma::sort(analytic_eigval);

    // Find the three lowest analytical eigenvalues
    double first_eigval_analytic = analytic_eigval_sorted(0);
    double second_eigval_analytic = analytic_eigval_sorted(1);
    double third_eigval_analytic = analytic_eigval_sorted(2);

    
    // Find the three lowest eigenvalues
    std::cout << "The three lowest analytical eigenvalues are: " << std::endl;
    std::cout << first_eigval_analytic << std::endl;
    std::cout << second_eigval_analytic << std::endl;
    std::cout << third_eigval_analytic << std::endl;

    // Find the index of the lowest eigenvalues in the original vector of eigenvalues
    arma::uvec index_first_eigval_analytic = find(analytic_eigval == first_eigval_analytic);
    arma::uvec index_second_eigval_analytic = find(analytic_eigval == second_eigval_analytic);
    arma::uvec index_third_eigval_analytic = find(analytic_eigval == third_eigval_analytic);

    // The three analytical eigenvectors corresponding to the three lowest eigenvalues
    arma::vec a1 = analytic_eigvec.col(index_first_eigval_analytic[0]);
    arma::vec a2 = analytic_eigvec.col(index_second_eigval_analytic[0]);
    arma::vec a3 = analytic_eigvec.col(index_third_eigval_analytic[0]);


    // Writing the eigenvectors to file

    //Discrete rep of x without BPs
    arma::vec x = arma::linspace(h,1-h,N);
    
    std::string filename = "task6Output"+std::to_string(n)+"_Jacobi.txt";

    std::ofstream ofile;
    ofile.open(filename);
    //Adds the first BP
    ofile <<"x" + std::to_string(0) + "  "<< std::setprecision(12) << std::scientific << 0 << "  j1(x"+std::to_string(0)+")  " << 0 << "  a2(x"+std::to_string(0)+")  " << 0 << "  a3(x"+std::to_string(0)+")  " << 0 << std::endl;
    //for loop for the vector points
    for (int i =0; i < N; i++)
      {
	ofile <<"x" + std::to_string(i+1) + "  "<< std::setprecision(12) << std::scientific << x[i]<< "  a1(x"+std::to_string(i+1)+")  " << a1[i] << "  a2(x"+std::to_string(i+1)+")  " << a2[i] << "  a3(x"+std::to_string(i+1)+")  " << a3[i] << std::endl;
      }
    // Adds the last BP
    ofile <<"x" + std::to_string(N+1) + "  "<< std::setprecision(12) << std::scientific << 1 << "  j1(x"+std::to_string(N+1)+")  " << 0 << "  j2(x"+std::to_string(N+1)+")  " << 0 << "  j3(x"+std::to_string(N+1)+")  " << 0 << std::endl;    
    ofile.close();

    std::string filename_analytic = "task6Output"+std::to_string(n)+"_analytic.txt";

    std::ofstream ofile_analytic;
    ofile_analytic.open(filename_analytic);
    //Adds the first BP
    ofile_analytic <<"x" + std::to_string(0) + "  "<< std::setprecision(12) << std::scientific << 0<< "  a1(x"+std::to_string(0)+")  " << 0 << "  a2(x"+std::to_string(0)+")  " << 0 << "  a3(x"+std::to_string(0)+")  " << 0 << std::endl;
    //for loop for the vector points
    for (int i = 0; i < N; i++)
      {
	ofile_analytic <<"x" + std::to_string(i+1) + "  "<< std::setprecision(12) << std::scientific << x[i]<< "  a1(x"+std::to_string(i+1)+")  " << a1[i] << "  a2(x"+std::to_string(i+1)+")  " << a2[i] << "  a3(x"+std::to_string(i+1)+")  " << a3[i] << std::endl;
      }
    // Adds the last BP
    ofile_analytic <<"x" + std::to_string(N+1) + "  "<< std::setprecision(12) << std::scientific << 1 << "  a1(x"+std::to_string(N+1)+")  " << 0 << "  a2(x"+std::to_string(N+1)+")  " << 0 << "  a3(x"+std::to_string(N+1)+")  " << 0 << std::endl;        
    ofile_analytic.close();
    
    return 0;
}
