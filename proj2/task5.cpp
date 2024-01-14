#include <iostream>
#include <math.h> 
#include <iomanip>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <armadillo>
#include "functions.hpp"
#include <iostream>
#include <string>
#include <fstream>

//Tol
double eps = 1e-70;

arma::vec vector = arma::vec("10 20 30 40 50 60 70 80 90 100");

// vector with the number of transformations for matrix A
arma::vec T_A = arma::vec(vector.size());
// vector with the size for matrix A
arma::vec N_A = arma::vec(vector.size());
  



int main(){

  for (int i = 0; i < vector.size(); i++){
    
    // Datapoints N
    int N = vector[i];
    
    //Steplength
    double h = 1.0/(N+1);

    //Super and subdiagonal
    double a = -1.0/(h*h);

    // Main diagonal
    double b = 2.0/(h*h);

    // Making the matrix
    arma::mat A = tridiagmat(N,N,a,b,a);
    
    //Setting up R
    arma::mat R = tridiagmat(N,N,0.0,1.0,0.0);

    int k;
    int l;

    // number of transformations
    int q = 0; 
    while (max_offdiag_symmetric(A,k,l)>eps){
        double Loffdiag = max_offdiag_symmetric(A,k,l);
        Jacobi_rotate(A,R,k,l);
	q++;
    }
    // vector for number of transforamtions and size of matrix
    T_A(i) = q;
    N_A(i) = N;
    std::cout << "Number of transformations for N = " << N << ": " << q << std::endl;
  }



  // make the file
  std::string filename = "task5Output.txt";
  
  std::ofstream ofile;
  ofile.open(filename);

  //for loop for the vector points
  for (int i = 0; i < vector.size(); i++)
    {
       
      ofile <<"N" + std::to_string(i+1) + "  "<< std::setprecision(12) << std::scientific << N_A[i]<< "  v(x"+std::to_string(i+1)+")  " << T_A[i] << std::endl;
    }
  ofile.close();  
  
  return 0;
}
