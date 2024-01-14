#include <iostream>
#include <math.h> 
#include <iomanip>
#include <armadillo>
#include <string>
#include <fstream>
#include "functions.hpp"

// Declares N
int N;



int main(int argc, char* argv[])
{

  // Declaring the size of the matrix and problem
  int N=atoi(argv[1]);

  //steplength
  double h = 1.0/(N+1);

  //Here we should define the super,sub and diagonal vectors of the matrix

  arma::vec btilde = arma::vec(N);
  arma::vec gtilde = arma::vec(N);
  arma::vec v = arma::vec(N);
  arma::vec g = arma::vec(N);


  
  //Discrete rep of x without BPs
  arma::vec x = arma::linspace(h,1-h,N);
  

  // Filling g with the correct values
  for (int i = 0; i< N; i++){
    g[i] = gfunc(x[i], h);
  }
  
  forwardback_special(g, btilde, gtilde, v);




  // make the file
  std::string filename = "task9Output"+std::to_string(N)+".txt";
  
  std::ofstream ofile;
  ofile.open(filename);
  //Adds the first BP
  ofile <<"x" + std::to_string(0) + "  "<< std::setprecision(12) << std::scientific << 0<< "  v(x"+std::to_string(0)+")  " << 0 << std::endl;
  //for loop for the vector points
  for (int i =0; i < N; i++)
    {
       
      ofile <<"x" + std::to_string(i+1) + "  "<< std::setprecision(12) << std::scientific << x[i]<< "  v(x"+std::to_string(i+1)+")  " << v[i] << std::endl;
    }
  // Adds the last BP
  ofile <<"x" + std::to_string(N+2) + "  "<< std::setprecision(12) << std::scientific << 1<< "  v(x"+std::to_string(N+2)+")  " << 0 << std::endl;
  ofile.close();

  return 0;
}
    
