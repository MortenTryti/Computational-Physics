#include <iostream>
#include <armadillo>
#include "functions.hpp"
#include <string>


arma::mat A = arma::mat(" 1.0 0.0 0.0 0.5; 0.0 1.0 -0.7 0.0; 0.0 -0.7 1.0 0.0; 0.5 0.0 0.0 1.0");

int i;
int j;

int main()
{   
    std::cout << "The largest modulus element was " << max_offdiag_symmetric(A,i,j)<<std::endl;
    
    std::cout<< "\n The largest element was found at \n -------------------------- \nColumn: " << i<< " |  Row: " << j << std::endl;

    
    A.print("The matrix is");

    return 0;
}




