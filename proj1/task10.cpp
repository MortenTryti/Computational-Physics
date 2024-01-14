#include <iostream>
#include <math.h> 
#include <iomanip>
#include <armadillo>
#include <string>
#include <fstream>
#include "functions.hpp"
#include <chrono>


// Declares N
int N;


//Takes the first argument to be number of datapoints
//The second argument is the number of times it should measure time to average over
int main(int argc, char* argv[])
{

    // Declaring the size of the matrix and problem
    int N=atoi(argv[1]);

    //steplength
    double h = 1.0/(N+1);

    //Here we should define the super,sub and diagonal vectors of the matrix
    arma::vec a = arma::vec(N).fill(-1.0);
    arma::vec b = arma::vec(N).fill(2.0);
    arma::vec c = arma::vec(N).fill(-1.0);

    //Other vectors we will use later
    arma::vec at = arma::vec(N);
    arma::vec bt = arma::vec(N);
    arma::vec ct = arma::vec(N);

    arma::vec btilde = arma::vec(N);
    arma::vec gtilde = arma::vec(N);
    arma::vec v = arma::vec(N);
    arma::vec g = arma::vec(N);


    
    //Discrete rep of x without BPs
    arma::vec x = arma::linspace(h,1-h,N);
    
    for (int i = 0; i< N; i++){
        g[i] = gfunc(x[i], h);
    }
    
    //Initialising the variable to take our timevalues
    double duration_seconds_g = 0;
    double duration_seconds_s = 0;
    //Number of iterations
    int I = atoi(argv[2]);
    //For loop to average over
    for (int i = 0; i<I;i++){
        // Start measuring time
        auto t1 = std::chrono::high_resolution_clock::now();

        forwardback( a,  b,  c,  g, btilde, gtilde, v);

        // Stop measuring time
        auto t2 = std::chrono::high_resolution_clock::now();


        // Calculate the elapsed time
        // We use chrono::duration<double>::count(), which by default returns duration in seconds
        duration_seconds_g += std::chrono::duration<double>(t2 - t1).count();
   }
    //Printing the results
  std::cout<< "General: <t>= "<< duration_seconds_g/I <<" |  n= "<< N << " | " << " I= "<<I <<"\n";

    //For loop for averaging over the time that the special algo uses
    for (int i= 0; i<I;i++){
        // Start measuring time
        auto t3 = std::chrono::high_resolution_clock::now();

        forwardback_special(g,btilde,gtilde,v);


        // Stop measuring time
        auto t4 = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time
        // We use chrono::duration<double>::count(), which by default returns duration in seconds
        duration_seconds_s += std::chrono::duration<double>(t4 - t3).count();
        
    }
    std::cout<< "Special: <t>= "<< duration_seconds_s/I <<" |  n= "<< N << " | " << " I= "<<I <<"\n\n";
    std::cout<< "Special time/general time : "<< duration_seconds_s/duration_seconds_g << "\n";


  return 0;
}
    