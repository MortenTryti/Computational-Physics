#include <armadillo>
#include <omp.h>
#include <random>
#include <algorithm>
#include <chrono>
#include "ising.hpp"


int n = 1000000;
int n_burn = n*0.1;
int seed = 1;


int main(int argc, char* argv[]){

  
  // Start measuring time
  auto t1 = std::chrono::high_resolution_clock::now();

  // determine if we zoom in on the temperature range 2.25-2.35 or 2.1-2.4
  int zoom = std::stod(argv[1]);
  // number of temperatures
  int T_points = std::stod(argv[2]);
  // lattice size
  int L = std::stod(argv[3]);

  double T_begin;
  double T_end;
  
  // if we do not zoom in
  if (zoom == 0){
    T_begin = 2.1;
    T_end = 2.4;
  }

  // if we zoom in
  else if (zoom == 1){
    T_begin = 2.2;
    T_end = 2.35;
  }

  std::cout << T_begin << " " << T_end << std::endl;

  arma::vec T_vec = arma::linspace(T_begin, T_end, T_points);
  arma::mat Tvec = arma::zeros(T_vec.size(), 1);
  
  for (int i = 0; i < T_vec.size(); i++){
    Tvec(i) = T_vec(i,0);
  }
  
  arma::mat variables = arma::zeros( T_points, 4);
  
  std::cout << "L = " << L << std::endl;    
  // parallelisation
  #pragma omp parallel for
  for (int j=0; j<T_vec.size(); j++){
    double T = T_vec[j];
    std::cout << "T = " << T << std::endl;
    // initialise the system
    Ising ising = Ising(L,T,n, seed + j,false); 
    // run the Monte Carlo simulation
    ising.MonteCarlo(n_burn);
    int N = L*L;
    // calculate the variables
    double E_mean = sum(ising.E_)/n;
    double E_square_mean = sum(square(ising.E_))/n;      
    double M_abs_mean = sum(abs(ising.M_))/n;
    double M_square_mean = sum(square(ising.M_))/n;
    double C_V = 1/(N*T*T)*(E_square_mean - pow(E_mean,2));
    double chi = 1/(N*T)*(M_square_mean - pow(M_abs_mean,2));
    std::cout << "E_mean = " << E_mean << std::endl;
    // save the variables in a matrix
    variables(j,0) = E_mean/N;
    variables(j,1) = M_abs_mean/N;
    variables(j,2) = C_V;
    variables(j,3) = chi;
      
  }
    
  // save the variables in a file 
  if (zoom == 0){
    Tvec.save("T_unzoom_"+std::to_string(T_points)+"_vec.bin");
    variables.save("variables_unzoom_"+std::to_string(T_points)+"_"+std::to_string(L)+".bin");
  }

  else if (zoom == 1){
    Tvec.save("T_zoom_"+std::to_string(T_points)+"_vec.bin");
    variables.save("variables_zoom_"+std::to_string(T_points)+"_"+std::to_string(L)+".bin");
  }
  
  // Stop measuring time
  auto t2 = std::chrono::high_resolution_clock::now();

  // Calculate the elapsed time
  // We use chrono::duration<double>::count(), which by default returns duration in seconds
  double duration_seconds = std::chrono::duration<double>(t2 - t1).count();

  // Print the elapsed time
  std::cout << "Elapsed time: " << duration_seconds << " s" << std::endl;
  
  return 0;
}
