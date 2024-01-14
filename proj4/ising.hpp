#include <math.h>
#include <stdlib.h>
#include <armadillo>
#include <random>
#ifndef __ising_hpp__
#define __ising_hpp__


// Should write the header for the theoretical functions here

//Theoretical mean energy
double theoretical_mean_energy(double T,int N);


//Theoretical mean energy square
double theoretical_mean_energy_square(double T,int N);

//Theoretical mean magnetisation
double theoretical_mean_magnetisation(double T,int N);

//Theoretical mean magnetisation square
double theoretical_mean_magnetisation_square(double T,int N);

//Theoretical mean heat cap
double theoretical_heat_cap(double T,int N);

//Theoretical suseptibility
double theoretical_suseptibility(double T,int N);



class Ising {
  public:
  //Define system variables here
  int L_;
  double T_;
  int n_;
  int seed_;
  std::vector<double> prob_frac;
  std::vector<int> energy_diff;
  // Spin matrix
  arma::mat A_;

  // Current energy and magnetization of the system
  double E_A;
  double M_A;

  // Vector of energy and magnetization values for each MC cycle
  arma::vec E_;
  arma::vec M_;

  // Generator
  std::mt19937 generator_;

  // Constructor, if the seed is not specified a random seed is used
  Ising(int L, double T, int n, int seed = -1, bool samespin = false);

  // Counts the number of parallel neighbouring spins to a given state in the spin matrix
  int count_parallel_spin(int row, int col);

  // Attepmts to flip one spin in the system via the Metropolis algorithm
  // Adds the energy and magnetization difference to E_A and M_A if a spin is flipped
  void one_attempted_spin_flip();

  // Does one cycle in the Markov chain Monte Carlo methon with L^2 appmpted spin flips
  void one_cycle();

  //Does n cycles
  void n_cycles();

  // The main function for initialize the Monte Carlo method over n cycles in the Ising model
  // Modifies the vectors E_ and M_ containing the energies and magnetizations of each cycle,
  // which can then be used further in the main program
  void MonteCarlo(int n_burn);

  // Calculate the energy of the spin configuration
  double energy();

  // Calculate the magnetization of the spin configuration
  double magnetization();

  // Calculate the specific hear capacity of the spin configuration
  double specific_heat_capacity();

  // Calculate the suceptibility of the spin configuration
  double susceptibility();


  ////////////////////////////////////////////////////////////////////////////////
  // old functions from the first implementation (slow) can remove later, maybe useful to have still
  ////////////////////////////////////////////////////////////////////////////////

  // Calculete the mean energy per particle
  double mean_energy();

  //Mean energy calculation w/ each iteration
  //stored in a vector. Goes until the diff between
  //the numerical and analytical is smaller than tol
  int mean_energy_convergence(double tol);

  // Calculete the mean energy squared per particle
  double mean_energy_square();

  // Calculate the absolute magnetization of the spin configuration
  double abs_magnetization();

  // Calculate the absolute magnetization of the spin configuration
  double mean_mag();


  //Mean magnetisation calculation w/ each iteration
  //stored in a vector. Goes until the diff between
  //the numerical and analytical is smaller than tol
  int mean_mag_convergence(double tol);

  int heatCap_convergence(double tol);

  int magSus_convergence(double tol);

  double mean_mag_square();






};


#endif
