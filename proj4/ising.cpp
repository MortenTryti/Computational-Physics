#include "ising.hpp"
#include <limits>
#include <algorithm>
#include <armadillo>
#include <random>
#include <ctime>
#include <math.h>
#include <cmath>
//Defining the theoretical variables here

double theoretical_mean_energy(double T,int N){
  // Returns the energy in units of [J]
  // T should be normalised as T/(k_B*J)
  double Ntihi = N;
 double eps =  -8./Ntihi * sinh (8./T)/(3.+ cosh (8./T));
 return eps;
}

double theoretical_mean_energy_square(double T, int N){
   // T should be normalised as T/(k_B*J)
  double Z = 4.*cosh (8./T) + 12.;
  double eps2 = (8.*8.)/(N*N)*(1.-3./(3.+cosh (8./T)));
  //double eps2 = (8*8)/(N*N*Z) *(2*exp(8/T)+exp(-8/T));
  return eps2;
}

double theoretical_mean_magnetisation(double T,int N){
   // T should be normalised as T/(k_B*J)
  double Z = 4.*cosh (8./T) + 12.;
  double Ntihi = N;
  double mean_m_abs = (8.*(exp(8./T) + 2.))/(Ntihi*Z);
  return mean_m_abs;

}

double theoretical_mean_magnetisation_square(double T,int N){
   // T should be normalised as T/(k_B*J)
  double Z = 4.*cosh (8./T) + 12.;
  double mean_m_abs2 = 32.*(exp(8./T) + 1.)/(N*N*Z);
  return mean_m_abs2;
}

double theoretical_heat_cap(double T, int N)
{
  double Z = 4.*cosh (8./T) + 12.;
  double cv = (8.*8.)/(N*T*T) * (1.+3*cosh(8./T))/( (3.+cosh (8./T))*(3.+cosh (8./T)) );
  return cv;
}

double theoretical_suseptibility(double T, int N){
  //T should be normalised to T/(k_B*J)
  double Z = 4.*cosh (8./T) + 12.;
  double Ntihi = N;
  //Defining <M^2>
  double a = 32.*(exp(8./T) + 1.)/Z;
  //Defining <M>
  double b =(8.*(exp(8./T) + 2.))/Z;
  //Making xai
  double xai = 1./(N*T) * (a-b*b) ;
  return xai;
}


// Define the constructor
Ising::Ising(int L, double T, int n, int seed, bool samespin){

    if (seed < 0){
      std::random_device random_seed;
      seed_ = random_seed();
      //std::cout<<"seed:"<<seed_<<std::endl;
    } else{
    seed_ = seed;
    //std::cout<<"seed:"<<seed_<<std::endl;
    }
    generator_.seed (seed_);
    std::uniform_int_distribution<int> distribution(1,2);
  

    // Length of lattice
    L_ = L;
    // Temperature of system
    T_ = T;
    // Number of cycles in the Markov chain Monte Carlo Method
    n_ = n;

    // Possible spin flip probs.
    prob_frac = {exp(8./T_),exp(4./T_),1,exp(-4./T_),exp(-8./T_)};
    // Possible energy differences of a single spin flip, index is number of parallel spins
    energy_diff = {-8, -4, 0, 4, 8};

    // Vector of energy and magnetization values for each MC cycle
    M_ = arma::vec(n_).fill(0);
    E_ = arma::vec(n_).fill(0);

    // Spin matrix for the lattice
    A_ = arma::mat(L_,L_).fill(1);
    int a = 1;
    int b = -1;
    if (samespin ==false){


    // Filling the spin matric randomlly
    for (int i=0; i<L;i++){
        for (int j=0;j<L;j++){
            int coinflip = distribution(generator_);
            if (coinflip ==1){
                A_(i,j) = b ;
            } else{
                A_(i,j) = a;
            }

        }
    }}
    M_A = magnetization();
}

// Counts the number of parallel neighbouring spins to a given state in the spin matrix
int Ising::count_parallel_spin(int row, int col){
    int spin_state = A_(row, col);
    // The spin of the neighbours with periodic boundary conditions (using modulo)
    // Need to add L to not get -1 index when col = 0 or row = 0 for the modulo operator
    // neighbours_spin = {up, down, left, right}
    std::vector<double> neighbours_spin = {A_((row+L_-1)%L_,col),A_((row+L_+1)%L_,col),A_(row,(col+L_-1)%L_),A_(row,(col+1)%L_)};
    // Count the number of states with parallel spins to the random state
    int parallel_spin = 0;
    for (int i=0;i < neighbours_spin.size();i++){
      if (neighbours_spin[i] == spin_state){
	parallel_spin += 1;
      }
    }
    return parallel_spin;
}



void Ising::one_attempted_spin_flip(){
    std::uniform_int_distribution<int> distribution(0,L_-1);
    // Generate random row and col for the spin to (possible) be flipped
    int row = distribution(generator_);
    int col = distribution(generator_);

    // The spin of the random state
    int random_state = A_(row,col);

    // Number of neighbouring states with parallel spin to the chosen state
    int parallel_spin = count_parallel_spin(row,col);

    // Acceptance probability
    double acceptance_prob = std::min(1.0,prob_frac[parallel_spin]);

    std::uniform_real_distribution<double> distribution2(0.0,1.0);
    double random_number = distribution2(generator_);

    // Flip the spin if the random number is less than the acceptance probability
    if (random_number <= acceptance_prob){
      // Flip the spin
      A_(row,col) = -random_state;
  
      // Calculate the energy of the new state by adding the energy difference from the previous state,
      // using the vector containting the 5 possible energy differences
      E_A += energy_diff[parallel_spin];
      
      // Calculate the magnetization of the new state by adding the magnetization difference from the previous state,
      // For each spin flip the magnetization difference is 2 times the original spin
      M_A += -2*random_state;
     
    }

}


// One cycle of spin flips consisting of L^2 attempted spin flips
void Ising::one_cycle(){
  // Run over one cycle corresponding to L^2 attempted spin flips
  for (int i=0; i < L_*L_; i++){
    one_attempted_spin_flip();
  }
}

// n cycles of spin flips consisting of L^2 attempted spin flips
void Ising::n_cycles(){
  // Run over one cycle corresponding to L^2 attempted spin flips
  for (int i=0; i < n_; i++){
    one_cycle();
  }
}

// The main function for initialize the Markov chain Monte Carlo method over n cycles in the Ising model
// Returns the vectors E_ and M_ containing the energies and magnetizations of each cycle,
// which can then be used further in the main program
void Ising::MonteCarlo(int n_burn){
  // Calculete the energy and magnetization of the initial state
  E_A = energy();
  M_A = magnetization();
  std::cout << "Initial energy " << E_A << std::endl;
  std::cout << "Initial mag " << M_A << std::endl;
  // Run over n_burn cycles before saving savint the energies and magnetizations
  if (n_burn >= 0){
    for (int i=0; i < n_burn; i++){
      one_cycle();
    }
  }
  for (int i=0; i<n_; i++){
    // Run one cycle and save the energy of the cycle in the energy and magnetication vector
    one_cycle();
    E_(i) = E_A;
    M_(i) = magnetization();
  }
}

// Finding the total energy of the system
double Ising::energy(){
  double E = 0;
  for (int i=0; i<L_;i++){
    for (int j=0;j<L_;j++){
      E += -A_(i,j)*(A_((i+1)%L_,j) + A_(i,(j+1)%L_));
    }
  }
  return E;
}

double Ising::magnetization(){
  double M = 0;
  for (int i=0; i<L_;i++){
    for (int j=0;j<L_;j++){
      M += A_(i,j);
    }
  }
  //  std::cout<<M<<std::endl;
  return M;
}

////////////////////////////////////////////////////////////////////////////////
// old functions from the first implementation (slow) can remove later, maybe useful to have still
////////////////////////////////////////////////////////////////////////////////

double Ising::abs_magnetization(){
  double M = 0;
  for (int i=0; i<L_;i++){
    for (int j=0;j<L_;j++){
      M += A_(i,j);
    }
  }
  return abs(M);
}

double Ising::mean_energy(){
  double total_energy_sample = 0;
  for (int i=0; i<n_; i++){
    one_cycle();
    total_energy_sample += energy();
  }
  double mean_total_energy_spin = total_energy_sample/(L_*L_*n_);
  return mean_total_energy_spin;
}

double Ising::mean_energy_square(){
  double total_energy_sample = 0;
  for (int i=0; i<n_; i++){
    one_cycle();
    total_energy_sample += energy()*energy();
  }
  double mean_total_energy_spin = total_energy_sample/(L_*L_*L_*L_*n_);
  return mean_total_energy_spin;

}



double Ising::mean_mag(){
  double total_mag_sample = 0;
  for (int i=0; i<n_; i++){
    one_cycle();
    total_mag_sample += abs_magnetization();
  }
  double mean_total_mag_spin = total_mag_sample/(L_*L_*n_);
  return mean_total_mag_spin;
}


double Ising::mean_mag_square(){
  double total_mag_sample = 0;
  for (int i=0; i<n_; i++){
    one_cycle();
    total_mag_sample += magnetization()*magnetization();
  }
  double denom = L_*L_*L_*L_*n_;
  double mean_total_mag_spin = total_mag_sample/denom;
  return mean_total_mag_spin;
}




double Ising::specific_heat_capacity(){
  //The calculated heat capasoty is calculated as C_V/k_B, so plot accordingly
int N = L_*L_;
double C_V = 0.;
double tot_eps = 0.;
double tot_eps2= 0.;
  for (int i=0; i<n_; i++){
    one_cycle();
    tot_eps += energy();
    tot_eps2 += energy()*energy();
  }
double norm_eps = tot_eps/(n_);
double norm_eps2 = tot_eps2/(n_);


C_V += 1./N *1/(T_*T_)  *(norm_eps2-norm_eps*norm_eps);
return C_V;
}


int Ising::mean_energy_convergence(double tol){
  // To count number of MC cycles used
  int numMC=1;
double Theo_value =theoretical_mean_energy(T_,L_*L_);
double total_energy_sample = energy()/(L_*L_);
  while (abs(Theo_value - total_energy_sample)>tol){
    //Does one cycle
    one_cycle();
    // Adds one to counter
    numMC +=1;
    //Updates value
    total_energy_sample = total_energy_sample*(numMC-1)/(numMC) +energy()/(L_*L_*numMC);
  }
  //std::cout<<"The theoretical value for the mean energy is "<<Theo_value<<std::endl;
  //std::cout<<"The last numerical value for avg energy is "<<total_energy_sample<<std::endl;
  //std::cout<< "#MC cycles = "<< numMC<< "| tol = "<< tol<<std::endl;

  return numMC;
}

int Ising::mean_mag_convergence(double tol){
  // To count number of MC cycles used
  int numMC=1;
double Theo_value =theoretical_mean_magnetisation(T_,L_*L_);
double total_mag_sample = abs_magnetization()/(L_*L_);
  while (abs(Theo_value - total_mag_sample)>tol){
    //Does one cycle
    one_cycle();
    // Adds one to counter
    numMC +=1;
    //Updates value
    total_mag_sample = total_mag_sample*(numMC-1)/(numMC) +abs_magnetization()/(L_*L_*numMC);
  }
  //std::cout<<"The theoretical value for the mean magnetisation is "<<Theo_value<<std::endl;
  //std::cout<<"The last numerical value for mean magnetisation is "<<total_mag_sample<<std::endl;
  //std::cout<< "#MC cycles = "<< numMC<< "| tol = "<< tol<<std::endl;

  return numMC;
}

int Ising::heatCap_convergence(double tol){
  // To count number of MC cycles used
  int numMC=1;
double Theo_value =theoretical_heat_cap(T_,L_*L_);
int N = L_*L_;
double C_V = 0.;
double tot_eps = 0.;
double tot_eps2= 0.;
    tot_eps += energy();
    tot_eps2 += energy()*energy();
double norm_eps = tot_eps;
double norm_eps2 = tot_eps2;
//std::cout<<"The theoretical value for C_V is "<<Theo_value<<std::endl;
C_V = 1./N *1/(T_*T_)  *(norm_eps2-norm_eps*norm_eps);
  while (abs(Theo_value - C_V)>tol){
    //Does one cycle
    one_cycle();
    // Adds one to counter
    numMC +=1;
    //Updates value
    tot_eps = tot_eps*(numMC-1)/(numMC) +energy()/(numMC);
    tot_eps2= tot_eps2*(numMC-1)/(numMC) +energy()*energy()/(numMC);
    C_V = 1./N *1./(T_*T_)  *(tot_eps2-tot_eps*tot_eps);


  }
  
  //std::cout<<"The last numerical value for C_V is "<<C_V<<std::endl;
  //std::cout<< "#MC cycles = "<< numMC<< "| tol = "<< tol<<std::endl;

  return numMC;
}

int Ising::magSus_convergence(double tol){
  // To count number of MC cycles used
  int numMC=1;
double Theo_value =theoretical_suseptibility(T_,L_*L_);
int N = L_*L_;
double xai = 0.;
double M = 0.;
double M2= 0.;
    M += abs_magnetization();
    M2 += magnetization()*magnetization();
//std::cout<<"The theoretical value for X is "<<Theo_value<<std::endl;
xai = 1./N *1./(T_)  *(M2-M*M);

  while (abs(Theo_value - xai)>tol){
    //Does one cycle
    one_cycle();
    // Adds one to counter
    numMC +=1;
    //Updates value
    M = M*(numMC-1)/(numMC) +abs_magnetization()/(numMC);
    M2= M2*(numMC-1)/(numMC) +magnetization()*magnetization()/(numMC);
    //New sus
    xai = 1./(N*T_)  *(M2-M*M);
  }
  
  //std::cout<<"The last numerical value for X is "<<xai<<std::endl;
  //std::cout<< "#MC cycles = "<< numMC<< "| tol = "<< tol<<std::endl;

  return numMC;
}
