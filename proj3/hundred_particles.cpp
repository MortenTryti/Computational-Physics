#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include "classes.hpp"
#include <math.h> 
#include <iomanip>

double b0 = 9.648e1;
double v0 = 2.41e6;
double d = 500;
double f = 0.7;
  
int q = 1;
double m = 40.078;


double t = 500;
double steps = 40000;

double dt = 500./steps;

int N_particles = 100;

int main(int argc, char* argv[])
{
  bool interaction = std::stoi(argv[1]);
  double f = std::stod(argv[2]);
  bool timedep = true;
  std::cout << interaction << std::endl;

  double stepsize = 0.005;
  double omega = 0.02;
  double omega_end = 2.5;
  int N_omegas = (omega_end-omega)/stepsize;
  std::cout << "Looping over " << N_omegas << " values for omega" << std::endl;
  arma::mat PT = arma::zeros(N_omegas, 2);
  arma::arma_rng::set_seed_random();

  for (int j = 0; j < N_omegas; j++) {
    PenningTrap trap = PenningTrap(b0,v0,d,f,omega);
    for (int i = 0; i < N_particles; i++) {
      arma::vec r = arma::vec(3).randn() * 0.1 *d;
      arma::vec v = arma::vec(3).randn() * 0.1 *d;    
      Particle p = Particle(q,m,r,v);
      trap.add_particle(p);
    }
    
    for (int i = 0; i < steps; i++)
      {
	double t = i*dt;
	trap.evolve_RK4(dt, t, timedep, interaction);    
      }    
    int number = trap.number_of_particles();
    PT(j,0) = omega;
    PT(j,1) = number;
    omega += stepsize;
    if (j % 10 == 0) {
      std::cout << "Loop " << j << " of " << N_omegas << std::endl;
    }  
  }


  // Create an output string for f
  std::ostringstream fstring;
  fstring << std::fixed;
  fstring << std::setprecision(2);
  fstring << f;
  std::string foutput = fstring.str();
  PT.save("hundred_particles_"+foutput+".bin");
  // PT.save("hundred_particles_"+std::to_string(f)+".bin");
  
  return 0;
}

