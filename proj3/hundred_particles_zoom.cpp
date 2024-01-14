#include <iostream>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <armadillo>
#include "classes.hpp"
#include <math.h> 
#include <iomanip>

double b0 = 9.648e1;
double v0 = 2.41e6;
double d = 500;
double f = 0.7;
double omega = 2e6;
  
int q = 1;
double m = 40.078;


double t = 500;
double steps = 80000;

double dt = t/steps;

int N_particles = 100;

int main(int argc, char* argv[])
{
  bool interaction = std::stoi(argv[1]);
  double f = std::stod(argv[2]);
  bool timedep = true;
  std::cout << interaction << std::endl;
  std::cout << f << std::endl;
  std::cout << timedep << std::endl;
  double stepsize = 0.00125;
  double omega = 1.25;
  double omega_end = 1.5;
  int N_omegas = (omega_end-omega)/stepsize;
  std::cout << "Looping over " << N_omegas << " values for omega" << std::endl;
  arma::mat PT = arma::zeros(N_omegas, 2);
  arma::arma_rng::set_seed_random();

  // #pragma omp parallel for
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
    PT(j,1) = number/N_particles;
    omega += stepsize;

    if (j % 10 == 0) {
      std::cout << "Loop " << j << " of " << N_omegas << std::endl;
    }
  }
  
  // Create an output string fro f
  std::ostringstream fstring;
  fstring << std::fixed;
  fstring << std::setprecision(2);
  fstring << f;
  std::string foutput = fstring.str();
  PT.save("hundred_particles_zoom_"+std::to_string(interaction)+"_"+foutput+".bin");
  // PT.save("hundred_particles_"+std::to_string(f)+".bin");
  
  return 0;
}
