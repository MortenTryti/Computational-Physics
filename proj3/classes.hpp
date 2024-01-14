#include <armadillo>
#ifndef __classes_hpp__
#define __classes_hpp__

class Particle {
  
  
  

  public:
  // Defining the constructor with the relevant arguments
    int q_;
    double m_;
    arma::vec r_;
    arma::vec v_;
    Particle(int q, double m, arma::vec r, arma::vec v);

    //New function to spit out charge
    int charge();

    //New function to spit out mass
    double mass();

    //New function to spit out position
    arma::vec pos();

    //New function to spit out velocity
    arma::vec vel();

};


class PenningTrap {

  public:
    // variables
    double B0;
    double V0;
    double d;
    double f;
    double omega;    
    std::vector<Particle> particles;

    // Constructor
    PenningTrap(double B0, double V0, double d, double f, double omega);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r, double t, bool timedep=false);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i, double t, bool timedep);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i, double t, bool timedep, bool interaction);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt, double t, bool timedep, bool interaction);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt, bool timedep, bool interaction);

    // Counting numbers of particles inside trap
    int number_of_particles();

		// Return the kinetic energy of the trap

};


#endif
