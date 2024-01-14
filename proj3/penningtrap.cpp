#include <iostream>
#include <armadillo>
#include "classes.hpp"
#include <omp.h>

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f_in, double omega_in)
{
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
  f = f_in;
  omega = omega_in;

  std::vector<Particle> particles;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}


// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r, double t, bool timedep)
{
    arma::vec pos_dep = arma::vec(3);
    pos_dep(0) = r(0);
    pos_dep(1) = r(1);
    pos_dep(2) = -2*r(2);
    if (arma::norm(r) > d){
      arma::vec Efield = arma::vec(3).fill(0.);
      return Efield;    
    } else{
      if (timedep){
	arma::vec Efield = V0*(1+f*std::cos(omega*t))/(d*d)* pos_dep;
	return Efield;
      }
      else{
	arma::vec Efield = V0/(d*d)* pos_dep;
	return Efield;}
    }
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
  //Magnetic field is constant in z-direction
  if (arma::norm(r) > d){
    arma::vec Bfield = arma::vec(3).fill(0.);
    return Bfield;
  } else{
    arma::vec Bfield = B0*arma::vec("0.0 0.0 1.0");
    return Bfield;
  }
}


// The total force on particle_i from the external fields
 arma::vec PenningTrap::total_force_external(int i, double t, bool timedep)
{
    Particle p_i = particles[i];
    arma::vec F_ext = (p_i.charge())*(external_E_field(p_i.pos(), t, timedep) + arma::cross(p_i.vel(),external_B_field(p_i.pos())) );
    return F_ext;
}

arma::vec PenningTrap::total_force_particles(int i)
{   double ke = 1.389e5;
//std::cout<<"T"<<std::endl;
    int qi = particles[i].charge();
    arma::vec ri = particles[i].pos();

    arma::vec F_p = arma::vec("0.0 0.0 0.0");

    if (particles.size()<2)
    {
     
      return F_p;
    }
    else
    {
      for (int j = 0; j < particles.size(); j++)
      {
        
        if (i!=j)
        {
          int qj = particles[j].charge();
          arma::vec rj = particles[j].pos();
	  
          F_p = F_p + (ke*qi*qj *  (ri-rj))/(std::pow(arma::norm(ri-rj),3));
        }
      }
      
      return F_p;
    }
    
    
}

arma::vec PenningTrap::total_force(int i, double t, bool timedep, bool interaction)
{
  //std::cout<<interaction<<std::endl;
  if (interaction){
    
    return total_force_external(i, t, timedep) + total_force_particles(i);
  }
  else{
    
    return total_force_external(i, t, timedep);
  }
}


// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, bool timedep, bool interaction)
{
  //Making copy of particle instances
  std::vector<Particle> particles_copy = particles;
  // For loop to run over all particles

  for (int i = 0; i < particles.size(); i++)
  {
    //Initialize where we store the different accelerations on each particle
    //Making new variables for readability and smaller notation
    Particle pc_i = particles_copy[i];
    Particle p_i = particles[i];

    //Acceleration from total force
    arma::vec a_tot = total_force(i, i*dt, timedep, interaction)/p_i.mass();
    //define new velocity and position
    arma::vec v_ip1 = arma::vec("0.0 0.0 0.0");
    arma::vec r_ip1 = arma::vec("0.0 0.0 0.0");
    
    //Does FE-algo in 3D
    
    v_ip1 = pc_i.vel() + dt*a_tot;
    r_ip1 = pc_i.pos() + dt*pc_i.vel();
    
    
    //updates particle with new particle instance
    particles[i] = Particle(p_i.charge(),p_i.mass(),r_ip1,v_ip1); 
  }


}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, double t, bool timedep, bool interaction){
     //Making copy of particle instances
     std::vector<Particle> particles_old = particles;
    
     std::vector <arma::vec> k1v(particles.size());
     std::vector <arma::vec> k1r(particles.size());
     std::vector <arma::vec> k2v(particles.size());
     std::vector <arma::vec> k2r(particles.size());
     std::vector <arma::vec> k3v(particles.size());
     std::vector <arma::vec> k3r(particles.size());
     std::vector <arma::vec> k4v(particles.size());
     std::vector <arma::vec> k4r(particles.size());
     
     #pragma omp parallel for
     //Find k1
     for (int j=0; j < particles.size();j++){
         Particle& pj = particles[j];
         k1v[j] = dt*total_force(j, t, timedep, interaction)/pj.mass();
         k1r[j] = dt*pj.vel();
     }
     #pragma omp parallel for
     //Update position so we can find k2
     for (int i=0; i < particles.size();i++){
         Particle& pi = particles[i];
         Particle& pi_old = particles_old[i];
         
         pi.r_ = pi_old.pos()+k1r[i]/2;
         pi.v_ = pi_old.vel()+k1v[i]/2;
         
         
     }
     #pragma omp parallel for
     //Find k2
     for (int j=0; j < particles.size();j++){
         Particle& pj = particles[j];
         k2v[j] = dt*total_force(j, t, timedep, interaction)/pj.mass();
         k2r[j] = dt*pj.vel();
     }
     #pragma omp parallel for
     //Update with k2 so we can find k3
     for (int i=0; i < particles.size();i++){
         Particle& pi = particles[i];
         Particle& pi_old = particles_old[i];
         
         pi.r_ = pi_old.pos()+k2r[i]/2;
         pi.v_ = pi_old.vel()+k2v[i]/2;
         
         
     }
     #pragma omp parallel for     
     //Find k3
     for (int j=0; j < particles.size();j++){
         Particle& pj = particles[j];
         k3v[j] = dt*total_force(j, t, timedep, interaction)/pj.mass();
         k3r[j] = dt*pj.vel();
     }
     #pragma omp parallel for     
     //Update with k3 so we can find k4
     for (int i=0; i < particles.size();i++){
         Particle& pi = particles[i];
         Particle& pi_old = particles_old[i];
         
         pi.r_ = pi_old.pos()+k3r[i];
         pi.v_ = pi_old.vel()+k3v[i];
         
     }
     #pragma omp parallel for     
     //Find k4
     for (int j=0; j < particles.size();j++){
         Particle& pj = particles[j];
         k4v[j] = dt*total_force(j, t, timedep, interaction)/pj.mass();
         k4r[j] = dt*pj.vel();
     }
     #pragma omp parallel for     
     //Finding the new position with the 4 k's. This is the final for loop
     for (int i=0; i < particles.size();i++){
         Particle& pi = particles[i];
         Particle& pi_old = particles_old[i];
        
        pi.r_ = pi_old.pos()+1./6 *(k1r[i]+2*k2r[i]+2*k3r[i]+k4r[i]);
        pi.v_ = pi_old.vel()+1./6 *(k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i]);
        
        
        
     }
 }

int PenningTrap::number_of_particles(){
  int count = 0;
  for (int i = 0; i < particles.size(); i++){
    double r = arma::norm(particles[i].pos());
    // std::cout << r << std::endl;
    if (r < d){
      count += 1;
    }
  }
  return count;
}


// }


// /*
// //Force on particle_i from particle_j
// arma::vec PenningTrap::force_particle(int i, int j){
//     //Find relevant constantvalue
//     double k;
//     //Define the relevant charges
//     int qi = particles[i].charge();
//     int qj = particles[j].charge();
//     //Position
//     arma::vec ri = particles[i].pos();
//     arma::vec rj = particles[j].pos();
//     //Defining the force on particle i from j in 3-space directions
//     arma::vec Fij = -k*qi*qj *  (ri-rj)/(arma::abs(ri-rj)*arma::abs(ri-rj)*arma::abs(ri-rj));
//     return Fij;
// }*/
