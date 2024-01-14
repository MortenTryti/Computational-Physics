#include <iostream>
#include <fstream>
#include <armadillo>
#include "classes.hpp"
#include <math.h> 
#include <iomanip>

double b0 = 9.648e1;
double v0 = 2.41e6;
double d = 500;
double f = 0.1;
double omega = 2;
  

int q1 = 1;
double m1 = 40.078;

int q2 = 1;
double m2 = 40.078;

arma::vec r1 = {20,0,20};
arma::vec v1 = {0,25,0};

arma::vec r2 = {25, 25,0};
arma::vec v2 = {0,40,5};

double t = 50;
double ft = 0.0005;


int N_particles = 2;

int main(int argc, char* argv[])
{
  double steps = std::stoi(argv[2]);
  double dt = 50./steps;
  bool interaction = std::stoi(argv[1]);
  bool timedep = false;
  std::cout << "Interaction: " << interaction << std::endl;
  PenningTrap trap = PenningTrap(b0,v0,d,f,omega);
  Particle p1 = Particle(q1,m1,r1,v1);
  Particle p2 = Particle(q2,m2,r2,v2);  
  trap.add_particle(p1);
  trap.add_particle(p2);
 


  arma::mat TXYZcord = arma::zeros(steps, 3*N_particles+1);
  TXYZcord(0,0) = 0;  
  for (int i=0; i<3;i++){
    TXYZcord(0,i+4) = r2(i);
    TXYZcord(0,i+1) = r1(i);
  }
  arma::mat XYZ_Vel = arma::zeros(steps, 3*N_particles);
  for (int i=0; i<3;i++){
    XYZ_Vel(0,i+3) = v2(i);
    XYZ_Vel(0,i) = v1(i);
  }
  
  
  
  
  
  

  for (int i = 1; i < steps; i++)
  {
    double t = 0;
    trap.evolve_forward_Euler(dt , timedep, interaction);    
    
    // Adding the i'th velocities and stuff to the matrix
    
    for (int j=0;j<3;j++){
        XYZ_Vel(i,j) = trap.particles[0].v_[j];
        XYZ_Vel(i,j+3)= trap.particles[1].v_[j];
        TXYZcord(i,j+1) = trap.particles[0].r_[j];
        TXYZcord(i,j+4)= trap.particles[1].r_[j];
    }


    
    TXYZcord(i,0) = i*dt;
    
  }
  bool ffs = true;
  std::string counter = "";
  counter += argv[2];

  if (interaction == true){
  TXYZcord.save("TXYZcordsEulertrue"+counter+".bin");    
  XYZ_Vel.save("XYZ_VelEulertrue"+counter+".bin"); 
} else if (interaction == false){
  TXYZcord.save("TXYZcordsEulerfalse"+counter+".bin");    
  XYZ_Vel.save("XYZ_VelEulerfalse"+counter+".bin"); 
}
  int number = trap.number_of_particles();
  std::cout << "number of particles " <<  number << std::endl;
  return 0;
}
