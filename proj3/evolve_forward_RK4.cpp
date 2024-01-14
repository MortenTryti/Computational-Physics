#include <iostream>
#include <fstream>
#include <armadillo>
#include "classes.hpp"
#include <math.h> 
#include <iomanip>

double b0 = 9.648e1;
double v0 = 2.41e6;
double d = 500;

int q1 = 1;
double m1 = 40.078;

int q2 = 1;
double m2 = 40.078;

arma::vec r1 = {20,0,20};
arma::vec v1 = {0,25,0};

arma::vec r2 = {25, 25,0};
arma::vec v2 = {0,40,5};

double t = 50;
double steps = 1000;

double dt = 50./steps;

double ft = 0.0005;



int main()
{

  PenningTrap trap = PenningTrap(b0,v0,d);
  Particle p1 = Particle(q1,m1,r1,v1);
  Particle p2 = Particle(q2,m2,r2,v2);  
  trap.add_particle(p1);
 


  arma::mat Zcord1 = arma::zeros(steps,2);
  Zcord1(0,0) = r1(2);
  Zcord1(0,1) = 0;  
  

  trap.evolve_RK4(dt,false);
  
  // make the file
  std::string filename = "evolve_RK4.txt";
  
  std::ofstream ofile;
  ofile.open(filename);

  for (int i = 0; i < steps; i++)
  {
    trap.evolve_RK4(dt, false);    
    
    Zcord1(i,0) = trap.particles[0].pos()(2);
    Zcord1(i,1) = i*dt;

    ofile <<"N" + std::to_string(i) + "  "<< std::setprecision(12) << std::scientific << Zcord1(i,0) << std::endl;
    
  }
  
  Zcord1.save("ZcordRK4.bin");    
  ofile.close();  

  return 0;
}