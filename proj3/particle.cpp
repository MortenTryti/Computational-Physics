#include <iostream>
#include <armadillo>
#include "classes.hpp"

// Actually defining the constructor
Particle::Particle(int q, double m, arma::vec r, arma::vec v)
{
  // Relevant variables
  q_ = q;
  m_ = m;
  r_ = r;
  v_ = v;
}

//Functions lets me get the relevant variables for particle

//New function to spit out charge
int Particle::charge(){
  return q_;
}

//New function to spit out mass
double Particle::mass(){
  return m_;
}

//New function to spit out position
arma::vec Particle::pos(){
  return r_;
}

//New function to spit out velocity
arma::vec Particle::vel(){
 return v_;
}
