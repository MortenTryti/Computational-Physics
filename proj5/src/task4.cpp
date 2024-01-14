#include <armadillo>
#include <iostream>
#include <complex>
#include <cmath>
#include "../include/header.hpp"


double h = 1e-2;

int main(){

    double stuff = 0.05;
    double momenta = 1;
    double center = 0.25;

    arma::cx_mat u0 = init_state(h,momenta,momenta,center,center,stuff,stuff);
    arma::mat u0_prob_sq = arma::real(u0)%arma::real(u0) + arma::imag(u0)%arma::imag(u0);
    //This is done to see if the wavefunc has been properly normalised
    double tot_prob = sqrt(arma::accu(u0_prob_sq));
    std::cout<<"The value of |psi|="<< std::sqrt(arma::cdot(u0,u0)) <<std::endl;
    //Save it so we can plot it with plot4.py
    u0_prob_sq.save("data/u0_prob_sq.bin");
    return 0;
}
