#include <armadillo>
#include <iostream>
#include <complex>
#include <cmath>
#include "../include/header.hpp"


int M = 5;
double h = 1;
double dt = 1;  

int dim = (M-2)*(M-2);

arma::cx_vec u1 = arma::cx_vec(dim).fill(1.0);


int main() {
    //Random V matrix
    arma::mat V = arma::mat("1,2,3; 6,7,8 ; 11,12,13");
    

    //Making the matrics
    arma::sp_cx_mat A = create_Amatrix(V,M,h,dt);
    arma::sp_cx_mat B = create_Bmatrix(V,M,h,dt);

    
    //This is solving the problem, note we check if the solver evolver does what we want it to/is consistent
    arma::cx_vec b = B*u1;
    arma::cx_vec u2 = arma::spsolve(A,b);
    (u2-solver_evolver(A,B,u1)).print("Consistency check of implimentation of solver_evolver");
    (A*u2-b).print("Consistency check to see if our solution actually solves the equation,\n if yes then this should be filled with elems close to 0");



    return 0;
}
