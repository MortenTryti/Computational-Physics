#include <armadillo>
#include <iostream>
#include <complex>
#include <math.h>

//This is the headerfile

//Index mapping
int k_index(int M,int i, int j);

//Finding tot prob of matrix
double tot_prob(arma::cx_mat U);


//A quick solver to find u^(n+1)
arma::cx_colvec solver_evolver(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec uold);


//For creating the A-matrix
arma::sp_cx_mat create_Amatrix(arma::mat V,int M, double h, double dt);

//For creating the B-matrix
arma::sp_cx_mat create_Bmatrix(arma::mat V,int M, double h, double dt);

arma::cx_mat create_nonsparse_Amatrix(arma::mat V,int M, double h, double dt);

arma::cx_mat create_nonsparse_Bmatrix(arma::mat V,int M, double h, double dt);


//initial state construction
arma::cx_mat init_state(double h, double px, double py , double xc, double yc, double sig_x, double sig_y);

// The potential set up for the 2 opening slit
arma::mat pot_V(int N,double v0,int nh = 2);

//Transforms U from matrix to vec u
arma::cx_vec from_mat_to_vec(arma::cx_mat U);

//Transforms vector u to matrix U
arma::cx_mat from_vec_to_mat(arma::cx_vec u);