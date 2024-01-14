#include <armadillo>
#include <iostream>
#include <complex>
#include <cmath>
#include "../include/header.hpp"


// Command line input should be in following order
// double h - spatial steplenght
// double dt - temporal steplenght
// double T - total runtime
// double x_c - center of initial wavepacket in x-dir
// double sig_x - root of gaussian variance in x-dir
// double p_x - momenta in x-dir
// double y_c - center of initial wavepacket in y-dir
// double sig_y - root of gaussian variance in y-dir
// double p_x - momenta in y-dir

int main(int argc, char* argv[]){
    

    //These are pretty self explanatory
    double h =  std::atof(argv[1]);
    double dt =  std::atof(argv[2]);
    double T =  std::atof(argv[3]);
    double x_c =  std::atof(argv[4]); 
    double sig_x =  std::atof(argv[5]);
    double p_x =  std::atof(argv[6]);
    double y_c =  std::atof(argv[7]); 
    double sig_y =  std::atof(argv[8]);
    double p_y =  std::atof(argv[9]);
    double v0 =  std::atof(argv[10]);
    //Used to save to different files, so no over writing
    //Declaring the number of datapoints in the spatial directions
    int M =1.0/h+1;
    std::cout<<"The value of M "<< M<<std::endl;
    //Declaring the number of datapoints in the temporal directions
    int N =T/dt+1;
    


    //Making the indices we are to save snapshots of 
    int i_t1 = 0.001/dt;
    int i_t2 = 0.002/dt;

    //Making the potentiale
    arma::mat V = pot_V(M,v0);
    //Making the initial state
    arma::cx_mat U0 = init_state(h,p_x,p_y,x_c,y_c,sig_x,sig_y);
    
    arma::cx_mat U0_internal = U0.submat(1,1,M-2,M-2);
    U0_internal = U0_internal.st();
    double U0_probsq = tot_prob(U0);

    std::cout<<"Initial norm if wavefunc is "<<arma::accu(U0_probsq)<<std::endl;
    std::cout<<"Initial norm if wavefunc is "<<arma::cdot(U0_internal,U0_internal)<<std::endl;
    
    //Making the A and B matrices
    std::cout<<"# Elems in V "<< V.n_elem <<std::endl;
    
    arma::mat V_int = V.submat(1,1,M-2,M-2);
    std::cout<<"# Elems in V_int "<< V_int.n_elem <<std::endl;
    arma::sp_cx_mat A = create_Amatrix(V_int,M,h,dt);
    arma::sp_cx_mat B = create_Bmatrix(V_int,M,h,dt);
    std::cout<<"# cols in A "<< A.n_cols <<std::endl;
    std::cout<<"# cols in B "<< B.n_cols <<std::endl;
    std::cout<<"# cols in U_int "<< U0_internal.n_cols <<std::endl;
    //making the matrix we stor things in
    arma::cx_cube the_imaginary_CUBE = arma::cx_cube(M-2,M-2,3) ;
    arma::cube the_real_CUBE = arma::cube(M-2,M-2,3) ;
    arma::cube the_imag_CUBE = arma::cube(M-2,M-2,3) ;
    arma::cube the_probability_CUBE = arma::cube(M-2,M-2,3) ;

    //The dimension of thee internal points
    int dim = (M-2)*(M-2);
    
    //the_imaginary_CUBE.slice(0) = U0_internal;
    
    arma::cx_colvec u0_vec = arma::cx_colvec(dim);
    arma::cx_colvec uNew_vec = arma::cx_colvec(dim);
    the_imaginary_CUBE.slice(0) = U0_internal;
    arma::vec dviance_probvec = arma::vec(N);
    //Running over timesteps
    for (int i = 0;i<N;i++ ){
        
        dviance_probvec.at(i) = abs(1- tot_prob(U0_internal));
        //Save U as it is now in a MxM slice at time pos i*dt
        //the_CUBE.slice(i) = U0;
        //Transform U to vec u, to use in iteration
        u0_vec = U0_internal.as_col();
        //Use solver to find new u
        uNew_vec = solver_evolver(A,B,u0_vec);
        //std::cout<<"The new norm of the new solution "<< arma::accu(arma::real(uNew_vec)%arma::real(uNew_vec) + arma::imag(uNew_vec)%arma::imag(uNew_vec)) <<std::endl;
        //Transforming the new vector to a mat so it can be savd and used in new iter
        std::cout<<"The deviance of the prob from 1 is "<<dviance_probvec.at(i)<<std::endl;
        
        U0_internal = arma::reshape(uNew_vec,M-2,M-2);
        
        if (i==i_t1){
            the_probability_CUBE.slice(1) =arma::reshape( arma::real(U0)*arma::real(U0) + arma::imag(U0)*arma::imag(U0),M-2,M-2);
            the_imaginary_CUBE.slice(1) = U0_internal;
            the_real_CUBE.slice(1) = arma::reshape(arma::real(U0_internal),M-2,M-2);
            the_imag_CUBE.slice(1) = arma::reshape(arma::imag(U0_internal),M-2,M-2);
        } else if (i==i_t2){
            the_imaginary_CUBE.slice(2) = U0_internal;
            the_probability_CUBE.slice(2) = arma::reshape(arma::real(U0)*arma::real(U0) + arma::imag(U0)*arma::imag(U0),M-2,M-2);
            the_real_CUBE.slice(2) = arma::reshape(arma::real(U0_internal),M-2,M-2);
            the_imag_CUBE.slice(2) = arma::reshape(arma::imag(U0_internal),M-2,M-2);
        }

        //Saving the new vector as the "old" for the next iteration
        u0_vec = uNew_vec;
        std::cout<<"----------ITERATION "<< i<<" ----------\n"<<std::endl;
        

        //repeat!

    }
    the_imaginary_CUBE.save("data/Task8_imag_cube.bin");
    the_probability_CUBE.save("data/Task8_prob_cube.bin");
    the_real_CUBE.save("data/Task8RE_cube.bin");
    the_imag_CUBE.save("data/Task8IM_cube.bin");
  

  return 0;
}
