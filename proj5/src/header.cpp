#include <armadillo>
#include <iostream>
#include <complex>
#include <math.h>
#include <cmath>

//This is where we actually write the code


//For finding total prob
double tot_prob(arma::cx_mat U){
    double tot_p = arma::accu(arma::real(U)%arma::real(U) + arma::imag(U)%arma::imag(U));
    return tot_p;

}


//For index mapping
int k_index(int M,int i, int j){
    //Used to pick out the (i,j) index of a matrix and convering it to an index k

    //Todo: Check if i and j are the correct indices.....
  return i + (M-2)*(j) ;
}


//Transforms U from matrix to vec u
arma::cx_vec from_mat_to_vec(arma::cx_mat U){
    

    int Mmin2 = U.n_cols;
    arma::cx_vec u = arma::cx_vec(Mmin2*Mmin2);
    for (int i=0;i<Mmin2;i++){
        for (int j=0;j<Mmin2;j++){
            int k = j+ i*Mmin2;
            u.at(k) = U.col(i)(j);
        }
    }
return u;
}

//Transforms vector u to matrix U
arma::cx_mat from_vec_to_mat(arma::cx_vec u){
    int N = u.n_elem/std::sqrt(u.n_elem);
    arma::cx_mat U = arma::cx_mat(N,N);
    for (int k = 0; k<u.n_elem;k++){
        int i = k%N;
        int j = floor(k/N);
        U(i,j) = u.at(k);
    }


    return U;
}


//Making the sparse A matrix
arma::sp_cx_mat create_Amatrix(arma::mat V,int M, double h, double dt){
    // set dimension of matrix and define constant r
    int dim = (M-2)*(M-2);
     arma::cx_double r = arma::cx_double(0,(dt)/(2*h*h));
    //arma::cx_double r(0,1);
    arma::cx_double coeff2 = arma::cx_double(0,dt/2);
    arma::cx_double coeff1 = 4.0*r;

    // create empty matrices
    arma::sp_cx_mat A=arma::sp_cx_mat(dim,dim);


    arma::cx_vec a = arma::cx_vec(dim);
    //Mapping the valus v_ij to a_k etc
    for (int j = 0; j<M-2;j++){
        for (int i = 0; i<M-2;i++){
            //finding k
            int k = k_index(M,i,j);
            //Using k to create a
            a.at(k) = 1.0+coeff1+coeff2*V(i,j);
        }
    }
    //Filling up the diagonal elements
    A.diag(0)=a;

    // off diagonal values
    A.diag((M-2)) -= r;
    A.diag(-((M-2))) -= r;

    //Filling up the periodic offdiagonal elements
    for (int i=1;i<(M-2)*(M-2)+1;i++){
    if (i%(M-2)!=0){
        A.diag(1).at(i-1) = -r;
        A.diag(-1).at(i-1) = -r;
    }
    }
    return A;
}


//Making the non-sparse A matrix, for better visual
arma::cx_mat create_nonsparse_Amatrix(arma::mat V,int M, double h, double dt){
    // set dimension of matrix and define constant r
    int dim = (M-2)*(M-2);
    arma::cx_double r = arma::cx_double(0,(dt)/(2*h*h));
    //arma::cx_double r(0,1);
    arma::cx_double coeff2 = arma::cx_double(0,dt/2);
    arma::cx_double coeff1 = 4.0*r;

    // create empty matrices
    arma::cx_mat A=arma::cx_mat(dim,dim);


    arma::cx_vec a = arma::cx_vec(dim);
    //Mapping the valus v_ij to a_k etc
    for (int j = 0; j<M-2;j++){
        for (int i = 0; i<M-2;i++){
            //finding k
            int k = k_index(M,i,j);
            //Using k to create a
            a.at(k) = 1.0+coeff1+coeff2*V(i,j);
        }
    }
    //Filling up the diagonal elements
    A.diag(0)=a;

    // off diagonal values
    A.diag((M-2)) -= r;
    A.diag(-((M-2))) -= r;

    //Filling up the periodic offdiagonal elements
    for (int i=1;i<(M-2)*(M-2)+1;i++){
    if (i%(M-2)!=0){
        A.diag(1).at(i-1) = -r;
        A.diag(-1).at(i-1) = -r;
    }
    }
    return A;
}






//Making the B matrix
arma::sp_cx_mat create_Bmatrix(arma::mat V,int M, double h, double dt){
    // set dimension of matrix and define constant r
    int dim = (M-2)*(M-2);
    //Defining the variables in making the diagonal
    arma::cx_double r = arma::cx_double(0,(dt)/(2*h*h));
    //arma::cx_double r(0,1);
    arma::cx_double coeff2 = arma::cx_double(0,dt/2);
    arma::cx_double coeff1 = 4.0*r;

    // create empty matrices
    arma::sp_cx_mat B = arma::sp_cx_mat(dim,dim);

    // off diagonal values

    arma::cx_vec b = arma::cx_vec(dim);

    //Mapping the valus v_ij to a_k etc
    for (int j = 0; j<M-2;j++){
        for (int i = 0; i<M-2;i++){
            //finding k
            int k = k_index(M,i,j);
            //Using k to create a
            b.at(k) = 1.0-coeff1 -coeff2 *V(i,j);
        }
    }


    B.diag(0)=b;

    B.diag((M-2)) += r;
    B.diag(-((M-2))) += r;


    for (int i=1;i<(M-2)*(M-2)+1;i++){
    if (i%(M-2)!=0){
        B.diag(1).at(i-1) = r;
        B.diag(-1).at(i-1) = r;
    }
    }



    return B;
}


//Making the non-sparse B matrix, for visualisation
arma::cx_mat create_nonsparse_Bmatrix(arma::mat V,int M, double h, double dt){
    // set dimension of matrix and define constant r
    int dim = (M-2)*(M-2);
     arma::cx_double r = arma::cx_double(0,(dt)/(2*h*h));
    //arma::cx_double r(0,1);
    arma::cx_double coeff2 = arma::cx_double(0,dt/2);
    arma::cx_double coeff1 = 4.0*r;

    // create empty matrices
    arma::cx_mat B = arma::cx_mat(dim,dim);

    // off diagonal values

    arma::cx_vec b = arma::cx_vec(dim);

    //Mapping the valus v_ij to a_k etc
    for (int j = 0; j<M-2;j++){
        for (int i = 0; i<M-2;i++){
            //finding k
            int k = k_index(M,i,j);
            //Using k to create a
            b.at(k) = 1.0-coeff1 -coeff2 *V(i,j);
        }
    }


    B.diag(0)=b;

    B.diag((M-2)) += r;
    B.diag(-((M-2))) += r;


    for (int i=1;i<(M-2)*(M-2)+1;i++){
    if (i%(M-2)!=0){
        B.diag(1).at(i-1) = r;
        B.diag(-1).at(i-1) = r;
    }
    }



    return B;
}



//initial state construction

arma::cx_mat init_state(double h, double px, double py , double xc, double yc, double sig_x, double sig_y){
    int N = (1.0/h)+1;
    std::cout<<"Number of points in U"<<N<<std::endl;
    double x0=0;
    double y0=0;
    
    arma::cx_mat u0 = arma::cx_mat(N,N);
    
    
    for (int i= 0;i<N;i++){
        for (int j= 0;j<N;j++){
            double gauss_exponent = -pow(x0+i*h - xc,2)/pow(2.0*sig_x,2)-pow(y0+j*h - yc,2)/pow(2.0*sig_y,2) ;
            double envelope = exp(gauss_exponent);
            double real = cos(px*i*h+py*j*h);
            double imag = sin(px*i*h+py*j*h);
            arma::cx_double r(real,imag);
            u0.at(i,j) = envelope * r;

        }
    }
    
    

    u0.row(0).zeros();
    u0.row(N-1).zeros();
    u0.col(0).zeros();
    u0.col(N-1).zeros();

    arma::mat abs_u0 =arma::real(u0)%arma::real(u0) + arma::imag(u0)%arma::imag(u0);
    double normconst = std::sqrt(arma::accu(abs_u0)) ;
    u0/=normconst;
    

    return u0;

}

//A quick solver to find u^(n+1)
arma::cx_colvec solver_evolver(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec uold){
    arma::cx_colvec b = B*uold;
    arma::cx_colvec u2 = arma::spsolve(A,b);

    return u2;
}


// The potential set up for the 2 opening slit
arma::mat pot_V(int N,double v0,int nh){

    // create potential matrix to fill
    arma::mat V = arma::mat(N,N,arma::fill::zeros);


    double h = 1/N;
    //Used to find the indices
    arma::vec xpos = arma::linspace(0,1,N);
    arma::vec ypos = arma::linspace(0,1,N);

    if (nh==2){
    //Setting up the first barrier, from y=0.0 to y=0.5-0.05-0.025
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
        //Checking if we are in the range we want to be
            if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)<=0.5-0.05-0.025){
            
            V(i,j)=v0;
            }
    }

    //Setting up the second barrir, from y = 0.5-0.025 to 0.5+0.025
    for (int j=0;j<N;j++){
        //Setting up the place between the splits
            if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)>=0.5-0.025 and ypos(i)<=0.5+0.025 ){
            V(i,j)=v0;
            }
    }

    //Setting up the third section of the barrier, from y=0.5+0.05+0.025 ot y=1.0
    for (int j=0;j<N;j++){
        //Setting up the last barrier, checking if we are in between the regions where it should be
            if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)>=0.5+0.05+0.025){
            
            V(i,j)=v0;
            }
    }


    } } else if (nh==1){
        //Setting up the first barrier, from y=0.0 to y=0.5-0.025
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
        //Checking if we are in the range we want to be
            if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)<=0.5-0.025){
            
            V(i,j)=v0;
            } else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)>=0.5+0.025){
                V(i,j) = v0;
            }
    }}} else if (nh==3){
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
        //Checking if we are in the range we want to be for the first wall
            //if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)<=0.5-0.025-0.05-0.05){
            if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and i<=75){
            
            V(i,j)=v0;
            } else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and i>=85 and i<=95 ){
                //else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)>=0.5-0.025-0.05 and ypos(i)<=0.5-0.025 )
                // In here we make the second wall segment
                V(i,j)=v0;
            } else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and i>=105 and i<=115 ){
                //else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)>=0.5+0.025 and ypos(i)<=0.5+0.025+0.05 ){
                // In here we make the third wall segment
                V(i,j)=v0;
            } else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and i>=125){
                //else if (0.5-0.01<xpos(j) and xpos(j)<0.5+0.01 and ypos(i)>=0.5+0.025+0.05+0.05){
                //Last wall segment
                V(i,j)=v0;
            }
    } }
    }


    // Boundary of the potential
    V.row(0).fill(v0);
    V.row(N-1).fill(v0);
    V.col(0).fill(v0);
    V.col(N-1).fill(v0);
    return V;
}
