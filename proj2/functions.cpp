#include <armadillo>
#include "functions.hpp"


arma::mat tridiagmat(int n, int m, double sub, double main, double super)
{
	// Create Matrix full of zeroes
	arma::mat A = arma::mat(n,m).fill(0.);
	// Set the subdiagonal, main diagonal and superdiagonal
	A.diag(-1).fill(sub);
	A.diag().fill(main);
	A.diag(1).fill(super);

	return A;

}


double max_offdiag_symmetric(arma::mat A, int &i, int &j)
{	//This function is to find the largest off diagonal element.
	// i is the column and j is the row
	A.diag().zeros();
	int N = A.n_cols;
	
	double maxval = arma::abs(A).max();
	int lin_index = arma::abs(A).index_max(); 
	j = int(lin_index) % N;
	
	i = int(lin_index) /N; 
    return maxval;

}


void Jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l){

	int N = A.n_cols;
	double tau,t,c,s;
        tau = (A(l,l)-A(k,k))/(2*A(k,l));


        //Defining the cosine, tangent and the sine
        if (A(k,l)==0) {
            s=0.0;
            c=1.0;
        } else{
            if (tau>0){
                t = 1/(tau + std::sqrt(1+tau*tau));
               
            } else{
                t = -1/(-tau + std::sqrt(1+tau*tau));
                }
            c = 1/std::sqrt(1+t*t);
            s = c*t;
        }
        

        //Updating values for m+1 iteration
        double a_kk = A(k,k);
        double a_ll = A(l,l);
        A(k,k) = a_kk*c*c -2.0*A(k,l)*c*s+a_ll*s*s;
        
        A(l,l) = a_ll *c*c + 2.0*A(k,l)*c*s+a_kk*s*s;
        
        A(l,k) = 0.0;
        
        A(k,l) = 0.0;
        

        //Now we change the remaining elements of A
        for (int i=0;i<N;i++){
            if (i != k && i != l){
                double a_ik = A(i,k);
                A(i,k) = A(i,k)*c - A(i,l)*s;
                A(k,i) = A(i,k);
                A(i,l)= A(i,l)*c + a_ik*s;
                A(l,i) = A(i,l);
            }
        }

    for (int i=0;i<N;i++){
        double r_ik = R(i,k);
        //Updating values in R to m+1 version
        R(i,k) = R(i,k)*c-R(i,l)*s;
        R(i,l) = R(i,l)*c + r_ik*s;
    }
    
    
    }

