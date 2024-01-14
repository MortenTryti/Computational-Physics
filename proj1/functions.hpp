#include <armadillo>



// declaring g
double gfunc(double x, double h);

// declaring the forward function
int forwardback(arma::vec a, arma::vec b, arma::vec c, arma::vec g,arma::vec btilde,arma::vec gtilde, arma::vec &v);

// declaring the specialiced forward function
int forwardback_special(arma::vec g, arma::vec btilde, arma::vec gtilde, arma::vec &v);




