#include <iostream>
#include <math.h> 
#include <iomanip>
#include <armadillo>
#include <string>
#include <fstream>




// define the u(x) function
double u(double x)
{
    double b = 1-(1-exp(-10))*x -   exp(-10 * x);
    return b;
}


int main(int argc, char* argv[])
{
    int N = atoi(argv[1]);
    
    // define vector
    arma::vec x = arma::linspace(0,1,N+2);


    // make the file
    std::string filename = "task2Output.txt";
    std::ofstream ofile;
    ofile.open(filename);

    //for loop for the vector points
    for (int i =0; i < N+2; i++)
        {ofile <<"x  "<< std::setprecision(4) << std::scientific << x[i]<< "  u(x)  " << u(x[i]) << std::endl;
        }
        ofile.close();

    return 0;
}