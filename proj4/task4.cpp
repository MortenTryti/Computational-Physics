#include <armadillo>
#include <random>
#include <algorithm>
#include "ising.hpp"
#include <iostream>
int L  = 2;


int seed = 42;
int N = L*L;

int main(int argc, char* argv[]){
    double T = std::stod(argv[1]);
    int n = std::stoi(argv[2]);
    std::cout <<"The temperature is "<< T<<std::endl;
//This is not that elegant, but oh well
Ising ising1 = Ising(L,T,n,seed, true);

ising1.MonteCarlo(n*0.1);

//Extracting the vectors
arma::vec E = ising1.E_;
arma::vec M = ising1.M_;

//Finding mean energies etc

double meanEnergy = arma::sum(E)/n;
double meanEnergy_sq = arma::sum(arma::square(E))/n;
double meanEnergy_perSpin = meanEnergy/(N);
double meanEnergy_square_perSpin = meanEnergy_sq/(N*N);
double meanmag = arma::sum(arma::abs(M))/(n);
double meanmag_perSpin = meanmag/N;
double meanmag_square = arma::sum(arma::square(M))/(n);
double meanmag_square_perspin = meanmag_square/(N*N);

// Getting the numerical quantities, names should be explanitory



double CV = 1./N * 1./(T*T) *(meanEnergy_sq - meanEnergy*meanEnergy);
/*double cv = ising1.specific_heat_capacity();
std::cout<<cv<<std::endl;*/
double xai = 1./N * 1./T * (meanmag_square - meanmag*meanmag);

std::cout<< "<e> = "<<meanEnergy_perSpin<<std::endl;
std::cout<< "<e_T> = "<<theoretical_mean_energy(T,N)<<std::endl;
std::cout<< "<e^2> = "<<meanEnergy_square_perSpin<<std::endl;
std::cout<< "<e_T^2> = "<<theoretical_mean_energy_square(T,N)<<std::endl;
std::cout<< "<|m|> = "<<meanmag_perSpin <<std::endl;
std::cout<< "<|m_T|> = "<<theoretical_mean_magnetisation(T,N) <<std::endl;
std::cout<< "<|m^2|> = "<<meanmag_square_perspin <<std::endl;
std::cout<< "<|m_T^2|> = "<<theoretical_mean_magnetisation_square(T,N) <<std::endl;
std::cout<< "C_V = "<<CV<<std::endl;
std::cout<< "C_V_T = "<<theoretical_heat_cap(T,N) << std::endl;
std::cout<<"X = "<<xai<<std::endl;
std::cout<<"X_T = "<<theoretical_suseptibility(T,N)<<std::endl;

//Setting tolerances and random seeds
std::vector<double> tolerances = {1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8};
std::vector<double> seedvec = {1,5,9076572,777,42};
//Saving them here
arma::mat MC_cycles = arma::mat(5,tolerances.size());

//Averaging over the number of MC cycles needed to get tol
for (int i =0;i<tolerances.size();i++){
    double tol = tolerances.at(i);
    MC_cycles(0,i) = tol;
    //std::cout<<tol<<std::endl;
    double v1=0;
    double v2=0;
    double v3=0;
    double v4=0;
    for (int k=0;k<seedvec.size();k++){
        Ising ising5 = Ising(L,T,n,seedvec.at(k));
        Ising ising6=  Ising(L,T,n,seedvec.at(k));
        Ising ising7=  Ising(L,T,n,seedvec.at(k));
        Ising ising8=  Ising(L,T,n,seedvec.at(k));

        v1 += ising5.mean_energy_convergence(tol);


        v2 += ising6.mean_mag_convergence(tol);

        v3 += ising7.heatCap_convergence(tol);
        v4 += ising8.magSus_convergence(tol);
    }
    double lenlen = seedvec.size();
    //Saving the average nr of MC cycles needed
    MC_cycles(1,i) = v1/lenlen;


    MC_cycles(2,i) = v2/lenlen;

    MC_cycles(3,i) = v3/lenlen;
    MC_cycles(4,i) = v4/lenlen;
}
//Used to print to check if ok
//MC_cycles.print();


MC_cycles.save("MC_Cycles.bin");
return 0;
}
