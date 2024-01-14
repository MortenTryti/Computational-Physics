[1mdiff --cc proj3/classes.hpp[m
[1mindex 3e770e9,7884dee..0000000[m
[1m--- a/proj3/classes.hpp[m
[1m+++ b/proj3/classes.hpp[m
[36m@@@ -61,13 -63,15 +63,25 @@@[m [mclass PenningTrap [m
      arma::vec total_force_particles(int i);[m
  [m
      // The total force on particle_i from both external fields and other particles[m
[32m++<<<<<<< HEAD[m
[32m +  arma::vec total_force(int i, bool interaction);[m
[32m +[m
[32m +    // Evolve the system one time step (dt) using Runge-Kutta 4th order[m
[32m +    void evolve_RK4(double dt, bool interaction);[m
[32m +[m
[32m +    // Evolve the system one time step (dt) using Forward Euler[m
[32m +  void evolve_forward_Euler(double dt, bool interaction);[m
[32m++=======[m
[32m+     arma::vec total_force(int i, double t, bool timedep=false, bool interaction=false);[m
[32m+ [m
[32m+     // Evolve the system one time step (dt) using Runge-Kutta 4th order[m
[32m+     void evolve_RK4(double dt, bool timedep=false, bool interaction=false);[m
[32m+ [m
[32m+     // Evolve the system one time step (dt) using Forward Euler[m
[32m+     void evolve_forward_Euler(double dt, bool timedep=false, bool interaction=false);[m
[32m+ [m
[32m+ 		// Return the kinetic energy of the trap[m
[32m++>>>>>>> 0f48c966a5d8250c447e7c2a0049cba620738fb3[m
  [m
  };[m
  [m
[1mdiff --cc proj3/makefile[m
[1mindex 7bf5584,4700cb4..0000000[m
[1m--- a/proj3/makefile[m
[1m+++ b/proj3/makefile[m
[36m@@@ -3,7 -3,5 +3,12 @@@[m [mcompileEuler[m
  [m
  compileRK4:[m
  	g++ evolve_forward_RK4.cpp penningtrap.cpp particle.cpp -o RK4 -larmadillo -std=c++11[m
[32m++<<<<<<< HEAD[m
[32m +[m
[32m +compileTwoP:[m
[32m +	g++ two_particles.cpp penningtrap.cpp particle.cpp -o Task8 -larmadillo -std=c++11[m
[32m +[m
[32m++=======[m
[32m+ compile_two_particles:[m
[32m+ 	g++ two_particles.cpp penningtrap.cpp particle.cpp -o two_particles -larmadillo -std=c++11[m
[32m++>>>>>>> 0f48c966a5d8250c447e7c2a0049cba620738fb3[m
[1mdiff --cc proj3/penningtrap.cpp[m
[1mindex ea16c04,7decb16..0000000[m
[1m--- a/proj3/penningtrap.cpp[m
[1m+++ b/proj3/penningtrap.cpp[m
[36m@@@ -183,6 -204,10 +204,13 @@@[m [mvoid PenningTrap::evolve_forward_Euler([m
       for (int i=0; i < particles.size();i++){[m
           Particle& pi = particles[i];[m
           Particle& pi_old = particles_old[i];[m
[32m++<<<<<<< HEAD[m
[32m++=======[m
[32m+         [m
[32m+         // (1./6*(k1r[i]+2*k2r[i]+2*k3r[i]+k4r[i])).print("The Ks");[m
[32m+         // pi_old.pos().print("old Pos");[m
[32m+         // (pi_old.pos()+1./6 *(k1r[i]+2*k2r[i]+2*k3r[i]+k4r[i])).print("Whole shabang");[m
[32m++>>>>>>> 0f48c966a5d8250c447e7c2a0049cba620738fb3[m
          pi.r_ = pi_old.pos()+1./6 *(k1r[i]+2*k2r[i]+2*k3r[i]+k4r[i]);[m
          pi.v_ = pi_old.vel()+1./6 *(k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i]);[m
          [m
[1mdiff --cc proj3/two_particles.cpp[m
[1mindex 0efe390,5fc454f..0000000[m
[1m--- a/proj3/two_particles.cpp[m
[1m+++ b/proj3/two_particles.cpp[m
[36m@@@ -22,7 -25,7 +25,11 @@@[m [marma::vec r2 = {25, 25,0}[m
  arma::vec v2 = {0,40,5};[m
  [m
  double t = 50;[m
[32m++<<<<<<< HEAD[m
[32m +double steps = 100000;[m
[32m++=======[m
[32m+ double steps = 10000;[m
[32m++>>>>>>> 0f48c966a5d8250c447e7c2a0049cba620738fb3[m
  [m
  double dt = 50./steps;[m
  [m
[36m@@@ -33,15 -36,10 +40,22 @@@[m [mint N_particles = 2[m
  [m
  int main(int argc, char* argv[])[m
  {[m
[32m++<<<<<<< HEAD[m
[32m +  std::string check = argv[1];[m
[32m +  bool inter;[m
[32m +  if (check == "true"){[m
[32m +     inter = true; [m
[32m +}  else if (check == "false"){[m
[32m +  inter = false;[m
[32m +}[m
[32m +  [m
[32m +  PenningTrap trap = PenningTrap(b0,v0,d);[m
[32m++=======[m
[32m+   bool interaction = std::stoi(argv[1]);[m
[32m+   bool timedep = false;[m
[32m+   std::cout << interaction << std::endl;[m
[32m+   PenningTrap trap = PenningTrap(b0,v0,d,f,omega);[m
[32m++>>>>>>> 0f48c966a5d8250c447e7c2a0049cba620738fb3[m
    Particle p1 = Particle(q1,m1,r1,v1);[m
    Particle p2 = Particle(q2,m2,r2,v2);  [m
    trap.add_particle(p1);[m
[36m@@@ -69,7 -67,7 +83,11 @@@[m
  [m
    for (int i = 0; i < steps; i++)[m
    {[m
[32m++<<<<<<< HEAD[m
[32m +    trap.evolve_RK4(dt, inter);    [m
[32m++=======[m
[32m+     trap.evolve_RK4(dt, timedep, interaction);    [m
[32m++>>>>>>> 0f48c966a5d8250c447e7c2a0049cba620738fb3[m
      [m
      // Adding the i'th velocities and stuff to the matrix[m
      [m
[36m@@@ -85,14 -83,9 +103,14 @@@[m
      TXYZcord(i,0) = i*dt;[m
      [m
    }[m
[31m -  [m
[31m -  TXYZcord.save("TXYZcordsRK4.bin");    [m
[31m -  XYZ_Vel.save("XYZ_VelRK4.bin"); [m
[32m +std::cout<<inter<<std::endl;[m
[32m +  if (check == "true"){[m
[32m +  TXYZcord.save("TXYZcordsRK4True.bin");    [m
[32m +  XYZ_Vel.save("XYZ_VelRK4True.bin"); [m
[32m +} else if (check == "false"){[m
[32m +  TXYZcord.save("TXYZcordsRK4false.bin");    [m
[32m +  XYZ_Vel.save("XYZ_VelRK4false.bin"); [m
[32m +}[m
  [m
    return 0;[m
[31m- }[m
[32m+ }[m
