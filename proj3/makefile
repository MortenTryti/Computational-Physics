compileEuler:
	g++ evolve_forward_Euler.cpp penningtrap.cpp particle.cpp -o Euler -larmadillo -std=c++11
compileRK4:
	g++ evolve_forward_RK4.cpp penningtrap.cpp particle.cpp -o RK4 -larmadillo -std=c++11
compile_two_particlesRK4:
	g++ two_particles.cpp penningtrap.cpp particle.cpp -o two_particles -larmadillo -std=c++11
compile_two_particlesEuler:
	g++ two_particlesEuler.cpp penningtrap.cpp particle.cpp -o two_particlesE -larmadillo -std=c++11
compile_hundred_particles:
	g++ -fopenmp hundred_particles.cpp penningtrap.cpp particle.cpp -o hundred_particles -larmadillo -std=c++11
compile_hundred_particles_O2:
	g++ -O2 -fopenmp hundred_particles.cpp penningtrap.cpp particle.cpp -o hundred_particles -larmadillo -std=c++11
compile_hundred_particles_O3:
	g++ -O3 -fopenmp hundred_particles.cpp penningtrap.cpp particle.cpp -o hundred_particles -larmadillo -std=c++11
compile_hundred_particles_zoom_O3:
	g++ -O3 -fopenmp hundred_particles_zoom.cpp penningtrap.cpp particle.cpp -o hundred_particles_zoom -larmadillo -std=c++11



