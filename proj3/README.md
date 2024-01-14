# README

### How to run the code for task 8

To use the program for making the datafiles for two particles in the Penning trap use the command
```
make compile_two_particlesRK4
```
this creates the machine code file two_particles, to run this file use 

```
./two_particles bool n
```
where bool is either 1 or 0, 1 means that we have interactions, 0 means no interaction. The second input n is the number of steps and should be an integer. The output of this file will be a set of .bin files, one containing the timesteps and positions of the two particles, the other contains the velocities of the two particles. 

To make the simulation files with the Euler method write 

```
make compile_two_particlesEuler
```
this gives out the machine code file ``` two_particlesE ``` which operates exactly like the one for RK4, just with FE as the numerical integration method.

Before running the plotting code we advice you to make all the simulation files beforehand, for all the different timesteps and interactions. To make the plot of the positions and phase space of the two particles simply run 
```
python3 plot_two_particles.py
```
this will generate many plots as pdf files. To run the error analysis write 
```
python3 two_particles_error_analysis.py
```
this outputs the convergence rate of both RK4 and FE as well as the relative error plots.

### 100 particles in the Penning trap

The file that is used in the case with 100 particles in the Penning trap is "hundred_particles.cpp". How to compile the file without any compiler optimization

```
make compile_hundred_particles
```

With compiler optimization -O2 or -O3 use
```
make compile_hundred_particles_O2
make compile_hundred_particles_O3
```

The machine code takes in two user inputs, the first determines if the Columb interaction between the particles is turned on of off (0 if off, 1 if on), while the second user option determines the ampliutude $f$ in the time dependent potential. In this task we have runned the program with the values 0.1, 0.4 and 0.7 for $f$ with the Columb interaction turned off. To reproduce these results the commands are
```
./hundred_particles 0 0.1
./hundred_particles 0 0.4
./hundred_particles 0 0.7
```
Each of these programs outputs a .bin file called "hundred_particles_f.bin", where "f" is the value of the amplitude given into the program (and is written to the filename with two decimals). A python program called "plot_hundred_particles.py" take the .bin file as input and makes the plots for the task. The pyhton file takes in the value for $f$ with two digits, the commands to make the plots are therefore
```
python3 plot_hundred_particles 0.10
python3 plot_hundred_particles 0.40
python3 plot_hundred_particles 0.70
```
The file that is used to zoom in on one of the resonances is called "hundred_particles_zoom.cpp", where the range and stepsize is set, and is compiled with the O3 option
```
make compile_hundred_particles_zoom_O3
```
The machine code takes in two user inputs, the first determines if the Columb interaction between the particles is turned on of off (0 if off, 1 if on), while the second user option determines the ampliutude $f$ in the time dependent potential. In this task we have runned the program with the values 0.1 with and without the Coulomb interaction. To reproduce these results the commands are
```
./hundred_particles_zoom 0 0.1
./hundred_particles_zoom 1 0.1
```
Each of these programs outputs a .bin file called "hundred_particles_zoom_interaction_f.bin", where "interaction" is 0 or 1 depending if the interaction is turned off or on and "f" is the value of the amplitude given into the program (and is written to the filename with two decimals). A python program called "plot_hundred_particles_zoom.py" take the .bin file as input and makes the plots for the task. The pyhton file takes in the value for $f$ with two digits, the commands to make the plots are therefore
```
python3 plot_hundred_particles_zoom 0 0.10
python3 plot_hundred_particles_zoom 1 0.10
```

