## README project 4

Project 4 was solved mainly with object oriented programming. The class Ising in the file ```ising.cpp``` and ```ising.hpp``` constructs the system and does most of the heavy lifting with the Monte Carlo method. It also includes a number of different functions with differing complexities where most of the mare explained in the docstrings.


### Task 4
To compile task 4 write 
```
make task4M
```
in the terminal. This will then make the machine code file, the machine code file takes the temperature $T$ and number of MC cycles $N$ as command line arguments, it can be run by typing the following in the terminal. 

```
./task4 T N
```
The output from this file will be the theoretical as well as the numerically obtained quantities. Further it also calculates an estimate for the average number of MC cycles a quantity need for the absolute difference between theoretical and numerical variables to be within a certain tolerance of each other. These quantities are then saved to a .bin file so they can be plotted in python.

To generate the plot write in the terminal .

```
python3 task4_plot.py
```

### Task 5

To compile task 5 
```
make task5M
```

and to run it 
```
./task5
```
this will generate some data which we need to plot. To plot this use
```
python3  plot_task5.py
```
### Task 6

To compile task 6 write
```
make task6M
```

to run
```
./task6
```
this saves the data to a .bin file which we then can use in the python script ```plot6.py```. To get the figures use

```
python3 plot6.py
```
this function also gives variances.

### Task 7 & 8 & 9
The data used in task 7,8 and 9 is made from task8.cpp which calculates $\langle \epsilon \rangle$, $\langle |m| \rangle$, $C_V$ and $\chi$ for different temperatures and can be compiled both with and without parallelisation. To compile task8.cpp without parallelisation run
```
make task8M
```
which produces the program ```task8```. To compile task8 with parallelisation run
```
make task8M-omp
```
which produces the program ```task8-omp```. Both task8 and ```task8-omp``` takes in three input parameters where the first is 0 or 1 that determines if the temperature region is [2.1,2.4] (called unzoom) or [2.2,2.35] (called zoom) respectively. The second argment is the number of temperature points $n$ used in the simulation and the third is the lattice size $L$. The general program can therefore be written as
```
./task8-omp (unzoom/zoom) n L
```
In the program the number of Monte Carlo cycles is set to one million with 10% of it as burn-in.
The output from ```task8``` and ```taks8-omp``` are for both a .bin file called ```variables_(zoom/unzoom)_n_L.bin``` depending on the values of the three variables described above. The variables 

### Task 7
To time the code with and without the parallelisation over the temeperature loop the following commands are used for the table in the report, where the elapsed time is printed in the end of the program

Without parallizarion
```
./task8 0 10 2
./task8 0 10 4
./task8 0 10 8
./task8 0 10 16
```
With parallizarion
```
./task8-omp 0 10 2
./task8-omp 0 10 4
./task8-omp 0 10 8
./task8-omp 0 10 16
```
### Task 8
In task 8 the .bin files from ```task8.cpp``` are read into a python script called ```plot8.py``` which is also used in taks 9. To produce the .bin files needed for task8 in the unzoomed range [2.1, 2.4] with 20 temereature points for $L=40,60,80$ and 100 run
```
./task8-omp 0 20 40
./task8-omp 0 20 60
./task8-omp 0 20 80
./task8-omp 0 20 100
```
After this run the python script ```plot8.py``` that reads in the .bin files, and takes also three inputrs where the two first are the same as ```task8.cpp```, so (unzoom/zoom) and the number of temperature points $n$. The third input is a string that is a list of all the lattice sizes one want to plot. The general program can therefore be written as
```
python3 plot8.py (zoom/unzoom) n "[L_1,L_2,...,L_N]
```
To make the plots of the variablas in task 8 one therefor needs to run
```
python3 plot8.py 0 20 "[40,60,80,100]"
```
which produces a pdf file for each variable.
### Task 9
In task 9 the new first argument is 1 (zoom) with temperature region [2.2,2.35] where 10 temperature points are used. The commands to produce the needed .bin files are therefore
```
./task8-omp 1 10 40
./task8-omp 1 10 60
./task8-omp 1 10 80
./task8-omp 1 10 100
```
The ```plot.py``` script takes agian in the .bin files from the commands and produce the plots for each variable. For zoom=1 ```plot8.py``` also makes a linear regression from the critical temperatures from obtained from the interpolation in the plots for $C_V$ and $\chi$. To produce all the plots for the variables and the linear regression one therefore needs to run
```
python3 plot8.py 1 10 "[40,60,80,100]"
```
whcih again saves the plots as pdf files.



