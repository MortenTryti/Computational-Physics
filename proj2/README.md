# README

Hello user! Welcome to the second project in FYS4150. This README is to better understand how to compile the code.


### General information about code structuring

So to compile all of this we have made a make file which is rather helpful. Moreover most if not all of the functions we have used have been written in the functions.cpp and functions.hpp files. 


### Task2 
How to compile the file
```
make compile2
```
To run the produced machine code type in the terminal 
```
./task2
```

and it should output the eigenvalues and eigenvectors found by Armadillo, then the analytical eigenvalues and eigenvectors. The eigenvalues are sorted, however some of the egienvalues will have to be matched by hand as they are not sorted, up to a multiplicative factor of -1 they should be the same.  


### Task3 
How to compile the file
```
make compile3
```
As you mught have guessed, the way to run the program is just 
```
./task3
```
where in the terminal you will be told what the largest modulus element is, which column and row to find it in and also a print of the matrix so the user can check by eye.


### Task4 
How to compile the file (I think you get the structure now)


```
make compile4
```
then to run the machine code 
```
./task4 N
```
where $N$ is the dimension of the $N\times N$ matrix we wish to solve for. The output here should be the diagonal form of the matrix A, just so the user can give it a quick check, the eigenvalues of the iterative method, the analytic method as well as the difference between the two. Naturally both are sorted. 

### Task5
How to compile the file
```
make compile5
```

Then run the machine code
```
./task5
```
The machine code calculates the number of transforamtions for different matrix sizes $N = 10, 20, ..., 100$ and writes it to a file called "task5Output.txt". The data is then taken as input into "task5Read.py" that makes a log-log plot of the data and another plot with a fit of the data. The plots are saved as "task5Plot_log.pdf" and "task5Plot_fit.pdf".
To make the plots run
```
python3 task5Read.py
```

### Task 6
How to compile the file

```
make compile6
```
Then to run the machine code 
```
./task4 N
```
where $N$ is the number of steps in the discretization of $x$. The machine code calculates the three lowest eigenvectors using the Jacobi rotation algorithm and the analytical solution and writes them to a file "task6OutputN_Jacobi.txt" and task6OutputN_analytic.txt", where $N$ is the number of steps. The data from these files are then taken as input into "task6Read.py" that makes a plot of the the three lowest eigenvectors for $N=10$ and $N=100$ as a function of $x$. The plots are saved in "task6Plot_10.pdf" and "task6Plot_100.pdf". To produce these plots one therefore has to run
```
./task 10
./task 100
python3 task6Read.py
```