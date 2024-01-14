# Project 1


## Information about the code used in the project
In this project the following list explains which code files are relevant for which task and how to run/compile them:


### The backbone, AKA functions.cpp and functions.hpp
Two files which are relevant for almost every task are functions.cpp and functions.hpp, these are the files in which we ahve declared and written some less important functions and more importantly the general and special algorithm. These are technically not to be compiled alone as they are essentially packages used by other scripts. 


### Task 2: 

The file Task2.cpp essentially just writes out the function u(x) for a given number of data-points. Note while you give it N datapoints it writes out N+2 as N just declares the "internal" datapoints. This is to keep consistency with later code and we do not need to think much about it.  


The .cpp file is compiled as 
```
g++ task2.cpp -o task2
```
and you run the machine code with


```
./task2 N
```
where N is the number of internal datapoints you want. It will then write u(x) for the linearly spaced number of points on the domain [0,1] to a file called "task2Output.txt" which is to be read by the file task2Read.py. The pythonfile task2Read.py creates a plot of the data from task2Output.txt. To run this file write the following in the terminal


```
python3 task2Read.py
```




### Task 7:

The relevant code in task 7 is found in the file task7.cpp and task7Read.py. Note task7.cpp uses functions.cpp, hence to compile it we must write

```
g++ task7.cpp functions.cpp -o task7
```
To run the compiled code 


```
./task7 N
```
The machine code takes in one argument N, which again is the number of data-points we wish to run over. The code the writes the relevant data given by the generalised algorithm to a file called "task7OutputN.txt" where N is the same as the agument which the code takes. The code also includes the external points which we enforce the value of. The file task7Read.py then takes in the files written by the machine code for N = 10,100,1000,10000 and graphs them agains the analytical solution. To run the python code 

```
python3 task7Read.py
```




### Task 8:

Task 8 was done with task7.cpp and task8.py, task7.cpp was used to write the relvant ".txt" files with data and the analysis of the data was completed in task8.py. Compiling task7.cpp was done as in the previous exercise. For 8.a and 8.b the code reads from the files "task7OutputN.txt" with N=10,100,1000,10000 and plots the figures. In 8.c it reads from the same files, but with N=10,100,1000,10000,1000000,10000000 and prints out a table with the relevant information to the terminal. To run the python code write in the temrinal  

```
python3 task8.py
```


### Task 9:

The code which is relvant for this task is task9.cpp and functions.cpp, task9.cpp is essentially a copy of task7.cpp, however now instead of using the general algorithm it uses the special, and the ".txt" files it writes out to are called "task9OutputN.txt" again with N being adaptable in the terminal. The special algorithm can be found in functions.cpp as the function "forwardback_special". To compile the code
```
g++ task9.cpp functions.cpp -o task9
```

```
./task9 N
```
with N being the number of internal points desired. 


### Task 10:


The only file relevant for solving task10 is task10.cpp. To compile this file write
```
g++ task10.cpp functions.cpp -o task10
```
and to run the compiled code
```
./task10 N I
```

where N is the number of internal datapoints and I is the number of iterations we are to average the measured time over, both should be integers. The output of this script to the terminal is the average time which the general algorithm used, the average time which the special algorithm used and the ratio between the two. An example of how the output should look in the terminal is 

```
General: <t>= 6.58e-06 |  n= 100 |  I= 10
Special: <t>= 3.85e-06 |  n= 100 |  I= 10

Special time/general time : 0.585106
```


