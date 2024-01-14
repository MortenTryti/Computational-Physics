
compile2:
	g++ ./src/task2.cpp ./src/header.cpp -o ./bin/task2 -larmadillo -std=c++11

run2:
	g++ ./src/task2.cpp ./src/header.cpp -o ./bin/task2 -larmadillo -std=c++11 
	./bin/task2

compile3:
	g++ ./src/task3.cpp ./src/header.cpp -o ./bin/task3 -larmadillo -std=c++11 

run3:
	g++ ./src/task3.cpp ./src/header.cpp -o ./bin/task3 -larmadillo -std=c++11 
	./bin/task3

compile4:
	g++ ./src/task4.cpp ./src/header.cpp -o ./bin/task4 -larmadillo -std=c++11 

run4:
	g++ ./src/task4.cpp ./src/header.cpp -o ./bin/task4 -larmadillo -std=c++11 
	./bin/task4

compile5:
	g++ ./src/task5.cpp ./src/header.cpp -o ./bin/task5 -larmadillo -std=c++11 


plot5_all:
	g++ ./src/task5.cpp ./src/header.cpp -o ./bin/task5 -larmadillo -std=c++11 
	./bin/task5 1
	./bin/task5 2
	./bin/task5 3
	python3 ./scripts/plot5.py


compile7:
	g++ ./src/task7.cpp ./src/header.cpp -o ./bin/task7 -larmadillo -std=c++11


run7Nopot:
	g++ ./src/task7.cpp ./src/header.cpp -o ./bin/task7 -larmadillo -std=c++11 
	./bin/task7 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.05 0.0 0.0 "noPot"


run7Wpot:
	g++ ./src/task7.cpp ./src/header.cpp -o ./bin/task7 -larmadillo -std=c++11 
	./bin/task7 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.10 0.0 1e10 "Wpot"

plot7_all:
	g++ ./src/task7.cpp ./src/header.cpp -o ./bin/task7 -larmadillo -std=c++11 
	./bin/task7 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.05 0.0 0.0 "noPot"
	g++ ./src/task7.cpp ./src/header.cpp -o ./bin/task7 -larmadillo -std=c++11 
	./bin/task7 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.10 0.0 1e10 "Wpot"
	python3 ./scripts/plot7.py

compile8:
	g++ ./src/task8.cpp ./src/header.cpp -o ./bin/task8 -larmadillo -std=c++11 


run8:
	g++ ./src/task8.cpp ./src/header.cpp -o ./bin/task8 -larmadillo -std=c++11 
	./bin/task8 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 



plot8:
	g++ ./src/task8.cpp ./src/header.cpp -o ./bin/task8 -larmadillo -std=c++11 
	./bin/task8 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 
	python3 ./scripts/plot8.py

compile9:
	g++ ./src/task9.cpp ./src/header.cpp -o ./bin/task9 -larmadillo -std=c++11 


run9_1s:
	g++ ./src/task9.cpp ./src/header.cpp -o ./bin/task9 -larmadillo -std=c++11 
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 1 "oneSlit"


run9_2s:
	g++ ./src/task9.cpp ./src/header.cpp -o ./bin/task9 -larmadillo -std=c++11 
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 2 "twoSlit"


run9_3s:
	g++ ./src/task9.cpp ./src/header.cpp -o ./bin/task9 -larmadillo -std=c++11 
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 3 "threeSlit"


run9_all:
	g++ ./src/task9.cpp ./src/header.cpp -o ./bin/task9 -larmadillo -std=c++11 
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 1 "oneSlit"
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 2 "twoSlit" 
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 3 "threeSlit"


plot9_all:
	g++ ./src/task9.cpp ./src/header.cpp -o ./bin/task9 -larmadillo -std=c++11 
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 1 "oneSlit"
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 2 "twoSlit"
	./bin/task9 0.005 0.000025 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 3 "threeSlit"
	python3 ./scripts/task9.py
