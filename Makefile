makeT: SO6.cpp Z2.cpp main.cpp
	g++ -pthread -Wall -ggdb main2.cpp SO6.cpp Z2.cpp T_Hist.cpp -O0 -o main.out
