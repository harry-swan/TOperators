makeT: SO6.cpp Z2.cpp main.cpp
	g++ -pthread -ggdb main2.cpp SO6.cpp Z2.cpp T_Hist.cpp -O3 -o main.out