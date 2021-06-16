makeT: SO6.cpp Z2.cpp main.cpp
	g++ -pthread main2.cpp SO6.cpp Z2.cpp T_Hist.cpp -O3 -o main.out