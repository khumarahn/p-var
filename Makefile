.PHONY: default short-test graph

default:
	g++ test.cpp p_var_real.cpp -Wall -O3 -march=native -o test.x

short-test:
	g++ short-test.cpp -Wall -O3 -march=native -o short-test.x

graph:
	g++ graph.cpp p_var_real.cpp -Wall -O3 -march=native -o graph.x
