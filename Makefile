.PHONY: default short-test graph

default:
	g++     test.cpp p_var_real.cpp -Wall -Wextra -pedantic -O3 -march=native -o test.gcc.x
	clang++ test.cpp p_var_real.cpp -Wall -Wextra -pedantic -O3 -march=native -o test.clang.x

short-test:
	g++ short-test.cpp -Wall -Wextra -pedantic -O3 -march=native -o short-test.x

graph:
	g++ graph.cpp p_var_real.cpp -Wall -Wextra -O3 -march=native -o graph.x
