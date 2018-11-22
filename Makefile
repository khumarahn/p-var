.PHONY: default

default:
	g++ test.cpp p_var_real.cpp -Wall -O3 -march=native -o test.x

test-dist:
	g++ test-dist.cpp -Wall -O3 -march=native -o test-dist.x
