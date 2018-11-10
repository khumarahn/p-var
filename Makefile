.PHONY: default

default:
	g++ test.cpp -Wall -O3 -march=native -o test.x
