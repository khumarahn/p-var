#include <iostream>
#include <vector>
#include <array>

// p_var is a one file header-only library
#include "p_var.h"
using p_var_ns::p_var;

// define a type for data points, in this case R^2
//typedef std::array<double, 2> vecR2;
typedef double vecR2[2];
// define a distance function, here L^1
double dist(const vecR2 & a, const vecR2 & b) {
	return std::abs(b[0]-a[0]) + std::abs(b[1]-a[1]);
}

int main () {
	// create a path
	//std::vector<vecR2> path {{0,0}, {1,1}};
	vecR2 path[2]{ {0,0}, {1,1} };

	std::cout << "3-variation wrt L^1 distance: " << p_var(path, 3, dist).value << std::endl;
	std::cout << "3-variation wrt Euclidean distance: " << p_var(path, 3).value << std::endl;

	return 0;
}
