// Copyright 2018 Alexey Korepanov & Terry Lyons

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <utility>
#include <vector>
#include <array>
#include <ctime>

#include "dist.h"

void test_dist_template()
{
	std::vector<float> X = { 1., 2., 3. };
	std::vector<float> Y = { 2., 4., 5. };
	std::vector<std::vector<double>> XX = { {1., 2.}, { 3.,0.} };
	std::vector<std::vector<double>> YY = { {2., 4.}, { 5.,0.} };
	std::complex<float> xcs(1);
	std::complex<float> ycs(4);
	std::complex<float> xc[]{ 1.,2.,3. };
	std::complex<float> yc[]{ 2.,4.,5. };
	double x[]{ 1.,2.,3. };
	double y[]{ 2.,4.,5. };
	int xi[]{ 1,2,3 };
	int yi[]{ 2,4,5 };
	std::cout <<
		dist(X, Y) << " " <<
		dist(XX, YY) << " " <<
		dist(xcs, ycs) << " " <<
		dist(xc, yc) << " " <<
		dist(x, y) << " " <<
		dist(xi, yi) << " " <<
		dist(0., 3.) << " " <<
		dist(0, 3);
}

// the test code
int main() {
	using std::cout;
	int test_no = 0;

	// run tests on the default generic distance function
	{
		cout << "*** TEST " << ++test_no << " ***\n";
		cout << "Testing dist(a,b), all answers should be three: ";
		test_dist_template();
		cout << "\n";
	}

	return 0;
}
