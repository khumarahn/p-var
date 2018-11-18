// Copyright 2018 Alexey Korepanov & Terry Lyons

#include <iostream>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <utility>
#include <vector>
#include <array>
#include <ctime>

#include "p_var.h"

// p_var requires a vector of points, a distance function
// and a scalar p specifying the p-variation to be computed

// if the points are double in R^1 we use the usual usual distance
double distR1(double a, double b) { return std::abs(b - a); };

// for R^d with the usual Euclidean distance we introduce a point type
// and metric
const size_t d = 2;
typedef std::array<double, d> vecRd;

double distRd(vecRd a, vecRd b) {
	double s = 0;
	for (size_t k = 0; k < d; k++) {
		double ds = std::abs(b[k] - a[k]);
		s += ds * ds;
	}
	return std::sqrt(s);
}

// reference implementation of p-variation for testing
double p_var_ref(const std::vector<double>& path, double p) {

	if (path.size() == 0) {
		return 	-std::numeric_limits<double>::infinity();
	}

	// running p-variation
	std::vector<double> run_p_var(path.size(), 0.0);

	for (size_t n = 1; n < path.size(); n++) {
		// compute run_p_var[n]
		double r = 0;
		for (size_t k = n; k-- > 0; ) { // for k in {n-1,...,0}
			// try last step [k,n]
			double prev = run_p_var[k];
			double step = std::pow(std::abs(path[n] - path[k]), p);
			r = std::max(r, prev + step);
		}
		run_p_var[n] = r;
	}

	return run_p_var.back();
}

// constructor for random Brownian paths
unsigned int random_seed(unsigned int & seed) {
	std::uniform_int_distribution<unsigned int> unif(0, std::numeric_limits<unsigned int>::max());
	std::random_device rnd_device;
	seed = unif(rnd_device);
	return seed;
}

std::vector<double> make_brownian_path(double sd, int steps) {
	unsigned int seed;
	std::default_random_engine generator(random_seed(seed));
	std::normal_distribution<double> gauss_dist(0.0, sd);
	std::vector<double> increments(steps), path(1, 0.);
	auto gauss_rv = [&]() { return gauss_dist(generator); }; // bind the generator
	std::generate(begin(increments), end(increments), gauss_rv);
	std::partial_sum(begin(increments), end(increments), back_inserter(path));
	return path;
}

void test_dist_template()
{
	std::vector<float> X = { 1., 2., 3. };
	std::vector<float> Y = { 2., 4., 5. };
	std::vector<std::vector<double>> XX = { {1., 2.}, { 3.,0.} };
	std::vector<std::vector<double>> YY = { {2., 4.}, { 5.,0.} };
	std::complex<float> xc[]{ 1.,2.,3. };
	std::complex<float> yc[]{ 2.,4.,5. };
	double x[]{ 1.,2.,3. };
	double y[]{ 2.,4.,5. };
	int xi[]{ 1,2,3 };
	int yi[]{ 2,4,5 };
	std::cout <<
		dist(X, Y) << " " <<
		dist(XX, YY) << " " <<
		dist(xc, yc) << " " <<
		dist(x, y) << " " <<
		dist(xi, yi) << " " <<
		dist(0., 3.) << " " <<
		dist(0, 3);
}

// the test code
int main() {
	using std::cout;

	// run tests on the default generic distance function
	test_dist_template(); // all answers should be three
	
	int test_no = 0;
	// check a long periodic path 0,1,4,0,1,4,...
	{
		cout << "*** TEST " << ++test_no << " ***\n";
		size_t rep = 100000;
		std::vector<double> path(rep * 4);
		for (size_t j = 0; j < path.size(); j++) {
			switch (j % 4) {
			case 0:
				path[j] = 0;
				break;
			case 1:
				path[j] = 1;
				break;
			case 2:
				path[j] = 1;
				break;
			case 3:
				path[j] = 4;
				break;
			}
		}
		cout << "Long periodic path: 0,1,1,4 repeated " << rep << " times\n";
		for (double p = 1.0; p < 3.01; p += 0.5) {
			//double pv = p_var(path, p, distR1);
			double pv = p_var(path, p);
			double pv_ref = std::pow(4, p) * (2 * rep - 1);
			cout << "  " << p << "-variation: " << pv
				<< ", error: " << pv - pv_ref
				<< "\n";
		}
	}

	// check a long monotone path from 0.0 to 1.0
	{
		cout << "\n*** TEST " << ++test_no << " ***\n";
		size_t steps = 10000;
		std::vector<double> path(steps + 1);
		for (size_t j = 0; j < path.size(); j++) {
			path[j] = double(j) / steps;
		}
		cout << "Long monotone path: "
			<< path[0] << ", " << path[1] << ", " << path[2]
			<< ", ... , " << path.back() << "\n";
		for (double p = 1.0; p < 3.01; p += 0.5) {
			//double pv = p_var(path, p, distR1);
			double pv = p_var(path, p);
			double pv_ref = 1.0;
			cout << "  " << p << "-variation: " << pv
				<< ", error: " << pv - pv_ref
				<< "\n";
		}
	}

	// path: unit square in R^2
	{
		cout << "\n*** TEST " << ++test_no << " ***\n";
		size_t steps = 4;
		std::vector<vecRd> path(steps + 1);
		for (size_t j = 0; j < path.size(); j++) {
			for (size_t k = 0; k < d; k++) {
				if (k == 0 && (j % 4 == 1 || j % 4 == 2)) {
					path[j][k] = 1;
				}
				else if (k == 1 && (j % 4 == 2 || j % 4 == 3)) {
					path[j][k] = 1;
				}
				else {
					path[j][k] = 0;
				}
			}
		}
		cout << "Simple path:";
		for (size_t j = 0; j < path.size(); j++) {
			cout << " (" << path[j][0];
			for (size_t k = 1; k < d; k++) {
				cout << ", " << path[j][k];
			}
			cout << ")";
		}
		cout << "\n";
		cout << "cumulative p-variation:\n";
		for (double p = 1.0; p < 3.01; p += 0.5) {
			cout << "  p=" << p << ": ";
			double pv = 0;
			std::vector<vecRd>::iterator it = path.begin();
			for (;;) {
				//pv = p_var(path.begin(), it, p, distRd);
				pv = p_var(path.begin(), it, p);
				cout << pv;
				if (it == path.end()) {
					break;
				}
				else {
					it++;
					cout << ", ";
				}
			}
			double pv_ref = (p > 2) ? std::pow(2.0, p*0.5 + 1.) : 4;
			cout << ";  error at the end: " << pv - pv_ref << "\n";
		}
	}

	// another short path in R^2
	{
		cout << "\n*** TEST " << ++test_no << " ***\n";
		double p = 3;
		size_t steps = 3;
		std::vector<vecRd> path(steps + 1);
		for (size_t j = 0; j < path.size(); j++) {
			for (size_t k = 0; k < d; k++) {
				if (k == 0 && (j % 4 == 2 || j % 4 == 3)) {
					path[j][k] = 1;
				}
				else if (k == 1 && (j % 4 == 0 || j % 4 == 1)) {
					path[j][k] = 1;
				}
				else {
					path[j][k] = 0;
				}
			}
		}
		cout << "Simple path:";
		for (size_t j = 0; j < path.size(); j++) {
			cout << " (" << path[j][0];
			for (size_t k = 1; k < d; k++) {
				cout << ", " << path[j][k];
			}
			cout << ")";
		}
		cout << "\n";
		//double pv = p_var(path, p, distRd);
		double pv = p_var(path, p);
		double pv_ref = pow(2.0, 0.5*p);
		cout << "  " << p << "-variation: " << pv
			<< ", error: " << pv - pv_ref
			<< "\n";
	}

	// compare to reference on random Brownian paths
	{
		cout << "\n*** TEST " << ++test_no << " ***\n";
		double p = 3;
		size_t count = 10;
		size_t steps = 10000;
		double max_err = 0.0;

		for (size_t c = 0; c < count; c++) {
			double sd = 1 / sqrt(double(steps));
			std::vector<double> path = make_brownian_path(sd, steps);
			//double pv = p_var(path, p, distR1);
			double pv = p_var(path, p);
			double pv_ref = p_var_ref(path, p);
			double err = std::abs(pv_ref - pv);
			max_err = std::max(max_err, err);
		}
		cout << count << " random Brownian paths of length " << steps
			<< " compared to reference\n";
		cout << "  max error: " << max_err << "\n";
	}

	// benchmark
	{
		cout << "\n*** TEST " << ++test_no << " (BENCHMARK) ***\n";
		double p = 3;
		for (size_t steps = 1; steps <= 10000000; steps *= 10) {
			double sd = 1 / sqrt(double(steps));
			std::vector<double> path = make_brownian_path(sd, steps);

			clock_t clock_begin = std::clock();
			//double pv = p_var(path, p, distR1);
			double pv = p_var(path, p);
			clock_t clock_end = std::clock();

			double elapsed_secs = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
			cout << "  a random Brownian path of length " << steps
				<< " has " << p << "-variation of " << pv << "\n"
				<< "  which took " << elapsed_secs << " seconds to compute\n";
		}
	}

	return 0;
}
