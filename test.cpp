// Copyright 2018 Alexey Korepanov & Terry Lyons

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <complex>
#include <numeric>
#include <algorithm>
#include <utility>
#include <vector>
#include <array>
#include <ctime>

#include "p_var.h"
#include "p_var_real.h"

using p_var_ns::p_var;

// usual distance in R^1
double distR1(double a, double b) { return std::abs(b - a); }

// R^d with L^1 distance
const size_t d = 2;
typedef std::array<double, d> vecRd;
double distRd(vecRd a, vecRd b) {
	double s = 0;
	for (size_t k = 0; k < d; k++) {
		double ds = std::abs(b[k] - a[k]);
		s += ds;
	}
	return s;
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
		for (size_t k = 0; k < n; k++) {
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

// path gelerated by an intermittent dynamical system, like Levy
// alpha should be in (1/2, 1)
std::vector<double> make_intermittent_path(int steps, double alpha = 0.6) {
	auto LSV = [alpha](double x) {
		if (x <= 0 || x >= 0.5) {
			return 0.0;
		}
		return x * (1. + std::pow(2*x, alpha));
	};
	auto SLSV = [&LSV](double x) {
		return (x <= 0.5) ? LSV(x) : (1. - LSV(1. - x));
	};
	auto phi = [](double x) {
		return x - 0.5;
	};
	
	unsigned int seed;
	std::default_random_engine generator(random_seed(seed));
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	double x = unif(generator);
	
	std::vector<double> path(steps + 1);
	path[0] = 0.;
	double norm = std::pow(steps + 1, -alpha);
	for (size_t k=1; k < path.size(); k++) {
		path[k] = path[k-1] + phi(SLSV(x)) * norm;
		x = SLSV(x);
	}
	return path;
}

void test_dist_template()
{
	using p_var_ns::p_var_dist_ns::dist;
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
	int test_no = 0;

	// run tests on the default generic distance function
	{
		cout << "*** TEST " << ++test_no << " ***\n";
		cout << "Testing dist(a,b), all answers should be three: ";
		test_dist_template();
		cout << "\n";
	}

	// check a long periodic path 0,1,4,0,1,4,...
	{
		cout << "\n*** TEST " << ++test_no << " ***\n";
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
		std::array<std::array<double, 2>, 5> path = {{{0,0}, {1,0}, {1,1}, {0,1}, {0,0}}};
		cout << "Simple path:";
		for (size_t j = 0; j < path.size(); j++) {
			cout << " (" << path[j][0];
			for (size_t k = 1; k < 2; k++) {
				cout << ", " << path[j][k];
			}
			cout << ")";
		}
		cout << "\n";
		cout << "cumulative p-variation:\n";
		for (double p = 1.0; p < 3.01; p += 0.5) {
			cout << "  p=" << p << ": ";
			double pv = 0;
			auto it = path.begin();
			for (;;) {
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
		double pv = p_var(path, p, distRd);
		double pv_ref = pow(2.0, p);
		cout << "  " << p << "-variation wrt L^1 distance: " << pv
			<< ", error: " << pv - pv_ref
			<< "\n";
		double pv2 = p_var(path, p);
		double pv_ref2 = pow(2.0, 0.5*p);
		cout << "  " << p << "-variation wrt Euclidean distance: " << pv2
			<< ", error: " << pv2 - pv_ref2
			<< "\n";
	}

	// compare to reference on random Brownian paths
	{
		cout << "\n*** TEST " << ++test_no << " ***\n";
		double p = 3;
		size_t count = 100;
		size_t steps = 10000;
		double max_err = 0.0;

		for (size_t c = 0; c < count; c++) {
			double sd = 1 / sqrt(double(steps));
			std::vector<double> path = make_brownian_path(sd, steps);
			//double pv = p_var(path, p, distR1);
			double pv = p_var(path, p);
			double pv_ref = p_var_real::pvar(path, p);
			double err = std::abs(pv_ref - pv);
			max_err = std::max(max_err, err);
		}
		cout << count << " random Brownian paths of length " << steps
			<< " compared to reference\n";
		cout << "  max error: " << max_err << "\n";
	}

	// benchmark
	{
		cout << "\n*** TEST " << ++test_no << ": BROWNIAN BENCHMARK ***\n";
		double p = 3;
		cout << "Computing p-variation with p=" << p << " for long random Brownian paths,\n"
			<< "benchmarking against the real line specific method\n"
			<< std::setw(15) << "Length"
			<< std::setw(15) << "p-variation"
			<< std::setw(15) << "Seconds"
			<< std::setw(15) << "R mthd secs"
			<< std::setw(15) << "Answers diff"
			<< "\n";
		for (size_t steps = 1; steps <= 10000000; steps *= 10) {
			double sd = 1 / sqrt(double(steps));
			std::vector<double> path = make_brownian_path(sd, steps);

			clock_t clock_begin = std::clock();
			double pv = p_var(path, p);
			clock_t clock_end = std::clock();

			clock_t ref_clock_begin = std::clock();
			double pv_ref = p_var_real::pvar(path, p);
			clock_t ref_clock_end = std::clock();

			double elapsed_secs = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
			double ref_elapsed_secs = double(ref_clock_end - ref_clock_begin) / CLOCKS_PER_SEC;

			cout	<< std::setw(15) << steps
				<< std::setw(15) << pv
				<< std::setw(15) << elapsed_secs
				<< std::setw(15) << ref_elapsed_secs
				<< std::setw(15) << pv - pv_ref
				<< "\n";
		}
	}
	
	// benchmark
	{
		cout << "\n*** TEST " << ++test_no << ": INTERMITTENT (LEVY) BENCHMARK ***\n";
		double alpha = 0.6;
		double p = 1. / alpha + 0.1;
		cout << "Computing p-variation with p=" << p << " for long paths generated by\n"
			<< "an intermittent dynamical system (LSV map in Levy regime),\n"
			<< "such paths have long monotone excursions;\n"
			<< "benchmarking against the real line specific method\n"
			<< std::setw(15) << "Length"
			<< std::setw(15) << "p-variation"
			<< std::setw(15) << "Seconds"
			<< std::setw(15) << "R mthd secs"
			<< std::setw(15) << "Answers diff"
			<< "\n";
		for (size_t steps = 1; steps <= 10000000; steps *= 10) {
			std::vector<double> path = make_intermittent_path(steps, alpha);

			clock_t clock_begin = std::clock();
			double pv = p_var(path, p);
			clock_t clock_end = std::clock();

			clock_t ref_clock_begin = std::clock();
			double pv_ref = p_var_real::pvar(path, p);
			clock_t ref_clock_end = std::clock();

			double elapsed_secs = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
			double ref_elapsed_secs = double(ref_clock_end - ref_clock_begin) / CLOCKS_PER_SEC;

			cout	<< std::setw(15) << steps
				<< std::setw(15) << pv
				<< std::setw(15) << elapsed_secs
				<< std::setw(15) << ref_elapsed_secs
				<< std::setw(15) << pv - pv_ref
				<< "\n";
		}
	}

	return 0;
}
