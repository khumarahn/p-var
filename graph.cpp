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

// constructor for random Brownian paths
std::vector<double> make_brownian_path(double sd, int steps) {
	std::default_random_engine generator(0);
	std::normal_distribution<double> gauss_dist(0.0, sd);
	std::vector<double> increments(steps), path(1, 0.);
	auto gauss_rv = [&]() { return gauss_dist(generator); }; // bind the generator
	std::generate(begin(increments), end(increments), gauss_rv);
	std::partial_sum(begin(increments), end(increments), back_inserter(path));
	return path;
}

// the test code
int main() {
	using std::cout;
	
	size_t steps = 1000*1000;
	double sd = 1 / sqrt(double(steps));
	std::vector<double> path = make_brownian_path(sd, steps);

	cout 	<< std::setw(15) << "p"
		<< std::setw(15) << "p-var"
		<< std::setw(15) << "p-var^{1/p}"
		<< std::setw(15) << "secs"
		<< std::setw(15) << "R mthd secs"
		<< std::setw(15) << "Answers diff"
		<< "\n";
	for (double p = 1.0; p < 8.01; p += 1./256) {
		clock_t clock_begin = std::clock();
		double pv = p_var(path, p);
		clock_t clock_end = std::clock();

		clock_t ref_clock_begin = std::clock();
		double pv_ref = p_var_real::pvar(path, p);
		clock_t ref_clock_end = std::clock();

		double elapsed_secs = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
		double ref_elapsed_secs = double(ref_clock_end - ref_clock_begin) / CLOCKS_PER_SEC;

		cout	<< std::setw(15) << p
			<< std::setw(15) << pv
			<< std::setw(15) << std::pow(pv, 1./p)
			<< std::setw(15) << elapsed_secs
			<< std::setw(15) << ref_elapsed_secs
			<< std::setw(15) << pv - pv_ref
			<< std::flush << "\n";
	}

	return 0;
}
