// Copyright 2018 Alexey Korepanov & Terry Lyons
#ifndef p_var_h__
#define p_var_h__

// This file contains p_var: a template function for computing p-variation
// in a metric space.
//************************************
// Method:    p_var
// FullName:  p_var
// Access:    public 
// Returns:   auto
// Qualifier: -> decltype(pow(dist(std::declval<point_t>(), std::declval<point_t>()), p))
// Parameter: const std::vector<point_t> & path
// Parameter: func_t dist
// Parameter: power_t p
//************************************

// * p_var returns float_t which inherits its definition from  pow(dist(point, point),power)
// * point_t: the data point type
// * dist(point_t, point_t): the distance function (assumed symmetric and
//   satisfying the triangle inequality), consumable by pow
// * float_t allowing the value -infinity (the p-variation of an empty interval)
// * power_t the type of p which must be acceptable to pow.  

#include <cmath>
#include <vector>
#include <numeric>
#include <type_traits>

template <typename power_t, typename point_t, typename  func_t>
auto p_var(const std::vector<point_t>& path, func_t dist, power_t p) -> decltype(pow(dist(std::declval<point_t>(), std::declval<point_t>()), p)) {
	typedef decltype(pow(dist(std::declval<point_t>(), std::declval<point_t>()), p)) float_t;
	if (path.size() == 0) {
		return -std::numeric_limits<float_t>::infinity();
	}
	if (path.size() <= 1) {
		return float_t(0);
	}

	// running p-variation in p-th power
	std::vector<float_t> run_p_var(path.size(), 0.0);

	// compute N = ceil(log2(path.size))
	size_t N = 1;
	for (size_t n = path.size(); n >>= 1; ) {
		N++;
	}

	// spatial index:
	// ind[n,k] = max { dist(path[k << n], path[(k << n) + m])  : 0 <= m < (1 << n) }
	std::vector < std::vector<float_t> > ind(N + 1);
	for (size_t n = 1; n <= N; n++) {
		ind[n].resize(1 << (N - n), 0.0);
	}

	float_t max_p_var = 0.0;

	for (size_t j = 1; j < path.size(); j++) {
		// update ind
		for (size_t n = 1; n <= N; n++) {
			size_t k = j >> n;
			ind[n][k] = std::max<float_t>(ind[n][k], dist(path[k << n], path[j]));
		}

		// compute run_p_var[j]: p-th power of p-variation of path[0..j]
		size_t m = j;
		float_t r = 0;
		do {
			// reduce m
			size_t n = 0;
			while (n < N) {
				size_t nn = n + 1,
					kk = m >> nn,
					mm = kk << nn;
				if (r < dist(path[mm], path[j]) + ind[nn][kk]) {
					break;
				}
				else {
					n++;
				}
			}
			m = (m >> n) << n;
			if (m > 0) {
				m--;
			}

			// check if we improve p-variation
			float_t d = dist(path[m], path[j]);
			if (d > r) {
				max_p_var = std::max<float_t>(max_p_var, run_p_var[m] + std::pow(d, p));
			}
			r = std::pow(max_p_var - run_p_var[m], 1. / p);

		} while (m != 0);

		run_p_var[j] = max_p_var;
	}

	return run_p_var.back();
}

#endif // p_var_h__
