// Copyright 2018 Alexey Korepanov & Terry Lyons
#ifndef p_var_h__
#define p_var_h__

/*
 * This file contains a template for p_var: a function for computing p-variation
 * in a metric space.
 *
 * Usage:
 * 1. p_var(path, dist, p)
 * where path is a vector, dist is a distance function and p is a positive real number.
 * 2. p_var(iterator_begin, iterator_end, dist, p)
 * if the path is in a container with random access iterators.
 * p_var(path, dist, p) is the same as p_var(path.begin(), path.end(), dist, p)
 * See test.cpp for examples and benchmarks.
 *
 * Notes:
 * 1. Return type is that of std::pow(dist(), p).
 * 2. We do not take p-th root of the sum of increments.
 * 3. p-variation of an empty path is negative infinity.
 * 4. p-variation of a path with one point is zero.
 */

#include <cmath>
#include <vector>
#include <numeric>
#include <type_traits>

template <typename power_t, typename func_t, typename iterator>
auto p_var(iterator path_begin, iterator path_end, func_t dist, power_t p) 
	-> decltype(pow(dist(typename std::iterator_traits<iterator>::value_type(), typename std::iterator_traits<iterator>::value_type()), p)) {

	typedef decltype(pow(dist(typename std::iterator_traits<iterator>::value_type(), typename std::iterator_traits<iterator>::value_type()), p)) float_t;

	// this computation uses only p and two other things:
	// path size
	size_t path_size = path_end - path_begin;
	// distances between path[a] and path[b]:
	auto path_dist = [&dist,path_begin](size_t a, size_t b) {
		return dist(*(path_begin + a), *(path_begin + b));
	};

	if (path_size == 0) {
		return -std::numeric_limits<float_t>::infinity();
	} else if (path_size == 1) {
		return float_t(0);
	}

	// running p-variation
	std::vector<float_t> run_p_var(path_size, 0);

	// compute N = ceil(log2(path_size))
	size_t N = 1;
	for (size_t n = path_size; n >>= 1; ) {
		N++;
	}

	// spatial index:
	// ind[n,k] = max { dist(path[k << n], path[(k << n) + m])  : 0 <= m < (1 << n) }
	std::vector < std::vector<float_t> > ind(N + 1);
	for (size_t n = 1; n <= N; n++) {
		ind[n].resize(1 << (N - n), float_t(0));
	}

	float_t max_p_var = float_t(0);

	for (size_t j = 1; j < path_size; j++) {
		// update ind
		for (size_t n = 1; n <= N; n++) {
			size_t k = j >> n;
			ind[n][k] = std::max<float_t>(ind[n][k], path_dist(k << n, j));
		}

		// compute run_p_var[j]: the p-variation of path[0..j]
		size_t m = j;
		float_t r = 0;
		do {
			// reduce m
			size_t n = 0;
			while (n < N) {
				size_t nn = n + 1,
					kk = m >> nn,
					mm = kk << nn;
				if (r < path_dist(mm, j) + ind[nn][kk]) {
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
			float_t d = path_dist(m, j);
			if (d > r) {
				max_p_var = std::max<float_t>(max_p_var, run_p_var[m] + std::pow(d, p));
			}
			r = std::pow(max_p_var - run_p_var[m], 1. / p);

		} while (m != 0);

		run_p_var[j] = max_p_var;
	}

	return run_p_var.back();
}

template <typename power_t, typename point_t, typename  func_t>
auto p_var(const std::vector<point_t>& path, func_t dist, power_t p) -> decltype(pow(dist(std::declval<point_t>(), std::declval<point_t>()), p)) {
	return p_var<power_t, func_t, typename std::vector<point_t>::const_iterator>(path.begin(), path.end(), dist, p);
}

#endif // p_var_h__
