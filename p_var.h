// Copyright 2018 Alexey Korepanov & Terry Lyons
#ifndef p_var_h__
#define p_var_h__

/*
 * This file contains a template for p_var: a function for computing p-variation
 * in a metric space.
 *
 * Usage:
 * 1. p_var(path, p, dist)
 * where path is a vector, dist is a distance function and p is a positive real number.
 * 2. p_var(iterator_begin, iterator_end, p, dist)
 * if the path is in a container with random access iterators.
 * p_var(path, p, dist) is the same as p_var(path.begin(), path.end(), p, dist)
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
#include <limits>
#include <type_traits>

template <typename power_t, typename func_t, typename const_iterator_t>
auto p_var(const_iterator_t path_begin, const_iterator_t path_end, power_t p, func_t dist)
-> decltype(pow(dist(typename std::iterator_traits<const_iterator_t>::value_type(), typename std::iterator_traits<const_iterator_t>::value_type()), p)) {
	typedef decltype(pow(dist(typename std::iterator_traits<const_iterator_t>::value_type(), typename std::iterator_traits<const_iterator_t>::value_type()), p)) float_t;
	typedef typename std::iterator_traits<const_iterator_t>::difference_type diff_t;
	// this computation uses only p and two other things:
	// path size
	diff_t path_size = path_end - path_begin;
	// distances between path[a] and path[b]:
	auto path_dist = [&dist, path_begin](diff_t a, diff_t b) {
		return dist(*(path_begin + a), *(path_begin + b));
	};

	if (path_size == 0) {
		return -std::numeric_limits<float_t>::infinity();
	}
	else if (path_size == 1) {
		return float_t(0);
	}

	// running p-variation
	std::vector<float_t> run_p_var(path_size, 0);

	// compute N = log2(path_size)
	diff_t N = 0;
	for (diff_t n = path_size; n >>= 1; ) {
		N++;
	}

	// spatial index:
	// for 0 <= j < path_size and 1 <= n <= N,
	//   ind(j, n) = max { path_dist(k, k + m)  :  k = (j << n) >> n  and  0 <= m < (1 << n) }
	// we store ind(j, n) in ind[ind_n(j,n)] with a suitable function ind_n
	std::vector<float_t> ind(path_size, 0);
	auto ind_n = [&N](diff_t j, diff_t n) {
		return (1 << (N-n)) - 1 + (j >> n);
	};

	float_t max_p_var = float_t(0);

	for (diff_t j = 1; j < path_size; j++) {
		// update ind
		for (diff_t n = 1; n <= N; n++) {
			diff_t k = (j >> n) << n;
			float_t &i = ind[ind_n(j, n)];
			i = std::max<float_t>(i, path_dist(k, j));
		}

		// compute max_p_var = p-variation of path[0..j] as
		//   max{run_p_var[m] + path_dist(m, j)^p}
		// as m goes through j-1,...,0
		diff_t m = j - 1;
		float_t delta = 0;
		diff_t delta_m = j;
		for (diff_t n=N;;) {
			diff_t k = (m >> n) << n;
			// using spatial index, we skip all m >= k when n > 0 and
			// ind[ind_n(m, n)] + path_dist(k, j) < (max_p_var - run_p_var[m])^(1/p)
			bool skip = false;
			if (n > 0) {
				float_t id = ind[ind_n(m, n)] + path_dist(k, j);
				if (delta >= id) {
					skip = true;
				}
				else if (m < delta_m) {
					delta = std::pow(max_p_var - run_p_var[m], 1. / p);
					delta_m = m;
					if (delta >= id) {
						skip = true;
					}
				}
			}

			if (skip) {
				if (k > 0) {
					m = k - 1;
					while (n < N && (k>>n) % 2 == 0) {
						n++;
					}
				}
				else {
					break;
				}
			}
			else {
				if (n > 1) {
					n--;
				}
				else {
					float_t d = path_dist(m, j);
					if (d > delta) {
						max_p_var = std::max<float_t>(max_p_var, run_p_var[m] + std::pow(d, p));
					}

					if (m > 0) {
						n = 0;
						while (n < N  &&  (m>>n) % 2 == 0) {
							n++;
						}
						m--;
					}
					else {
						break;
					}
				}
			}
		}

		run_p_var[j] = max_p_var;
	}

	return run_p_var.back();
}

template <typename power_t, typename point_t, typename  func_t>
auto p_var(const std::vector<point_t>& path, power_t p, func_t dist)
-> decltype(pow(dist(std::declval<point_t>(), std::declval<point_t>()), p)) {
	return p_var(std::cbegin(path), std::cend(path), p, dist);
}

#endif // p_var_h__
