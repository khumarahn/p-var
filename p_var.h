// Copyright 2018 Alexey Korepanov & Terry Lyons
#pragma once

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
 * 5. If dist is omitted, we try a default function for Euclidean distance.
 *    It works on scalar types and multidimensional arrays of scalars.
 */

#include <cmath>
#include <vector>
#include <limits>
#include <numeric>
#include <iterator>

namespace p_var_ns {

namespace internal {
	// forward declaration of the default distance function
	template <class point_t>
	auto dist(const point_t & a, const point_t & b);
	// alias for the type of distance function
	template <typename point_t>
	using dist_func_t = decltype(dist(std::declval<point_t>(), std::declval<point_t>()))(*) (const point_t &, const point_t &);

	// helper type aliases
	template <typename iterator_t>
	using iterator_value_t = typename std::iterator_traits<iterator_t>::value_type;
	template <typename container_t>
	using container_iterator_value_t = iterator_value_t<typename container_t::iterator>;
}

// forward declaration of the p-variation backbone computation
template <typename func_t, typename power_t>
auto p_var_backbone(size_t path_size, power_t p, func_t path_dist);


// *** INTERFACE ***
// iterators
template <typename power_t, typename const_iterator_t,
	 typename func_t = internal::dist_func_t<internal::iterator_value_t<const_iterator_t> > >
auto p_var(const_iterator_t path_begin, const_iterator_t path_end, power_t p, func_t dist = internal::dist) {
	auto path_dist = [&path_begin,&dist](size_t a, size_t b) {
		return dist(*(path_begin + a), *(path_begin + b));
	};
	return p_var_backbone(path_end - path_begin, p, path_dist);
}
// vector or array or else
template <typename power_t, typename vector_t, typename func_t = internal::dist_func_t<internal::container_iterator_value_t<vector_t> > >
auto p_var(vector_t path, power_t p, func_t dist = internal::dist) {
	return p_var(std::cbegin(path), std::cend(path), p, dist);
}


// *** BACKBONE ***
// Input:
// * path_size >= 0 integer
// * p >= 1 real
// * path_dist(a, b), a function defined for all 0 <= a,b < path_size such that
//   a)  path_dist(a,b) + path_dist(b,c) <= path_dist(a,c)
//   b)  path_dist(a,b) = path_dist(b,a)
// Output:
// \max \sum_k path_dist(a_k, a_{k+1})^p
// over all increasing subsequences a_k of 0,...,path_size-1
template <typename func_t, typename power_t>
auto p_var_backbone(size_t path_size, power_t p, func_t path_dist)
{
        typedef decltype(path_dist(0, 0)) dist_t;
        typedef decltype(std::pow(path_dist(0, 0), p)) real_t;

	if (path_size == 0) {
		return -std::numeric_limits<real_t>::infinity();
	}
	else if (path_size == 1) {
		return real_t(0);
	}

	// running p-variation
	std::vector<real_t> run_p_var(path_size, 0);

	size_t s = path_size - 1;

	size_t N = 1;
	while (s >> N) {
		N++;
	}

	// spatial index:
	// for 0 <= j < path_size and 1 <= n <= N,
	// * let  a = (j << n) >> n  and  b = min{a + (1 >> n), path_size}
	//   think of [a,b) is a "level n diadic interval" containing j
	// * choose k = ind_k(j, n) so that it is somewhere in the middle of [a,b)
	// * compute ind(j, n) = max { path_dist(k, m) : a <= m < b}
	// * store ind(j, n) in a flat array ind[] at position ind_n(j,n) with a suitable function ind_n
	std::vector<dist_t> ind(s, 0.0);
	auto ind_n = [s](size_t j, size_t n) {
		return (s >> n) + (j >> n);
	};
	auto ind_k = [s](size_t j, size_t n) {
		return std::min<size_t>(((j >> n) << n) + (1 << (n-1)), s);
	};

	real_t max_p_var = real_t(0);

	for (size_t j = 0; j < path_size; j++) {
		// update ind
		for (size_t n = 1; n <= N; n++) {
			if (!(j >> n == s >> n && (s >> (n-1)) % 2 == 0)) {
				dist_t &i = ind[ind_n(j, n)];
				i = std::max<dist_t>(i, path_dist(ind_k(j, n), j));
			}
		}
		if (j == 0) {
			continue;
		}

		// compute max_p_var = p-variation of path[0..j] as
		//   max{run_p_var[m] + path_dist(m, j)^p}
		// as m goes through j-1,...,0
		size_t m = j - 1;
		real_t delta = 0;
		size_t delta_m = j;
		for (size_t n=0;;) {
			while (n > 0 && m >> n == s >> n && (s >> (n-1)) % 2 == 0) {
				n--;
			}

			// using spatial index, we skip all m >= ((m << n) >> n) when n > 0 and
			// ind[ind_n(m, n)] + path_dist(ind_k(m, n), j) < (max_p_var - run_p_var[m])^(1/p)
			bool skip = false;
			if (n > 0) {
				dist_t id = ind[ind_n(m, n)] + path_dist(ind_k(m, n), j);
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
				size_t k = (m >> n) << n;
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
					dist_t d = path_dist(m, j);
					if (d > delta) {
						max_p_var = std::max<real_t>(max_p_var, run_p_var[m] + std::pow(d, p));
					}

					if (m > 0) {
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

namespace internal {
	// We define a template function dist(a,b) for computing euclidean distance:
	// for arithmetic and complex types : std::abs(b - a)
	// for sequential container types and arrays it is defined recursively:
	// sqrt(sum(dist(x_i,y_i)^2))
	// and is terminated by reaching an arithmetic or complex type.
	//
	// dist(a,b) returns the the distance in the type returned by sqrt
	// for the underlying arithmetic or complex type

	template<int I> struct rank : rank<I-1> { static_assert(I > 0, ""); };
	template<> struct rank<0> {};

	// integer type, distance is double
	template <typename point_t>
	auto euclidean_distance(point_t a, point_t b, rank<1>) ->
	typename std::enable_if<std::is_integral<point_t>::value, double>::type
	{
		return decltype(std::sqrt(std::abs(b - a)))(std::abs(b - a));
	}

	// type without iterators, on which std::abs works, hopefully scalar
	template <typename number_t>
	auto euclidean_distance(const number_t & a, const number_t & b, rank<0>) ->
	decltype(std::abs(b - a))
	{
		return std::abs(b - a);
	}

	// (nested) vector or array type
	template<typename point_t>
	auto euclidean_distance(const point_t & a, const point_t & b, rank<2>) ->
	decltype(std::sqrt(dist(*std::cbegin(a), *std::cbegin(a))))
	{
		typedef decltype(dist(*std::begin(a), *std::begin(a))) return_t;
		return std::sqrt(std::inner_product(std::begin(a), std::end(a), std::begin(b), return_t(0)
					, std::plus<>()
					, [](auto a, auto b) { auto ds = dist(a,b); return ds * ds; }
					));
	}

	// the main function template for dist(a,b)
	template <class point_t>
	auto dist(const point_t & a, const point_t & b) { return euclidean_distance(a, b, rank<2>()); }

} // namespace internal

} // namespace p_var_ns
