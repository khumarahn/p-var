#pragma once
#include <type_traits>
#include <numeric>
#include <complex>

// defines a template function dist(a,b) for computing euclidean distance:
// for arithmetic and complex types : std::abs(b - a)
// for sequential container types and arrays it is defined recursively:
// sqrt(sum(d(x_i,y_i)^2)) 
// and is terminated by reaching an arithmetic or complex type.
//
// dist(a,b) returns the the distance in the type returned by sqrt 
// for the underlying arithmetic or complex type
//
template <class point_t>
auto dist(const point_t & a, const point_t & b);

namespace kdfghaglhk
{
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
	

} //namespace kdfghaglhk

// the main function template for dist(a,b)
template <class point_t>
auto dist(const point_t & a, const point_t & b) { return kdfghaglhk::euclidean_distance(a, b, kdfghaglhk::rank<2>()); }

