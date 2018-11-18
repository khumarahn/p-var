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

	// integer type
	template <typename point_t>
	auto euclidean_distance(point_t a, point_t b, rank<10>) ->
	typename std::enable_if<std::is_integral<point_t>::value, double>::type
	{
		return decltype(sqrt(std::abs(b - a)))(std::abs(b - a));
	}
	
	// scalar type
	template <typename point_t>
	auto euclidean_distance(point_t a, point_t b, rank<5>) ->
	typename std::enable_if<std::is_arithmetic<point_t>::value, point_t>::type
	{
		return std::abs(b - a);
	}
	
	// (nested) vector or array type
	template<typename point_t>
	auto euclidean_distance(const point_t & a, const point_t & b, rank<0>) ->
	decltype(sqrt(dist(*std::cbegin(a), *std::cbegin(a))))
	{
		typedef decltype(dist(*std::begin(a), *std::begin(a))) return_t;
		return std::sqrt(std::inner_product(std::begin(a), std::end(a), std::begin(b), return_t(0)
			, std::plus<>()
//			, [](auto a, auto b) { auto ds = std::abs(b - a); return ds * ds ;}
            , [](auto a, auto b) { auto ds = dist(a,b); return ds * ds; }
		));
	}
	

} //namespace kdfghaglhk

// new basic types with distance operators should be added here
// complex types
template<typename coord_t>
coord_t dist(const std::complex<coord_t>& a, const std::complex<coord_t>& b)
{
	return std::abs(b - a);
}

// the main function template for dist(a,b)
template <class point_t>
auto dist(const point_t & a, const point_t & b) { return kdfghaglhk::euclidean_distance(a, b, kdfghaglhk::rank<20>()); }

