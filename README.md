# P-var
Efficient computation of p-variation in metric spaces.

Let p≥1. For a discrete time path X<sub>0</sub>,...,X<sub>N</sub> in a metric space with distance d, its p-variation is
<p align="center">
  sup Σ<sub>k</sub> d(X<sub>n<sub>k</sub></sub>, X<sub>n<sub>k-1</sub></sub>)<sup>p</sup>,
</p>
where the supremum is taken over all increasing subsequences n<sub>k</sub> of 0,...,N.
Sometimes one takes p-th root of the sum, but here we don't.

There is a very efficient algorithm of computing p-variation for real-valued processes, see
[paper by Vygantas Butkus and Rimas Norvaiša](https://link.springer.com/article/10.1007/s10986-018-9414-3).
But it does not work for example in R<sup>2</sup>. Here we rectify this, at least partially.
We provide a short C++ function which computes p-variation in a general metric space,
and is sufficiently fast for paths with millions of points.

In addition, we improve the speed of the one-dimensional algorithm by a factor of 2 to 3.

## Usage
Below is a minimal working example. Details and important notes are in [`p_var.h`](p_var.h).
```
#include <iostream>
#include <vector>
#include <array>

// p_var is a one file header-only library
#include "p_var.h"
using p_var_ns::p_var;

// define a type for data points, in this case R^2
typedef std::array<double, 2> vecR2;

// define a distance function, here L^1
double dist(vecR2 a, vecR2 b) {
	return std::abs(b[0]-a[0]) + std::abs(b[1]-a[1]);
}

int main () {
	// create a path
	std::vector<vecR2> path({{0,0}, {1,1}});

	auto pv = p_var(path, 3, dist);
	// pv.value is the p-variation, and
	// pv.points is vector<size_t> with the maximising subsequence

	std::cout << "3-variation wrt L^1 distance: " << pv.value << std::endl;
	std::cout << "3-variation wrt Euclidean distance: " << p_var(path, 3).value << std::endl;

	return 0;
}
```

## Limitations
Our method is fast on data such as simulated Brownian paths, with complexity of
perhaps N log(N). But its worst case complexity is N<sup>2</sup>.
This happens in pathological situations such as monotone paths in R.
The [one-dimensional method](https://link.springer.com/article/10.1007/s10986-018-9414-3)
works with such paths much faster.

## Authors and contributors (in alphabetic order)
* Alexey Korepanov
* Terry Lyons
* Pavel Zorin-Kranich
