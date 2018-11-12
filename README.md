# P-var
Efficient computation of p-variation in metric spaces.

Let p≥1. For a discrete time path X<sub>0</sub>,...,X<sub>N</sub> with in a metric space with distance d, its p-variation is
<p align="center">
  sup Σ<sub>k</sub> d(X<sub>n<sub>k</sub></sub>, X<sub>n<sub>k-1</sub></sub>)<sup>p</sup>,
</p>
where the supremum is taken over all increasing subsequences n<sub>k</sub> of 0,...,N.
Sometimes one takes p-th root of the sum, but here we don't.

There is a [very efficient way](https://link.springer.com/article/10.1007/s10986-018-9414-3)
of computing p-variation for real-valued processes, but it does not work
for example in R<sup>2</sup>. Here we rectify this, at least partially. We provide a short C++ function
which computes p-variation in a general metric space, and is sufficiently fast to
work with paths with millions of points.

## Usage
```
#include <iostream>
#include <vector>
#include <array>

// p_var is a template, so it is enough to include the header
#include "p_var.h"

// type for data points, in this case planar vectors
typedef std::array<double, 2> vecR2;

// distance function
double dist(vecR2 a, vecR2 b) {
	return std::abs(b[0]-a[0]) + std::abs(b[1]-a[1]);
}

int main () {
	// create and initialize a path
	std::vector<vecR2> path(2);
	path[0][0] = path[0][1] = 0;
	path[1][0] = path[1][1] = 1;
	
	// get its 3-variation
	double pv = p_var(path, 3, dist); 
	std::cout << pv;
}

```

# Limitations
Our method is fast on data such as simulated Brownian paths, with complexity of
perhaps N log(N). But its worst case complexity is N<sup>2</sup>.
This happens in pathological situations such as monotone paths in R.
The [one-dimensional method](https://link.springer.com/article/10.1007/s10986-018-9414-3)
works with monotone paths much faster.
