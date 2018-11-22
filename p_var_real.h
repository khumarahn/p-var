#pragma once

#include <numeric>
#include <cmath>
#include <queue>
#include <vector>
#include <list>

namespace p_var_real {
	typedef std::vector<double> NumericVector;
	double pvar(const NumericVector& x, double p);
}
