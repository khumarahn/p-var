#pragma once

#include <vector>

namespace p_var_real {
	typedef std::vector<double> NumericVector;
	double pvar(const NumericVector& x, double p);
}
