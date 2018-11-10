// Copyright 2018 Alexey Korepanov & Terry Lyons

// This file contains p_var_: a template function for computing p-variation
// in a metric space.
// It requires as arguments:
// * float_t: base type for distances and p (use float, double or long double)
// * point_t: data point type
// * float_t dist(point_t, point_t): distance function (assumed symmetric and
//   satisfying the triangle inequality)

#include <cmath>
#include <vector>

template <typename float_t, typename point_t, float_t (*dist)(point_t, point_t)>
float_t p_var_(const std::vector<point_t>& path, float_t p) {

    if (path.size() <= 1) {
        return 0.0;
    }

    // running p-variation in p-th power
    std::vector<float_t> run_p_var(path.size(), 0.0);

    // compute N = ceil(log2(path.size))
    size_t N = 1;
    for (size_t n = path.size(); n >>= 1; ) {
        N++;
    }

    // spatial index:
    // ind[n,k] = max { dist(path[k << n], path[(k << n) + m]  : 0 <= m < (1 << n) }
    std::vector < std::vector<float_t> > ind(N + 1);
    for (size_t n=1; n<=N; n++) {
        ind[n].resize( 1 << (N - n), 0.0);
    }

    float_t max_p_var = 0.0;

    for (size_t j=1; j<path.size(); j++) {
        // update ind
        for (size_t n=1; n<=N; n++) {
            size_t k = j >> n;
            ind[n][k] = std::max(ind[n][k], dist(path[k << n], path[j]));
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
                if (r  <  dist(path[mm], path[j]) + ind[nn][kk]) {
                    break;
                } else {
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
                max_p_var = std::max(max_p_var, run_p_var[m] + std::pow(d, p));
            }
            r = std::pow(max_p_var - run_p_var[m], 1./p);

        } while (m != 0);

        run_p_var[j] = max_p_var;
    }

    return std::pow(run_p_var.back(), 1./p);
}
