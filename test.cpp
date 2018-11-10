// Copyright 2018 Alexey Korepanov & Terry Lyons

#include <iostream>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <utility>
#include <vector>
#include <array>
#include <ctime>

#include "p_var.h"

bool debug = false;

// specialize p_var_ to R^1 with usual distance
double dist(double a, double b) {return std::abs(b-a);};
double p_var(const std::vector<double>& path, double d) {
    return p_var_<double, double, dist> (path,d);
}

// specialize p_var_ to R^d with Euclidean distance
const size_t d = 2;
typedef std::array<double, d> vec_Xd;

double dist_Xd(vec_Xd a, vec_Xd b) {
    double s = 0;
    for (size_t k=0; k < d; k++) {
        double ds = std::abs(b[k] - a[k]);
        s += ds * ds;
    }
    return std::sqrt(s);
}
double p_var_Xd(const std::vector<vec_Xd>& path, double d) {
    return p_var_<double, vec_Xd, dist_Xd> (path,d);
}

// reference implementation of p-variation
double p_var_ref(const std::vector<double>& path, double p) {

    if (path.size() == 0) {
        return 0;
    }

    // running p-variation in p-th power
    std::vector<double> run_p_var(path.size(), 0.0);

    for (size_t n=1; n<path.size(); n++) {
        // compute run_p_var[n]
        double r = 0;
        for (size_t k = n; k-- > 0; ) { // for k in {n-1,...,0}
            // try last step [k,n]
            double prev = run_p_var[k];
            double step = std::pow(std::abs(path[n] - path[k]), p);
            r = std::max(r, prev + step);
        }
        run_p_var[n] = r;
    }

    return std::pow(run_p_var.back(), 1./p);
}

std::vector<double> make_test_path(int steps) {
    std::vector<double> path(steps+1);
    for (size_t n=0; n<path.size(); n++) {
        path[n] = n % 3;
    }
    return path;
}

unsigned int random_seed(unsigned int & seed) {
    std::uniform_int_distribution<unsigned int> unif(0, std::numeric_limits<unsigned int>::max());
    std::random_device rnd_device;
    seed = unif(rnd_device);
    return seed;
}

std::vector<double> make_brownian_path(double sd, int steps ) {
    unsigned int seed;
    std::default_random_engine generator(random_seed(seed));
    std::normal_distribution<double> gauss_dist(0.0, sd);
    std::vector<double> increments(steps), path(1,0.);
    auto gauss_rv = [&](){ return gauss_dist(generator); }; // bind the generator
    std::generate(begin(increments), end(increments), gauss_rv);
    std::partial_sum(begin(increments), end(increments), back_inserter(path));
    return path;
}


int main() {
    using std::cout;

    int test_no = 0;

    // check a long periodic path 0,1,4,0,1,4,...
    {
        cout << "*** TEST " << ++test_no << " ***\n";
        size_t rep = 100000;
        std::vector<double> path(rep*4);
        for (size_t j=0; j<path.size(); j++) {
            switch(j % 4) {
                case 0:
                    path[j] = 0;
                    break;
                case 1:
                    path[j] = 1;
                    break;
                case 2:
                    path[j] = 1;
                    break;
                case 3:
                    path[j] = 4;
                    break;
            }
        }
        cout << "Long periodic path: 0,1,1,4 repeated " << rep << " times\n";
        for (double p = 1.0; p<3.01; p+=0.5) {
            double pv = p_var(path,p);
            double pv_ref = 4 * std::pow(2 * rep - 1, 1./p);
            cout << "  " << p << "-variation: " << pv
                << ", error: " << pv - pv_ref
                << "\n";
        }
    }

    // check a long monotone path from 0.0 to 1.0
    {
        cout << "\n*** TEST " << ++test_no << " ***\n";
        size_t steps = 10000;
        std::vector<double> path(steps+1);
        for (size_t j=0; j<path.size(); j++) {
            path[j] = double(j) / steps;
        }
        cout << "Long monotone path: "
            << path[0] << ", " << path[1] << ", " << path[2]
            << ", ... , " << path.back() << "\n";
        for (double p = 1.0; p<3.01; p+=0.5) {
            double pv = p_var(path,p);
            double pv_ref = 1.0;
            cout << "  " << p << "-variation: " << pv
                << ", error: " << pv - pv_ref
                << "\n";
        }
    }

    // check a simple path
    {
        cout << "\n*** TEST " << ++test_no << " ***\n";
        size_t steps = 4;
        std::vector<vec_Xd> path(steps + 1);
        for (size_t j=0; j<path.size(); j++) {
            for (size_t k=0; k<d; k++) {
                if (k==0 && (j%4==1 || j%4==2)) {
                    path[j][k] = 1;
                } else if (k==1 && (j%4==2 || j%4==3)) {
                    path[j][k] = 1;
                } else {
                    path[j][k] = 0;
                }
            }
        }
        cout << "Simple path:";
        for (size_t j=0; j<path.size(); j++) {
            cout << " (" << path[j][0];
            for (size_t k=1; k<d; k++) {
                cout << ", " << path[j][k];
            }
            cout << ")";
        }
        cout << "\n";
        for (double p = 1.0; p<3.01; p+=0.5) {
            double pv = p_var_Xd(path,p);
            double pv_ref = (p > 2) ? std::pow(2.0, 0.5 + 1./p) : std::pow(4, 1./p);
            cout << "  " << p << "-variation: " << pv
                << ", error: " << pv - pv_ref
                << "\n";
        }
    }

    // another simple path
    {
        cout << "\n*** TEST " << ++test_no << " ***\n";
        double p = 3;
        size_t steps = 3;
        std::vector<vec_Xd> path(steps + 1);
        for (size_t j=0; j<path.size(); j++) {
            for (size_t k=0; k<d; k++) {
                if (k==0 && (j%4==2 || j%4==3)) {
                    path[j][k] = 1;
                } else if (k==1 && (j%4==0 || j%4==1)) {
                    path[j][k] = 1;
                } else {
                    path[j][k] = 0;
                }
            }
        }
        cout << "Simple path:";
        for (size_t j=0; j<path.size(); j++) {
            cout << " (" << path[j][0];
            for (size_t k=1; k<d; k++) {
                cout << ", " << path[j][k];
            }
            cout << ")";
        }
        cout << "\n";
        double pv = p_var_Xd(path,p);
        double pv_ref = std::pow(2.0, 0.5);
        cout << "  " << p << "-variation: " << pv
            << ", error: " << pv - pv_ref
            << "\n";
    }

    // compare to reference on random Brownian paths
    {
        cout << "\n*** TEST " << ++test_no << " ***\n";
        double p = 3;
        size_t count = 10;
        size_t steps = 10000;
        double max_err = 0.0;

        for (size_t c=0; c < count; c++) {
            double sd = 1 / sqrt(double(steps));
            std::vector<double> path = make_brownian_path(sd, steps);
            double pv = p_var(path,p);
            double pv_ref = p_var_ref(path,p);
            double err = std::abs(pv_ref - pv);
            max_err = std::max(max_err, err);
        }
        cout << count << " random Brownian paths of length " << steps
            << " compared to reference\n";
        cout << "  max error: " << max_err << "\n";
    }

    // benchmark
    {
        cout << "\n*** TEST " << ++test_no << " (BENCHMARK) ***\n";
        double p = 3;
        for (size_t steps=1; steps <= 10000000; steps*=10) {
            double sd = 1 / sqrt(double(steps));
            std::vector<double> path = make_brownian_path(sd, steps);

            clock_t clock_begin = std::clock();
            double pv = p_var(path,p);
            clock_t clock_end = std::clock();

            double elapsed_secs = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
            cout << "  a random Brownian path of length " << steps
                << " has " << p << "-variation of " << pv << "\n"
                << "  which took " << elapsed_secs << " seconds to compute\n";
        }
    }
}
