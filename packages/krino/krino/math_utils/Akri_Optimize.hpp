#ifndef KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_
#define KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_

#include <functional>
#include <tuple>

namespace krino {

std::tuple<double,double> line_search_armijo_1d(const std::function<double(double)> & fn, const double f0, const double dirDeriv0, const double xmin = 1.e-8);

template<typename VEC>
std::tuple<VEC,double> line_search_armijo_vector(const std::function<double(const VEC&)> & fn,
    const VEC& x0,
    const VEC& dir,
    const double f0,
    const VEC& gradf0,
    const double xmin = 1.e-8);

template<typename VEC>
VEC bfgs(const std::function<double(const VEC&)> & fn,
    const std::function<VEC(const VEC&)> & gradient,
    const VEC& x0,
    const double tol = 1e-6,
    const int maxIter = 50);

}

#endif /* KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_ */
