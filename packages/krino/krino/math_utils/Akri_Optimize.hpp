#ifndef KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_
#define KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_

#include <functional>
#include <stk_math/StkVector.hpp>
#include <Akri_DistributedVector.hpp>

namespace krino {

using Vector3dObjectiveFn = std::function<double(const stk::math::Vector3d&)>;
using Vector3dObjectiveSensFn = std::function<void(const stk::math::Vector3d&, stk::math::Vector3d&)>;
using DistributedVectorObjectiveFn = std::function<double(const DistributedVector&)>;
using DistributedVectorObjectiveSensFn = std::function<void(const DistributedVector&, DistributedVector&)>;

template<typename VEC>
void steepest_descent(const std::function<double(const VEC&)> & calc_objective,
    const std::function<void(const VEC&, VEC&)> & fill_gradient,
    VEC& x,
    const double xTol = 1e-6,
    const double gradTol = 1e-6,
    const unsigned maxIter = 200);

template<typename VEC>
void lbfgs(const std::function<double(const VEC&)> & calc_objective,
    const std::function<void(const VEC&, VEC&)> & fill_gradient,
    VEC& x,
    const double xTol = 1e-6,
    const double gradTol = 1e-6,
    const unsigned maxIter = 200,
    const unsigned maxLevels = 10);
}

#endif /* KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_ */
