#ifndef KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_
#define KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_

#include <functional>
#include <stk_math/StkVector.hpp>

namespace krino {

template<typename Objective, typename VEC>
void steepest_descent(const Objective & objFn,
    VEC& x,
    const double xTol = 1e-6,
    const double gradTol = 1e-6,
    const double objectiveScale = 1.0,
    const unsigned maxIter = 200,
    const bool doOutput = false);

template<typename Objective, typename VEC>
void lbfgs(const Objective & objFn,
    VEC& x,
    const double xTol = 1e-6,
    const double gradTol = 1e-6,
    const double objectiveScale = 1.0,
    const unsigned maxIter = 200,
    const unsigned maxLevels = 10,
    const bool doOutput = false);

template<typename Objective, typename VEC>
std::function<void(const Objective &objFn, VEC &x)>
build_steepest_descent_optimizer_with_parameters(const double xTol = 1e-6,
    const double gradTol = 1e-6,
    const double objectiveScale = 1.0,
    const unsigned maxIter = 200,
    const bool doOutput = false)
{
  auto fn = [=](const Objective &objFn, VEC &x) { steepest_descent(objFn, x, xTol, gradTol, objectiveScale, maxIter, doOutput); };
  return fn;
}

template<typename Objective, typename VEC>
std::function<void(const Objective &objFn, VEC &x)>
build_lbgfs_optimizer_with_parameters(const double xTol = 1e-6,
    const double gradTol = 1e-6,
    const double objectiveScale = 1.0,
    const unsigned maxIter = 200,
    const unsigned maxLevels = 10,
    const bool doOutput = false)
{
  auto fn = [=](const Objective &objFn, VEC &x) { lbfgs(objFn, x, xTol, gradTol, objectiveScale, maxIter, maxLevels, doOutput); };
  return fn;
}

}

#endif /* KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_ */
