#ifndef KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_
#define KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_

#include <functional>
#include <stk_math/StkVector.hpp>
#include <Akri_OptimizeTypes.hpp>

namespace krino {

template<typename Objective, typename VEC>
void steepest_descent(const Objective & objFn,
    VEC& x,
    const double xTol = OptimizationDefaults::stepTol,
    const double gradTol = OptimizationDefaults::gradTol,
    const double objectiveScale = OptimizationDefaults::objectiveScale,
    const unsigned maxIter = OptimizationDefaults::maxIter,
    const bool doOutput = OptimizationDefaults::doOutput);

template<typename Objective, typename VEC>
void lbfgs(const Objective & objFn,
    VEC& x,
    const double xTol = OptimizationDefaults::stepTol,
    const double gradTol = OptimizationDefaults::gradTol,
    const double objectiveScale = OptimizationDefaults::objectiveScale,
    const unsigned maxIter = OptimizationDefaults::maxIter,
    const unsigned maxLevels = OptimizationDefaults::maxLevelsBFGS,
    const bool doOutput = OptimizationDefaults::doOutput);

template<typename OptimizerFunction>
OptimizerFunction
build_steepest_descent_optimizer_with_parameters(const double xTol = OptimizationDefaults::stepTol,
    const double gradTol = OptimizationDefaults::gradTol,
    const double objectiveScale = OptimizationDefaults::objectiveScale,
    const unsigned maxIter = OptimizationDefaults::maxIter,
    const bool doOutput = OptimizationDefaults::doOutput)
{
  using Objective = typename OptimizerFunctionTraits<OptimizerFunction>::objective_type;
  using VEC = typename OptimizerFunctionTraits<OptimizerFunction>::vector_type;

  auto fn = [=](const Objective &objFn, VEC &x) { steepest_descent(objFn, x, xTol, gradTol, objectiveScale, maxIter, doOutput); };
  return fn;
}

template<typename OptimizerFunction>
OptimizerFunction
build_lbgfs_optimizer_with_parameters(const double xTol = OptimizationDefaults::stepTol,
    const double gradTol = OptimizationDefaults::gradTol,
    const double objectiveScale = OptimizationDefaults::objectiveScale,
    const unsigned maxIter = OptimizationDefaults::maxIter,
    const unsigned maxLevels = OptimizationDefaults::maxLevelsBFGS,
    const bool doOutput = OptimizationDefaults::doOutput)
{
  using Objective = typename OptimizerFunctionTraits<OptimizerFunction>::objective_type;
  using VEC = typename OptimizerFunctionTraits<OptimizerFunction>::vector_type;

  auto fn = [=](const Objective &objFn, VEC &x) { lbfgs(objFn, x, xTol, gradTol, objectiveScale, maxIter, maxLevels, doOutput); };
  return fn;
}

}

#endif /* KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZE_HPP_ */
