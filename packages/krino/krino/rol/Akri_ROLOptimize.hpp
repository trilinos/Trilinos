#ifndef KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_
#define KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_

#include <Akri_OptimizeTypes.hpp>

namespace krino {

template<typename ObjectiveType, typename VEC>
void rol_optimize(const ObjectiveType & objFn,
    VEC& x,
    const double xTol = OptimizationDefaults::stepTol,
    const double gradTol = OptimizationDefaults::gradTol,
    const double objectiveScale = OptimizationDefaults::objectiveScale,
    const int maxIter = OptimizationDefaults::maxIter,
    const bool doOutput = OptimizationDefaults::doOutput);

template<typename OptimizerFunction>
OptimizerFunction
build_rol_optimizer_with_parameters(const double xTol = OptimizationDefaults::stepTol,
    const double gradTol = OptimizationDefaults::gradTol,
    const double objectiveScale = OptimizationDefaults::objectiveScale,
    const unsigned maxIter = OptimizationDefaults::maxIter,
    const bool doOutput = OptimizationDefaults::doOutput)
{
  using Objective = typename OptimizerFunctionTraits<OptimizerFunction>::objective_type;
  using VEC = typename OptimizerFunctionTraits<OptimizerFunction>::vector_type;

  auto fn = [=](const Objective &objFn, VEC &x) { rol_optimize(objFn, x, xTol, gradTol, objectiveScale, maxIter, doOutput); };
  return fn;
}

}

#endif /* KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_ */
