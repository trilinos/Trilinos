#ifndef KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_
#define KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_

#include <functional>
#include <tuple>
#include <Akri_OptimizeTypes.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_ObjectiveInterface.hpp>

namespace krino {

template<typename ObjectiveType, typename VEC>
void rol_optimize(const ObjectiveType & objFn,
    VEC& x,
    const double xTol = 1.e-6,
    const double gradTol = 1e-6,
    const double objectiveScale = 1.0,
    const int maxIter = 200,
    const bool doOutput = false);

template<typename Objective, typename VEC>
std::function<void(const Objective &objFn, VEC &x)>
build_rol_optimizer_with_parameters(const double xTol = 1.e-6,
    const double gradTol = 1e-6,
    const double objectiveScale = 1.0,
    const unsigned maxIter = 200,
    const bool doOutput = false)
{
  auto fn = [=](const Objective &objFn, VEC &x) { rol_optimize(objFn, x, xTol, gradTol, objectiveScale, maxIter, doOutput); };
  return fn;
}

}

#endif /* KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_ */
