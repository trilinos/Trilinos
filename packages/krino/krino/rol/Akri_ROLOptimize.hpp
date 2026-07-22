#ifndef KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_
#define KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_

#include <functional>
#include <tuple>

namespace krino {

template<typename VEC>
void rol_optimize(const std::function<double(const VEC&)> & calc_objective,
    const std::function<void(const VEC&, VEC&)> & fill_gradient,
    VEC& x0,
    const double xTol = 1.e-6,
    const double gradTol = 1e-6,
    const int maxIter = 50);

}

#endif /* KRINO_KRINO_ROL_AKRI_ROLOPTIMIZE_HPP_ */
