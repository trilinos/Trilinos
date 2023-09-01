#ifndef PATCH_ENERGY_OBJECTIVE_H
#define PATCH_ENERGY_OBJECTIVE_H

#include "patch_objective.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

// objective function for point i = \sum_j (x_i - x_j)^2 + (y_i - y_j)^2

class PatchEnergyObjective : public PatchObjective
{
  public:
    double compute_quality(const utils::Point& ptIn) override;

    utils::Point compute_quality_rev(const utils::Point& ptIn, double q3Bar = 1) override;

    utils::impl::Mat2x2<double> compute_hessian(const utils::Point& ptIn) override;
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
