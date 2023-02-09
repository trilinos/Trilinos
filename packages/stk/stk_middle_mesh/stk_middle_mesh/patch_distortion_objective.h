#ifndef PATCH_DISTORTION_OBJECTIVE_H
#define PATCH_DISTORTION_OBJECTIVE_H

#include <vector>

#include "active_vert_data.h"
#include "distortion_metric.h"
#include "matrix.h"
#include "patch_objective.h"
#include "taylor_patch.h"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class PatchDistortionObjective : public opt::impl::PatchObjective
{
  public:
    explicit PatchDistortionObjective(std::shared_ptr<DistortionMetric> metric)
      : m_metric(metric)
    {}

    // computes quality when the active vert is at the given (parametric)
    // point
    double compute_quality(const utils::Point& ptIn) override;

    // computes derivative of computeQuality wrt the active vert only
    // Computes derivative wrt the parameterized coords
    // pt is the current position of the active vert
    utils::Point compute_quality_rev(const utils::Point& ptIn, double q3Bar = 1) override;

    // computes the Hessian of the quality metric wrt the parameterized
    // coordinates of the active vertex
    // pt_in is the current position of the active vert
    utils::impl::Mat2x2<double> compute_hessian(const utils::Point& ptIn) override;

  private:
    std::shared_ptr<DistortionMetric> m_metric;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
