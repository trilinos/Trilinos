#ifndef PATCH_DISTORTION_OBJECTIVE_H
#define PATCH_DISTORTION_OBJECTIVE_H

#include <vector>

#include "active_vert_data.hpp"
#include "distortion_metric.hpp"
#include "matrix.hpp"
#include "patch_objective.hpp"
#include "taylor_patch.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class PatchDistortionObjective : public opt::impl::PatchObjective
{
  public:
    explicit PatchDistortionObjective(std::shared_ptr<DistortionMetric<double>> metric)
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

    template <typename T>
    T compute_qualityT(std::shared_ptr<DistortionMetric<T>> metric, const utils::PointT<T>& ptIn);

    template <typename T>
    utils::PointT<T> compute_quality_revT(std::shared_ptr<DistortionMetric<T>> metric, const utils::PointT<T>& ptIn, T q3Bar = 1);

    template <typename T>
    utils::impl::Mat2x2<T> compute_hessianT(std::shared_ptr<DistortionMetric<T>> metric, const utils::PointT<T>& ptIn);

    void set_min_denominator(const utils::Point& ptIn);

    std::shared_ptr<DistortionMetric<double>> m_metric;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
