#ifndef CREATE_MESH_QUALITY_IMPROVER_H
#define CREATE_MESH_QUALITY_IMPROVER_H

#include "mesh_quality_improver.h"
#include "mesh_quality_improver_opts.h"

#include "optimization_step.h"
#include "optimization_step_hessian.h"
#include "optimization_step_laplacian.h"
#include "patch_distortion_objective.h"
#include "regularized_distortion_metric.h"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

template <typename Tfunc>
std::shared_ptr<MeshQualityImprover>
make_standard_improver(std::shared_ptr<Mesh> mesh, Tfunc filter,
                       const MeshQualityImproverOpts& opts = MeshQualityImproverOpts())
{
  // std::shared_ptr<RegularizedDistortionMetric> metric(delta);
  auto metric     = std::make_shared<RegularizedDistortionMetric>(opts.delta);
  auto qualityObj = std::make_shared<PatchDistortionObjective>(metric);
  auto qualityOpt = std::make_shared<opt::impl::OptimizationStepHessian>(qualityObj);

  return std::make_shared<MeshQualityImprover>(mesh, filter, opts.nlayers, qualityOpt, opts.maxDeltaX, opts.itermax);
}

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
