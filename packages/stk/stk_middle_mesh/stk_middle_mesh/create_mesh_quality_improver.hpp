#ifndef CREATE_MESH_QUALITY_IMPROVER_H
#define CREATE_MESH_QUALITY_IMPROVER_H

#include "mesh_quality_improver.hpp"
#include "mesh_quality_improver_opts.hpp"

#include "optimization_step.hpp"
#include "optimization_step_hessian.hpp"
#include "optimization_step_laplacian.hpp"
#include "patch_distortion_objective.hpp"
#include "regularized_distortion_metric.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

template <typename Tfunc>
std::shared_ptr<MeshQualityImprover>
make_standard_improver(std::shared_ptr<Mesh> mesh, Tfunc filter,
                       const MeshQualityImproverOpts& opts = MeshQualityImproverOpts())
{
  auto metric     = std::make_shared<RegularizedDistortionMetric<double>>(opts.delta);
  auto qualityObj = std::make_shared<PatchDistortionObjective>(metric);
  auto qualityOpt = std::make_shared<opt::impl::OptimizationStepHessian>(qualityObj);

  return std::make_shared<MeshQualityImprover>(mesh, filter, opts.nlayers, qualityOpt, opts.maxDeltaX, opts.itermax, opts.verboseOutput);
}

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
