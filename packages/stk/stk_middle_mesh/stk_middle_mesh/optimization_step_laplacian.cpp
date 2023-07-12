#include "optimization_step_laplacian.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

utils::Point OptimizationStepLaplacian::compute_updated_point(ActiveVertData& active)
{
  m_obj.set_active_patch(active);
  const std::vector<utils::Point>& verts = active.get_unique_verts();
  // auto& pts_orig = active.getPointsOrig();

  utils::Point ptOut;

  // auto pt_in = m_obj.computeParameterization(vert->get_point_orig(0));
  // utils::Point pt_in0 = m_obj.computeParameterization(pts_orig[0]);
  double weightSum = 0;
  for (unsigned int i = 1; i < verts.size(); ++i)
  {
    // computes weights based on the original vertex position, but apply them to the
    // current position
    // hopefully this makes the new and old meshes similar
    auto ptI = m_obj.compute_parameterization(verts[i]);
    // auto pt_i0 = m_obj.computeParameterization(pts_orig[i]);
    // auto disp = pt_i0 - pt_in0;
    // auto weight_i = 1.0/std::sqrt(dot(disp, disp));
    double weightI = 1;

    ptOut += weightI * ptI;
    weightSum += weightI;
  }

  return ptOut / weightSum;
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
