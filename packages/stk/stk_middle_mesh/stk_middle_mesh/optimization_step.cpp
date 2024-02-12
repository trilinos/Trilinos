#include "optimization_step.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

double OptimizationStep::improve_quality_single(ActiveVertData& active)
{
  PatchObjective& obj = get_objective();
  obj.set_active_patch(active);

  mesh::MeshEntityPtr vert = active.get_current_vert();
  auto ptPrime0            = obj.compute_parameterization(vert->get_point_orig(0));

  // take step
  utils::Point ptPrime = compute_updated_point(active);

  // update vertex
  auto pt = obj.compute_inverse_parameterization(ptPrime);

  vert->set_point_orig(0, pt);

  // update observables
  auto deltaXPt    = ptPrime - ptPrime0;
  //std::cout << "vert " << active.get_current_vert()->get_id() << " delta_x = " << deltaXPt << std::endl;
  deltaXPt.z       = 0;
  double deltaXMag = std::sqrt(dot(deltaXPt, deltaXPt));

  return deltaXMag;
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
