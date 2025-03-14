#include "patch_energy_objective.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

double PatchEnergyObjective::compute_quality(const utils::Point& ptIn)
{
  double q          = 0;
  const std::vector<utils::Point>& verts = m_data->get_unique_verts();
  for (unsigned int i = 1; i < verts.size(); ++i)
  {
    auto pt = compute_parameterization(verts[i]);
    // auto pt = applyutils::impl::PlaneProjection(m_proj, verts[i]->get_point_orig(0));

    double dx = ptIn.x - pt.x;
    double dy = ptIn.y - pt.y;
    q += dx * dx + dy * dy;
  }

  return q;
}

utils::Point PatchEnergyObjective::compute_quality_rev(const utils::Point& ptIn, double qBar)
{
  // double q          = 0;
  const std::vector<utils::Point>& verts = m_data->get_unique_verts();
  utils::Point ptInBar;
  for (unsigned int i = 1; i < verts.size(); ++i)
  {
    auto pt = compute_parameterization(verts[i]);
    // auto pt = applyutils::impl::PlaneProjection(m_proj, verts[i]->get_point_orig(0));

    double dx = ptIn.x - pt.x;
    double dy = ptIn.y - pt.y;
    // q += dx * dx + dy * dy;

    double dxBar = 2 * dx * qBar;
    double dyBar = 2 * dy * qBar;

    ptInBar.x += dxBar;
    ptInBar.y += dyBar;
  }

  return ptInBar;
}

utils::impl::Mat2x2<double> PatchEnergyObjective::compute_hessian(const utils::Point& /*ptIn*/)
{
  auto npts = m_data->get_num_verts() - 1;
  utils::impl::Mat2x2<double> hess;
  hess(0, 0) = 2 * npts;
  hess(1, 1) = 2 * npts;

  return hess;
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
