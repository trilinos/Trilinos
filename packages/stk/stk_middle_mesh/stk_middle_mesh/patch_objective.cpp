#include "patch_objective.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

utils::Point PatchObjective::compute_inverse_parameterization(const utils::Point& pt)
{
  auto pt1 = m_patch.eval_point(pt.x, pt.y);
  auto pt2 = m_proj.project_back(pt1);

  return pt2;
}

utils::Point PatchObjective::compute_parameterization(const utils::Point& pt)
{

  return m_proj.project_forward(pt);
}

bool PatchObjective::has_invalid()
{
  auto vertIn = m_data->get_current_vert();
  auto ptIn   = compute_parameterization(vertIn->get_point_orig(0));

  utils::impl::Matrix<utils::Point> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  utils::impl::Mat2x2<double> w;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    mesh::impl::TriPts<double> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};
    w           = mesh::impl::DistortionMetric<double>::compute_w(ptsI);
    double detW = det2x2(w);
    if (detW < 0)
      return true;
  }

  return false;
}

void PatchObjective::set_parameterization(ActiveVertData& data)
{
  // compute the average normal vector (taking into account inverted elements)
  utils::Point nrm;
  for (int i = 0; i < data.get_num_elements(); ++i)
  {
    auto pts  = data.get_element_verts_orig(i);
    auto b1   = pts[1] - pts[0];
    auto b2   = pts[2] - pts[0];
    auto nrmI = cross(b1, b2);

    // if (computeTriArea(pts[0], pts[1], pts[2]) < 0)
    //   nrm_i = -nrm_i;

    nrm = nrm + nrmI / data.get_num_elements();
  }

  auto basis = utils::impl::compute_basis(nrm);
  m_proj     = utils::impl::ChangeOfBasis(basis);

  /*
    // arbitrarilly chose the first element to decide the projection
    // plane
    auto pts = data.getElementVerts(0);
    m_proj = getutils::impl::PlaneProjectionEnum(pts[0]->get_point_orig(0),
                                    pts[1]->get_point_orig(0),
                                    pts[2]->get_point_orig(0));

    for (int i=0; i < data.get_num_elements(); ++i)
    {
      auto pts = data.getElementVerts(i);
      std::cout << "for points " << std::endl;
      for (int j=0; j < 3; ++j)
        std::cout << pts[j]->get_point_orig(0) << std::endl;

      auto proj = getutils::impl::PlaneProjectionEnum(pts[0]->get_point_orig(0),
                                         pts[1]->get_point_orig(0),
                                         pts[2]->get_point_orig(0));
      std::cout << "best plane = " << proj << std::endl;



    }
  */
}

void PatchObjective::set_patch(ActiveVertData& data)
{
  /*
  const auto& verts = data.getUniqueVerts();
  std::vector<utils::Point> pts(verts.size());
  for (unsigned int i=0; i < verts.size(); ++i)
    pts[i] = m_proj.projectForward(verts[i]->get_point_orig(0));

  m_patch.constructPatch(pts);
  */

  // TODO: when using the original points, there is no need to recompute
  //       the taylor approximation each time, we can just keep the old result
  //       (although that would require the utils::TaylorPatch living on the ActiveVertData)
  const auto& pts = data.get_points_orig();
  std::vector<utils::Point> ptsProjected(pts.size());
  for (unsigned int i = 0; i < pts.size(); ++i)
  {
    ptsProjected[i] = m_proj.project_forward(pts[i]);
  }

  m_patch.construct_patch(ptsProjected);
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
