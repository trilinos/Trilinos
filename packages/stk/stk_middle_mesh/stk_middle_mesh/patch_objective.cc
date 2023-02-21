#include "patch_objective.h"
#include "distortion_metric.h"
#include "mesh_io.h"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

utils::Point PatchObjective::compute_inverse_parameterization(const utils::Point& pt)
{
  auto pt1 = m_patch.eval_point(pt.x, pt.y);
  auto pt2 = m_proj.project_back(pt1);
  /*
  if (output)
  {
    std::cout << "pt = " << pt2 << ", r = " << std::sqrt(dot(pt2, pt2)) << std::endl;

    double max_u=0.0, max_v=0.0, min_u=0.0, min_v=0.0;
    auto pt0 = m_proj.projectForward(m_data->getCurrentVert()->get_point_orig(0));
    std::cout << "patch is centered on " << m_data->getCurrentVert()->get_point_orig(0) << ", node " << pt0 <<
  std::endl; for (auto vert : m_data->getUniqueVerts())
    {
      auto pt = m_proj.projectForward(vert->get_point_orig(0));
      max_u = std::max(max_u, pt.x - pt0.x);
      max_v = std::max(max_v, pt.y - pt0.y);

      min_u = std::min(min_u, pt.x - pt0.x);
      min_v = std::min(min_v, pt.y - pt0.y);
    }

    std::cout << "patch was formed with u in range " << min_u << ", " << max_u
              << ", v range " << min_v << ", " << max_v << std::endl;

    std::cout << "check current vert locations" << std::endl;
    int idx = 0;
    for (auto vert: m_data->getUniqueVerts())
    {
      auto pt_xyz = vert->get_point_orig(0);
      auto pt_prime = m_proj.projectForward(pt_xyz);
      auto pt_prime_evaluated = m_patch.evalutils::Point(pt_prime.x, pt_prime.y, false);
      auto pt_xyz_evaluated = m_proj.projectBack(pt_prime_evaluated);
      std::cout << "patch vert " << idx << ", r_original = " << std::sqrt(dot(pt_xyz, pt_xyz))
                << ", r_evaluated = " << std::sqrt(dot(pt_xyz_evaluated, pt_xyz_evaluated)) << std::endl;
      idx++;
    }

    std::cout << "check original vert locations" << std::endl;
    idx = 0;
    for (auto pt_xyz: m_data->getPointsOrig())
    {
      auto pt_prime = m_proj.projectForward(pt_xyz);
      auto pt_prime_evaluated = m_patch.evalutils::Point(pt_prime.x, pt_prime.y, true);
      auto pt_xyz_evaluated = m_proj.projectBack(pt_prime_evaluated);
      std::cout << "patch vert " << idx << ", r_original = " << std::sqrt(dot(pt_xyz, pt_xyz))
                << ", r_evaluated = " << std::sqrt(dot(pt_xyz_evaluated, pt_xyz_evaluated)) << std::endl;
      idx++;
    }
  }
  */

  return pt2;
}

utils::Point PatchObjective::compute_parameterization(const utils::Point& pt)
{
  // if (m_data->getCurrentVert()->get_id() == 267)
  //   std::cout << "projecting point to u, v, pt = " << pt << ", r = " << std::sqrt(dot(pt, pt))<< std::endl;
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
    mesh::impl::TriPts ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};
    w           = mesh::impl::DistortionMetric::compute_w(ptsI);
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
    ptsProjected[i] = m_proj.project_forward(pts[i]);

  m_patch.construct_patch(ptsProjected);
}

void PatchObjective::compute_parameterization(const utils::Point& pt, utils::impl::Matrix<utils::Point>& triPts,
                                              std::vector<utils::impl::Mat2x2<double>>& wMats)
{
  triPts.resize(m_data->get_num_elements(), 3);
  wMats.resize(m_data->get_num_elements());

  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    auto pts     = m_data->get_element_verts(i);
    auto ptsOrig = m_data->get_element_verts_orig(i);
    for (int j = 0; j < 3; ++j)
    {
      triPts(i, j) = m_proj.project_forward(pts[j]->get_point_orig(0));
      ptsOrig[j]   = m_proj.project_forward(ptsOrig[j]);
    }
    triPts(i, m_data->get_current_vert_idx(i)) = pt;

    wMats[i] = mesh::impl::DistortionMetric::compute_w(ptsOrig);
    // W_mats[i] = m_data->getW(i, m_proj);
  }
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
