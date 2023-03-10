#include "stk_middle_mesh/predicates/point_classifier_normal_interpolation.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

PointRecordForTriangle PointClassifierNormalInterpolation::classify(mesh::MeshEntityPtr el, utils::Point pt,
                                                                    utils::Point normal)
{
  assert(el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);

  std::array<mesh::MeshEntityPtr, 3> verts;
  get_downward(el, 0, verts.data());

  utils::Point rhs = pt - verts[0]->get_point_orig(0);

  utils::Point dx1 = verts[1]->get_point_orig(0) - verts[0]->get_point_orig(0);
  utils::Point dx2 = verts[2]->get_point_orig(0) - verts[0]->get_point_orig(0);
  auto& mat        = m_mat3x3;
  mat(0, 0)        = dx1.x;
  mat(0, 1)        = dx2.x;
  mat(0, 2)        = -normal.x;
  mat(1, 0)        = dx1.y;
  mat(1, 1)        = dx2.y;
  mat(1, 2)        = -normal.y;
  mat(2, 0)        = dx1.z;
  mat(2, 1)        = dx2.z;
  mat(2, 2)        = -normal.z;

  std::array<int, 3> ipiv;
  solve_linear_system(mat, ipiv.data(), &(rhs.x));

  auto r = m_xiClassifier.classify_triangle(rhs);
  r.el   = el;

  return r;
}

PointRecordForTriangle PointClassifierNormalInterpolation::classify_reverse(mesh::MeshEntityPtr el, utils::Point pt,
                                                                            std::array<utils::Point, 4>& normals,
                                                                            bool logicalResultOnly)
{
  assert(el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);

  std::array<mesh::MeshEntityPtr, 3> verts;
  get_downward(el, 0, verts.data());
  std::array<utils::Point, 3> vertPts = {verts[0]->get_point_orig(0), verts[1]->get_point_orig(0),
                                         verts[2]->get_point_orig(0)};

  ReverseClassificationGammaSolver solver(vertPts, pt, normals, m_tolerances.newtonTol, m_tolerances.newtonItermax);
  bool didConverge        = false;
  std::array<double, 3> x = solver.solve(didConverge);

  if (!logicalResultOnly)
  {
    if (!didConverge)
      x = solve_reverse3x3(verts, pt, normals);

    check_solution_reverse(verts, pt, normals, x);
  }

  auto r = m_xiClassifier.classify_triangle(utils::Point(x[0], x[1]));
  r.el   = el;
  return r;
}

std::array<double, 3>
PointClassifierNormalInterpolation::solve_reverse3x3(const std::array<mesh::MeshEntityPtr, 3>& verts,
                                                     const utils::Point& pt, const std::array<utils::Point, 4>& normals)
{
  utils::Point rhs = -pt + verts[0]->get_point_orig(0);
  utils::Point dx1 = verts[1]->get_point_orig(0) - verts[0]->get_point_orig(0);
  utils::Point dx2 = verts[2]->get_point_orig(0) - verts[0]->get_point_orig(0);
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  auto func = [&](const std::vector<double>& x, std::vector<double>& f) {
    double lambda1 = x[0], lambda2 = x[1], gamma = x[2];
    for (int i = 0; i < 3; ++i)
      f[i] =
          dx1[i] * lambda1 + dx2[i] * lambda2 - gamma * (normals[0][i] + dn1[i] * lambda1 + dn2[i] * lambda2) + rhs[i];
  };

  auto jac = [&](const std::vector<double>& x, utils::impl::Matrix<double>& mat) {
    double lambda1 = x[0], lambda2 = x[1], gamma = x[2];
    for (int i = 0; i < 3; ++i)
    {
      mat(i, 0) = dx1[i] - gamma * dn1[i];
      mat(i, 1) = dx2[i] - gamma * dn2[i];
      mat(i, 2) = -(normals[i][i] + dn1[i] * lambda1 + dn2[i] * lambda2);
    }
  };

  m_x0[0] = 0.5;
  m_x0[1] = 0.5;
  m_x0[2] = 0; // TODO: maybe there is a better way to do this:
               //       compute vector from pt to centroid of
               //       element, project onto normal vector
               //       at the centroid, compute gamma

  int retVal = m_newton.solve(func, jac, m_x0);
  if (retVal > 0)
    throw std::runtime_error("Newton failed to converge");

  return {m_x0[0], m_x0[1], m_x0[2]};
}

void PointClassifierNormalInterpolation::check_solution_reverse(const std::array<mesh::MeshEntityPtr, 3>& verts,
                                                                const utils::Point& pt,
                                                                const std::array<utils::Point, 4>& normals,
                                                                const std::array<double, 3>& solution)
{
  double lambda1 = solution[0];
  double lambda2 = solution[1];
  double lambda0 = 1 - lambda1 - lambda2;
  double gamma   = solution[2];

  // std::cout << "lambdas = " << utils::Point(lambda0, lambda1, lambda2) << ", gamma = " << gamma << std::endl;

  utils::Point p1 = pt + gamma * (lambda0 * normals[0] + lambda1 * normals[1] + lambda2 * normals[2]);
  utils::Point p2 = lambda0 * verts[0]->get_point_orig(0) + lambda1 * verts[1]->get_point_orig(0) +
                    lambda2 * verts[2]->get_point_orig(0);
  utils::Point disp = p2 - p1;
  double dist       = std::sqrt(dot(disp, disp));
  // std::cout << "dist = " << dist << std::endl;
  if (dist > 1e-11 || std::isnan(dist))
    throw std::runtime_error("Newton produced incorrect solution");
}

} // namespace impl
} // namespace predicates
} // namespace middle_mesh
} // namespace stk
