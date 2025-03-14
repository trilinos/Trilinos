#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/patch_distortion_objective.hpp"
#include "stk_middle_mesh/patch_energy_objective.hpp"
#include "stk_middle_mesh/regularized_distortion_metric.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

namespace {

MeshEntityPtr find_closest_vert(std::shared_ptr<Mesh> mesh, const utils::Point& pt)
{
  double minDist        = std::numeric_limits<double>::max();
  MeshEntityPtr minVert = nullptr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      auto pt2    = vert->get_point_orig(0);
      auto disp   = pt - pt2;
      double dist = dot(disp, disp);
      if (dist < minDist)
      {
        minDist = dist;
        minVert = vert;
      }
    }

  return minVert;
}

void expect_float_eq2(const utils::Point& pt1, const utils::Point& pt2)
{
  EXPECT_FLOAT_EQ(pt1.x, pt2.x);
  EXPECT_FLOAT_EQ(pt1.y, pt2.y);
}

void test_gradient(opt::impl::PatchObjective& obj, opt::impl::ActiveVertData& active, utils::Point& pt,
                   const double tol)
{
  obj.set_active_patch(active);
  utils::Point grad = obj.compute_quality_rev(pt);

  // check finite difference
  double eps = 1e-7;
  double q0  = obj.compute_quality(pt);
  // std::cout << "q0 = " << q0 << std::endl;

  pt.x += eps;
  double q1 = obj.compute_quality(pt);
  // std::cout << "q1 = " << q1 << std::endl;
  double dqDx1 = (q1 - q0) / eps;
  pt.x -= eps;

  pt.y += eps;
  double q2    = obj.compute_quality(pt);
  double dqDx2 = (q2 - q0) / eps;
  pt.y -= eps;

  // std::cout << "deriv_fd = " << utils::Point(dq_dx1, dq_dx2) << std::endl;
  // std::cout << "deriv    = " << pt << std::endl;
  EXPECT_NEAR(dqDx1, grad.x, tol);
  EXPECT_NEAR(dqDx2, grad.y, tol);
}

void test_hessian(opt::impl::PatchObjective& obj, opt::impl::ActiveVertData& active, utils::Point& pt, const double /*tol*/)

{
  obj.set_active_patch(active);
  auto hessian = obj.compute_hessian(pt);

  utils::impl::Mat2x2<double> hessian2;
  double eps = 1e-7;
  auto col0  = obj.compute_quality_rev(pt, 1);

  pt.x += eps;
  auto col1 = obj.compute_quality_rev(pt, 1);
  pt.x -= eps;

  pt.y += eps;
  auto col2 = obj.compute_quality_rev(pt, 1);
  pt.y -= eps;

  hessian2(0, 0) = (col1.x - col0.x) / eps;
  hessian2(1, 0) = (col1.y - col0.y) / eps;
  hessian2(0, 1) = (col2.x - col0.x) / eps;
  hessian2(1, 1) = (col2.y - col0.y) / eps;

  std::cout << "hessian = \n" << hessian << std::endl;
  // std::cout << "hessian2 = \n" << hessian2 << std::endl;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      // assert(std::abs(hessian(i, j) - hessian2(i, j)) < 5e-5);
      EXPECT_NEAR(hessian(i, j), hessian2(i, j), 5e-5);
}

} // namespace

TEST(PatchObjective, computeDeriv)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 2;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);
  MeshEntityPtr vert         = find_closest_vert(mesh, utils::Point(0.5, 0.5));
  opt::impl::ActiveVertData active(mesh, vert);

  auto metric = std::make_shared<RegularizedDistortionMetric<double>>();
  PatchDistortionObjective obj(metric);
  opt::impl::PatchEnergyObjective obj2;
  // evaluate somewhere other than the initial position
  // (where the W matrix was calculated -> quality = 1 -> gradient = 0)
  // vert->set_point_orig(0, utils::Point(0.55, 0.65));
  utils::Point pt = utils::Point(0.55, 0.65);

  test_gradient(obj, active, pt, 2e-6);
  test_gradient(obj2, active, pt, 2e-6);
}

TEST(PatchObjective, computeHessian)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 2;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);
  MeshEntityPtr vert         = find_closest_vert(mesh, utils::Point(0.5, 0.5));
  opt::impl::ActiveVertData active(mesh, vert);

  auto metric = std::make_shared<RegularizedDistortionMetric<double>>();
  PatchDistortionObjective obj(metric);

  opt::impl::PatchEnergyObjective obj2;
  // evaluate somewhere other than the initial position
  // (where the W matrix was calculated -> quality = 1 -> gradient = 0)
  utils::Point pt = utils::Point(0.55, 0.65);

  test_hessian(obj, active, pt, 5e-5);
  test_hessian(obj2, active, pt, 5e-5);
}

TEST(PatchObjective, Parameterization)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 2;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);
  MeshEntityPtr vert         = find_closest_vert(mesh, utils::Point(0.5, 0.5));
  opt::impl::ActiveVertData active(mesh, vert);

  auto metric = std::make_shared<RegularizedDistortionMetric<double>>();
  // auto f = [](MeshEntityPtr) { return true; };
  PatchDistortionObjective obj(metric);
  obj.set_active_patch(active);

  expect_float_eq2(obj.compute_parameterization(utils::Point(0.5, 0.5)), utils::Point(0.5, 0.5));
  expect_float_eq2(obj.compute_parameterization(utils::Point(0.6, 0.7)), utils::Point(0.6, 0.7));

  auto pt1 = obj.compute_inverse_parameterization(utils::Point(0.5, 0.5));
  expect_float_eq2(pt1, utils::Point(0.5, 0.5));
  expect_float_eq2(obj.compute_parameterization(utils::Point(0.6, 0.7)), utils::Point(0.6, 0.7));
}

// TODO: test compute_parameterization and compute_inverse_parameterization

} // namespace impl
} // namespace middle_mesh
} // namespace stk
