#include "gtest/gtest.h"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/quad_metrics.hpp"
#include "stk_middle_mesh/utils.hpp"

namespace {

using namespace stk::middle_mesh;

template <typename T>
std::array<double, 4> makeFunction(mesh::MeshEntityPtr quad, T func)
{
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  mesh::get_downward(quad, 0, verts.data());

  std::array<double, 4> vals;
  for (int i=0; i < 4; ++i)
  {
    utils::Point pt = verts[i]->get_point_orig(0);
    vals[i] = func(pt);
  }

  return vals;
}

const double one_over_rad3 = 1.0/std::sqrt(3.0);

const std::vector<utils::Point> quadPts = { mesh::convert_xi_coords_from_range(-1, 1, {-one_over_rad3, -one_over_rad3, 0}),
                                             mesh::convert_xi_coords_from_range(-1, 1, { one_over_rad3, -one_over_rad3, 0}),
                                             mesh::convert_xi_coords_from_range(-1, 1, { one_over_rad3,  one_over_rad3, 0}),
                                             mesh::convert_xi_coords_from_range(-1, 1, {-one_over_rad3,  one_over_rad3, 0})};
const std::vector<double> quadWeights = {0.25, 0.25, 0.25, 0.25};
}

TEST(QuadMetrics, Standard)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto mesh = mesh::make_empty_mesh();
  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(2, 0, 0);
  auto v3 = mesh->create_vertex(2, 2, 0);
  auto v4 = mesh->create_vertex(0, 2, 0);
  auto quad = mesh->create_quad_from_verts(v1, v2, v3, v4);

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y + 3*pt.z; };
  auto vals = makeFunction(quad, f);

  std::vector<double> detJ;
  mesh::impl::QuadMetrics metrics;
  metrics.compute_det_jacobian(quad, quadPts, detJ);

  double val = 0;
  for (size_t i=0; i < quadPts.size(); ++i)
  {
    val += mesh::interpolate(quad->get_type(), vals.data(), quadPts[i]) * quadWeights[i] * detJ[i];
  }

  EXPECT_NEAR(val, 12.0, 1e-13);
}

TEST(QuadMetrics, Skewed)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto mesh = mesh::make_empty_mesh();
  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(2, 0, 0);
  auto v3 = mesh->create_vertex(2, 4, 0);
  auto v4 = mesh->create_vertex(0, 2, 0);
  auto quad = mesh->create_quad_from_verts(v1, v2, v3, v4);

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y + 3*pt.z; };
  auto vals = makeFunction(quad, f);

  std::vector<double> detJ;
  mesh::impl::QuadMetrics metrics;
  metrics.compute_det_jacobian(quad, quadPts, detJ);

  double val = 0;
  for (size_t i=0; i < quadPts.size(); ++i)
  {
    val += mesh::interpolate(quad->get_type(), vals.data(), quadPts[i]) * quadWeights[i] * detJ[i];
  }

  std::cout << "val = " << val << std::endl;

  EXPECT_NEAR(val, 12.0 + 8*5.0/3, 1e-13);
}