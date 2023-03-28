#include "gtest/gtest.h"
#include "stk_middle_mesh/quad_point_finder.hpp"

using namespace stk::middle_mesh;

TEST(QuadPointFinder, Exactness)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(2, 0, 0);
  auto v3 = mesh->create_vertex(2, 4, 0);
  auto v4 = mesh->create_vertex(0, 2, 0);
  auto quad = mesh->create_quad_from_verts(v1, v2, v3, v4);

  mesh::impl::QuadPointFinder finder;
  {
    utils::Point ptXyz(0.25, 0.25), ptXiGuess(0, 0, 0);
    utils::Point ptXi = finder.compute_nearest_point_on_quad(quad, ptXyz, ptXiGuess);
    utils::Point ptXyzFound = mesh::compute_coords_from_xi_3d(quad, ptXi);

    for (int i=0; i < 3; ++i)
      EXPECT_NEAR(ptXyz[i], ptXyzFound[i], 1e-13);
  }
}