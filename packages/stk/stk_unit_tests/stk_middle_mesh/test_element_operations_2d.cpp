#include "stk_middle_mesh/element_operations_2d.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class ElementOperations2DTester : public ::testing::Test
{
  protected:
    ElementOperations2DTester()
      : proj(utils::Point(0, 0, 0), utils::Point(1, 0, 0), utils::Point(0, 1, 0))
    {
      elemOps.set_projection(&proj);
    }

    void setup(const utils::Point& pt1, const utils::Point& pt2, const utils::Point& pt3, const utils::Point& pt4)
    {
      mesh = mesh::make_empty_mesh();

      v1 = mesh->create_vertex(pt1);
      v2 = mesh->create_vertex(pt2);
      v3 = mesh->create_vertex(pt3);
      v4 = mesh->create_vertex(pt4);

      el  = mesh->create_quad_from_verts(v1, v2, v3, v4);
      tri = mesh->create_triangle_from_verts(v1, v2, v3);
    }

    std::shared_ptr<mesh::Mesh> mesh;
    mesh::MeshEntityPtr v1;
    mesh::MeshEntityPtr v2;
    mesh::MeshEntityPtr v3;
    mesh::MeshEntityPtr v4;
    mesh::MeshEntityPtr el;
    mesh::MeshEntityPtr tri;
    mesh::impl::ElementOperations2D elemOps;
    utils::impl::Projection proj;
};

} // namespace

TEST_F(ElementOperations2DTester, compute_edge_coords)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  setup({0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0});

  mesh::MeshEntityPtr edges[mesh::MAX_DOWN];
  get_downward(el, 1, edges);

  // compute_edge_coords
  utils::Point pt = elemOps.compute_edge_coords(edges[0], 0);
  EXPECT_FLOAT_EQ(pt.get_x(), 0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = elemOps.compute_edge_coords(edges[0], 0.5);
  EXPECT_FLOAT_EQ(pt.get_x(), 1);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = elemOps.compute_edge_coords(edges[0], 1);
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);
}

TEST_F(ElementOperations2DTester, compute_quad_xi_coords)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point pt;
  setup({0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0});

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(0, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(2, 2));
  EXPECT_FLOAT_EQ(pt.get_x(), 1);
  EXPECT_FLOAT_EQ(pt.get_y(), 1);

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(1, 1));
  EXPECT_FLOAT_EQ(pt.get_x(), 0.5);
  EXPECT_FLOAT_EQ(pt.get_y(), 0.5);
}

TEST_F(ElementOperations2DTester, compute_quad_xi_coordsShifted)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point pt;
  setup({1, 1, 0}, {2, 1, 0}, {2, 3, 0}, {1, 2, 0});

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(1, 1));
  EXPECT_NEAR(pt.get_x(), 0, 1e-13);
  EXPECT_NEAR(pt.get_y(), 0, 1e-13);

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(2, 1));
  EXPECT_NEAR(pt.get_x(), 1, 1e-13);
  EXPECT_NEAR(pt.get_y(), 0, 1e-13);

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(2, 3));
  EXPECT_NEAR(pt.get_x(), 1, 1e-13);
  EXPECT_NEAR(pt.get_y(), 1, 1e-13);

  pt = elemOps.compute_quad_xi_coords(el, utils::Point(1, 2));
  EXPECT_NEAR(pt.get_x(), 0, 1e-13);
  EXPECT_NEAR(pt.get_y(), 1, 1e-13);
}

TEST_F(ElementOperations2DTester, compute_tri_xi_coords)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point pt;
  setup({0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0});

  // Triangle xi coordinates
  pt = elemOps.compute_tri_xi_coords(tri, v1->get_point_orig(0));
  EXPECT_FLOAT_EQ(pt.get_x(), 0.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0.0);

  pt = elemOps.compute_tri_xi_coords(tri, v2->get_point_orig(0));
  EXPECT_FLOAT_EQ(pt.get_x(), 1.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0.0);

  pt = elemOps.compute_tri_xi_coords(tri, v3->get_point_orig(0));
  EXPECT_FLOAT_EQ(pt.get_x(), 0.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 1.0);

  pt = elemOps.compute_tri_xi_coords(tri, utils::Point(4.0 / 3.0, 2.0 / 3.0));
  EXPECT_FLOAT_EQ(pt.get_x(), 1.0 / 3.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 1.0 / 3.0);
}

TEST_F(ElementOperations2DTester, Area)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  setup({0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0});

  EXPECT_FLOAT_EQ(elemOps.compute_area(el), 4.0);
  EXPECT_FLOAT_EQ(elemOps.compute_area(tri), 2.0);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
