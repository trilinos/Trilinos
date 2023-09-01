#include "gtest/gtest.h"

#include "stk_middle_mesh/predicates/intersection_common.hpp"
#include "stk_middle_mesh/predicates/triangle_coord_utils.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

namespace {

using utils::Point;

class TriangleCoordUtilsTester : public ::testing::Test
{
  protected:
    TriangleCoordUtilsTester() { setup(); }

    void setup() { setup(Point(0, 0, 0), Point(1, 0, 1), Point(0, 1, 1)); }

    void setup(const Point& p1, const Point& p2, const Point& p3)
    {
      mesh = mesh::make_empty_mesh();
      v1   = mesh->create_vertex(p1);
      v2   = mesh->create_vertex(p2);
      v3   = mesh->create_vertex(p3);
      tri  = mesh->create_triangle_from_verts(v1, v2, v3);
      e1   = tri->get_down(0);
      e2   = tri->get_down(1);
      e3   = tri->get_down(2);

      mesh::MeshEntityPtr v1Reversed = mesh->create_vertex(p1);
      mesh::MeshEntityPtr v2Reversed = mesh->create_vertex(p2);
      mesh::MeshEntityPtr v3Reversed = mesh->create_vertex(p3);
      e1Reversed                     = mesh->create_edge(v2Reversed, v1Reversed);
      e2Reversed                     = mesh->create_edge(v3Reversed, v2Reversed);
      e3Reversed                     = mesh->create_edge(v1Reversed, v3Reversed);
      triReversed                    = mesh->create_triangle_from_verts(v1Reversed, v2Reversed, v3Reversed);
    }

    std::shared_ptr<mesh::Mesh> mesh;
    mesh::MeshEntityPtr v1, v2, v3;
    mesh::MeshEntityPtr e1, e2, e3;
    mesh::MeshEntityPtr tri;
    mesh::MeshEntityPtr e1Reversed, e2Reversed, e3Reversed;
    mesh::MeshEntityPtr triReversed;
    predicates::impl::TriangleCoordUtils triCoordUtils;
    // utils::impl::Projection proj_xy;
    // utils::impl::Projection proj_yz;
    // utils::impl::Projection proj_zx;
};

void expect_near(const Point& pt1, const Point& pt2, double tol)
{
  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(pt1[i], pt2[i], tol);
}

} // namespace

TEST_F(TriangleCoordUtilsTester, compute_xyz_coords_Verts)
{
  PointRecordForTriangle r1(PointClassification::Vert, 0, tri);
  PointRecordForTriangle r2(PointClassification::Vert, 1, tri);
  PointRecordForTriangle r3(PointClassification::Vert, 2, tri);

  expect_near(triCoordUtils.compute_xyz_coords(r1), v1->get_point_orig(0), 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r2), v2->get_point_orig(0), 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r3), v3->get_point_orig(0), 1e-13);
}

TEST_F(TriangleCoordUtilsTester, compute_xyz_coords_Edges)
{
  PointRecordForTriangle r1(PointClassification::Edge, 0, tri, Point(0.5, 0.0));
  PointRecordForTriangle r2(PointClassification::Edge, 1, tri, Point(0.5, 0.5));
  PointRecordForTriangle r3(PointClassification::Edge, 2, tri, Point(0.0, 0.5));

  expect_near(triCoordUtils.compute_xyz_coords(r1), (v1->get_point_orig(0) + v2->get_point_orig(0)) / 2, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r2), (v2->get_point_orig(0) + v3->get_point_orig(0)) / 2, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r3), (v3->get_point_orig(0) + v1->get_point_orig(0)) / 2, 1e-13);
}

TEST_F(TriangleCoordUtilsTester, compute_xyz_coords_Interior)
{
  PointRecordForTriangle r1(PointClassification::Interior, -1, tri, Point(1.0 / 3.0, 1.0 / 3.0));

  Point centroid = (v1->get_point_orig(0) + v2->get_point_orig(0) + v3->get_point_orig(0)) / 3;
  expect_near(triCoordUtils.compute_xyz_coords(r1), centroid, 1e-13);
}

TEST_F(TriangleCoordUtilsTester, get_edge_xi_Standard)
{
  PointRecordForTriangle r1(PointClassification::Edge, 0, tri, Point(1.0 / 3.0, 0));
  EXPECT_NEAR(triCoordUtils.get_edge_xi(r1), 1.0 / 3.0, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r1), mesh::compute_edge_coords_orig(e1, triCoordUtils.get_edge_xi(r1)),
              1e-13);

  PointRecordForTriangle r2(PointClassification::Edge, 1, tri, Point(2.0 / 3.0, 1.0 / 3.0));
  EXPECT_NEAR(triCoordUtils.get_edge_xi(r2), 1.0 / 3.0, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r2), mesh::compute_edge_coords_orig(e2, triCoordUtils.get_edge_xi(r2)),
              1e-13);

  PointRecordForTriangle r3(PointClassification::Edge, 2, tri, Point(0, 2.0 / 3.0));
  EXPECT_NEAR(triCoordUtils.get_edge_xi(r3), 1.0 / 3.0, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r3), mesh::compute_edge_coords_orig(e3, triCoordUtils.get_edge_xi(r3)),
              1e-13);
}

TEST_F(TriangleCoordUtilsTester, get_edge_xi_Reversed)
{
  PointRecordForTriangle r1(PointClassification::Edge, 0, triReversed, Point(1.0 / 3.0, 0));
  EXPECT_NEAR(triCoordUtils.get_edge_xi(r1), 2.0 / 3.0, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r1),
              mesh::compute_edge_coords_orig(e1Reversed, triCoordUtils.get_edge_xi(r1)), 1e-13);

  PointRecordForTriangle r2(PointClassification::Edge, 1, triReversed, Point(2.0 / 3.0, 1.0 / 3.0));
  EXPECT_NEAR(triCoordUtils.get_edge_xi(r2), 2.0 / 3.0, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r2),
              mesh::compute_edge_coords_orig(e2Reversed, triCoordUtils.get_edge_xi(r2)), 1e-13);

  PointRecordForTriangle r3(PointClassification::Edge, 2, triReversed, Point(0, 2.0 / 3.0));
  EXPECT_NEAR(triCoordUtils.get_edge_xi(r3), 2.0 / 3.0, 1e-13);
  expect_near(triCoordUtils.compute_xyz_coords(r3),
              mesh::compute_edge_coords_orig(e3Reversed, triCoordUtils.get_edge_xi(r3)), 1e-13);
}

TEST_F(TriangleCoordUtilsTester, compute_orthogonal_dist)
{
  EXPECT_NEAR(triCoordUtils.compute_orthogonal_dist(mesh::MeshEntityType::Triangle, Point(0.5, -0.5), 0), 0.5, 1e-13);
  EXPECT_NEAR(triCoordUtils.compute_orthogonal_dist(mesh::MeshEntityType::Triangle, Point(-0.5, 0.5), 2), 0.5, 1e-13);
  EXPECT_NEAR(triCoordUtils.compute_orthogonal_dist(mesh::MeshEntityType::Triangle, Point(1, 1), 1), 1.0 / std::sqrt(2),
              1e-13);
}

TEST_F(TriangleCoordUtilsTester, orthogonal_proj_triangle)
{
  expect_near(triCoordUtils.orthogonal_proj_triangle(Point(0.5, -0.5), 0), Point(0.5, 0), 1e-13);
  expect_near(triCoordUtils.orthogonal_proj_triangle(Point(-0.5, 0.5), 2), Point(0, 0.5), 1e-13);
  expect_near(triCoordUtils.orthogonal_proj_triangle(Point(1, 1), 1), Point(-0.5, 0.5), 1e-13);
}

TEST_F(TriangleCoordUtilsTester, create_tri_record_edge)
{
  setup({0, 0, 0}, {1, 0, 0}, {0, 1, 0});

  // edge 0
  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 0, 0);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 0);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 0, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0.0, 0, 0), 1e-13);
  }

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 0, 0.25);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 0);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0.25, 0, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0.25, 0, 0), 1e-13);
  }  

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 0, 1);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 0);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(1, 0, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(1, 0, 0), 1e-13);
  }  

  // edge 2
  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 1, 0);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 1);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(1, 0, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(1, 0, 0), 1e-13);
  }

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 1, 0.25);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 1);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0.75, 0.25, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0.75, 0.25, 0), 1e-13);
  }  

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 1, 1);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 1);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 1, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0, 1, 0), 1e-13);

  }  

  // edge 3
  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 2, 0);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 2);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 1, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0, 1, 0), 1e-13);
  }

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 2, 0.25);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 2);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 0.75, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0, 0.75, 0), 1e-13);
  }  

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 2, 1);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Edge);
    EXPECT_EQ(record.id, 2);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 0, 0), 1e-13);
    expect_near(triCoordUtils.compute_xyz_coords(record), utils::Point(0, 0, 0), 1e-13);

  }    
}

TEST_F(TriangleCoordUtilsTester, create_tri_record_vert)
{
  setup({0, 0, 0}, {1, 0, 0}, {0, 1, 0});

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 0);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Vert);
    EXPECT_EQ(record.id, 0);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 0, 0), 1e-13);
  }

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 1);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Vert);
    EXPECT_EQ(record.id, 1);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(1, 0, 0), 1e-13);
  }  

  {
    PointRecordForTriangle record = triCoordUtils.create_record(tri, 2);
    EXPECT_EQ(record.el, tri);
    EXPECT_EQ(record.type, PointClassification::Vert);
    EXPECT_EQ(record.id, 2);
    expect_near(triCoordUtils.compute_xi_coords(record), utils::Point(0, 1, 0), 1e-13);
  }  
}

TEST(TriangleCoordUtils, classify_onto)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex({0, 0, 0});
  auto v2   = mesh->create_vertex({1, 0, 0});
  auto v3   = mesh->create_vertex({0, 1, 0});
  auto v4   = mesh->create_vertex({1, 1, 0});
  auto v5   = mesh->create_vertex({-1, 0, 0});
  auto v6   = mesh->create_vertex({0, -1, 0});
  auto tri1  = mesh->create_triangle_from_verts(v1, v2, v3);  
  auto tri2  = mesh->create_triangle_from_verts(v1, v4, v2);
  auto tri3  = mesh->create_triangle_from_verts(v2, v5, v3);  
  auto tri4  = mesh->create_triangle_from_verts(v1, v3, v6);
  predicates::impl::TriangleCoordUtils triCoordUtils;

  // vertex 1
  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 0);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri2);

    EXPECT_EQ(record2.el, tri2);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0, 0, 0}, 1e-13);
  }

  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 0);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri4);

    EXPECT_EQ(record2.el, tri4);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0, 0, 0}, 1e-13);
  }  

  // vertex 2
  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 1);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri2);

    EXPECT_EQ(record2.el, tri2);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 2);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {1, 0, 0}, 1e-13);
  }

  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 1);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri3);

    EXPECT_EQ(record2.el, tri3);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {1, 0, 0}, 1e-13);
  }

  // vertex 3   
  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 2);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri3);

    EXPECT_EQ(record2.el, tri3);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 2);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0, 1, 0}, 1e-13);
  }

  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 2);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri4);

    EXPECT_EQ(record2.el, tri4);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 1);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0, 1, 0}, 1e-13);
  }   

  // edge 1
  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 0, 0.25);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri2);

    EXPECT_EQ(record2.el, tri2);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 2);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0.25, 0, 0}, 1e-13);
  }

  // edge 2
  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 1, 0.25);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri3);

    EXPECT_EQ(record2.el, tri3);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 2);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0.75, 0.25, 0}, 1e-13);
  }  

  // edge 3
  {
    PointRecordForTriangle record1 = triCoordUtils.create_record(tri1, 2, 0.25);
    PointRecordForTriangle record2 = triCoordUtils.classify_onto(record1, tri4);

    EXPECT_EQ(record2.el, tri4);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 0);
    expect_near(triCoordUtils.compute_xyz_coords(record2), {0, 0.75, 0}, 1e-13);
  }  
}

TEST_F(TriangleCoordUtilsTester, ExteriorDeviation)
{
  setup({0, 0, 0}, {1, 0, 0}, {0, 1, 0});

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, tri, Point(0.5, -0.5));
    EXPECT_NEAR(triCoordUtils.compute_exterior_deviation(r1), 0.5, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, tri, Point(-0.5, 0.5));
    EXPECT_NEAR(triCoordUtils.compute_exterior_deviation(r1), 0.5, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, tri, Point(0.75, 0.75));
    EXPECT_NEAR(triCoordUtils.compute_exterior_deviation(r1), 0.5, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, tri, Point(1.5, 0.5));
    EXPECT_NEAR(triCoordUtils.compute_exterior_deviation(r1), 1.5, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, tri, Point(0.5, 1.5));
    EXPECT_NEAR(triCoordUtils.compute_exterior_deviation(r1), 1.5, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, tri, Point(0.5, -1));
    EXPECT_NEAR(triCoordUtils.compute_exterior_deviation(r1), 1.5, 1e-13);
  }  

}


} // namespace impl
} // namespace middle_mesh
} // namespace stk
