#include "gtest/gtest.h"

#include <cmath>

#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/predicates/intersection.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

namespace {

void expect_near(const utils::Point& pt1, const utils::Point& pt2, double tol)
{
  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(pt1[i], pt2[i], tol);
}

} // namespace

TEST(PointClassifier, Triangle)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  PointClassifier c;
  utils::Point pt;
  PointRecord r;

  // verts
  pt = utils::Point(0, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  // edges
  pt = utils::Point(0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(0.5, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  // interior
  pt = utils::Point(0.25, 0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.5, 0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, -0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.5, 2.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.25, -0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.25, 0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.25, 0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST(PointClassifier, Quad)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  PointClassifier c;
  utils::Point pt;
  PointRecord r;

  // verts
  pt = utils::Point(0, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(1, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  pt = utils::Point(0, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 3);

  // edges
  pt = utils::Point(0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0.5, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 3);

  // interior
  pt = utils::Point(0.25, 0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.75, 0.25);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.75, 0.75);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.75);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1, 1.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1, -0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST(PointClassifier, EdgeXiQuad)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  PointClassifier c;
  utils::Point pt;
  PointRecord r;

  pt = utils::Point(0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);

  pt = utils::Point(1, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);

  pt = utils::Point(0.5, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);
}

TEST(PointClassifier, EdgeXiTriangle)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  PointClassifier c;
  utils::Point pt;
  PointRecord r;

  pt = utils::Point(0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);

  pt = utils::Point(0.5, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_FLOAT_EQ(c.get_edge_xi(r), 0.5);
}

TEST(PointClassifier, LineEdgeXiQuad)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  const double eps           = 1e-8; // set eps large to test fallback cases
  const double smallSlopeEps = 1e-6;
  PointClassifier c(eps, smallSlopeEps);
  utils::Point pt1, pt2;
  PointRecord r1, r2;
  PossibleEdgeIntersection p;

  pt1 = utils::Point(0.5, 0.5);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Interior);

  // bisect edges
  pt2 = utils::Point(0.5, -1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt2 = utils::Point(1.5, 0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 1);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt2 = utils::Point(0.5, 2);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 2);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt2 = utils::Point(-0.5, 0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 3);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  // diagonal lines
  pt2 = utils::Point(0.25, -0.25);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 1.0 / 3.0);

  pt2 = utils::Point(1.25, 0.25);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 1);
  EXPECT_FLOAT_EQ(p.edgeXi, 1.0 / 3.0);

  // fallback cases
  pt1 = utils::Point(2 * eps, 0.25);
  pt2 = utils::Point(-2 * eps, 0.75);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Interior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 3);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt1 = utils::Point(0.25, 2 * eps);
  pt2 = utils::Point(0.75, -2 * eps);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Interior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  // test with point inside the lower left triangle
  pt1 = utils::Point(0.25, 0.5);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Interior);

  // edge intersections
  pt2 = utils::Point(0.25, -1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.25);

  pt2 = utils::Point(1.5, 0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 1);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.50);

  pt2 = utils::Point(0.25, 2);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 2);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.75);

  pt2 = utils::Point(-0.5, 0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 3);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  // vertex intersections
  // (there are two possible solutions for each of these, which  result is
  //  returned is implementation defined)
  pt2 = utils::Point(-1, -2);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 0);
  // EXPECT_FLOAT_EQ(p.edgeXi, 0);

  pt2 = utils::Point(1.75, -0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 1);
  // EXPECT_FLOAT_EQ(p.edgeXi, 1);

  // TODO: missing vertex 2 case

  pt2 = utils::Point(-1, 3);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 3);
  // EXPECT_FLOAT_EQ(p.edgeXi, 1);

  // test with point inside the upper right triangle
  pt1 = utils::Point(0.75, 0.5);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Interior);

  // edge intersections
  pt2 = utils::Point(-1, -3);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt2 = utils::Point(2, 0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 1);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt2 = utils::Point(0, 2);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 2);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt2 = utils::Point(-1, 0.5);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 3);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  // vertex intersections
  pt2 = utils::Point(-1.5, -1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 0);
  // EXPECT_FLOAT_EQ(p.edgeXi, 0);

  pt2 = utils::Point(2, -2);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 1);
  // EXPECT_FLOAT_EQ(p.edgeXi, 1);

  pt2 = utils::Point(2, 3);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 2);
  // EXPECT_FLOAT_EQ(p.edgeXi, 0);

  pt2 = utils::Point(-1.5, 2);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Vert);
  EXPECT_EQ(p.record.id, 3);
  // EXPECT_FLOAT_EQ(p.edgeXi, 1);
}

TEST(PointClassifier, LineEdgeXiTriangle)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  const double eps           = 1e-8; // set eps large to test fallback cases
  const double smallSlopeEps = 1e-6;
  PointClassifier c(eps, smallSlopeEps);
  utils::Point pt1, pt2;
  PointRecord r1, r2;
  PossibleEdgeIntersection p;

  pt1 = utils::Point(0.25, 0.25);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Interior);

  pt2 = utils::Point(0.25, -0.25);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.25);

  pt2 = utils::Point(-0.25, 0.25);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 2);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.75);

  pt2 = utils::Point(1.25, 0.25);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 1);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.25);

  pt2 = utils::Point(0.25, 1.25);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 1);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.75);

  // test fallback cases
  pt1 = utils::Point(2 * eps, 0.25);
  pt2 = utils::Point(-2 * eps, 0.75);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Interior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 2);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

  pt1 = utils::Point(0.25, 2 * eps);
  pt2 = utils::Point(0.75, -2 * eps);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Interior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  p = c.get_edge_xi(r1, r2);
  EXPECT_EQ(p.record.type, PointClassification::Edge);
  EXPECT_EQ(p.record.id, 0);
  EXPECT_FLOAT_EQ(p.edgeXi, 0.5);
}

TEST(PointClassifier, EdgeCorner)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  const double eps           = 1e-8; // set eps large to test fallback cases
  const double smallSlopeEps = 1e-6;
  PointClassifier c(eps, smallSlopeEps);
  utils::Point pt1, pt2;
  PointRecord r1, r2;
  CornerRecord r3;

  // horizontal line
  pt1 = utils::Point(0.25, -0.25);
  pt2 = utils::Point(0.25, 1.25);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_TRUE(r3.hasIntersection);
  EXPECT_EQ(r3.record1.type, PointClassification::Edge);
  EXPECT_EQ(r3.record2.type, PointClassification::Edge);
  EXPECT_EQ(r3.record1.id, 0);
  EXPECT_EQ(r3.record2.id, 2);
  EXPECT_FLOAT_EQ(r3.xi1, 0.25);
  EXPECT_FLOAT_EQ(r3.xi2, 0.75);

  // vertical line
  pt1 = utils::Point(-0.25, 0.25);
  pt2 = utils::Point(1.25, 0.25);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_TRUE(r3.hasIntersection);
  EXPECT_EQ(r3.record1.type, PointClassification::Edge);
  EXPECT_EQ(r3.record2.type, PointClassification::Edge);
  EXPECT_EQ(r3.record1.id, 3);
  EXPECT_EQ(r3.record2.id, 1);
  EXPECT_FLOAT_EQ(r3.xi1, 0.75);
  EXPECT_FLOAT_EQ(r3.xi2, 0.25);

  // general line
  pt1 = utils::Point(-0.5, 0);
  pt2 = utils::Point(1, 1.5);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_TRUE(r3.hasIntersection);
  EXPECT_EQ(r3.record1.type, PointClassification::Edge);
  EXPECT_EQ(r3.record2.type, PointClassification::Edge);
  EXPECT_EQ(r3.record1.id, 3);
  EXPECT_EQ(r3.record2.id, 2);
  EXPECT_FLOAT_EQ(r3.xi1, 0.5);
  EXPECT_FLOAT_EQ(r3.xi2, 0.5);

  pt1 = utils::Point(-1.5, 3);
  pt2 = utils::Point(2, -0.5);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_TRUE(r3.hasIntersection);
  EXPECT_EQ(r3.record1.type, PointClassification::Edge);
  EXPECT_EQ(r3.record2.type, PointClassification::Edge);
  EXPECT_EQ(r3.record1.id, 2);
  EXPECT_EQ(r3.record2.id, 1);
  EXPECT_FLOAT_EQ(r3.xi1, 0.5);
  EXPECT_FLOAT_EQ(r3.xi2, 0.5);

  pt1 = utils::Point(-1.5, 2);
  pt2 = utils::Point(2, -1.5);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_TRUE(r3.hasIntersection);
  EXPECT_EQ(r3.record1.type, PointClassification::Edge);
  EXPECT_EQ(r3.record2.type, PointClassification::Edge);
  EXPECT_EQ(r3.record1.id, 0);
  EXPECT_EQ(r3.record2.id, 3);
  EXPECT_FLOAT_EQ(r3.xi1, 0.5);
  EXPECT_FLOAT_EQ(r3.xi2, 0.5);

  // non intersection
  pt1 = utils::Point(-0.25, -0.25);
  pt2 = utils::Point(-0.25, 1.25);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_FALSE(r3.hasIntersection);

  pt1 = utils::Point(-0.5, 0);
  pt2 = utils::Point(-0.25, 1.5);
  r1  = c.classify(el, pt1);
  r2  = c.classify(el, pt2);
  EXPECT_EQ(r1.type, PointClassification::Exterior);
  EXPECT_EQ(r2.type, PointClassification::Exterior);
  r3 = c.get_edge_xi_corner(r1, r2);
  EXPECT_FALSE(r3.hasIntersection);
}

TEST(PointClassifier, EdgeIntersectionXiQuad)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  const double eps           = 1e-13;
  const double smallSlopeEps = 1e-13;
  PointClassifier c(eps, smallSlopeEps);
  utils::Point pt1, pt2;
  PointRecord r1, r2;
  PossibleEdgeIntersection p;
  ///*
  // test vertex to external
  pt1 = utils::Point(0, 0);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Vert);
  EXPECT_EQ(r1.id, 0);

  // intersection points are excluded verts: should be excluded
  {
    pt2 = utils::Point(2, 0);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 1);

    pt2 = utils::Point(0, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 3);
  }

  // no intersections at all
  {
    pt2 = utils::Point(0.5, -0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);

    pt2 = utils::Point(-0.5, 0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);
  }

  // edge intersections
  {
    pt2 = utils::Point(2, 1);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 1);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

    pt2 = utils::Point(2, 0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 1);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.25);

    pt2 = utils::Point(0.5, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.75);
  }

  // vertex intersection
  {
    pt2 = utils::Point(2, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 2);
  }

  // test shared vertex to external
  pt1 = utils::Point(1, 0);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Vert);
  EXPECT_EQ(r1.id, 1);

  // intersection points are excluded verts: should be excluded
  {
    pt2 = utils::Point(-1, 0);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(r1.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 0);

    pt2 = utils::Point(1, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(r1.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 2);
  }

  // no intersections at all
  {
    pt2 = utils::Point(0.5, -0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);

    pt2 = utils::Point(1.5, 0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);
  }

  // edge intersections
  {
    pt2 = utils::Point(0, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

    pt2 = utils::Point(-1, 1);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 3);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);
  }

  // vertex intersection
  {
    pt2 = utils::Point(-1, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 3);
  }
  // test edge to external
  pt1 = utils::Point(0.5, 0);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Edge);
  EXPECT_EQ(r1.id, 0);

  // intersection points are excluded verts: should be excluded
  {
    pt2 = utils::Point(-1, 0);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 0);

    pt2 = utils::Point(2, 0);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 1);
  }

  // no intersections at all
  {
    pt2 = utils::Point(0.25, -0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);

    pt2 = utils::Point(0.75, -0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);

    pt2 = utils::Point(0.50, -0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);
  }

  // edge intersections
  {
    pt2 = utils::Point(2, 1.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 1);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

    pt2 = utils::Point(2, 6);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.25);

    pt2 = utils::Point(0, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.75);

    pt2 = utils::Point(0.5, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.50);

    pt2 = utils::Point(-1, 1.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 3);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);
  }

  // vertex intersection
  {
    pt2 = utils::Point(2, 3);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 2);

    pt2 = utils::Point(-1, 3);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 3);
  }
}

TEST(PointClassifier, EdgeIntersectionXiQuadOverlappedEdges)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  const double eps           = 1e-13;
  const double smallSlopeEps = 1e-13;
  PointClassifier c(eps, smallSlopeEps);
  PointRecord r1, r2;
  PossibleEdgeIntersection p, p2;

  {
    // edge bisects element
    r1 = c.classify(el, {0.5, 0});
    r2 = c.classify(el, {0, 1.5});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
  }

  {
    // edges overlap
    r1 = c.classify(el, {0.5, 0});
    r2 = c.classify(el, {1.5, 0});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 1);
  }

  {
    // edge passes through non-adjacent vertex
    r1 = c.classify(el, {0.5, 0});
    r2 = c.classify(el, {1.5, 2});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 2);
  }

  {
    // edge passes through non-adjacent vertex
    r1 = c.classify(el, {0.5, 0});
    r2 = c.classify(el, {-1.5, 4});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 3);
  }

  {
    // edges has endpoint on edge 1, passes through another edge
    r1 = c.classify(el, {1, 0.5});
    r2 = c.classify(el, {-1, 0.5});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 3);
  }

  {
    // edges has endpoint on edge 1, passes through another edge
    r1 = c.classify(el, {1, 0.5});
    r2 = c.classify(el, {-1, -1});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 0);
  }

  {
    // edge passes through adjacent vert, starting on edge 1
    r1 = c.classify(el, {1, 0.5});
    r2 = c.classify(el, {1, 2});
    p  = c.get_edge_intersection_xi(r1, r2);

    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 2);
  }
}

TEST(PointClassifier, EdgeIntersectionXiTriangle)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  const double eps           = 1e-13;
  const double smallSlopeEps = 1e-13;
  PointClassifier c(eps, smallSlopeEps);
  utils::Point pt1, pt2;
  PointRecord r1, r2;
  PossibleEdgeIntersection p;

  // test vertex to external
  pt1 = utils::Point(0, 0);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Vert);
  EXPECT_EQ(r1.id, 0);

  // intersection points are excluded verts: should be excluded
  {
    pt2 = utils::Point(2, 0);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);

    pt2 = utils::Point(0, 2);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);
  }

  // no intersections at all
  {
    pt2 = utils::Point(0.5, -0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);

    pt2 = utils::Point(-0.5, 0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.id, -1);
  }

  // edge intersections
  {
    pt2 = utils::Point(1, 1);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 1);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);
  }

  // test with internal point on edge
  pt1 = utils::Point(0.5, 0);
  r1  = c.classify(el, pt1);
  EXPECT_EQ(r1.type, PointClassification::Edge);
  EXPECT_EQ(r1.id, 0);

  {
    pt2 = utils::Point(-0.5, 1);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 2);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.5);

    pt2 = utils::Point(1, 0.5);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Edge);
    EXPECT_EQ(p.record.id, 1);
    EXPECT_FLOAT_EQ(p.edgeXi, 0.25);

    // vertex intersection
    pt2 = utils::Point(-1, 3);
    r2  = c.classify(el, pt2);
    EXPECT_EQ(r2.type, PointClassification::Exterior);
    p = c.get_edge_intersection_xi(r1, r2);
    EXPECT_EQ(p.record.type, PointClassification::Vert);
    EXPECT_EQ(p.record.id, 2);
  }
}

TEST(PointClassifier, OrthogonalProj)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(1, 1);
  auto v4   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_quad_from_verts(v1, v2, v3, v4);

  PointClassifier c;
  utils::Point pt;
  PointRecord r;

  pt = utils::Point(0.25, 0.35);
  r  = c.classify(el, pt);
  EXPECT_FLOAT_EQ(c.compute_orthogonal_dist(r, 0), 0.35);
  EXPECT_FLOAT_EQ(c.compute_orthogonal_dist(r, 1), 0.75);
  EXPECT_FLOAT_EQ(c.compute_orthogonal_dist(r, 2), 0.65);
  EXPECT_FLOAT_EQ(c.compute_orthogonal_dist(r, 3), 0.25);

  EXPECT_FLOAT_EQ(c.get_edge_xi_orthogonal(r, 0), 0.25);
  EXPECT_FLOAT_EQ(c.get_edge_xi_orthogonal(r, 1), 0.35);
  EXPECT_FLOAT_EQ(c.get_edge_xi_orthogonal(r, 2), 0.75);
  EXPECT_FLOAT_EQ(c.get_edge_xi_orthogonal(r, 3), 0.65);
}

TEST(PointClassifier, XYZCoordsTriangle)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  PointClassifier c;
  utils::Point pt;
  PointRecord r;

  pt = utils::Point(0.25, 0.35);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);

  pt = utils::Point(0, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);

  pt = utils::Point(1, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);

  pt = utils::Point(0, 1);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);

  pt = utils::Point(0.5, 0);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);

  pt = utils::Point(0.5, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  expect_near(c.compute_xyz_coords(r), pt, 1e-13);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
