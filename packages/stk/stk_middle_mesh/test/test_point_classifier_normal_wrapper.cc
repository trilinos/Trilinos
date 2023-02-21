#include "predicates/point_classifier_normal_wrapper.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

namespace testing {

class PointClassifierNormalWrapperTester : public ::testing::Test
{
  protected:
    PointClassifierNormalWrapperTester() { mesh = mesh::make_empty_mesh(); }

    utils::Point compute_normal_vector(predicates::impl::PointClassifierNormalWrapper& c, mesh::MeshEntityPtr face,
                                       utils::Point pt)
    {
      return c.compute_normal_vector(face, pt);
    }

    utils::Point compute_triangle_xi_coords(predicates::impl::PointClassifierNormalWrapper& c,
                                            std::array<utils::Point, 3>& verts, const utils::Point& pt)
    {
      return c.compute_triangle_xi_coords(verts, pt);
    }

    std::shared_ptr<mesh::Mesh> mesh;
};
} // namespace testing

namespace {
void expect_near(const utils::Point& pt1, const utils::Point& pt2)
{
  EXPECT_NEAR(pt1.x, pt2.x, 1e-13);
  EXPECT_NEAR(pt1.y, pt2.y, 1e-13);
  EXPECT_NEAR(pt1.z, pt2.z, 1e-13);
}

void expect_vectors_parallel(const utils::Point& v1, const utils::Point& v2)
{
  double v1Len    = std::sqrt(dot(v1, v1));
  double v2Len    = std::sqrt(dot(v2, v2));
  double cosTheta = dot(v1, v2) / (v1Len * v2Len);
  EXPECT_NEAR(cosTheta, 1.0, 1e-13);
}
} // namespace

using testing::PointClassifierNormalWrapperTester;

TEST_F(PointClassifierNormalWrapperTester, Triangle)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  auto meshSrc = mesh::make_empty_mesh();
  auto v5      = meshSrc->create_vertex(0, 0);
  auto v6      = meshSrc->create_vertex(1, 0);
  auto v7      = meshSrc->create_vertex(0, 1);
  auto elSrc   = meshSrc->create_triangle_from_verts(v5, v6, v7);

  predicates::impl::PointClassifierNormalWrapper c(meshSrc);
  utils::Point pt;
  predicates::impl::PointRecord r;

  // verts
  pt = utils::Point(0, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 1);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  // edges
  pt = utils::Point(0.5, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(0.5, 0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  // interior
  pt = utils::Point(0.25, 0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.5, 0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, -0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.5, 2.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.25, -0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.25, 0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.25, 0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST_F(PointClassifierNormalWrapperTester, Quad)
{
  auto v1 = mesh->create_vertex(0, 0);
  auto v2 = mesh->create_vertex(1, 0);
  auto v3 = mesh->create_vertex(1, 1);
  auto v4 = mesh->create_vertex(0, 1);
  auto el = mesh->create_quad_from_verts(v1, v2, v3, v4);

  auto meshSrc = mesh::make_empty_mesh();
  auto v5      = meshSrc->create_vertex(0, 0);
  auto v6      = meshSrc->create_vertex(1, 0);
  auto v7      = meshSrc->create_vertex(0, 1);
  auto elSrc   = meshSrc->create_triangle_from_verts(v5, v6, v7);

  predicates::impl::PointClassifierNormalWrapper c(meshSrc);
  utils::Point pt;
  predicates::impl::PointRecord r;

  // verts
  pt = utils::Point(0, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(1, 1);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  pt = utils::Point(0, 1);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 3);

  // edges
  pt = utils::Point(0.5, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0.5, 1);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  pt = utils::Point(0, 0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 3);

  // interior
  pt = utils::Point(0.25, 0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.75, 0.25);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.75, 0.75);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.75);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, 1);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 1);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1, 1.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1, -0.5);
  r  = c.classify(el, elSrc, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST_F(PointClassifierNormalWrapperTester, TriangleReverse)
{
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0, 0);
  auto v2   = mesh->create_vertex(1, 0);
  auto v3   = mesh->create_vertex(0, 1);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  predicates::impl::PointClassifierNormalWrapper c(mesh);
  utils::Point pt;
  predicates::impl::PointRecord r;

  // verts
  pt = utils::Point(0, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 1);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  // edges
  pt = utils::Point(0.5, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(0.5, 0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  // interior
  pt = utils::Point(0.25, 0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.5, 0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, -0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.5, 2.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.25, -0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.25, 0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.25, 0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST_F(PointClassifierNormalWrapperTester, QuadReverse)
{
  auto v1 = mesh->create_vertex(0, 0);
  auto v2 = mesh->create_vertex(1, 0);
  auto v3 = mesh->create_vertex(1, 1);
  auto v4 = mesh->create_vertex(0, 1);
  auto el = mesh->create_quad_from_verts(v1, v2, v3, v4);

  predicates::impl::PointClassifierNormalWrapper c(mesh);
  utils::Point pt, normal(0, 0, 1);
  predicates::impl::PointRecord r;

  // verts
  pt = utils::Point(0, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(1, 1);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  pt = utils::Point(0, 1);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 3);

  // edges
  pt = utils::Point(0.5, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0.5, 1);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  pt = utils::Point(0, 0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 3);

  // interior
  pt = utils::Point(0.25, 0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.75, 0.25);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.75, 0.75);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.75);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, 1);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 1);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1, 1.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1, -0.5);
  r  = c.classify_reverse(el, pt);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST_F(PointClassifierNormalWrapperTester, TriangleNormal_xyplane)
{
  auto v1  = mesh->create_vertex(0, 0, 0);
  auto v2  = mesh->create_vertex(1, 0, 0);
  auto v3  = mesh->create_vertex(0, 1, 0);
  auto tri = mesh->create_triangle_from_verts(v1, v2, v3);
  predicates::impl::PointClassifierNormalWrapper c(mesh);

  auto normal = compute_normal_vector(c, tri, utils::Point(0.25, 0.25, 0));
  std::cout << "normal = " << normal << std::endl;
  expect_vectors_parallel(normal, utils::Point(0, 0, 1));
}

TEST_F(PointClassifierNormalWrapperTester, TriangleNormal_tilted)
{
  auto v1  = mesh->create_vertex(0, 0, 0);
  auto v2  = mesh->create_vertex(1, 0, 2);
  auto v3  = mesh->create_vertex(0, 1, 2);
  auto tri = mesh->create_triangle_from_verts(v1, v2, v3);
  predicates::impl::PointClassifierNormalWrapper c(mesh);

  auto normal = compute_normal_vector(c, tri, compute_tri_coords_from_xi(tri, utils::Point(0.25, 0.35)));

  utils::Point b1              = v2->get_point_orig(0) - v1->get_point_orig(0);
  utils::Point b2              = v3->get_point_orig(0) - v1->get_point_orig(0);
  utils::Point normalDirection = cross(b1, b2);
  expect_vectors_parallel(normal, normalDirection);
}

TEST_F(PointClassifierNormalWrapperTester, QuadNormal)
{
  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);
  auto v3 = mesh->create_vertex(1, 1, 1); // make the second triangle in a different plane
  auto v4 = mesh->create_vertex(0, 1, 0);

  auto quad = mesh->create_quad_from_verts(v1, v2, v3, v4);
  predicates::impl::PointClassifierNormalWrapper c(mesh);

  utils::Point b1      = v2->get_point_orig(0) - v1->get_point_orig(0);
  utils::Point b2      = v4->get_point_orig(0) - v1->get_point_orig(0);
  utils::Point b3      = v4->get_point_orig(0) - v3->get_point_orig(0);
  utils::Point b4      = v2->get_point_orig(0) - v3->get_point_orig(0);
  utils::Point normal1 = cross(b1, b2);
  utils::Point normal2 = cross(b3, b4);

  utils::Point pt1 = v1->get_point_orig(0);
  auto normal      = compute_normal_vector(c, quad, pt1);
  expect_vectors_parallel(normal, normal1);

  utils::Point pt2 = v3->get_point_orig(0);
  normal           = compute_normal_vector(c, quad, pt2);
  expect_vectors_parallel(normal, normal2);
}

TEST_F(PointClassifierNormalWrapperTester, TriangleXiCoords)
{
  predicates::impl::PointClassifierNormalWrapper c(mesh);

  {
    // xy plane
    utils::Point pt1(0, 0, 0), pt2(1, 0, 0), pt3(0, 1, 0);
    std::array<utils::Point, 3> pts = {pt1, pt2, pt3};
    utils::Point centroid           = (pt1 + pt2 + pt3) / 3;
    utils::Point ptXi               = compute_triangle_xi_coords(c, pts, centroid);
    expect_near(ptXi, utils::Point(1.0 / 3.0, 1.0 / 3.0, 0));
  }

  {
    // yz plane
    utils::Point pt1(0, 0, 0), pt2(0, 1, 0), pt3(0, 0, 1);
    std::array<utils::Point, 3> pts = {pt1, pt2, pt3};
    utils::Point centroid           = (pt1 + pt2 + pt3) / 3;
    utils::Point ptXi               = compute_triangle_xi_coords(c, pts, centroid);
    expect_near(ptXi, utils::Point(1.0 / 3.0, 1.0 / 3.0, 0));
  }

  {
    // xz plane
    utils::Point pt1(0, 0, 0), pt2(1, 0, 0), pt3(0, 0, 1);
    std::array<utils::Point, 3> pts = {pt1, pt2, pt3};
    utils::Point centroid           = (pt1 + pt2 + pt3) / 3;
    utils::Point ptXi               = compute_triangle_xi_coords(c, pts, centroid);
    expect_near(ptXi, utils::Point(1.0 / 3.0, 1.0 / 3.0, 0));
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
