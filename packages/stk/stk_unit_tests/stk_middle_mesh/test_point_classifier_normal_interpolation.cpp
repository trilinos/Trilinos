#include "stk_middle_mesh/predicates/point_classifier_normal_interpolation.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

TEST(PointClassifierNormalInterpolation, Triangle)
{
  utils::Point offset = utils::Point(1, 2);
  auto mesh           = mesh::make_empty_mesh();
  auto v1             = mesh->create_vertex(0 + offset.x, 0 + offset.y);
  auto v2             = mesh->create_vertex(1 + offset.x, 0 + offset.y);
  auto v3             = mesh->create_vertex(0 + offset.x, 1 + offset.y);
  auto el             = mesh->create_triangle_from_verts(v1, v2, v3);
  utils::Point normal = {0, 0, 1};

  predicates::impl::PointClassifierNormalInterpolation c;
  utils::Point pt;
  predicates::impl::PointRecordForTriangle r;

  // verts
  pt = utils::Point(0, 0) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 1) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  // edges
  pt = utils::Point(0.5, 0) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(0.5, 0.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 0.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  // interior
  pt = utils::Point(0.25, 0.25) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.5, 0.25) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, -0.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.5, 2.5) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.25, -0.25) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.25, 0.25) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.25, 0.25) + offset;
  r  = c.classify(el, pt, normal);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST(PointClassifierNormalInterpolation, TriangleReverse)
{
  utils::Point offset                 = utils::Point(1, 2);
  auto mesh                           = mesh::make_empty_mesh();
  auto v1                             = mesh->create_vertex(0 + offset.x, 0 + offset.y);
  auto v2                             = mesh->create_vertex(1 + offset.x, 0 + offset.y);
  auto v3                             = mesh->create_vertex(0 + offset.x, 1 + offset.y);
  auto el                             = mesh->create_triangle_from_verts(v1, v2, v3);
  utils::Point normal                 = {0, 0, 1};
  std::array<utils::Point, 4> normals = {normal, normal, normal, normal};

  predicates::impl::PointClassifierNormalInterpolation c;
  utils::Point pt;
  predicates::impl::PointRecordForTriangle r;

  // verts
  pt = utils::Point(0, 0) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(1, 0) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 1) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Vert);
  EXPECT_EQ(r.id, 2);

  // edges
  pt = utils::Point(0.5, 0) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 0);

  pt = utils::Point(0.5, 0.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 1);

  pt = utils::Point(0, 0.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Edge);
  EXPECT_EQ(r.id, 2);

  // interior
  pt = utils::Point(0.25, 0.25) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.5, 0.25) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Interior);

  pt = utils::Point(0.25, 0.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Interior);

  // exterior
  pt = utils::Point(1.5, 0) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.5, 0) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, 1.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0, -0.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.5, -0.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.5, 2.5) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(0.25, -0.25) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(-0.25, 0.25) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);

  pt = utils::Point(1.25, 0.25) + offset;
  r  = c.classify_reverse(el, pt, normals);
  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST(PointClassifierNormalInterpolation, ClassifyReverseRegression)
{
  // For this case, solving the 3x3 system for (lambda1, lambda2, gamma) fails
  // to converge, but solving the reduced system f(gamma, lambda(gamma)) = 0
  // converges easily
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(0.4635254915624213, 1.42658477444273, 2.107342425544701e-08);
  auto v2   = mesh->create_vertex(0.3428565549738922, 1.2741065117333, 0.7163027771654847);
  auto v3   = mesh->create_vertex(0.8816778784387098, 1.213525491562421, 0);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  utils::Point pt                     = utils::Point(0.2163118960624633, 1.365739561406607, 0.3741657439457501);
  std::array<utils::Point, 4> normals = {utils::Point(0.1026012399719781, 0.3157741471555242, 0.08450264273327771),
                                         utils::Point(0.04164258531798171, 0.1609265953454822, 0.06026809628897314),
                                         utils::Point(0.1324859315710616, 0.2041829292500395, 0.05967711046837183)};

  predicates::impl::PointClassifierNormalInterpolation c;
  predicates::impl::PointRecordForTriangle r = c.classify_reverse(el, pt, normals);

  EXPECT_EQ(r.type, PointClassification::Exterior);
}

TEST(PointClassifierNormalInterpolation, ClassifyReverseRegression2)
{
  // For this case, solving the 3x3 system for (lambda1, lambda2, gamma)
  // converges, but solving the reduced system f(gamma, lambda(gamma)) = 0
  // converges does not when using a change of basis of (delta_x1, delta_x2, cross(delta_x1, delta_x2))
  // Using a change of basis where the normal vector is included in the first two directions (the ones
  // used for the lambda solve) converges easily
  auto mesh = mesh::make_empty_mesh();
  auto v1   = mesh->create_vertex(1.198232336377407, 0, 0.8864670496909353);
  auto v2   = mesh->create_vertex(1.259782861986013, 0, 0.8029206973440202);
  auto v3   = mesh->create_vertex(1.19999381338826, 0.11076823646181, 0.8931680573159947);
  auto el   = mesh->create_triangle_from_verts(v1, v2, v3);

  utils::Point pt                     = utils::Point(1.365739561406607, 0.2163118960624631, 0.3741657447130542);
  std::array<utils::Point, 4> normals = {
      utils::Point(0.008856188164528336, -0.0005355423887644713, 0.006524558177620529),
      utils::Point(0.007914267788425789, -0.0001960064088777902, 0.00583062370170439),
      utils::Point(0.008563251912758543, 0.0003354223762162611, 0.006316635925236954)};

  predicates::impl::PointClassifierNormalInterpolation c;
  predicates::impl::PointRecordForTriangle r =
      c.classify_reverse(el, pt, normals, true); // TODO: can't figure out how to make this converge

  EXPECT_EQ(r.type, PointClassification::Exterior);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
