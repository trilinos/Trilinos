#include "stk_middle_mesh/predicates/element_xi_classifier.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {
using namespace predicates::impl;
using namespace utils;

void test_record(const predicates::impl::PointRecordForTriangle& r, PointClassification type, int id)
{
  EXPECT_EQ(r.type, type);
  EXPECT_EQ(r.id, id);
}

} // namespace

TEST(ElementXiClassifier, TriangleVerts)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_triangle(Point(0, 0)), PointClassification::Vert, 0);
  test_record(classifier.classify_triangle(Point(eps / 2, 0)), PointClassification::Vert, 0);
  test_record(classifier.classify_triangle(Point(0, eps / 2)), PointClassification::Vert, 0);

  test_record(classifier.classify_triangle(Point(1, 0)), PointClassification::Vert, 1);
  test_record(classifier.classify_triangle(Point(1 - eps / 2, 0)), PointClassification::Vert, 1);
  test_record(classifier.classify_triangle(Point(1, eps / 2)), PointClassification::Vert, 1);

  test_record(classifier.classify_triangle(Point(0, 1)), PointClassification::Vert, 2);
  test_record(classifier.classify_triangle(Point(eps / 2, 1)), PointClassification::Vert, 2);
  test_record(classifier.classify_triangle(Point(0, 1 + eps / 2)), PointClassification::Vert, 2);
}

TEST(ElementXiClassifier, TriangleEdges)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_triangle(Point(0.5, 0)), PointClassification::Edge, 0);
  test_record(classifier.classify_triangle(Point(0.5, eps / 2)), PointClassification::Edge, 0);
  test_record(classifier.classify_triangle(Point(2 * eps, 0)), PointClassification::Edge, 0);
  test_record(classifier.classify_triangle(Point(1 - 2 * eps, 0)), PointClassification::Edge, 0);

  Point offset = 2 * eps * Point(-1, 1) / std::sqrt(2);
  test_record(classifier.classify_triangle(Point(0.5, 0.5)), PointClassification::Edge, 1);
  test_record(classifier.classify_triangle(Point(0.5, 0.5) + offset), PointClassification::Edge, 1);
  test_record(classifier.classify_triangle(Point(1, 0) + offset), PointClassification::Edge, 1);
  test_record(classifier.classify_triangle(Point(0, 1) - offset), PointClassification::Edge, 1);
}

TEST(ElementXiClassifier, TriangleInterior)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_triangle(Point(2 * eps, 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_triangle(Point(1 - 4 * eps, 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_triangle(Point(2 * eps, 1 - 4 * eps)), PointClassification::Interior, 0);

  test_record(classifier.classify_triangle(Point(0.5, 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_triangle(Point(0.5, 0.5 - 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_triangle(Point(0.5 - 2 * eps, 0.5)), PointClassification::Interior, 0);
  test_record(classifier.classify_triangle(Point(0 + 2 * eps, 0.5)), PointClassification::Interior, 0);
}

TEST(ElementXiClassifier, TriangleExterior)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_triangle(Point(-2 * eps, -2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_triangle(Point(1 + 4 * eps, -2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_triangle(Point(-2 * eps, 1 + 4 * eps)), PointClassification::Exterior, -1);

  test_record(classifier.classify_triangle(Point(0.5, -2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_triangle(Point(0.5, 0.5 + 2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_triangle(Point(0.5 + 2 * eps, 0.5)), PointClassification::Exterior, -1);
  test_record(classifier.classify_triangle(Point(0 - 2 * eps, 0.5)), PointClassification::Exterior, -1);
}

TEST(ElementXiClassifier, QuadVerts)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_quad(Point(0, 0)), PointClassification::Vert, 0);
  test_record(classifier.classify_quad(Point(eps / 2, 0)), PointClassification::Vert, 0);
  test_record(classifier.classify_quad(Point(-eps / 2, 0)), PointClassification::Vert, 0);
  test_record(classifier.classify_quad(Point(0, eps / 2)), PointClassification::Vert, 0);
  test_record(classifier.classify_quad(Point(0, -eps / 2)), PointClassification::Vert, 0);

  test_record(classifier.classify_quad(Point(1, 0)), PointClassification::Vert, 1);
  test_record(classifier.classify_quad(Point(1 - eps / 2, 0)), PointClassification::Vert, 1);
  test_record(classifier.classify_quad(Point(1 + eps / 2, 0)), PointClassification::Vert, 1);
  test_record(classifier.classify_quad(Point(1, eps / 2)), PointClassification::Vert, 1);
  test_record(classifier.classify_quad(Point(1, -eps / 2)), PointClassification::Vert, 1);

  test_record(classifier.classify_quad(Point(1, 1)), PointClassification::Vert, 2);
  test_record(classifier.classify_quad(Point(1 - eps / 2, 1)), PointClassification::Vert, 2);
  test_record(classifier.classify_quad(Point(1 + eps / 2, 1)), PointClassification::Vert, 2);
  test_record(classifier.classify_quad(Point(1, 1 - eps / 2)), PointClassification::Vert, 2);
  test_record(classifier.classify_quad(Point(1, 1 + eps / 2)), PointClassification::Vert, 2);

  test_record(classifier.classify_quad(Point(0, 1)), PointClassification::Vert, 3);
  test_record(classifier.classify_quad(Point(0 - eps / 2, 1)), PointClassification::Vert, 3);
  test_record(classifier.classify_quad(Point(0 + eps / 2, 1)), PointClassification::Vert, 3);
  test_record(classifier.classify_quad(Point(0, 1 - eps / 2)), PointClassification::Vert, 3);
  test_record(classifier.classify_quad(Point(0, 1 + eps / 2)), PointClassification::Vert, 3);
}

TEST(ElementXiClassifier, QuadEdges)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_quad(Point(0.5, 0)), PointClassification::Edge, 0);
  test_record(classifier.classify_quad(Point(0.5, eps / 2)), PointClassification::Edge, 0);
  test_record(classifier.classify_quad(Point(0.5, -eps / 2)), PointClassification::Edge, 0);
  test_record(classifier.classify_quad(Point(0 + 2 * eps, 0)), PointClassification::Edge, 0);
  test_record(classifier.classify_quad(Point(1 - 2 * eps, 0)), PointClassification::Edge, 0);

  test_record(classifier.classify_quad(Point(1, 0.5)), PointClassification::Edge, 1);
  test_record(classifier.classify_quad(Point(1 + eps / 2, 0.5)), PointClassification::Edge, 1);
  test_record(classifier.classify_quad(Point(1 - eps / 2, 0.5)), PointClassification::Edge, 1);
  test_record(classifier.classify_quad(Point(1, 2 * eps)), PointClassification::Edge, 1);
  test_record(classifier.classify_quad(Point(1, 1 - 2 * eps)), PointClassification::Edge, 1);

  test_record(classifier.classify_quad(Point(0.5, 1)), PointClassification::Edge, 2);
  test_record(classifier.classify_quad(Point(0.5, 1 + eps / 2)), PointClassification::Edge, 2);
  test_record(classifier.classify_quad(Point(0.5, 1 - eps / 2)), PointClassification::Edge, 2);
  test_record(classifier.classify_quad(Point(0 + 2 * eps, 1)), PointClassification::Edge, 2);
  test_record(classifier.classify_quad(Point(1 - 2 * eps, 1)), PointClassification::Edge, 2);

  test_record(classifier.classify_quad(Point(0, 0.5)), PointClassification::Edge, 3);
  test_record(classifier.classify_quad(Point(0 + eps / 2, 0.5)), PointClassification::Edge, 3);
  test_record(classifier.classify_quad(Point(0 - eps / 2, 0.5)), PointClassification::Edge, 3);
  test_record(classifier.classify_quad(Point(0, 2 * eps)), PointClassification::Edge, 3);
  test_record(classifier.classify_quad(Point(0, 1 - 2 * eps)), PointClassification::Edge, 3);
}

TEST(ElementXiClassifier, QuadInterior)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_quad(Point(2 * eps, 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_quad(Point(1 - 2 * eps, 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_quad(Point(1 - 2 * eps, 1 - 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_quad(Point(2 * eps, 1 - 2 * eps)), PointClassification::Interior, 0);

  test_record(classifier.classify_quad(Point(0.5, 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_quad(Point(1 - 2 * eps, 0.5)), PointClassification::Interior, 0);
  test_record(classifier.classify_quad(Point(0.5, 1 - 2 * eps)), PointClassification::Interior, 0);
  test_record(classifier.classify_quad(Point(2 * eps, 0.5)), PointClassification::Interior, 0);
}

TEST(ElementXiClassifier, QuadExterior)
{
  double eps = 0.1; // make this large so testing is easier
  predicates::impl::ElementXiClassifier classifier(eps);

  test_record(classifier.classify_quad(Point(-2 * eps, -2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_quad(Point(1 + 2 * eps, -2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_quad(Point(1 + 2 * eps, 1 + 2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_quad(Point(-2 * eps, 1 + 2 * eps)), PointClassification::Exterior, -1);

  test_record(classifier.classify_quad(Point(0.5, -2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_quad(Point(1 + 2 * eps, 0.5)), PointClassification::Exterior, -1);
  test_record(classifier.classify_quad(Point(0.5, 1 + 2 * eps)), PointClassification::Exterior, -1);
  test_record(classifier.classify_quad(Point(-2 * eps, 0.5)), PointClassification::Exterior, -1);
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
