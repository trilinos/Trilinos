#include "gtest/gtest.h"

#include <cmath>

#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/predicates/intersection.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

namespace {

class QuadIntersectionTester : public ::testing::Test
{
  protected:
    QuadIntersectionTester()
    {
      make_standard_mesh();
      make_reversed_mesh();
    }

    void make_standard_mesh()
    {
      m_mesh  = mesh::make_empty_mesh();
      auto v1 = m_mesh->create_vertex(1, 1);
      auto v2 = m_mesh->create_vertex(2, 1);
      auto v3 = m_mesh->create_vertex(2, 2);
      auto v4 = m_mesh->create_vertex(1, 2);
      m_el    = m_mesh->create_quad_from_verts(v1, v2, v3, v4);
    }

    void make_reversed_mesh()
    {
      m_meshReversed         = mesh::make_empty_mesh();
      auto v1                = m_mesh->create_vertex(1, 1);
      auto v2                = m_mesh->create_vertex(2, 1);
      auto v3                = m_mesh->create_vertex(2, 2);
      auto v4                = m_mesh->create_vertex(1, 2);
      mesh::MeshEntityPtr e1 = m_mesh->create_edge(v2, v1);
      mesh::MeshEntityPtr e2 = m_mesh->create_edge(v3, v2);
      mesh::MeshEntityPtr e3 = m_mesh->create_edge(v4, v3);
      mesh::MeshEntityPtr e4 = m_mesh->create_edge(v1, v4);

      m_elReversed = m_mesh->create_quad(e1, e2, e3, e4, mesh::EntityOrientation::Reversed);
    }

    void test_get_edge_xi(const utils::Point& ptInterior, const utils::Point& ptExterior, PointClassification type,
                          int id, double xiForStandardEdge = -1)
    {
      test_get_edge_xi(m_el, ptInterior, ptExterior, type, id, xiForStandardEdge, xiForStandardEdge != -1);
      test_get_edge_xi(m_elReversed, ptInterior, ptExterior, type, id, 1 - xiForStandardEdge, xiForStandardEdge != -1);
    }

    void test_get_edge_xi_corner_vert(const utils::Point& pt1, const utils::Point& pt2, int id)
    {
      test_get_edge_xi_corner_vert(m_el, pt1, pt2, id);
      test_get_edge_xi_corner_vert(m_elReversed, pt1, pt2, id);
    }

    void test_get_edge_xi_corner_vert_vert(const utils::Point& pt1, const utils::Point& pt2, int id1, int id2)
    {
      test_get_edge_xi_corner_vert_vert(m_el, pt1, pt2, id1, id2);
      test_get_edge_xi_corner_vert_vert(m_elReversed, pt1, pt2, id1, id2);
    }

    void test_get_edge_xi_corner_vert_edge(const utils::Point& pt1, const utils::Point& pt2, int idVert, int idEdge,
                                           double xiEdge)
    {
      test_get_edge_xi_corner_vert_edge(m_el, pt1, pt2, idVert, idEdge, xiEdge);
      test_get_edge_xi_corner_vert_edge(m_elReversed, pt1, pt2, idVert, idEdge, 1 - xiEdge);
    }

    void test_get_edge_xi_corner_edge_edge(const utils::Point& pt1, const utils::Point& pt2, int id1, int id2,
                                           double edgeXi1, double edgeXi2)
    {
      test_get_edge_xi_corner_edge_edge(m_el, pt1, pt2, id1, id2, edgeXi1, edgeXi2);
      test_get_edge_xi_corner_edge_edge(m_elReversed, pt1, pt2, id1, id2, 1 - edgeXi1, 1 - edgeXi2);
    }

    void test_get_edge_intersection_xi_vert_vert(const utils::Point& ptVert, const utils::Point& ptExterior,
                                                 int otherVertId)
    {
      test_get_edge_intersection_xi_vert_vert(m_el, ptVert, ptExterior, otherVertId);
      test_get_edge_intersection_xi_vert_vert(m_elReversed, ptVert, ptExterior, otherVertId);
    }

    void test_get_edge_intersection_xi_vert_edge(const utils::Point& ptVert, const utils::Point& ptExterior,
                                                 int otherVertId)
    {
      test_get_edge_intersection_xi_vert_edge(m_el, ptVert, ptExterior, otherVertId);
      test_get_edge_intersection_xi_vert_edge(m_elReversed, ptVert, ptExterior, otherVertId);
    }

    void test_get_edge_intersection_xi_edge_edge(const utils::Point& ptVert, const utils::Point& ptExterior,
                                                 int otherEdgeId, double otherEdgeXi)
    {
      test_get_edge_intersection_xi_edge_edge(m_el, ptVert, ptExterior, otherEdgeId, otherEdgeXi);
      test_get_edge_intersection_xi_edge_edge(m_elReversed, ptVert, ptExterior, otherEdgeId, 1 - otherEdgeXi);
    }

  private:
    void test_get_edge_xi(mesh::MeshEntityPtr el, const utils::Point& ptInterior, const utils::Point& ptExterior,
                          PointClassification type, int id, double edgeXi, bool testXi)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, ptInterior);
      EXPECT_EQ(r1.type, PointClassification::Interior);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, ptExterior);
      EXPECT_EQ(r2.type, PointClassification::Exterior);
      PossibleEdgeIntersection p = m_classifier.get_edge_xi(r1, r2);
      EXPECT_EQ(p.record.type, type);
      EXPECT_EQ(p.record.id, id);
      if (testXi)
      {
        EXPECT_NEAR(p.edgeXi, edgeXi, 1e-13);
      }
    }

    void test_get_edge_xi_corner_vert(mesh::MeshEntityPtr el, const utils::Point& pt1, const utils::Point& pt2, int id)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, pt1);
      EXPECT_EQ(r1.type, PointClassification::Exterior);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, pt2);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      CornerRecord corner = m_classifier.get_edge_xi_corner(r1, r2);

      EXPECT_TRUE(corner.hasIntersection);
      EXPECT_EQ(corner.record1.type, PointClassification::Vert);
      EXPECT_EQ(corner.record1.id, id);
    }

    void test_get_edge_xi_corner_vert_vert(mesh::MeshEntityPtr el, const utils::Point& pt1, const utils::Point& pt2,
                                           int id1, int id2)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, pt1);
      EXPECT_EQ(r1.type, PointClassification::Exterior);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, pt2);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      CornerRecord corner = m_classifier.get_edge_xi_corner(r1, r2);

      predicates::impl::PointRecord& rFirst  = corner.record1.id == id1 ? corner.record1 : corner.record2;
      predicates::impl::PointRecord& rSecond = corner.record1.id == id1 ? corner.record2 : corner.record1;

      EXPECT_TRUE(corner.hasIntersection);
      EXPECT_EQ(rFirst.type, PointClassification::Vert);
      EXPECT_EQ(rFirst.id, id1);
      EXPECT_EQ(rSecond.type, PointClassification::Vert);
      EXPECT_EQ(rSecond.id, id2);
    }

    void test_get_edge_xi_corner_edge_edge(mesh::MeshEntityPtr el, const utils::Point& pt1, const utils::Point& pt2,
                                           int id1, int id2, double edgeXi1, double edgeXi2)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, pt1);
      EXPECT_EQ(r1.type, PointClassification::Exterior);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, pt2);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      CornerRecord corner = m_classifier.get_edge_xi_corner(r1, r2);

      predicates::impl::PointRecord& rFirst  = corner.record1.id == id1 ? corner.record1 : corner.record2;
      predicates::impl::PointRecord& rSecond = corner.record1.id == id1 ? corner.record2 : corner.record1;
      double xi1                             = corner.record1.id == id1 ? corner.xi1 : corner.xi2;
      double xi2                             = corner.record1.id == id1 ? corner.xi2 : corner.xi1;

      EXPECT_TRUE(corner.hasIntersection);
      EXPECT_EQ(rFirst.type, PointClassification::Edge);
      EXPECT_EQ(rFirst.id, id1);
      EXPECT_NEAR(xi1, edgeXi1, 1e-13);
      EXPECT_EQ(rSecond.type, PointClassification::Edge);
      EXPECT_EQ(rSecond.id, id2);
      EXPECT_NEAR(xi2, edgeXi2, 1e-13);
    }

    void test_get_edge_xi_corner_vert_edge(mesh::MeshEntityPtr el, const utils::Point& pt1, const utils::Point& pt2,
                                           int idVert, int idEdge, double xiEdge)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, pt1);
      EXPECT_EQ(r1.type, PointClassification::Exterior);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, pt2);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      CornerRecord corner = m_classifier.get_edge_xi_corner(r1, r2);

      predicates::impl::PointRecord& rVert =
          corner.record1.type == PointClassification::Vert ? corner.record1 : corner.record2;
      predicates::impl::PointRecord& rEdge =
          corner.record1.type == PointClassification::Vert ? corner.record2 : corner.record1;
      double xi = corner.record1.type == PointClassification::Vert ? corner.xi2 : corner.xi1;

      // std::cout << "r_vert = " << r_vert << std::endl;
      // std::cout << "r_edge = " << r_edge << std::endl;

      EXPECT_TRUE(corner.hasIntersection);
      EXPECT_EQ(rVert.type, PointClassification::Vert);
      EXPECT_EQ(rVert.id, idVert);
      EXPECT_EQ(rEdge.type, PointClassification::Edge);
      EXPECT_EQ(rEdge.id, idEdge);
      EXPECT_NEAR(xi, xiEdge, 1e-13);
    }

    void test_get_edge_intersection_xi_vert_vert(mesh::MeshEntityPtr el, const utils::Point& ptVert,
                                                 const utils::Point& ptExterior, int otherVertId)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, ptVert);
      EXPECT_EQ(r1.type, PointClassification::Vert);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, ptExterior);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      PossibleEdgeIntersection p = m_classifier.get_edge_intersection_xi(r1, r2);

      EXPECT_EQ(p.record.type, PointClassification::Vert);
      EXPECT_EQ(p.record.id, otherVertId);
    }

    void test_get_edge_intersection_xi_vert_edge(mesh::MeshEntityPtr el, const utils::Point& ptVert,
                                                 const utils::Point& ptExterior, int otherVertId)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, ptVert);
      EXPECT_EQ(r1.type, PointClassification::Edge);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, ptExterior);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      PossibleEdgeIntersection p = m_classifier.get_edge_intersection_xi(r1, r2);

      EXPECT_EQ(p.record.type, PointClassification::Vert);
      EXPECT_EQ(p.record.id, otherVertId);
    }

    void test_get_edge_intersection_xi_edge_edge(mesh::MeshEntityPtr el, const utils::Point& ptVert,
                                                 const utils::Point& ptExterior, int otherEdgeId, double otherEdgeXi)
    {
      predicates::impl::PointRecord r1 = m_classifier.classify(el, ptVert);
      EXPECT_EQ(r1.type, PointClassification::Edge);

      predicates::impl::PointRecord r2 = m_classifier.classify(el, ptExterior);
      EXPECT_EQ(r2.type, PointClassification::Exterior);

      PossibleEdgeIntersection p = m_classifier.get_edge_intersection_xi(r1, r2);

      EXPECT_EQ(p.record.type, PointClassification::Edge);
      EXPECT_EQ(p.record.id, otherEdgeId);
      EXPECT_NEAR(p.edgeXi, otherEdgeXi, 1e-13);
    }

    std::shared_ptr<mesh::Mesh> m_mesh;
    mesh::MeshEntityPtr m_el;

    std::shared_ptr<mesh::Mesh> m_meshReversed;
    mesh::MeshEntityPtr m_elReversed;
    PointClassifier m_classifier;
};

} // namespace

TEST_F(QuadIntersectionTester, get_edge_xiTriangle1)
{
  // pass through verts 1 to 4
  test_get_edge_xi({1.25, 1.5}, {0, -1}, PointClassification::Vert, 0);
  test_get_edge_xi({1.25, 1.5}, {3, 1.0 / 3.0}, PointClassification::Vert, 1);
  test_get_edge_xi({1.25, 1.5}, {5, 4}, PointClassification::Vert, 2);
  test_get_edge_xi({1.25, 1.5}, {0, 4}, PointClassification::Vert, 3);

  // pass through edges
  test_get_edge_xi({1.25, 1.5}, {1.25, -1}, PointClassification::Edge, 0, 0.25);
  test_get_edge_xi({1.25, 1.5}, {5, 0.25}, PointClassification::Edge, 1, 0.25);
  test_get_edge_xi({1.25, 1.5}, {3, 3.25}, PointClassification::Edge, 2, 0.25);
  test_get_edge_xi({1.25, 1.5}, {0, 2.75}, PointClassification::Edge, 3, 0.25);
}

TEST_F(QuadIntersectionTester, get_edge_xiTriangle2)
{
  // pass through verts 1 to 4
  test_get_edge_xi({1.75, 1.5}, {-0.5, 0}, PointClassification::Vert, 0);
  test_get_edge_xi({1.75, 1.5}, {3, -1}, PointClassification::Vert, 1);
  test_get_edge_xi({1.75, 1.5}, {3, 4}, PointClassification::Vert, 2);
  test_get_edge_xi({1.75, 1.5}, {-2, 4}, PointClassification::Vert, 3);

  // pass through edges
  test_get_edge_xi({1.75, 1.5}, {-1, -1.25}, PointClassification::Edge, 0, 0.25);
  test_get_edge_xi({1.75, 1.5}, {3, 0.25}, PointClassification::Edge, 1, 0.25);
  test_get_edge_xi({1.75, 1.5}, {1.75, 3}, PointClassification::Edge, 2, 0.25);
  test_get_edge_xi({1.75, 1.5}, {-1.0625, 2.4375}, PointClassification::Edge, 3, 0.25);
}

TEST_F(QuadIntersectionTester, get_edge_xiDiagonal)
{
  // pass through verts 1 to 4
  test_get_edge_xi({1.5, 1.5}, {0, 0}, PointClassification::Vert, 0);
  test_get_edge_xi({1.5, 1.5}, {3, 0}, PointClassification::Vert, 1);
  test_get_edge_xi({1.5, 1.5}, {3, 3}, PointClassification::Vert, 2);
  test_get_edge_xi({1.5, 1.5}, {0, 3}, PointClassification::Vert, 3);

  // pass through edges
  test_get_edge_xi({1.5, 1.5}, {0, -1.5}, PointClassification::Edge, 0, 0.25);
  test_get_edge_xi({1.5, 1.5}, {3, 0.75}, PointClassification::Edge, 1, 0.25);
  test_get_edge_xi({1.5, 1.5}, {3, 4.5}, PointClassification::Edge, 2, 0.25);
  test_get_edge_xi({1.5, 1.5}, {0, 2.25}, PointClassification::Edge, 3, 0.25);
}

TEST_F(QuadIntersectionTester, get_edge_xi_cornerVertIntersection)
{
  // test cases where the line intersects exactly one vertex on the element

  test_get_edge_xi_corner_vert({1.5, 0}, {0, 3}, 0);
  test_get_edge_xi_corner_vert({1.5, 0}, {3, 3}, 1);

  test_get_edge_xi_corner_vert({3, 1.5}, {1, 0.5}, 1);
  test_get_edge_xi_corner_vert({3, 1.5}, {1, 2.5}, 2);

  test_get_edge_xi_corner_vert({1.5, 3}, {3, 0}, 2);
  test_get_edge_xi_corner_vert({1.5, 3}, {0, 0}, 3);

  test_get_edge_xi_corner_vert({0.5, 1.5}, {3, 4}, 3);
  test_get_edge_xi_corner_vert({0.5, 1.5}, {3, -1}, 0);
}

TEST_F(QuadIntersectionTester, get_edge_xi_cornerVertVertIntersection)
{
  // test cases where the line passes through two verts

  // edge overlaps
  test_get_edge_xi_corner_vert_vert({0, 1}, {3, 1}, 0, 1);
  test_get_edge_xi_corner_vert_vert({2, 0}, {2, 3}, 1, 2);
  test_get_edge_xi_corner_vert_vert({0, 2}, {3, 2}, 2, 3);
  test_get_edge_xi_corner_vert_vert({1, 0}, {1, 3}, 0, 3);

  // diagonals
  test_get_edge_xi_corner_vert_vert({0, 0}, {3, 3}, 0, 2);
  test_get_edge_xi_corner_vert_vert({0, 3}, {3, 0}, 1, 3);
}

TEST_F(QuadIntersectionTester, get_edge_xi_cornerVertEdgeIntersection)
{
  test_get_edge_xi_corner_vert_edge({1.5, 0}, {3, 6}, 2, 0, 0.75);
  test_get_edge_xi_corner_vert_edge({1.5, 0}, {0, 6}, 3, 0, 0.25);

  test_get_edge_xi_corner_vert_edge({3, 1.5}, {0, 0.75}, 0, 1, 0.25);
  test_get_edge_xi_corner_vert_edge({3, 1.5}, {0, 2.25}, 3, 1, 0.75);

  test_get_edge_xi_corner_vert_edge({1.5, 3}, {0, -3}, 0, 2, 0.75);
  test_get_edge_xi_corner_vert_edge({1.5, 3}, {3, -3}, 1, 2, 0.25);

  test_get_edge_xi_corner_vert_edge({0, 1.5}, {3, 0.75}, 1, 3, 0.75);
  test_get_edge_xi_corner_vert_edge({0, 1.5}, {3, 2.25}, 2, 3, 0.25);
}

TEST_F(QuadIntersectionTester, get_edge_xi_cornerEdgeEdgeIntesection)
{
  test_get_edge_xi_corner_edge_edge({0, 4.75}, {3, -4.25}, 0, 3, 0.25, 0.25);
  test_get_edge_xi_corner_edge_edge({1.25, 0}, {1.25, 3}, 0, 2, 0.25, 0.75);
  test_get_edge_xi_corner_edge_edge({0, -4.25}, {3, 4.75}, 0, 1, 0.75, 0.75);

  test_get_edge_xi_corner_edge_edge({0, -0.25}, {3, 2.75}, 1, 0, 0.75, 0.25);
  test_get_edge_xi_corner_edge_edge({0, 1.25}, {3, 1.25}, 1, 3, 0.25, 0.75);
  test_get_edge_xi_corner_edge_edge({0, 3.25}, {3, 0.25}, 1, 2, 0.25, 0.75);

  test_get_edge_xi_corner_edge_edge({0, -1.75}, {3, 7.25}, 2, 3, 0.75, 0.75);
  test_get_edge_xi_corner_edge_edge({1.25, 0}, {1.25, 3}, 2, 0, 0.75, 0.25);
  test_get_edge_xi_corner_edge_edge({0, 3.75}, {3, 0.75}, 2, 1, 0.25, 0.75);

  test_get_edge_xi_corner_edge_edge({0, 2.75}, {3, -0.25}, 3, 0, 0.25, 0.75);
  test_get_edge_xi_corner_edge_edge({0, 1.25}, {3, 1.25}, 3, 1, 0.75, 0.25);
  test_get_edge_xi_corner_edge_edge({0, 0.75}, {3, 3.75}, 3, 2, 0.25, 0.75);
}

TEST_F(QuadIntersectionTester, get_edge_xiIntersectionVertVert)
{
  test_get_edge_intersection_xi_vert_vert({1, 1}, {3, 1}, 1);
  test_get_edge_intersection_xi_vert_vert({1, 1}, {3, 3}, 2);
  test_get_edge_intersection_xi_vert_vert({1, 1}, {1, 3}, 3);

  test_get_edge_intersection_xi_vert_vert({2, 1}, {0, 1}, 0);
  test_get_edge_intersection_xi_vert_vert({2, 1}, {2, 3}, 2);
  test_get_edge_intersection_xi_vert_vert({2, 1}, {0, 3}, 3);

  test_get_edge_intersection_xi_vert_vert({2, 2}, {0, 0}, 0);
  test_get_edge_intersection_xi_vert_vert({2, 2}, {2, 0}, 1);
  test_get_edge_intersection_xi_vert_vert({2, 2}, {0, 2}, 3);

  test_get_edge_intersection_xi_vert_vert({1, 2}, {1, 0}, 0);
  test_get_edge_intersection_xi_vert_vert({1, 2}, {3, 0}, 1);
  test_get_edge_intersection_xi_vert_vert({1, 2}, {3, 2}, 2);
}

TEST_F(QuadIntersectionTester, get_edge_xiIntersectionVertEdge)
{
  test_get_edge_intersection_xi_vert_edge({1.25, 1}, {0, 1}, 0);
  test_get_edge_intersection_xi_vert_edge({1.25, 1}, {3, 1}, 1);
  test_get_edge_intersection_xi_vert_edge({1.25, 1}, {5, 6}, 2);
  test_get_edge_intersection_xi_vert_edge({1.25, 1}, {0, 6}, 3);

  test_get_edge_intersection_xi_vert_edge({2, 1.25}, {0, 0.75}, 0);
  test_get_edge_intersection_xi_vert_edge({2, 1.25}, {2, 0}, 1);
  test_get_edge_intersection_xi_vert_edge({2, 1.25}, {2, 3}, 2);
  test_get_edge_intersection_xi_vert_edge({2, 1.25}, {0, 2.75}, 3);

  test_get_edge_intersection_xi_vert_edge({1.25, 2}, {0, -3}, 0);
  test_get_edge_intersection_xi_vert_edge({1.25, 2}, {5, -3}, 1);
  test_get_edge_intersection_xi_vert_edge({1.25, 2}, {3, 2}, 2);
  test_get_edge_intersection_xi_vert_edge({1.25, 2}, {0, 2}, 3);

  test_get_edge_intersection_xi_vert_edge({1, 1.25}, {1, 0}, 0);
  test_get_edge_intersection_xi_vert_edge({1, 1.25}, {3, 0.75}, 1);
  test_get_edge_intersection_xi_vert_edge({1, 1.25}, {3, 2.75}, 2);
  test_get_edge_intersection_xi_vert_edge({1, 1.25}, {1, 3}, 3);
}

TEST_F(QuadIntersectionTester, get_edge_xiIntersectionEdgeEdge)
{
  test_get_edge_intersection_xi_edge_edge({1.25, 1}, {3, 2.75}, 1, 0.75);
  test_get_edge_intersection_xi_edge_edge({1.25, 1}, {1.25, 3}, 2, 0.75);
  test_get_edge_intersection_xi_edge_edge({1.25, 1}, {0, 4.75}, 3, 0.25);

  test_get_edge_intersection_xi_edge_edge({2, 1.25}, {0, -0.75}, 0, 0.75);
  test_get_edge_intersection_xi_edge_edge({2, 1.25}, {0, 3.25}, 2, 0.75);
  test_get_edge_intersection_xi_edge_edge({2, 1.25}, {0, 1.25}, 3, 0.75);

  test_get_edge_intersection_xi_edge_edge({1.25, 2}, {1.25, 0}, 0, 0.25);
  test_get_edge_intersection_xi_edge_edge({1.25, 2}, {3, 0.25}, 1, 0.25);
  test_get_edge_intersection_xi_edge_edge({1.25, 2}, {0, -1.75}, 3, 0.75);

  test_get_edge_intersection_xi_edge_edge({1, 1.25}, {3, -0.75}, 0, 0.25);
  test_get_edge_intersection_xi_edge_edge({1, 1.25}, {3, 1.25}, 1, 0.25);
  test_get_edge_intersection_xi_edge_edge({1, 1.25}, {3, 7.25}, 2, 0.75);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
