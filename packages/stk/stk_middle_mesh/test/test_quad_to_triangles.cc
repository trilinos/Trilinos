#include "mesh.h"
#include "predicates/quad_to_triangles.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

namespace {

class QuadToTrianglesTester : public ::testing::Test
{
  protected:
    QuadToTrianglesTester()
    {
      mesh     = mesh::make_empty_mesh();
      verts[0] = mesh->create_vertex(0, 0, 0);
      verts[1] = mesh->create_vertex(1, 0, 0);
      verts[2] = mesh->create_vertex(1, 1, 0);
      verts[3] = mesh->create_vertex(0, 1, 0);
      el       = mesh->create_quad_from_verts(verts[0], verts[1], verts[2], verts[3]);
      quadToTriangles.set_triangles(el);
    }

    std::shared_ptr<mesh::Mesh> mesh;
    mesh::MeshEntityPtr el;
    std::array<mesh::MeshEntityPtr, 4> verts;
    QuadToTriangles quadToTriangles;
};

void expect_near(const utils::Point& pt1, const utils::Point& pt2, double tol)
{
  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(pt1[i], pt2[i], tol);
}
} // namespace

TEST_F(QuadToTrianglesTester, Setup)
{
  for (int i = 0; i < 4; ++i)
    expect_near(verts[0]->get_point_orig(0), quadToTriangles.verts[0]->get_point_orig(0), 1e-13);

  std::array<mesh::MeshEntityPtr, 3> el1Verts, el2Verts;
  get_downward(quadToTriangles.el1, 0, el1Verts.data());
  get_downward(quadToTriangles.el2, 0, el2Verts.data());
  for (int i = 0; i < 3; ++i)
  {
    mesh::MeshEntityPtr vertFromMap = quadToTriangles.verts[quadToTriangles.VERTMAP_TRI1_TO_QUAD[i]];
    mesh::MeshEntityPtr vertFromEl  = el1Verts[i];
    EXPECT_EQ(vertFromMap, vertFromEl);

    vertFromMap = quadToTriangles.verts[quadToTriangles.VERTMAP_TRI2_TO_QUAD[i]];
    vertFromEl  = el2Verts[i];
    EXPECT_EQ(vertFromMap, vertFromEl);
  }
}

TEST_F(QuadToTrianglesTester, Vert1)
{
  PointRecordForTriangle r1(PointClassification::Vert, 0, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Vert);
  EXPECT_EQ(r3.id, 0);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Vert2FirstEl)
{
  PointRecordForTriangle r1(PointClassification::Vert, 1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Vert);
  EXPECT_EQ(r3.id, 1);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Vert2SecondEl)
{
  PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Vert, 2, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Vert);
  EXPECT_EQ(r3.id, 1);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Vert3)
{
  PointRecordForTriangle r1(PointClassification::Exterior, 0, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Vert, 0, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Vert);
  EXPECT_EQ(r3.id, 2);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Vert4FirstEl)
{
  PointRecordForTriangle r1(PointClassification::Vert, 2, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Vert);
  EXPECT_EQ(r3.id, 3);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Vert4SecondEl)
{
  PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Vert, 1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Vert);
  EXPECT_EQ(r3.id, 3);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Edge1)
{
  PointRecordForTriangle r1(PointClassification::Edge, 0, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Edge);
  EXPECT_EQ(r3.id, 0);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Edge2)
{
  PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Edge, 2, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Edge);
  EXPECT_EQ(r3.id, 1);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Edge3)
{
  PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Edge, 0, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Edge);
  EXPECT_EQ(r3.id, 2);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Edge4)
{
  PointRecordForTriangle r1(PointClassification::Edge, 2, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Edge);
  EXPECT_EQ(r3.id, 3);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, InteriorEdge1)
{
  PointRecordForTriangle r1(PointClassification::Edge, 1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Interior);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, InteriorEdge2)
{
  PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Edge, 1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Interior);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Interior1)
{
  PointRecordForTriangle r1(PointClassification::Interior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Interior);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, Interior2)
{
  PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1);
  PointRecordForTriangle r2(PointClassification::Interior, -1, quadToTriangles.el2);
  PointRecord r3 = quadToTriangles.get_quad_record(el, r1, r2);

  EXPECT_EQ(r3.type, PointClassification::Interior);
  EXPECT_EQ(r3.el, el);
}

TEST_F(QuadToTrianglesTester, RecordConsistency_El1Vert1)
{
  PointRecordForTriangle r1(PointClassification::Vert, 1, quadToTriangles.el1);

  std::vector<PointClassification> types = {PointClassification::Interior, PointClassification::Exterior,
                                            PointClassification::Edge, PointClassification::Edge};
  std::vector<int> ids                   = {-1, -1, 1, 2};

  for (size_t i = 0; i < types.size(); ++i)
  {
    PointRecordForTriangle r2(types[i], ids[i], quadToTriangles.el2);
    quadToTriangles.enforce_record_consistency(r1, r2);

    EXPECT_EQ(r1.type, PointClassification::Vert);
    EXPECT_EQ(r1.id, 1);
    EXPECT_EQ(r1.el, quadToTriangles.el1);

    EXPECT_EQ(r2.type, PointClassification::Vert);
    EXPECT_EQ(r2.id, 2);
    EXPECT_EQ(r2.el, quadToTriangles.el2);
  }
}

TEST_F(QuadToTrianglesTester, RecordConsistency_El1Vert2)
{
  PointRecordForTriangle r1(PointClassification::Vert, 2, quadToTriangles.el1);

  std::vector<PointClassification> types = {PointClassification::Interior, PointClassification::Exterior,
                                            PointClassification::Edge, PointClassification::Edge};
  std::vector<int> ids                   = {-1, -1, 0, 1};

  for (size_t i = 0; i < types.size(); ++i)
  {
    PointRecordForTriangle r2(types[i], ids[i], quadToTriangles.el2);
    quadToTriangles.enforce_record_consistency(r1, r2);

    EXPECT_EQ(r1.type, PointClassification::Vert);
    EXPECT_EQ(r1.id, 2);
    EXPECT_EQ(r1.el, quadToTriangles.el1);

    EXPECT_EQ(r2.type, PointClassification::Vert);
    EXPECT_EQ(r2.id, 1);
    EXPECT_EQ(r2.el, quadToTriangles.el2);
  }
}

TEST_F(QuadToTrianglesTester, RecordConsistency_El1Edge1)
{
  PointRecordForTriangle r1(PointClassification::Edge, 1, quadToTriangles.el1);

  std::vector<PointClassification> types = {PointClassification::Interior, PointClassification::Exterior};
  std::vector<int> ids                   = {-1, -1};

  for (size_t i = 0; i < types.size(); ++i)
  {
    PointRecordForTriangle r2(types[i], ids[i], quadToTriangles.el2);
    quadToTriangles.enforce_record_consistency(r1, r2);

    EXPECT_EQ(r1.type, PointClassification::Edge);
    EXPECT_EQ(r1.id, 1);
    EXPECT_EQ(r1.el, quadToTriangles.el1);

    EXPECT_EQ(r2.type, PointClassification::Edge);
    EXPECT_EQ(r2.id, 1);
    EXPECT_EQ(r2.el, quadToTriangles.el2);
  }
}

TEST_F(QuadToTrianglesTester, RecordConsistency_El2Vert1)
{
  PointRecordForTriangle r2(PointClassification::Vert, 1, quadToTriangles.el2);

  std::vector<PointClassification> types = {PointClassification::Interior, PointClassification::Exterior,
                                            PointClassification::Edge, PointClassification::Edge};
  std::vector<int> ids                   = {-1, -1, 1, 2};

  for (size_t i = 0; i < types.size(); ++i)
  {
    PointRecordForTriangle r1(types[i], ids[i], quadToTriangles.el1);
    quadToTriangles.enforce_record_consistency(r1, r2);

    EXPECT_EQ(r1.type, PointClassification::Vert);
    EXPECT_EQ(r1.id, 2);
    EXPECT_EQ(r1.el, quadToTriangles.el1);

    EXPECT_EQ(r2.type, PointClassification::Vert);
    EXPECT_EQ(r2.id, 1);
    EXPECT_EQ(r2.el, quadToTriangles.el2);
  }
}

TEST_F(QuadToTrianglesTester, RecordConsistency_El2Vert2)
{
  PointRecordForTriangle r2(PointClassification::Vert, 2, quadToTriangles.el2);

  std::vector<PointClassification> types = {PointClassification::Interior, PointClassification::Exterior,
                                            PointClassification::Edge, PointClassification::Edge};
  std::vector<int> ids                   = {-1, -1, 0, 1};

  for (size_t i = 0; i < types.size(); ++i)
  {
    PointRecordForTriangle r1(types[i], ids[i], quadToTriangles.el1);
    quadToTriangles.enforce_record_consistency(r1, r2);

    EXPECT_EQ(r1.type, PointClassification::Vert);
    EXPECT_EQ(r1.id, 1);
    EXPECT_EQ(r1.el, quadToTriangles.el1);

    EXPECT_EQ(r2.type, PointClassification::Vert);
    EXPECT_EQ(r2.id, 2);
    EXPECT_EQ(r2.el, quadToTriangles.el2);
  }
}

TEST_F(QuadToTrianglesTester, RecordConsistency_El2Edge1)
{
  PointRecordForTriangle r2(PointClassification::Edge, 1, quadToTriangles.el2);

  std::vector<PointClassification> types = {PointClassification::Interior, PointClassification::Exterior};
  std::vector<int> ids                   = {-1, -1};

  for (size_t i = 0; i < types.size(); ++i)
  {
    PointRecordForTriangle r1(types[i], ids[i], quadToTriangles.el1);
    quadToTriangles.enforce_record_consistency(r1, r2);

    EXPECT_EQ(r1.type, PointClassification::Edge);
    EXPECT_EQ(r1.id, 1);
    EXPECT_EQ(r1.el, quadToTriangles.el1);

    EXPECT_EQ(r2.type, PointClassification::Edge);
    EXPECT_EQ(r2.id, 1);
    EXPECT_EQ(r2.el, quadToTriangles.el2);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
