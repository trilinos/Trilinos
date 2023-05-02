#include "stk_middle_mesh/predicates/intersection_common.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/predicates/quad_to_triangles.hpp"
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


TEST_F(QuadToTrianglesTester, QuadXiCoordinatesVert)
{
  {
    PointRecordForTriangle r1(PointClassification::Vert, 0, quadToTriangles.el1),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0, 0}, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Vert, 1, quadToTriangles.el1),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {1, 0}, 1e-13);
  }  

  {
    PointRecordForTriangle r1(PointClassification::Exterior, 0, quadToTriangles.el1),
                           r2(PointClassification::Vert, 2, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {1, 0}, 1e-13);
  }   

  {
    PointRecordForTriangle r1(PointClassification::Exterior, 0, quadToTriangles.el1),
                           r2(PointClassification::Vert, 0, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {1, 1}, 1e-13);
  } 

  {
    PointRecordForTriangle r1(PointClassification::Vert, 2, quadToTriangles.el1),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0, 1}, 1e-13);
  }  

  {
    PointRecordForTriangle r1(PointClassification::Exterior, 0, quadToTriangles.el1),
                           r2(PointClassification::Vert, 1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0, 1}, 1e-13);
  }     
}

TEST_F(QuadToTrianglesTester, QuadXiCoordinatesEdges)
{
  {
    PointRecordForTriangle r1(PointClassification::Edge, 0, quadToTriangles.el1, {0.5, 0}),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.5, 0}, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1),
                           r2(PointClassification::Edge, 2, quadToTriangles.el2, {0, 0.5});
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {1, 0.5}, 1e-13);
  } 

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1),
                           r2(PointClassification::Edge, 0, quadToTriangles.el2, {0.5, 0});
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.5, 1}, 1e-13);
  }  

  {
    PointRecordForTriangle r1(PointClassification::Edge, 2, quadToTriangles.el1, {0, 0.5}),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.0, 0.5}, 1e-13);
  }    
}

TEST_F(QuadToTrianglesTester, QuadXiCoordinatesInteriorEdge)
{
  {
    PointRecordForTriangle r1(PointClassification::Edge, 1, quadToTriangles.el1, {0.5, 0.5}),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.5, 0.5}, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1),
                           r2(PointClassification::Edge, 1, quadToTriangles.el2, {0.5, 0.5});
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.5, 0.5}, 1e-13);
  }  
}

TEST_F(QuadToTrianglesTester, QuadXiCoordinatesInterior)
{
  {
    PointRecordForTriangle r1(PointClassification::Interior, -1, quadToTriangles.el1, {0.25, 0.25}),
                           r2(PointClassification::Exterior, -1, quadToTriangles.el2);
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.25, 0.25}, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1),
                           r2(PointClassification::Interior, -1, quadToTriangles.el2, {0.25, 0.25});
    PointRecord quadRecord = quadToTriangles.get_quad_record(el, r1, r2);
    expect_near(quadToTriangles.get_quad_xi_coords(quadRecord), {0.75, 0.75}, 1e-13);
  }  
}

TEST_F(QuadToTrianglesTester, create_record_edge)
{
  PointRecord record = quadToTriangles.create_record(el, 0, 0);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 0);
  expect_near(quadToTriangles.get_quad_xi_coords(record), {0, 0}, 1e-13);

  record = quadToTriangles.create_record(el, 0, 1);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 0);  
  expect_near(quadToTriangles.get_quad_xi_coords(record), {1, 0}, 1e-13);

  record = quadToTriangles.create_record(el, 1, 0);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 1);  
  expect_near(quadToTriangles.get_quad_xi_coords(record), {1, 0}, 1e-13);  

  record = quadToTriangles.create_record(el, 1, 1);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 1);  
  expect_near(quadToTriangles.get_quad_xi_coords(record), {1, 1}, 1e-13);  

  
  record = quadToTriangles.create_record(el, 2, 0);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 2);  
  expect_near(quadToTriangles.get_quad_xi_coords(record), {1, 1}, 1e-13);  

  record = quadToTriangles.create_record(el, 2, 1);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 2);
  expect_near(quadToTriangles.get_quad_xi_coords(record), {0, 1}, 1e-13);   

  record = quadToTriangles.create_record(el, 3, 0);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 3);  
  expect_near(quadToTriangles.get_quad_xi_coords(record), {0, 1}, 1e-13);  

  record = quadToTriangles.create_record(el, 3, 1);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Edge);
  EXPECT_EQ(record.id, 3);  
  expect_near(quadToTriangles.get_quad_xi_coords(record), {0, 0}, 1e-13);   

}

TEST_F(QuadToTrianglesTester, create_record_vert)
{
  PointRecord record = quadToTriangles.create_record(el, 0);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Vert);
  EXPECT_EQ(record.id, 0);
  expect_near(quadToTriangles.get_quad_xi_coords(record), {0, 0}, 1e-13);  

  record = quadToTriangles.create_record(el, 1);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Vert);
  EXPECT_EQ(record.id, 1);
  expect_near(quadToTriangles.get_quad_xi_coords(record), {1, 0}, 1e-13);    

  record = quadToTriangles.create_record(el, 2);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Vert);
  EXPECT_EQ(record.id, 2);
  expect_near(quadToTriangles.get_quad_xi_coords(record), {1, 1}, 1e-13);   

  record = quadToTriangles.create_record(el, 3);
  EXPECT_EQ(record.el, el);
  EXPECT_EQ(record.type, PointClassification::Vert);
  EXPECT_EQ(record.id, 3);
  expect_near(quadToTriangles.get_quad_xi_coords(record), {0, 1}, 1e-13);      
}


TEST(QuadToTriangles, classify_onto)
{
  auto mesh     = mesh::make_empty_mesh();
  auto vert0 = mesh->create_vertex(0, 0, 0);
  auto vert1 = mesh->create_vertex(1, 0, 0);
  auto vert2 = mesh->create_vertex(1, 1, 0);
  auto vert3 = mesh->create_vertex(0, 1, 0);
  auto vert4 = mesh->create_vertex(0, -1, 0);
  auto vert5 = mesh->create_vertex(1, -1, 0);
  auto vert6 = mesh->create_vertex(2, 0, 0);
  auto vert7 = mesh->create_vertex(2, 1, 0);
  auto vert8 = mesh->create_vertex(1, 2, 0);
  auto vert9 = mesh->create_vertex(0, 2, 0);
  auto vert10 = mesh->create_vertex(-1, 1, 0);
  auto vert11 = mesh->create_vertex(-1, 0, 0);
  auto el0       = mesh->create_quad_from_verts(vert0, vert1, vert2, vert3);
  auto el1       = mesh->create_quad_from_verts(vert0, vert4, vert5, vert1);
  auto el2       = mesh->create_quad_from_verts(vert1, vert6, vert7, vert2);
  auto el3       = mesh->create_quad_from_verts(vert3, vert2, vert8, vert9);
  auto el4       = mesh->create_quad_from_verts(vert0, vert3, vert10, vert11);
  QuadToTriangles quadToTriangles;
  quadToTriangles.set_triangles(el0);

  // vertex 1
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 0);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el1);

    EXPECT_EQ(record2.el, el1);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    quadToTriangles.set_triangles(el1);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0, 0, 0}, 1e-13);
  }

  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 0);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el4);

    EXPECT_EQ(record2.el, el4);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    quadToTriangles.set_triangles(el4);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0, 0, 0}, 1e-13);
  }

  // vertex 2
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 1);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el1);

    EXPECT_EQ(record2.el, el1);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 3);
    //quadToTriangles.compute_xyz_coords(record2)  
    quadToTriangles.set_triangles(el1);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {1, 0, 0}, 1e-13);
  }

  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 1);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el2);

    EXPECT_EQ(record2.el, el2);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    quadToTriangles.set_triangles(el2);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {1, 0, 0}, 1e-13);
  }  

  // vertex 3
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 2);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el2);

    EXPECT_EQ(record2.el, el2);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 3);
    quadToTriangles.set_triangles(el2);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {1, 1, 0}, 1e-13);
  }


  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 2);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el3);

    EXPECT_EQ(record2.el, el3);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 1);
    quadToTriangles.set_triangles(el3);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {1, 1, 0}, 1e-13);
  }  

  // vertex 4
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 3);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el3);

    EXPECT_EQ(record2.el, el3);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 0);
    quadToTriangles.set_triangles(el3);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0, 1, 0}, 1e-13);
  }

  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 3);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el4);

    EXPECT_EQ(record2.el, el4);
    EXPECT_EQ(record2.type, PointClassification::Vert);
    EXPECT_EQ(record2.id, 1);
    quadToTriangles.set_triangles(el4);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0, 1, 0}, 1e-13);
  } 

  // edge 0
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 0, 0.25);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el1);

    EXPECT_EQ(record2.el, el1);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 3);
    quadToTriangles.set_triangles(el1);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0.25, 0, 0}, 1e-13);
  }

  // edge 1
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 1, 0.25);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el2);

    EXPECT_EQ(record2.el, el2);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 3);
    quadToTriangles.set_triangles(el2);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {1, 0.25, 0}, 1e-13);
  }  

  // edge 2
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 2, 0.25);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el3);

    EXPECT_EQ(record2.el, el3);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 0);
    quadToTriangles.set_triangles(el3);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0.75, 1, 0}, 1e-13);
  }   

  // edge 3
  {
    quadToTriangles.set_triangles(el0);  
    PointRecord record1 = quadToTriangles.create_record(el0, 3, 0.25);
    PointRecord record2 = quadToTriangles.classify_onto(record1, el4);

    EXPECT_EQ(record2.el, el4);
    EXPECT_EQ(record2.type, PointClassification::Edge);
    EXPECT_EQ(record2.id, 0);
    quadToTriangles.set_triangles(el4);  
    expect_near(quadToTriangles.compute_xyz_coords(record2), {0, 0.75, 0}, 1e-13);
  }      

}

TEST_F(QuadToTrianglesTester, compute_xyz_coords_exterior)
{
  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1, utils::Point(0.5, -0.5));
    PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2, utils::Point(0.5, -1.5));
    PointRecord record(PointClassification::Exterior, -1, el, r1, r2);
    expect_near(quadToTriangles.compute_xyz_coords(record, true), {0.5, -0.5}, 1e-13);
  }

  {
    PointRecordForTriangle r1(PointClassification::Exterior, -1, quadToTriangles.el1, utils::Point(0.5,  1.5));
    PointRecordForTriangle r2(PointClassification::Exterior, -1, quadToTriangles.el2, utils::Point(0.5, -0.5));
    PointRecord record(PointClassification::Exterior, -1, el, r1, r2);
    expect_near(quadToTriangles.compute_xyz_coords(record, true), {0.5, 1.5}, 1e-13);
  }
}



} // namespace impl
} // namespace middle_mesh
} // namespace stk
