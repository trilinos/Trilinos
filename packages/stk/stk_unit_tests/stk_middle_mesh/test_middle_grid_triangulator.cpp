#include "element_operations_2d.hpp"
#include "mesh_io.hpp"
#include "nonconformal4.hpp"
#include "util/nonconformal_interface_helpers.hpp"
#include "gtest/gtest.h"

namespace {

using smmPoint = stk::middle_mesh::utils::Point;
using smmMesh = stk::middle_mesh::mesh::Mesh;
using smmMeshEntityPtr = stk::middle_mesh::mesh::MeshEntityPtr;
template<typename T>
using smmFieldPtr = stk::middle_mesh::mesh::FieldPtr<T>;
template<typename T>
using smmVariableSizeFieldPtr = stk::middle_mesh::mesh::VariableSizeFieldPtr<T>;
using smmNonconformal4 = stk::middle_mesh::nonconformal4::impl::Nonconformal4;

class Nonconformal4Tester : public ::testing::Test
{
  protected:
    Nonconformal4Tester()
      : mMesh1(stk::middle_mesh::mesh::make_empty_mesh())
      , mMesh2(stk::middle_mesh::mesh::make_empty_mesh())
      , mMeshIn(stk::middle_mesh::mesh::make_empty_mesh())
    {}

    void create_mesh1(const std::vector<smmPoint>& verts, const std::vector<std::array<int, 4>>& elementVerts)
    {
      create_mesh(mMesh1, verts, elementVerts);
    }

    void create_mesh2(const std::vector<smmPoint>& verts, const std::vector<std::array<int, 4>>& elementVerts)
    {
      create_mesh(mMesh2, verts, elementVerts);
    }

    void run()
    {
      double eps = 1e-13;
      stk::middle_mesh::NormalProjectionOpts opts;
      opts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
      opts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);

      mNonconformalGenerator      = std::make_shared<smmNonconformal4>(mMesh1, mMesh2, opts);
      mMeshIn                     = mNonconformalGenerator->create();
      mMeshInElementsToMesh1      = mNonconformalGenerator->get_mesh1_classification();
      mMeshInElementsToMesh2      = mNonconformalGenerator->get_mesh2_classification();
      mMesh1InverseClassification = mNonconformalGenerator->compute_mesh1_inverse_classification();
      mMesh2InverseClassification = mNonconformalGenerator->compute_mesh2_inverse_classification();

      stk::middle_mesh::test_util::test_areas_positive(mMeshIn);
      stk::middle_mesh::test_util::test_total_areas_same(mMesh1, mMesh2, mMeshIn);
      stk::middle_mesh::test_util::test_every_element_classified(mMeshIn, mMeshInElementsToMesh1);
      stk::middle_mesh::test_util::test_every_element_classified(mMeshIn, mMeshInElementsToMesh2);
      stk::middle_mesh::test_util::test_every_element_classified_inverse(mMesh1, mMesh1InverseClassification);
      stk::middle_mesh::test_util::test_every_element_classified_inverse(mMesh2, mMesh2InverseClassification);
      stk::middle_mesh::test_util::test_area_per_element(mMesh1, mMesh1InverseClassification);
      stk::middle_mesh::test_util::test_area_per_element(mMesh2, mMesh2InverseClassification);
    }

    std::shared_ptr<smmMesh> mMesh1;
    std::shared_ptr<smmMesh> mMesh2;
    std::shared_ptr<smmNonconformal4> mNonconformalGenerator;
    std::shared_ptr<smmMesh> mMeshIn;
    smmFieldPtr<smmMeshEntityPtr> mMeshInElementsToMesh1;
    smmFieldPtr<smmMeshEntityPtr> mMeshInElementsToMesh2;
    smmVariableSizeFieldPtr<smmMeshEntityPtr> mMesh1InverseClassification;
    smmVariableSizeFieldPtr<smmMeshEntityPtr> mMesh2InverseClassification;

  private:
    void create_mesh(std::shared_ptr<smmMesh> mesh, const std::vector<smmPoint>& vertCoords,
                     const std::vector<std::array<int, 4>>& elementVerts)
    {
      for (auto& vert : vertCoords)
        mesh->create_vertex(vert);

      auto& verts = mesh->get_vertices();
      for (auto& elVertIds : elementVerts)
      {
        if (elVertIds[3] == -1)
        {
          auto v0 = verts[elVertIds[0]];
          auto v1 = verts[elVertIds[1]];
          auto v2 = verts[elVertIds[2]];
          mesh->create_triangle_from_verts(v0, v1, v2);
        } else
        {
          auto v0 = verts[elVertIds[0]];
          auto v1 = verts[elVertIds[1]];
          auto v2 = verts[elVertIds[2]];
          auto v3 = verts[elVertIds[3]];
          mesh->create_quad_from_verts(v0, v1, v2, v3);
        }
      }
    }
};

} // namespace

TEST_F(Nonconformal4Tester, El1ContainedInEl2)
{
  // No intersections.  Element 1 is completely inside element 2
  std::vector<smmPoint> mesh1Verts = {
      smmPoint(0, 0),         smmPoint(1, 0),        smmPoint(1, 1),       smmPoint(0, 1),
      smmPoint(-0.25, -0.25), smmPoint(1.25, -0.25), smmPoint(1.25, 1.25), smmPoint(-0.25, 1.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {3, 2, 6, 7}, {0, 3, 7, 4}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {mesh1Verts[4], mesh1Verts[5], mesh1Verts[6], mesh1Verts[7]};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  create_mesh2(mesh2Verts, mesh2Els);

  stk::middle_mesh::mesh::impl::print_vert_edges("mesh1", mMesh1);
  stk::middle_mesh::mesh::impl::print_vert_edges("mesh2", mMesh2);

  run();
  stk::middle_mesh::mesh::impl::print_vert_edges("mesh_in", mMeshIn);

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

TEST_F(Nonconformal4Tester, SingleEdgeOverlap)
{
  // The top edge of el1 overlaps with the top edge of el2.  No other
  // edges itnersect
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),    smmPoint(1, 0),         smmPoint(1, 1),
                                              smmPoint(0, 1),    smmPoint(-0.25, -0.25), smmPoint(1.25, -0.25),
                                              smmPoint(1.25, 1), smmPoint(-0.25, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {0, 3, 7, 4}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts = {smmPoint(-0.25, -0.25), smmPoint(1.25, -0.25), smmPoint(1.25, 1),
                                          smmPoint(-0.25, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 8);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapBisection)
{
  // One edge of mesh2 cuts el1 in half.
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0), smmPoint(1, 0),   smmPoint(1, 1),
                                              smmPoint(0, 1), smmPoint(-1, -1), smmPoint(2, -1),
                                              smmPoint(2, 2), smmPoint(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {3, 2, 6, 7}, {4, 0, 3, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(-1, -1),  smmPoint(2, -1), smmPoint(2, 0.5),
                                              smmPoint(-1, 0.5), smmPoint(2, 2),  smmPoint(-1, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {3, 2, 4, 5}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 16);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapSharedCorner)
{
  // The top left corner of el1 overlaps the top left corner of el2, but el2
  // is larger than el1
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, -1), smmPoint(1, -1), smmPoint(2, -1),
                                              smmPoint(0, 0),  smmPoint(1, 0),  smmPoint(2, 0),
                                              smmPoint(0, 1),  smmPoint(1, 1),  smmPoint(2, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 4, 3}, {1, 2, 5, 4}, {3, 4, 7, 6}, {4, 5, 8, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(0, -1), smmPoint(2, -1), smmPoint(2, 1),
                                              smmPoint(0, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 8);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapCutCorner)
{
  // An edge of the first element on mesh2 cuts off the top left corner of el1.
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),    smmPoint(1, 0),     smmPoint(1, 1),
                                              smmPoint(0, 1),    smmPoint(0.75, 2),  smmPoint(1, 3),
                                              smmPoint(-0.5, 3), smmPoint(-0.5, -1), smmPoint(2, -1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},  {7, 8, 1, 0},  {1, 8, 2, -1}, {2, 8, 4, -1},
                                              {2, 4, 3, -1}, {3, 4, 5, -1}, {3, 5, 6, -1}, {0, 3, 6, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(-0.5, 0), smmPoint(0.75, 2),  smmPoint(1, 3),
                                              smmPoint(-0.5, 3), smmPoint(-0.5, -1), smmPoint(2, -1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {4, 5, 1, 0}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 16);
}

TEST_F(Nonconformal4Tester, ThreeEdgeOverlapCutCorners)
{
  // edges of an element on mesh2 cut off the top left and top right corners
  // of el1.  No other edges are intersected
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0), smmPoint(1, 0),     smmPoint(1, 1),
                                              smmPoint(0, 1), smmPoint(-0.5, -1), smmPoint(1.25, -0.5),
                                              smmPoint(2, 3), smmPoint(-0.5, 3)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7}, {0, 3, 7, 4}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(-0.5, 0), smmPoint(0.75, 2),  smmPoint(1, 3),
                                              smmPoint(-0.5, 3), smmPoint(-0.5, -1), smmPoint(1.25, -.5),
                                              smmPoint(2, 3)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 22);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapBisection2)
{
  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<smmPoint> mesh1Verts     = {smmPoint(-1, 0), smmPoint(0, 0),  smmPoint(1, 0),
                                              smmPoint(3, 0),  smmPoint(-1, 1), smmPoint(0, 1),
                                              smmPoint(1, 1),  smmPoint(3, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(-1, 0), smmPoint(0.5, 0), smmPoint(0.5, 1),
                                              smmPoint(-1, 1), smmPoint(3, 0),   smmPoint(3, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 8);
}

TEST_F(Nonconformal4Tester, TwoVertexOverlapDiagonal)
{
  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<smmPoint> mesh1Verts = {
      smmPoint(0, 0),  smmPoint(1, 0), smmPoint(1, 1),   smmPoint(0, 1),  smmPoint(-1, -1),
      smmPoint(2, -1), smmPoint(3, 0), smmPoint(3, 1.5), smmPoint(-1, 2), smmPoint(-2, 0.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, -1}, {1, 6, 7, 2},
                                              {2, 7, 8, 3}, {0, 3, 8, 9}, {4, 0, 9, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(-1, -1),  smmPoint(2, -1), smmPoint(-1, 2),
                                              smmPoint(-2, 0.5), smmPoint(3, 0),  smmPoint(3, 1.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 12);
}

TEST_F(Nonconformal4Tester, FourVertexOverlapDiagonals)
{
  // mesh2 edges cross both diagonals of element 1, dividing it into 4 triangles
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),   smmPoint(1, 0),    smmPoint(1, 1),
                                              smmPoint(0, 1),   smmPoint(0.5, -3), smmPoint(2, -1),
                                              smmPoint(3, 0.5), smmPoint(2, 2),    smmPoint(0.5, 3),
                                              smmPoint(-1, 2),  smmPoint(-3, 0.5), smmPoint(-1, -1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(0.5, 0.5), smmPoint(-1, -1), smmPoint(2, -1),
                                              smmPoint(2, 2),     smmPoint(-1, 2),  smmPoint(0.5, -3),
                                              smmPoint(3, 0.5),   smmPoint(0.5, 3), smmPoint(-3, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 16);
}

TEST_F(Nonconformal4Tester, ThreeVertexOverlapDiagonals)
{
  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),   smmPoint(1, 0),    smmPoint(1, 1),
                                              smmPoint(0, 1),   smmPoint(0.5, -3), smmPoint(2, -1),
                                              smmPoint(3, 0.5), smmPoint(2, 2),    smmPoint(0.5, 3),
                                              smmPoint(-1, 2),  smmPoint(-3, 0.5), smmPoint(-1, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(0.5, 0.5), smmPoint(-1, 0.25), smmPoint(2, -1),
                                              smmPoint(2, 2),     smmPoint(-1, 2),    smmPoint(0.5, -3),
                                              smmPoint(3, 0.5),   smmPoint(0.5, 3),   smmPoint(-3, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 18);
}

TEST_F(Nonconformal4Tester, ThreeVertexNonCentroid)
{
  // mesh2 edges cross 3 vertices and one edge of el1, but the interior point is not
  // on the centroid of el1
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),   smmPoint(1, 0),    smmPoint(1, 1),
                                              smmPoint(0, 1),   smmPoint(0.5, -3), smmPoint(3, -2.0 / 3.0),
                                              smmPoint(3, 0.5), smmPoint(2, 2),    smmPoint(0.5, 3.5),
                                              smmPoint(-1, 4),  smmPoint(-3, 2),   smmPoint(-1, 0.125)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts = {
      smmPoint(0.25, 0.25), smmPoint(-1, 0.125), smmPoint(3, -2.0 / 3.0),
      smmPoint(2, 2),       smmPoint(-1, 4),     smmPoint(0.5, -3),
      smmPoint(3, 0.5),     smmPoint(0.5, 3.5),  smmPoint(-3, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 18);
}

TEST_F(Nonconformal4Tester, DoubleEdgeIntersection)
{
  // Two edges of a mesh2 element intersect the same edge of the mesh1 element
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),     smmPoint(1, 0),       smmPoint(1, 1),
                                              smmPoint(0, 1),     smmPoint(-0.25, -1),  smmPoint(1.5, -1),
                                              smmPoint(1.5, 2.0), smmPoint(0.75, 2),    smmPoint(0.25, 2.5),
                                              smmPoint(-1, 0.75), smmPoint(-1.5, 0.25), smmPoint(-1, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},   {4, 5, 1, 0},  {5, 6, 2, 1},   {2, 6, 7, -1},
                                              {2, 7, 8, 3},   {3, 8, 9, -1}, {3, 9, 11, -1}, {9, 10, 11, -1},
                                              {0, 3, 11, -1}, {0, 11, 4, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts = {
      smmPoint(0.5, 0.5),  smmPoint(-0.25, -1), smmPoint(0.5, -1),    smmPoint(-1, 0.25),
      smmPoint(1.5, -1),   smmPoint(1.5, 0.5),  smmPoint(1.5, 2),     smmPoint(0.75, 2),
      smmPoint(0.25, 2.5), smmPoint(-1, 0.75),  smmPoint(-1.5, 0.25),
  };
  std::vector<std::array<int, 4>> mesh2Els = {{1, 2, 0, 3}, {2, 4, 5, 0}, {0, 5, 6, 7}, {0, 7, 8, 9}, {3, 0, 9, 10}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 30);
}

TEST_F(Nonconformal4Tester, EdgeCrossThroughCornerVertex)
{
  // Edge starts on midpoint of one edge of el1 and passes through a vertex that
  // is not adjacent to the edge
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),  smmPoint(1, 0),  smmPoint(1, 1),
                                              smmPoint(0, 1),  smmPoint(0, -1), smmPoint(2, -2),
                                              smmPoint(3, -1), smmPoint(3, 2),  smmPoint(0, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 1, -1}, {1, 6, 7, 2}, {2, 7, 8, 3}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts     = {smmPoint(0, -1), smmPoint(2, -2), smmPoint(0.5, 1),
                                              smmPoint(0, 1),  smmPoint(3, -1), smmPoint(3, 2),
                                              smmPoint(0, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}, {3, 2, 5, 6}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 11);
}

TEST_F(Nonconformal4Tester, ThreeElementsWithCutCorner)
{
  // two elements cover most of el1, but a third element cuts off a corner
  std::vector<smmPoint> mesh1Verts     = {smmPoint(0, 0),       smmPoint(1, 0),      smmPoint(1, 1),
                                              smmPoint(0, 1),       smmPoint(-1, 0),     smmPoint(0.5, -0.25),
                                              smmPoint(1.5, -0.25), smmPoint(1.5, 1.25), smmPoint(-0.25, 1.25),
                                              smmPoint(-0.25, 0.5), smmPoint(-0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {5, 6, 1, 0},  {1, 6, 7, 2}, {2, 7, 8, 3},
                                              {0, 3, 8, 9}, {0, 9, 10, 4}, {0, 4, 5, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<smmPoint> mesh2Verts = {smmPoint(0.5, -0.25), smmPoint(-0.25, 0.5), smmPoint(-0.5, 0.5),
                                          smmPoint(-1, 0),      smmPoint(1.5, -0.25), smmPoint(1.5, 1.25),
                                          smmPoint(0.5, 1.25),  smmPoint(-0.25, 1.25)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {0, 4, 5, 6}, {0, 6, 7, 1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 22);
}

