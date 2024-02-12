#include "stk_middle_mesh/element_operations_2d.hpp"
#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "util/nonconformal_interface_helpers.hpp"
#include "gtest/gtest.h"

namespace {

using namespace stk::middle_mesh;
/*
using utils::Point = stk::middle_mesh::utils::Point;
using mesh::Mesh = stk::middle_mesh::mesh::Mesh;
using mesh::MeshEntityPtr = stk::middle_mesh::mesh::MeshEntityPtr;
template<typename T>
using mesh::FieldPtr = stk::middle_mesh::mesh::FieldPtr<T>;
template<typename T>
using mesh::VariableSizeFieldPtr = stk::middle_mesh::mesh::VariableSizeFieldPtr<T>;
using nonconformal4::impl::Nonconformal4 = stk::middle_mesh::nonconformal4::impl::Nonconformal4;
*/
class Nonconformal4Tester : public ::testing::Test
{
  protected:
    Nonconformal4Tester()
      : mMesh1(stk::middle_mesh::mesh::make_empty_mesh())
      , mMesh2(stk::middle_mesh::mesh::make_empty_mesh())
      , mMeshIn(stk::middle_mesh::mesh::make_empty_mesh())
    {}

    void create_mesh1(const std::vector<utils::Point>& verts, const std::vector<std::array<int, 4>>& elementVerts)
    {
      create_mesh(mMesh1, verts, elementVerts);
    }

    void create_mesh2(const std::vector<utils::Point>& verts, const std::vector<std::array<int, 4>>& elementVerts)
    {
      create_mesh(mMesh2, verts, elementVerts);
    }

    void run()
    {
      double eps = 1e-13;
      stk::middle_mesh::NormalProjectionOpts opts;
      opts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
      opts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);

      mNonconformalGenerator      = std::make_shared<nonconformal4::impl::Nonconformal4>(mMesh1, mMesh2, opts);
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

    std::shared_ptr<mesh::Mesh> mMesh1;
    std::shared_ptr<mesh::Mesh> mMesh2;
    std::shared_ptr<nonconformal4::impl::Nonconformal4> mNonconformalGenerator;
    std::shared_ptr<mesh::Mesh> mMeshIn;
    mesh::FieldPtr<mesh::MeshEntityPtr> mMeshInElementsToMesh1;
    mesh::FieldPtr<mesh::MeshEntityPtr> mMeshInElementsToMesh2;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> mMesh1InverseClassification;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> mMesh2InverseClassification;

  private:
    void create_mesh(std::shared_ptr<mesh::Mesh> mesh, const std::vector<utils::Point>& vertCoords,
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
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),         utils::Point(1, 0),        utils::Point(1, 1),       utils::Point(0, 1),
      utils::Point(-0.25, -0.25), utils::Point(1.25, -0.25), utils::Point(1.25, 1.25), utils::Point(-0.25, 1.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {3, 2, 6, 7}, {0, 3, 7, 4}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {mesh1Verts[4], mesh1Verts[5], mesh1Verts[6], mesh1Verts[7]};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

TEST_F(Nonconformal4Tester, El1ContainedInEl2Tri)
{
  // No intersections.  Element 1 is completely inside element 2
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),         utils::Point(1, 0),    utils::Point(0, 1),
      utils::Point(-0.25, -0.5),  utils::Point(2, -0.5), utils::Point(-0.25, 1.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1},
                                              {3, 4, 0, -1}, {0, 4, 1, -1},
                                              {1, 4, 2, -1}, {2, 4, 5, -1},
                                              {3, 0, 5, -1}, {0, 2, 5, -1}
                                              };
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {mesh1Verts[3], mesh1Verts[4], mesh1Verts[5]};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 7);
}

TEST_F(Nonconformal4Tester, SingleEdgeOverlap)
{
  // The top edge of el1 overlaps with the top edge of el2.  No other
  // edges itnersect
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0),         utils::Point(1, 1),
                                              utils::Point(0, 1),    utils::Point(-0.25, -0.25), utils::Point(1.25, -0.25),
                                              utils::Point(1.25, 1), utils::Point(-0.25, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {0, 3, 7, 4}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {utils::Point(-0.25, -0.25), utils::Point(1.25, -0.25), utils::Point(1.25, 1),
                                          utils::Point(-0.25, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 8);
}

TEST_F(Nonconformal4Tester, SingleEdgeOverlapTri)
{
  // The left edge of el1 overlaps with the top edge of el2.  No other
  // edges itnersect
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0), utils::Point(0, 1),    
                                              utils::Point(0.25, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 3, 2, -1}, {0, 1, 3, -1}, {3, 1, 2, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {mesh1Verts[0], mesh1Verts[1], mesh1Verts[2]};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 3);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapBisection)
{
  // One edge of mesh2 cuts el1 in half.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0),   utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {3, 2, 6, 7}, {4, 0, 3, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1),  utils::Point(2, -1), utils::Point(2, 0.5),
                                              utils::Point(-1, 0.5), utils::Point(2, 2),  utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {3, 2, 4, 5}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 16);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapBisectionTri)
{
  // One edge of mesh2 cuts el1 in half.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1), utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {0, 3, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 2);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapSharedCorner)
{
  // The top left corner of el1 overlaps the top left corner of el2, but el2
  // is larger than el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, -1), utils::Point(1, -1), utils::Point(2, -1),
                                              utils::Point(0, 0),  utils::Point(1, 0),  utils::Point(2, 0),
                                              utils::Point(0, 1),  utils::Point(1, 1),  utils::Point(2, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 4, 3}, {1, 2, 5, 4}, {3, 4, 7, 6}, {4, 5, 8, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, -1), utils::Point(2, -1), utils::Point(2, 1),
                                              utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 8);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapCutCorner)
{
  // An edge of the first element on mesh2 cuts off the top left corner of el1.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(1, 1),
                                              utils::Point(0, 1),    utils::Point(0.75, 2),  utils::Point(1, 3),
                                              utils::Point(-0.5, 3), utils::Point(-0.5, -1), utils::Point(2, -1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},  {7, 8, 1, 0},  {1, 8, 2, -1}, {2, 8, 4, -1},
                                              {2, 4, 3, -1}, {3, 4, 5, -1}, {3, 5, 6, -1}, {0, 3, 6, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-0.5, 0), utils::Point(0.75, 2),  utils::Point(1, 3),
                                              utils::Point(-0.5, 3), utils::Point(-0.5, -1), utils::Point(2, -1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {4, 5, 1, 0}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 16);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapCutCornerTri)
{
  // An edge of the first element on mesh2 cuts off the top left corner of el1.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.75, 0.75),  utils::Point(-0.25, 0.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1},  {1, 3, 2, -1},  {0, 2, 4, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.75, 0.75),  utils::Point(-0.25, 0.75)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {0, 3, 4, -1}, {4, 3, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 9);
}

TEST_F(Nonconformal4Tester, ThreeEdgeOverlapCutCorners)
{
  // edges of an element on mesh2 cut off the top left and top right corners
  // of el1.  No other edges are intersected
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0),     utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(-0.5, -1), utils::Point(1.25, -0.5),
                                              utils::Point(2, 3), utils::Point(-0.5, 3)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7}, {0, 3, 7, 4}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-0.5, 0), utils::Point(0.75, 2),  utils::Point(1, 3),
                                              utils::Point(-0.5, 3), utils::Point(-0.5, -1), utils::Point(1.25, -.5),
                                              utils::Point(2, 3)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 22);
}

TEST_F(Nonconformal4Tester, ThreeEdgeOverlapCutCornersTri)
{
  // edges of an element on mesh2 cut off the top left and top right corners
  // of el1.  No other edges are intersected
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),     utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.2, 0.9), utils::Point(0.9, 0.2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1}, {1, 4, 2, -1}, {4, 3, 2, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0),     utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.2, 0.9), utils::Point(0.9, 0.2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, -1}, {0, 4, 3, -1}, {0, 3, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 9);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapBisection2)
{
  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts     = {utils::Point(-1, 0), utils::Point(0, 0),  utils::Point(1, 0),
                                              utils::Point(3, 0),  utils::Point(-1, 1), utils::Point(0, 1),
                                              utils::Point(1, 1),  utils::Point(3, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, 0), utils::Point(0.5, 0), utils::Point(0.5, 1),
                                              utils::Point(-1, 1), utils::Point(3, 0),   utils::Point(3, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 8);
}

TEST_F(Nonconformal4Tester, TwoEdgeOverlapBisection2Tri)
{
  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0), utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(0.5, 0), utils::Point(0.5, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 4, 3, -1}, {4, 5, 3, -1}, {4, 2, 5, -1}, {4, 1, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 4);
}

TEST_F(Nonconformal4Tester, TwoVertexOverlapDiagonal)
{
  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),  utils::Point(1, 0), utils::Point(1, 1),   utils::Point(0, 1),  utils::Point(-1, -1),
      utils::Point(2, -1), utils::Point(3, 0), utils::Point(3, 1.5), utils::Point(-1, 2), utils::Point(-2, 0.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, -1}, {1, 6, 7, 2},
                                              {2, 7, 8, 3}, {0, 3, 8, 9}, {4, 0, 9, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1),  utils::Point(2, -1), utils::Point(-1, 2),
                                              utils::Point(-2, 0.5), utils::Point(3, 0),  utils::Point(3, 1.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 12);
}

TEST_F(Nonconformal4Tester, TwoVertexOverlapDiagonalTri)
{
  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),   utils::Point(1, 0),  utils::Point(1, 1), utils::Point(0, 1),
      utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7},
                                              {4, 0, 3, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {1, 2, 3, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

TEST_F(Nonconformal4Tester, FourVertexOverlapDiagonals)
{
  // mesh2 edges cross both diagonals of element 1, dividing it into 4 triangles
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(0.5, -3), utils::Point(2, -1),
                                              utils::Point(3, 0.5), utils::Point(2, 2),    utils::Point(0.5, 3),
                                              utils::Point(-1, 2),  utils::Point(-3, 0.5), utils::Point(-1, -1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0.5, 0.5), utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(2, 2),     utils::Point(-1, 2),  utils::Point(0.5, -3),
                                              utils::Point(3, 0.5),   utils::Point(0.5, 3), utils::Point(-3, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 16);
}

TEST_F(Nonconformal4Tester, FourVertexOverlapDiagonalsTri)
{
  // mesh2 edges cross both diagonals of element 1, dividing it into 4 triangles
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),   utils::Point(1, 0),  utils::Point(1, 1), utils::Point(0, 1),
      utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7},
                                              {4, 0, 3, 7}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2),
                                              utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 4, 3, -1}, {0, 1, 4, -1}, {1, 2, 4, -1}, {2, 3, 4, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 12);
}

TEST_F(Nonconformal4Tester, ThreeVertexOverlapDiagonals)
{
  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(0.5, -3), utils::Point(2, -1),
                                              utils::Point(3, 0.5), utils::Point(2, 2),    utils::Point(0.5, 3),
                                              utils::Point(-1, 2),  utils::Point(-3, 0.5), utils::Point(-1, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0.5, 0.5), utils::Point(-1, 0.25), utils::Point(2, -1),
                                              utils::Point(2, 2),     utils::Point(-1, 2),    utils::Point(0.5, -3),
                                              utils::Point(3, 0.5),   utils::Point(0.5, 3),   utils::Point(-3, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 18);
}

TEST_F(Nonconformal4Tester, ThreeVertexOverlapDiagonalsTri)
{
  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 2, -1},
                                              {3, 2, 6, -1}, {4, 0, 3, 6}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(2, -1), utils::Point(1, 1),
                                              utils::Point(-1, 2),  utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, -1}, {1, 2, 4, -1}, {4, 2, 3, -1}, {0, 4, 3, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

TEST_F(Nonconformal4Tester, ThreeVertexNonCentroid)
{
  // mesh2 edges cross 3 vertices and one edge of el1, but the interior point is not
  // on the centroid of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(0.5, -3), utils::Point(3, -2.0 / 3.0),
                                              utils::Point(3, 0.5), utils::Point(2, 2),    utils::Point(0.5, 3.5),
                                              utils::Point(-1, 4),  utils::Point(-3, 2),   utils::Point(-1, 0.125)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {
      utils::Point(0.25, 0.25), utils::Point(-1, 0.125), utils::Point(3, -2.0 / 3.0),
      utils::Point(2, 2),       utils::Point(-1, 4),     utils::Point(0.5, -3),
      utils::Point(3, 0.5),     utils::Point(0.5, 3.5),  utils::Point(-3, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 18);
}

TEST_F(Nonconformal4Tester, ThreeVertexNonCentroidTri)
{
  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(-1, -1), utils::Point(1.75, -0.25),
                                              utils::Point(-0.25, 1.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 2, -1},
                                              {3, 2, 6, -1}, {4, 0, 3, 6}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(1.75, -0.25), utils::Point(1, 1),
                                              utils::Point(-0.25, 1.75),  utils::Point(0.25, 0.25)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, -1}, {1, 2, 4, -1}, {4, 2, 3, -1}, {0, 4, 3, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

TEST_F(Nonconformal4Tester, DoubleEdgeIntersection)
{
  // Two edges of a mesh2 element intersect the same edge of the mesh1 element
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),     utils::Point(1, 0),       utils::Point(1, 1),
                                              utils::Point(0, 1),     utils::Point(-0.25, -1),  utils::Point(1.5, -1),
                                              utils::Point(1.5, 2.0), utils::Point(0.75, 2),    utils::Point(0.25, 2.5),
                                              utils::Point(-1, 0.75), utils::Point(-1.5, 0.25), utils::Point(-1, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},   {4, 5, 1, 0},  {5, 6, 2, 1},   {2, 6, 7, -1},
                                              {2, 7, 8, 3},   {3, 8, 9, -1}, {3, 9, 11, -1}, {9, 10, 11, -1},
                                              {0, 3, 11, -1}, {0, 11, 4, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {
      utils::Point(0.5, 0.5),  utils::Point(-0.25, -1), utils::Point(0.5, -1),    utils::Point(-1, 0.25),
      utils::Point(1.5, -1),   utils::Point(1.5, 0.5),  utils::Point(1.5, 2),     utils::Point(0.75, 2),
      utils::Point(0.25, 2.5), utils::Point(-1, 0.75),  utils::Point(-1.5, 0.25),
  };
  std::vector<std::array<int, 4>> mesh2Els = {{1, 2, 0, 3}, {2, 4, 5, 0}, {0, 5, 6, 7}, {0, 7, 8, 9}, {3, 0, 9, 10}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 30);
}

TEST_F(Nonconformal4Tester, DoubleEdgeIntersectionTri)
{
  // Two edges of a mesh2 element intersect the same edge of the mesh1 element
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),     utils::Point(1, 0),         utils::Point(1, 1),
                                              utils::Point(0, 1),     utils::Point(-0.25, 0.25),  utils::Point(-0.25, 0.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},   {4, 0, 3, 5}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {utils::Point(0, 0),     utils::Point(1, 0),         utils::Point(1, 1),
                                          utils::Point(0, 1),     utils::Point(-0.25, 0.25),  utils::Point(-0.25, 0.75),
                                          utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 6, 4, -1}, {0, 1, 6, -1}, {1, 2, 6, -1}, {6, 2, 3, -1}, {5, 6, 3, -1}, {4, 6, 5, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

TEST_F(Nonconformal4Tester, EdgeCrossThroughCornerVertex)
{
  // Edge starts on midpoint of one edge of el1 and passes through a vertex that
  // is not adjacent to the edge
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),  utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1),  utils::Point(0, -1), utils::Point(2, -2),
                                              utils::Point(3, -1), utils::Point(3, 2),  utils::Point(0, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 1, -1}, {1, 6, 7, 2}, {2, 7, 8, 3}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, -1), utils::Point(2, -2), utils::Point(0.5, 1),
                                              utils::Point(0, 1),  utils::Point(3, -1), utils::Point(3, 2),
                                              utils::Point(0, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}, {3, 2, 5, 6}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 11);
}

TEST_F(Nonconformal4Tester, EdgeCrossThroughCornerVertexTri)
{
  // Edge starts on midpoint of one edge of el1 and passes through a vertex that
  // is not adjacent to the edge
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),  utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1),  utils::Point(0, -1), utils::Point(1.5, -1),
                                              utils::Point(1.5, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, -1}, {1, 6, 2, -1}, {3, 2, 6, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, -1), utils::Point(1.5, -1), utils::Point(1.5, 2),
                                              utils::Point(0, 1),  utils::Point(0.5, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, 3}, {1, 2, 4, -1}, {3, 4, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 9);
}

TEST_F(Nonconformal4Tester, ThreeElementsWithCutCorner)
{
  // two elements cover most of el1, but a third element cuts off a corner
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),       utils::Point(1, 0),      utils::Point(1, 1),
                                              utils::Point(0, 1),       utils::Point(-1, 0),     utils::Point(0.5, -0.25),
                                              utils::Point(1.5, -0.25), utils::Point(1.5, 1.25), utils::Point(-0.25, 1.25),
                                              utils::Point(-0.25, 0.5), utils::Point(-0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {5, 6, 1, 0},  {1, 6, 7, 2}, {2, 7, 8, 3},
                                              {0, 3, 8, 9}, {0, 9, 10, 4}, {0, 4, 5, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {utils::Point(0.5, -0.25), utils::Point(-0.25, 0.5), utils::Point(-0.5, 0.5),
                                          utils::Point(-1, 0),      utils::Point(1.5, -0.25), utils::Point(1.5, 1.25),
                                          utils::Point(0.5, 1.25),  utils::Point(-0.25, 1.25)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {0, 4, 5, 6}, {0, 6, 7, 1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 22);
}

TEST_F(Nonconformal4Tester, ThreeElementsWithCutCornerTri)
{
  // two elements cover most of el1, but a third element cuts off a corner
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),       utils::Point(1, 0),      utils::Point(1, 1),
                                              utils::Point(0, 1),       utils::Point(-1, -1),    utils::Point(0.75, -0.5),
                                              utils::Point(-0.5, 0.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0},  {4, 0, 3, 6}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts = {utils::Point(0.75, -0.5), utils::Point(1, 0),   utils::Point(0, 1),
                                          utils::Point(-0.5, 0.75), utils::Point(-1, -1), utils::Point(1, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 5, 2, -1}, {4, 0, 3, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  stk::middle_mesh::test_util::test_number_of_elements(mMeshIn, 10);
}

