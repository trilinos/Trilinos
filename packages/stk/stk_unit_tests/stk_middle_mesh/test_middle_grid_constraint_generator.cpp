#include "stk_middle_mesh/mesh_projection_calculator.hpp"
#include "stk_middle_mesh/middle_grid_constraint_generator.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh::impl;

namespace {

class MiddleGridConstraintGeneratorTester : public ::testing::Test
{
  protected:
    MiddleGridConstraintGeneratorTester()
      : mMesh1(mesh::make_empty_mesh())
      , mMesh2(mesh::make_empty_mesh())
      , mMeshIn(mesh::make_empty_mesh())
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
      double eps  = 1e-13;
      mClassifier = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(
          mMesh2, PointClassifierNormalWrapperTolerances(eps));
      mRelationalData = std::make_shared<nonconformal4::impl::MeshRelationalData>(mMesh1, mMesh2, mMeshIn);
      mProjector      = std::make_shared<nonconformal4::impl::MeshProjectionCalculator>(
          mMesh1, mMesh2, mRelationalData, mClassifier, middle_mesh::impl::EdgeTracerTolerances(eps));
      mConstraintGenerator = std::make_shared<nonconformal4::impl::MiddleGridConstraintGenerator>(
          mMesh1, mMesh2, mMeshIn, mRelationalData, mClassifier);
      mProjector->project();
      mConstraintGenerator->generate();
    }

    std::shared_ptr<mesh::Mesh> mMesh1;
    std::shared_ptr<mesh::Mesh> mMesh2;
    std::shared_ptr<mesh::Mesh> mMeshIn;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> mClassifier;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> mRelationalData;
    std::shared_ptr<nonconformal4::impl::MeshProjectionCalculator> mProjector;
    std::shared_ptr<nonconformal4::impl::MiddleGridConstraintGenerator> mConstraintGenerator;

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

void check_entity_counts(std::shared_ptr<mesh::Mesh> mesh, int numVerts, int numEdges, int numEls)
{
  EXPECT_EQ(count_valid(mesh->get_vertices()), numVerts);
  EXPECT_EQ(count_valid(mesh->get_edges()), numEdges);
  EXPECT_EQ(count_valid(mesh->get_elements()), numEls);
}
} // namespace

TEST_F(MiddleGridConstraintGeneratorTester, El1ContainedInEl2)
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

  check_entity_counts(mMeshIn, 8, 12, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, El1ContainedInEl2Tri)
{
  // No intersections.  Element 1 is completely inside element 2
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),         utils::Point(1, 0),        utils::Point(0, 1),
      utils::Point(-0.25, -0.5), utils::Point(2, -0.5), utils::Point(-0.25, 1.5)};
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

  check_entity_counts(mMeshIn, 6, 12, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, SingleEdgeOverlap)
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

  check_entity_counts(mMeshIn, 8, 11, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, SingleEdgeOverlapTri)
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

  check_entity_counts(mMeshIn, 4, 6, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapBisection)
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

  check_entity_counts(mMeshIn, 12, 19, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapBisectionTri)
{
  // One edge of mesh2 cuts el1 in half.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1}};
  create_mesh1(mesh1Verts, mesh1Els);

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1), utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {0, 3, 2, -1}};
  create_mesh2(mesh2Verts, mesh2Els);

  run();

  check_entity_counts(mMeshIn, 4, 5, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapSharedCorner)
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

  check_entity_counts(mMeshIn, 9, 12, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapCutCorner)
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

  check_entity_counts(mMeshIn, 12, 22, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapCutCornerTri)
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

  check_entity_counts(mMeshIn, 8, 15, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeEdgeOverlapCutCorners)
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

  check_entity_counts(mMeshIn, 15, 25, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeEdgeOverlapCutCornersTri)
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

  check_entity_counts(mMeshIn, 8, 15, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapBisection2)
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

  check_entity_counts(mMeshIn, 10, 13, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoEdgeOverlapBisection2Tri)
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

  check_entity_counts(mMeshIn, 6, 9, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoVertexOverlapDiagonal)
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

  check_entity_counts(mMeshIn, 10, 17, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, TwoVertexOverlapDiagonalTri)
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

  check_entity_counts(mMeshIn, 8, 13, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, FourVertexOverlapDiagonals)
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

  check_entity_counts(mMeshIn, 13, 24, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, FourVertexOverlapDiagonalsTri)
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

  check_entity_counts(mMeshIn, 9, 16, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeVertexOverlapDiagonals)
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
  check_entity_counts(mMeshIn, 14, 26, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeVertexOverlapDiagonalsTri)
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
  check_entity_counts(mMeshIn, 8, 15, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeVertexNonCentroid)
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
  check_entity_counts(mMeshIn, 14, 26, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeVertexNonCentroidTri)
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
  check_entity_counts(mMeshIn, 8, 15, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, DoubleEdgeIntersection)
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
  check_entity_counts(mMeshIn, 21, 40, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, DoubleEdgeIntersectionTri)
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
  check_entity_counts(mMeshIn, 9, 17, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, EdgeCrossThroughCornerVertex)
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
  check_entity_counts(mMeshIn, 10, 16, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, EdgeCrossThroughCornerVertexTri)
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
  check_entity_counts(mMeshIn, 8, 14, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeElementsWithCutCorner)
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
  check_entity_counts(mMeshIn, 16, 28, 0);
}

TEST_F(MiddleGridConstraintGeneratorTester, ThreeElementsWithCutCornerTri)
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
  check_entity_counts(mMeshIn, 9, 15, 0);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
