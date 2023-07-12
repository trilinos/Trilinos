#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/element_operations_2d.hpp"
#include "stk_middle_mesh/mesh_boundary_snapper.hpp"
#include "stk_middle_mesh/mesh_io.hpp" //TODO: DEBUGGING
#include "stk_middle_mesh/field.hpp"

#include "util/meshes.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace stk::middle_mesh::mesh;
using namespace stk::middle_mesh::mesh::impl;

namespace {

std::shared_ptr<Mesh> make_square_boundary(std::vector<double>& xCoords)
{
  std::sort(xCoords.begin(), xCoords.end());

  auto mesh = make_empty_mesh();
  std::vector<MeshEntityPtr> vertsBottom, vertsTop;
  for (auto& xCoord : xCoords)
    vertsBottom.push_back(mesh->create_vertex(xCoord, 0, 0));

  for (auto& xCoord : xCoords)
    vertsTop.push_back(mesh->create_vertex(xCoord, 1, 0));

  for (size_t i = 1; i < vertsBottom.size(); ++i)
    mesh->create_quad_from_verts(vertsBottom[i - 1], vertsBottom[i], vertsTop[i], vertsTop[i - 1]);

  return mesh;
}

void test_edge_lengths_positive(std::shared_ptr<Mesh> mesh, int nverts)
{
  auto& verts = mesh->get_vertices();
  for (int i = 1; i < nverts; ++i)
    EXPECT_GE(verts[i]->get_point_orig(0).x, verts[i - 1]->get_point_orig(0).x);
}

void test_edges_only(std::vector<double> verts1, std::vector<double> verts2)
{
  auto mesh1 = make_square_boundary(verts1);
  auto mesh2 = make_square_boundary(verts2);

  MeshBoundarySnapper snapper;
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_edge_lengths_positive(mesh1, verts1.size());
  test_edge_lengths_positive(mesh2, verts2.size());

  mesh1 = make_square_boundary(verts2);
  mesh2 = make_square_boundary(verts1);
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_edge_lengths_positive(mesh1, verts2.size());
  test_edge_lengths_positive(mesh2, verts1.size());
}

void test_areas_positive(std::shared_ptr<Mesh> mesh)
{
  mesh::impl::ElementOperations2D elemOps;
  for (auto& e : mesh->get_elements())
    if (e)
    {
      EXPECT_GE(elemOps.compute_area(e), 0);
    }
}

void test_areas_equal(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2, MPI_Comm unionComm = MPI_COMM_WORLD)
{
  mesh::impl::ElementOperations2D elemOps;
  double area1Local = 0.0, area2Local = 0.0;

  if (mesh1)
    area1Local = elemOps.compute_area(mesh1);

  if (mesh2)
    area2Local = elemOps.compute_area(mesh2);

  const int root = 0;
  int nprocs = utils::impl::comm_size(unionComm);
  std::vector<double> areas1(nprocs), areas2(nprocs);
  MPI_Gather(&area1Local, 1, MPI_DOUBLE, areas1.data(), 1, MPI_DOUBLE, root, unionComm);
  MPI_Gather(&area2Local, 1, MPI_DOUBLE, areas2.data(), 1, MPI_DOUBLE, root, unionComm);

  double area1Sum = 0.0, area2Sum = 0.0;
  for (int i=0; i < nprocs; ++i)
  {
    area1Sum += areas1[i];
    area2Sum += areas2[i];
  }

  EXPECT_FLOAT_EQ(area1Sum, area2Sum);
}

} // namespace

TEST(MeshBoundarySnapper, Identical1D)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  test_edges_only({0, 0.2, 0.4, 1}, {0, 0.2, 0.4, 1});
}

TEST(MeshBoundarySnapper, Subdivision1D)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  test_edges_only({0, 0.2, 0.4, 1}, {0, 0.1, 0.2, 0.3, 0.4, 0.7, 1});
}

TEST(MeshBoundarySnapper, Dense1D)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  test_edges_only({0, 0.2, 0.4, 1}, {0, 0.01, 0.02, 0.03, 0.04, 0.1, 0.2, 0.3, 0.4, 0.7, 1});
  // testEdgesOnly({0, 0.01, 0.02, 0.03, 0.04, 0.1, 0.2, 0.3, 0.4, 0.7, 1}, {0, 0.2, 0.4, 1});
  // testEdgesOnly({0, 0.03, 0.06, 0.09, 0.012, 0.3, 0.6, 0.9, 1.2, 2.1, 3}, {0, 0.6, 1.2, 3});
}

TEST(MeshBoundarySnapper, Identical)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);

  MeshBoundarySnapper snapper;
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_areas_equal(mesh1, mesh2);
}

TEST(MeshBoundarySnapper, Identical2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  MeshSpec spec2;
  spec2.numelX = 3;
  spec2.numelY = 4;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

  MeshBoundarySnapper snapper;
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_areas_equal(mesh1, mesh2);
}

TEST(MeshBoundarySnapper, Subdivision)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MeshSpec spec, spec2;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 6;
  spec2.numelY = 8;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

  MeshBoundarySnapper snapper;
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_areas_equal(mesh1, mesh2);
}

TEST(MeshBoundarySnapper, QuarterAnnulus)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // std::cout << std::setprecision(16) << std::endl;
  double pi   = std::atan(1) * 4;
  int nmeshes = 30;

  for (int i = 0; i < nmeshes; ++i)
  {
    // std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;
    MeshSpec spec, spec2;
    spec.numelX = 5;
    spec.numelY = 5;
    spec.xmin   = 0.5;
    spec.xmax   = 1.5;
    spec.ymin   = 0;
    spec.ymax   = pi / 2;

    spec2.numelX = 5 + i;
    spec2.numelY = 5 + i;
    spec2.xmin   = 0.5;
    spec2.xmax   = 1.5;
    spec2.ymin   = 0;
    spec2.ymax   = pi / 2;

    // remap coordinates to annulus

    //    auto func = [&](const utils::Point& pt) { return pt; };

    auto func = [&](const utils::Point& pt) {
      // interpret x and y as r and theta
      double r     = pt.x;
      double theta = pt.y;

      double x = r * std::cos(theta);
      double y = r * std::sin(theta);
      double z = 0.0;
      utils::Point pt2(x, y, z);
      return pt2;
    };

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

    // printVertEdges("mesh1_initial", mesh1);
    // printVertEdges("mesh2_initial", mesh2);

    MeshBoundarySnapper snapper;
    snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

    // printVertEdges("mesh1_final", mesh1);
    // printVertEdges("mesh2_final", mesh2);

    test_areas_positive(mesh1);
    test_areas_positive(mesh2);
    test_areas_equal(mesh1, mesh2);
  }
}


TEST(MeshBoundarySnapper, PeriodicAnnulus)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // std::cout << std::setprecision(16) << std::endl;
  double pi     = std::atan(1) * 4;
  int nmeshes   = 640;
  double dtheta = pi / (16 * nmeshes);

  for (int i = 0; i < nmeshes; ++i)
  {
    // std::cout << "mesh " << i + 1 << " / " << nmeshes << std::endl;

    // std::cout << "creating mesh1" << std::endl;
    std::shared_ptr<Mesh> mesh1 = make_annulus_mesh(10, 10, 0.5, 1.5, 0);

    // std::cout << "\ncreating mesh2" << std::endl;
    std::shared_ptr<Mesh> mesh2 = make_annulus_mesh(13, 13, 0.5, 1.5, i * dtheta);

    /*
        printVertEdges("mesh1_initial", mesh1);
        printVertEdges("mesh2_initial", mesh2);
    */

    MeshBoundarySnapper snapper;
    snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);
    /*
        printVertEdges("mesh1_final", mesh1);
        printVertEdges("mesh2_final", mesh2);
    */
    test_areas_positive(mesh1);
    test_areas_positive(mesh2);
    test_areas_equal(mesh1, mesh2);
  }
}

TEST(MeshBoundarySnapper, PeriodicAnnulusMPMD)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 8 || utils::impl::comm_size(MPI_COMM_WORLD) <= 1)
    GTEST_SKIP();
  
  int myRank   = utils::impl::comm_rank(MPI_COMM_WORLD);
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  int commSizeUsed           = commSize - (commSize % 2);
  int color;
  if (myRank < commSizeUsed)
    color = myRank % 2 == 0;
  else
    color = MPI_UNDEFINED;

  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

  double pi     = std::atan(1) * 4;
  int nmeshes   = 640;
  double dtheta = pi / (16 * nmeshes);

  for (int i = 0; i < nmeshes; ++i)
  {
    std::shared_ptr<Mesh> mesh1, mesh2;

    if (color == 0)
      mesh1 = make_annulus_mesh(10, 10, 0.5, 1.5, 0, meshComm);

    if (color == 1)
      mesh2 = make_annulus_mesh(13, 13, 0.5, 1.5, i * dtheta, meshComm);

    /*
        printVertEdges("mesh1_initial", mesh1);
        printVertEdges("mesh2_initial", mesh2);
    */

    MeshBoundarySnapper snapper;
    snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);
    /*
        printVertEdges("mesh1_final", mesh1);
        printVertEdges("mesh2_final", mesh2);
    */

    if (mesh1 != nullptr)
      test_areas_positive(mesh1);
    if (mesh2 != nullptr)
      test_areas_positive(mesh2);

    test_areas_equal(mesh1, mesh2, MPI_COMM_WORLD);
  }

  if (meshComm != MPI_COMM_NULL)
    MPI_Comm_free(&meshComm);
}

TEST(MeshBoundarySnapper, TestMeshExchangeBoundaryEdgesGather)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) <= 1)
    GTEST_SKIP();

  int myRank = utils::impl::comm_rank(MPI_COMM_WORLD);

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  MeshExchangeBoundaryEdges exchange(mesh, mesh->get_comm(), 0);
  auto boundaryEdgeMesh = exchange.get_boundary_edge_mesh();

  if (myRank == 0) {
    const unsigned numExpectedGatheredEdges = 14;
    const unsigned numExpectedGatheredVertices = 14;

    EXPECT_EQ(numExpectedGatheredEdges, boundaryEdgeMesh->get_edges().size());
    EXPECT_EQ(numExpectedGatheredVertices, boundaryEdgeMesh->get_vertices().size());
  }
}

TEST(MeshBoundarySnapper, TestMeshExchangeBoundaryEdgesRemoteUpdate)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) <= 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);
  auto myRank = utils::impl::comm_rank(MPI_COMM_WORLD);

  auto origVertices = mesh->get_vertices();
  auto field = create_field<utils::Point>(mesh, FieldShape(1, 0, 0), 1);
  auto& fieldRef = *field;

  for (auto vertex : origVertices) {
    fieldRef(vertex, 0, 0) = vertex->get_point_orig(0);
  }

  MeshExchangeBoundaryEdges exchange(mesh, mesh->get_comm(), 0);
  auto boundaryEdgeMesh = exchange.get_boundary_edge_mesh();

  if (myRank == 0) {
    auto vertices = boundaryEdgeMesh->get_vertices();
    for (auto vertex : vertices) {
      utils::Point pt = vertex->get_point_orig(0);
      pt = pt * 2;
      vertex->set_point_orig(0, pt);
    }
  }
  exchange.update_remote_vertices();

  for (auto vertex : origVertices) {
    std::vector<MeshEntityPtr> edges;
    get_upward(vertex, 1, edges);

    for (auto edge : edges)
      if (edge->count_up() == 1 && edge->count_remote_shared_entities() == 0) {
        for (auto i = 0; i < 3; ++i)
          EXPECT_DOUBLE_EQ(vertex->get_point_orig(0)[i], fieldRef(vertex, 0, 0)[i] * 2);
      }
  }
}

TEST(MeshBoundarySnapper, SimpleMeshesFinerMesh1)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MeshSpec spec, spec2;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 4;
  spec2.numelY = 4;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

  MeshBoundarySnapper snapper;
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_areas_positive(mesh1);
  test_areas_positive(mesh2);
  test_areas_equal(mesh1, mesh2);
}

TEST(MeshBoundarySnapper, SimpleMeshesFinerMesh2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MeshSpec spec, spec2;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 5;
  spec2.numelY = 5;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

  MeshBoundarySnapper snapper;
  snapper.snap(mesh1, mesh2, MPI_COMM_WORLD);

  test_areas_positive(mesh1);
  test_areas_positive(mesh2);
  test_areas_equal(mesh1, mesh2);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
