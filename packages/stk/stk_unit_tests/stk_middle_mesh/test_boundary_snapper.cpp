#include "gtest/gtest.h"

#include "create_mesh.hpp"
#include "element_operations_2d.hpp"
#include "mesh_boundary_snapper.hpp"
#include "mesh_io.hpp" //TODO: DEBUGGING

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
  snapper.snap(mesh1, mesh2);

  test_edge_lengths_positive(mesh1, verts1.size());
  test_edge_lengths_positive(mesh2, verts2.size());

  mesh1 = make_square_boundary(verts2);
  mesh2 = make_square_boundary(verts1);
  snapper.snap(mesh1, mesh2);

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

void test_areas_equal(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2)
{
  mesh::impl::ElementOperations2D elemOps;
  EXPECT_FLOAT_EQ(elemOps.compute_area(mesh1), elemOps.compute_area(mesh2));
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
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
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
  snapper.snap(mesh1, mesh2);

  test_areas_equal(mesh1, mesh2);
}

TEST(MeshBoundarySnapper, Subdivision)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
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
  snapper.snap(mesh1, mesh2);

  test_areas_equal(mesh1, mesh2);
}

TEST(MeshBoundarySnapper, QuarterAnnulus)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
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
    snapper.snap(mesh1, mesh2);

    // printVertEdges("mesh1_final", mesh1);
    // printVertEdges("mesh2_final", mesh2);

    test_areas_positive(mesh1);
    test_areas_positive(mesh2);
    test_areas_equal(mesh1, mesh2);
  }
}

TEST(MeshBoundarySnapper, PeriodicAnnulus)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
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
    snapper.snap(mesh1, mesh2);
    /*
        printVertEdges("mesh1_final", mesh1);
        printVertEdges("mesh2_final", mesh2);
    */
    test_areas_positive(mesh1);
    test_areas_positive(mesh2);
    test_areas_equal(mesh1, mesh2);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
