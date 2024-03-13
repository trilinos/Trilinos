#include "gtest/gtest.h"

#include "stk_middle_mesh/boundary_fixture.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/create_mesh_quality_improver.hpp"
#include "stk_middle_mesh/mesh_boundary_snapper.hpp"
#include "stk_middle_mesh/mesh_quality_improver.hpp"
#include "stk_middle_mesh/regularized_distortion_metric.hpp"
#include "util/meshes.hpp"

#include "stk_middle_mesh/mesh_io.hpp" //TODO: DEBUGGING

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

namespace {

MeshEntityPtr find_closest_vert(std::shared_ptr<Mesh> mesh, const utils::Point& pt, double tol=10000)
{
  double minDist        = std::numeric_limits<double>::max();
  MeshEntityPtr minVert = nullptr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      auto pt2    = vert->get_point_orig(0);
      auto disp   = pt - pt2;
      double dist = dot(disp, disp);
      if (dist < minDist && dist < tol)
      {
        minDist = dist;
        minVert = vert;
      }
    }

  return minVert;
}

void test_meshes_identical(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2, const double tol)
{
  auto& verts1 = mesh1->get_vertices();
  auto& verts2 = mesh2->get_vertices();

  EXPECT_EQ(verts1.size(), verts2.size());

  for (unsigned int i = 0; i < verts1.size(); ++i)
  {
    auto v1 = verts1[i];
    auto v2 = verts2[i];

    EXPECT_NEAR(v1->get_point_orig(0).x, v2->get_point_orig(0).x, tol);
    EXPECT_NEAR(v1->get_point_orig(0).y, v2->get_point_orig(0).y, tol);
  }
}

} // namespace

TEST(MeshQualityImprover, Ideal)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh  = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);

  BoundaryFixture filter(mesh);

  auto fixer = make_standard_improver(mesh, filter, {.nlayers = -1});
  EXPECT_EQ(fixer->count_invalid_points(), 0);
  fixer->run();

  test_meshes_identical(mesh, mesh2, 1e-13);
}

TEST(MeshQualityImprover, SinglePoint)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4 || utils::impl::comm_size(MPI_COMM_WORLD) == 3)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 2;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh  = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);
  MeshEntityPtr vert          = find_closest_vert(mesh, utils::Point(0.5, 0.5));

  auto metric = std::make_shared<RegularizedDistortionMetric<double>>();

  BoundaryFixture filter(mesh);

  auto fixer = make_standard_improver(mesh, filter, {.nlayers = -1, .maxDeltaX = 1e-13});
  EXPECT_EQ(fixer->count_invalid_points(), 0);

  // move the vertex
  vert->set_point_orig(0, utils::Point(0.55, 0.65));

  // the fixer should return the vertex to its original location
  fixer->run();

  test_meshes_identical(mesh, mesh2, 1e-13);
}

TEST(MeshQualityImprover, FourPoint)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // in this case, there are 4 dofs, but only one of them is moved from
  // the initial position
  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh  = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);
  MeshEntityPtr vert          = find_closest_vert(mesh, utils::Point(1.0 / 3.0, 1.0 / 3.0), 0.1);

  auto metric = std::make_shared<RegularizedDistortionMetric<double>>();

  BoundaryFixture filter(mesh);

  auto fixer = make_standard_improver(mesh, filter, {.nlayers = -1, .maxDeltaX = 1e-15, .itermax = 100});

  // move the vertex
  if (vert)
    vert->set_point_orig(0, utils::Point(0.4, 0.45));
  EXPECT_EQ(fixer->count_invalid_points(), 0);

  // the fixer should return the vertex to its original location
  fixer->run();

  test_meshes_identical(mesh, mesh2, 1e-11);
}

TEST(MeshQualityImprover, FourPointInvalid)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // in this case, there are 4 dofs, but only one of them is moved from
  // the initial position
  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh  = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);
  MeshEntityPtr vert          = find_closest_vert(mesh, utils::Point(1.0 / 3.0, 1.0 / 3.0), 0.1);

  auto metric = std::make_shared<RegularizedDistortionMetric<double>>(1e-3);

  BoundaryFixture filter(mesh);

  auto fixer =
      make_standard_improver(mesh, filter, {.nlayers = -1, .maxDeltaX = 1e-14, .itermax = 100, .delta = 1e-3});

  // move the vertex to an invalid location
  if (vert)
    vert->set_point_orig(0, utils::Point(0.75, 0.75));
  EXPECT_EQ(fixer->count_invalid_points(), 4);

  // the fixer should return the vertex to its original location
  fixer->run();

  EXPECT_EQ(fixer->count_invalid_points(), 0);
  test_meshes_identical(mesh, mesh2, 1e-10);
}

TEST(MeshQualityImprover, FourPointAllInvalid)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // in this case, there are 4 dofs, but only one of them is moved from
  // the initial position
  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return utils::Point(pt.x, 0, pt.y); };

  std::shared_ptr<Mesh> mesh  = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);
  MeshEntityPtr vert1         = find_closest_vert(mesh, utils::Point(1.0 / 3.0, 0, 1.0 / 3.0), 0.1);
  MeshEntityPtr vert2         = find_closest_vert(mesh, utils::Point(2.0 / 3.0, 0, 1.0 / 3.0), 0.1);
  MeshEntityPtr vert3         = find_closest_vert(mesh, utils::Point(2.0 / 3.0, 0, 2.0 / 3.0), 0.1);
  MeshEntityPtr vert4         = find_closest_vert(mesh, utils::Point(1.0 / 3.0, 0, 2.0 / 3.0), 0.1);

  BoundaryFixture filter(mesh);

  auto fixer =
      make_standard_improver(mesh, filter, {.nlayers = -1, .maxDeltaX = 1e-14, .itermax = 100, .delta = 1e-3});

  // move all vertices to an invalid location
  if (vert1)
    vert1->set_point_orig(0, utils::Point(0.75, 0, 0.75));
  if (vert2)
    vert2->set_point_orig(0, utils::Point(0.25, 0, 0.75));
  if (vert3)
    vert3->set_point_orig(0, utils::Point(0.25, 0, 0.25));
  if (vert4)
    vert4->set_point_orig(0, utils::Point(0.75, 0, 0.25));

  EXPECT_EQ(fixer->count_invalid_points(), 4);

  // the fixer should return the vertex to its original location
  fixer->run();

  EXPECT_EQ(fixer->count_invalid_points(), 0);
  test_meshes_identical(mesh, mesh2, 1e-8);
}

TEST(MeshQualityImprover, 4By4)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // in this case, there are 4 dofs, but only one of them is moved from
  // the initial position
  MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh  = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec, func);

  BoundaryFixture filter(mesh);

  auto fixer =
      make_standard_improver(mesh, filter, {.nlayers = -1, .maxDeltaX = 1e-10, .itermax = 1000, .delta = 1e-3});

  double delta = 0.25;
  for (int i=1; i <= 3; ++i)
    for (int j=1; j <= 3; ++j)
    {
      MeshEntityPtr vert = find_closest_vert(mesh, utils::Point(i*delta, j*delta), 0.01);
      if (vert)
      {
        utils::Point pt = vert->get_point_orig(0);
        std::cout << "perturbing vert at " << pt << std::endl;
        vert->set_point_orig(0, pt + utils::Point(0.125, 0.125));
      }
    }

  // the fixer should return the vertex to its original location
  fixer->run();

  EXPECT_EQ(fixer->count_invalid_points(), 0);
  test_meshes_identical(mesh, mesh2, 1e-8);
}

// TODO: do a rotating partial sphere

} // namespace impl
} // namespace middle_mesh
} // namespace stk
