#include "gtest/gtest.h"
#include "stk_middle_mesh/remote_coordinate_updator.hpp"
#include "stk_middle_mesh/create_mesh.hpp"

namespace {

using namespace stk::middle_mesh;

mesh::MeshEntityPtr get_closest_vert(std::shared_ptr<stk::middle_mesh::mesh::Mesh> mesh, const utils::Point& pt, double tol=1e-12)
{
  double minDist = tol;
  mesh::MeshEntityPtr minVert = nullptr;
  for (auto vert : mesh->get_vertices())
  {
    auto disp = vert->get_point_orig(0) - pt;
    double dist = std::sqrt(dot(disp, disp));
    if (dist < minDist)
    {
      minDist = dist;
      minVert = vert;
    }
  }

  return minVert;
}

void increment_vert(stk::middle_mesh::mesh::MeshEntityPtr vert, const utils::Point& pt)
{
  auto vertPt = vert->get_point_orig(0);
  vert->set_point_orig(0, vertPt + pt);
}

}

TEST(RemoteCoordinateUpdator, 4Procs)
{
  if (stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  stk::middle_mesh::mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };
  auto mesh = create_mesh(spec, func);

  int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
  if (myrank == 0)
  {
    increment_vert(get_closest_vert(mesh, {0.5,  0.5}), {0.01, 0.0});
    increment_vert(get_closest_vert(mesh, {0.25, 0.5}), {0.0, 0.01});
  } else if (myrank == 1)
  {
    increment_vert(get_closest_vert(mesh, {0.75, 0.5}), {0.0, 0.01});
  }

  mesh::impl::RemoteCoordinateUpdator updator(mesh);
  updator.start_update();
  updator.finish_update();

  if (myrank == 0)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.51, 0.5}));
    EXPECT_TRUE(get_closest_vert(mesh, {0.25, 0.51}));
  } else if (myrank == 1)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.51, 0.5}));
  } else if (myrank == 2)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.51, 0.50}));
    EXPECT_TRUE(get_closest_vert(mesh, {0.25, 0.51}));
  } else if (myrank == 3)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.51, 0.50}));
    EXPECT_TRUE(get_closest_vert(mesh, {0.75, 0.51}));
  }


  if (myrank == 0)
  {
    increment_vert(get_closest_vert(mesh, {0.51,  0.5}), {0.01, 0.0});
    increment_vert(get_closest_vert(mesh, {0.25, 0.51}), {0.0, 0.01});
  } else if (myrank == 1)
  {
    increment_vert(get_closest_vert(mesh, {0.75, 0.51}), {0.0, 0.01});
  }

  updator.start_update();
  updator.finish_update();

  if (myrank == 0)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.52, 0.5}));
    EXPECT_TRUE(get_closest_vert(mesh, {0.25, 0.52}));
  } else if (myrank == 1)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.52, 0.5}));
  } else if (myrank == 2)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.52, 0.50}));
    EXPECT_TRUE(get_closest_vert(mesh, {0.25, 0.52}));
  } else if (myrank == 3)
  {
    EXPECT_TRUE(get_closest_vert(mesh, {0.52, 0.50}));
    EXPECT_TRUE(get_closest_vert(mesh, {0.75, 0.52}));
  }  

}