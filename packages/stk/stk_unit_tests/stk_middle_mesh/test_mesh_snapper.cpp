#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_snapper.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

namespace {

void test_meshes_identical(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2)
{
  auto& verts1 = mesh1->get_vertices();
  auto& verts2 = mesh2->get_vertices();

  EXPECT_EQ(verts1.size(), verts2.size());

  for (unsigned int i = 0; i < verts1.size(); ++i)
  {
    auto v1 = verts1[i];
    auto v2 = verts2[i];

    EXPECT_FLOAT_EQ(v1->get_point_orig(0).x, v2->get_point_orig(0).x);
    EXPECT_FLOAT_EQ(v1->get_point_orig(0).y, v2->get_point_orig(0).y);
  }
}

/*
void expect_float_eq(const utils::Point& pt1, const utils::Point& pt2)
{
  EXPECT_FLOAT_EQ(pt1.x, pt2.x);
  EXPECT_FLOAT_EQ(pt1.y, pt2.y);
  EXPECT_FLOAT_EQ(pt1.z, pt2.z);
}
*/

} // namespace

TEST(MeshSnapper, Identical)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec, spec2;
  spec.numelX = 5;
  spec.numelY = 5;
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

  //  auto func = [&](const utils::Point& pt) { return pt;};
  // rotate to x-z plane
  auto func = [&](const utils::Point& pt) { return utils::Point(pt.x, 0, pt.y); };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

  MeshSnapper snapper({1e-8});
  snapper.snap(mesh1, mesh2);

  test_meshes_identical(mesh1, mesh2);
}

TEST(MeshSnapper, Offset)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  double tol = 1e-8;
  MeshSpec spec, spec2;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 5;
  spec2.numelY = 5;
  spec2.xmin   = 0 + tol / (spec.numelX * 2);
  spec2.xmax   = 1 + tol / (spec.numelX * 2);
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  //  auto func = [&](const utils::Point& pt) { return pt;};
  // rotate to x-z plane
  auto func = [&](const utils::Point& pt) { return utils::Point(pt.x, 0, pt.y); };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<Mesh> mesh2 = create_mesh(spec2, func);

  MeshSnapper snapper({tol});
  snapper.snap(mesh1, mesh2);

  test_meshes_identical(mesh1, mesh2);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
