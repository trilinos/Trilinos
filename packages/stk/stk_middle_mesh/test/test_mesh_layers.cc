#include "gtest/gtest.h"

#include "create_mesh.h"
#include "mesh_layers.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

TEST(MeshLayers, NoFilter)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  MeshLayers layers;

  std::vector<MeshEntityPtr> roots, output;
  for (auto& v : mesh->get_vertices())
    if (v)
    {
      auto pt = v->get_point_orig(0);
      if (std::abs(pt.x) < 1e-13)
        roots.push_back(v);
    }

  auto filter = [](MeshEntityPtr) { return true; };

  // test one layer
  layers.get_layers(mesh, filter, roots, 1, output);

  EXPECT_EQ(output.size(), static_cast<unsigned int>(8));
  for (auto v : output)
  {
    double xCoord  = v->get_point_orig(0).x;
    bool isCorrect = std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13;
    EXPECT_TRUE(isCorrect);
  }

  // test two layer
  layers.get_layers(mesh, filter, roots, 2, output);

  EXPECT_EQ(output.size(), static_cast<unsigned int>(12));
  for (auto v : output)
  {
    double xCoord = v->get_point_orig(0).x;
    bool isCorrect =
        std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13 || std::abs(xCoord - 2.0 / 3.0) < 1e-13;
    EXPECT_TRUE(isCorrect);
  }

  // test three layers
  layers.get_layers(mesh, filter, roots, 3, output);
  EXPECT_EQ(output.size(), static_cast<unsigned int>(16));
  EXPECT_TRUE(is_unique(output));
}

TEST(MeshLayers, Filter)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // mesh1 and mesh2 are identical
  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  MeshLayers layers;

  std::vector<MeshEntityPtr> roots, output;
  MeshEntityPtr vExclude = nullptr;
  for (auto& v : mesh->get_vertices())
    if (v)
    {
      auto pt = v->get_point_orig(0);
      if (std::abs(pt.x) < 1e-13)
        roots.push_back(v);

      if (std::abs(pt.x - 1.0 / 3.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13)
        vExclude = v;
    }

  assert(vExclude);

  auto filter = [vExclude](MeshEntityPtr v) { return v != vExclude; };

  // test one layer
  layers.get_layers(mesh, filter, roots, 1, output);

  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(7));
  for (auto v : output)
  {
    double xCoord  = v->get_point_orig(0).x;
    double yCoord  = v->get_point_orig(0).y;
    bool isCorrect = std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13;

    if (std::abs(xCoord - 1.0 / 3.0) < 1e-13 && std::abs(yCoord - 1.0 / 3.0) < 1e-13)
      EXPECT_FALSE(true);
    else
      EXPECT_TRUE(isCorrect);
  }

  // test two layers
  layers.get_layers(mesh, filter, roots, 2, output);

  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(10));
  for (auto v : output)
  {
    double xCoord = v->get_point_orig(0).x;
    double yCoord = v->get_point_orig(0).y;
    bool isCorrect =
        std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13 || std::abs(xCoord - 2.0 / 3.0) < 1e-13;

    bool isExcluded = (std::abs(xCoord - 1.0 / 3.0) < 1e-13 && std::abs(yCoord - 1.0 / 3.0) < 1e-13) ||
                      (std::abs(xCoord - 2.0 / 3.0) < 1e-13 && std::abs(yCoord - 1.0 / 3.0) < 1e-13);

    if (isExcluded)
      EXPECT_FALSE(true);
    else
      EXPECT_TRUE(isCorrect);
  }

  // test three layers
  layers.get_layers(mesh, filter, roots, 3, output);

  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(14));

  std::set<MeshEntityPtr> verts;
  for (auto v : mesh->get_vertices())
  {
    auto pt         = v->get_point_orig(0);
    bool isExcluded = (std::abs(pt.x - 1.0 / 3.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13) ||
                      (std::abs(pt.x - 1.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13);
    if (!isExcluded)
      verts.insert(v);
  }

  for (auto v : output)
  {
    EXPECT_TRUE(verts.count(v) == 1);
  }

  // test four layers
  layers.get_layers(mesh, filter, roots, 4, output);
  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(15));

  verts.clear();
  for (auto v : mesh->get_vertices())
  {
    auto pt         = v->get_point_orig(0);
    bool isExcluded = (std::abs(pt.x - 1.0 / 3.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13);
    if (!isExcluded)
      verts.insert(v);
  }

  for (auto v : output)
  {
    EXPECT_TRUE(verts.count(v) == 1);
  }

  // test all layers
  layers.get_all_layers(mesh, filter, roots, output);
  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(15));

  for (auto v : output)
  {
    EXPECT_TRUE(verts.count(v) == 1);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
