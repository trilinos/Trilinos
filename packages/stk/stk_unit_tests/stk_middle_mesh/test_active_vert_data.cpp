#include "active_vert_data.hpp"
#include "create_mesh.hpp"
#include "mesh.hpp"
#include "gtest/gtest.h"

#include <limits>
#include <unordered_map>

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

mesh::MeshEntityPtr find_closest_vert(std::shared_ptr<mesh::Mesh> mesh, const utils::Point& pt)
{
  double minDist              = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minVert = nullptr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      auto pt2    = vert->get_point_orig(0);
      auto disp   = pt - pt2;
      double dist = dot(disp, disp);
      if (dist < minDist)
      {
        minDist = dist;
        minVert = vert;
      }
    }

  return minVert;
}

} // namespace

TEST(ActiveVertData, Mesh2x2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 2;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);

  mesh::MeshEntityPtr vert = find_closest_vert(mesh1, utils::Point(0.5, 0.5));

  opt::impl::ActiveVertData active(vert);

  EXPECT_EQ(active.get_num_verts(), 9);
  EXPECT_EQ(active.get_num_elements(), 12);

  // check unique verts
  std::unordered_map<mesh::MeshEntityPtr, int> vertCounts;
  for (auto& v : mesh1->get_vertices())
    vertCounts[v] = 0;

  /*
  vert_counts[vert] = 0;
  for (int i=0; i < vert->count_up(); ++i)
  {
    mesh::MeshEntityPtr edge = vert->get_up(i);
    mesh::MeshEntityPtr other_vert = edge->get_down(0) == vert ? edge->get_down(1) : edge->get_down(0);
    vert_counts[other_vert] = 0;
  }
  // add corner verts
  auto v = findClosestVert(mesh1, utils::Point(1, 0));
  vert_counts[v] = 0;

  v = findClosestVert(mesh1, utils::Point(0, 1));
  vert_counts[v] = 0;
  */

  for (auto& v : active.get_unique_verts())
  {
    EXPECT_EQ(vertCounts.count(v), static_cast<unsigned int>(1));
    vertCounts[v] += 1;
  }
  for (auto& p : vertCounts)
    EXPECT_EQ(p.second, 1);

  // check all triangle verts are in unique verts
  for (int i = 0; i < active.get_num_elements(); ++i)
  {
    auto vertsI = active.get_element_verts(i);
    for (int j = 0; j < 3; ++j)
      EXPECT_EQ(vertCounts.count(vertsI[j]), static_cast<unsigned int>(1));
  }

  // check current vertex
  EXPECT_EQ(active.get_current_vert(), vert);
  for (int i = 0; i < active.get_num_elements(); ++i)
  {
    auto vertsI = active.get_element_verts(i);
    int idx     = active.get_current_vert_idx(i);
    EXPECT_EQ(vertsI[idx], vert);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
