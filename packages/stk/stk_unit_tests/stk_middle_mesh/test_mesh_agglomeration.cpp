#include "gtest/gtest.h"

#include "create_mesh.hpp"
#include "mesh_agglomerator.hpp"
#include <unordered_map>
#include <unordered_set>

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

TEST(MeshAgglomerator, All)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);

  // get the partition edges
  double xval = 0.5;
  double yval = 0.5;
  std::unordered_set<MeshEntityPtr> partitionEdges;
  for (auto& edge : mesh1->get_edges())
  {
    auto pt1 = edge->get_down(0)->get_point_orig(0);
    auto pt2 = edge->get_down(1)->get_point_orig(0);
    if ((std::abs(pt1.x - xval) < 1e-13 && std::abs(pt2.x - xval) < 1e-13) ||
        (std::abs(pt1.y - yval) < 1e-13 && std::abs(pt2.y - yval) < 1e-13))
      partitionEdges.insert(edge);
  }

  for (auto& edge : partitionEdges)
    std::cout << "partition edge = " << edge << std::endl;

  auto funcE = [partitionEdges](MeshEntityPtr edge) { return partitionEdges.count(edge) > 0; };

  MeshAgglomerator agg(mesh1, funcE);

  EXPECT_EQ(agg.get_num_groups(), 4);
  for (int i = 0; i < agg.get_num_groups(); ++i)
  {
    auto elsI = agg.get_group_elements(i);
    // check that each partition edge only appears once
    std::unordered_map<MeshEntityPtr, int> edgeCounts;
    for (auto& el : elsI)
      for (int j = 0; j < el->count_down(); ++j)
        edgeCounts[el->get_down(j)] += 1;

    EXPECT_EQ(elsI.size(), static_cast<unsigned int>(4));
    for (auto& p : edgeCounts)
      if (partitionEdges.count(p.first) > 0)
      {
        EXPECT_TRUE(edgeCounts[p.first] <= 1);
      }
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
