#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_agglomerator.hpp"
#include <unordered_map>
#include <unordered_set>

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {
using namespace mesh;
using namespace mesh::impl;

class MeshAgglomeratorTester : public ::testing::Test
{
  protected:
    MeshAgglomeratorTester()
    {
      MeshSpec spec;
      spec.numelX = 4;
      spec.numelY = 4;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;

      auto func = [&](const utils::Point& pt) { return pt; };

      mesh1 = create_mesh(spec, func);

      // get the partition edges
      double xval = 0.5;
      double yval = 0.5;
      for (auto& edge : mesh1->get_edges())
      {
        auto pt1 = edge->get_down(0)->get_point_orig(0);
        auto pt2 = edge->get_down(1)->get_point_orig(0);
        if ((std::abs(pt1.x - xval) < 1e-13 && std::abs(pt2.x - xval) < 1e-13) ||
            (std::abs(pt1.y - yval) < 1e-13 && std::abs(pt2.y - yval) < 1e-13))
          partitionEdges.insert(edge);
      }

      auto funcE = [this](MeshEntityPtr edge) { return partitionEdges.count(edge) > 0; };

      agg = std::make_shared<MeshAgglomerator>(mesh1, funcE);
    }

    utils::Point get_lower_left_corner(MeshEntityPtr el)
    {
      utils::Point centroid = mesh::compute_centroid(el);

      if (centroid.x < 0.5 && centroid.y < 0.5)
      {
        return {0, 0};
      } else if (centroid.x > 0.5 && centroid.y < 0.5)
      {
        return {0.5, 0};
      } else if (centroid.x < 0.5 && centroid.y > 0.5)
      {
        return {0, 0.5};
      } else
      {
        return {0.5, 0.5};
      }
    }

    std::shared_ptr<Mesh> mesh1;
    std::unordered_set<MeshEntityPtr> partitionEdges;
    std::shared_ptr<MeshAgglomerator> agg;
};


MeshEntityPtr get_entity(std::shared_ptr<Mesh> mesh, int dim, const utils::Point& pt, double tol=1e-13)
{
  double closestDist = tol;
  MeshEntityPtr closestEntity = nullptr;
  for (auto entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      utils::Point centroid = mesh::compute_centroid(entity);
      double dist = std::sqrt(dot(centroid - pt, centroid - pt));
      if (dist <= closestDist)
      {
        closestDist = dist;
        closestEntity = entity;
      }
    }

  return closestEntity;
}


}

TEST_F(MeshAgglomeratorTester, GroupCreation)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  EXPECT_EQ(agg->get_num_groups(), 4);
  for (int i = 0; i < agg->get_num_groups(); ++i)
  {
    auto elsI = agg->get_group_elements(i);
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

  for (int i=0; i < agg->get_num_groups(); ++i)
  {
    auto els_i = agg->get_group_elements(i);
    EXPECT_EQ(els_i.size(), 4u);

    utils::Point lowerLeftCorner = get_lower_left_corner(*(els_i.begin()));
    std::vector<utils::Point> centroidOffsets = { {0.125, 0.125}, {0.375, 0.125}, {0.125, 0.375}, {0.375, 0.375}};
    for (int j=0; j < 4; ++j)
    {
      auto centroid = lowerLeftCorner + centroidOffsets[j];
      auto pred = [&](const MeshEntityPtr& entity)
      {
        utils::Point elCentroid = mesh::compute_centroid(entity);
        double dist = std::sqrt(dot(elCentroid - centroid, elCentroid - centroid));
        return dist < 1e-13;
      };

      auto it = std::find_if(els_i.begin(), els_i.end(), pred);
      EXPECT_NE(it, els_i.end());
    }
  }
}

TEST_F(MeshAgglomeratorTester, Verts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();
  
  EXPECT_EQ(agg->get_num_groups(), 4);

  for (int i=0; i < agg->get_num_groups(); ++i)
  {
    auto& groupEls = agg->get_group_elements(i);
    auto& groupVerts = agg->get_group_verts(i);
    std::set<MeshEntityPtr, MeshEntityCompare> expectedVerts, actualVerts(groupVerts.begin(), groupVerts.end());

    for (auto el : groupEls)
    {
      std::array<MeshEntityPtr, mesh::MAX_DOWN> verts;
      int ndown = mesh::get_downward(el, 0, verts.data());
      for (int j=0; j < ndown; ++j)
      {
        expectedVerts.insert(verts[j]);
      }
    }

    EXPECT_EQ(expectedVerts.size(), actualVerts.size());
    for (auto vert : expectedVerts)
    {
      EXPECT_EQ(actualVerts.count(vert), 1u);
    }
  }
}

TEST_F(MeshAgglomeratorTester, GetGroupFromVerts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshEntityPtr v1 = get_entity(mesh1, 0, {0, 0});
  MeshEntityPtr v2 = get_entity(mesh1, 0, {0.5, 0});
  MeshEntityPtr v3 = get_entity(mesh1, 0, {0, 0.5});
  MeshAgglomerator::SetType<MeshEntityPtr> verts = {v1, v2, v3};
  std::vector<int> groupIdxs = agg->get_group_idxs(verts, 2);

  EXPECT_EQ(groupIdxs.size(), 1u);
  utils::Point lowerLeftCorner = get_lower_left_corner(*(agg->get_group_elements(groupIdxs[0]).begin()));
  EXPECT_FLOAT_EQ(lowerLeftCorner.x, 0.0);
  EXPECT_FLOAT_EQ(lowerLeftCorner.y, 0.0);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
