#include "stk_middle_mesh/active_vert_data.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh.hpp"
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

void expect_near(const utils::Point& lhs, const utils::Point& rhs, double tol)
{
  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(lhs[i], rhs[i], tol);
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

  opt::impl::ActiveVertData active(mesh1, vert);

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

  const std::vector<utils::Point>& verts = active.get_unique_verts();
  for (auto& pt : verts)
  {
    mesh::MeshEntityPtr v = find_closest_vert(mesh1, pt);
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
    {
      EXPECT_EQ(vertCounts.count(find_closest_vert(mesh1, vertsI[j])), static_cast<unsigned int>(1));
    }
  }

  // check current vertex
  EXPECT_EQ(active.get_current_vert(), vert);
  for (int i = 0; i < active.get_num_elements(); ++i)
  {
    auto vertsI = active.get_element_verts(i);
    int idx     = active.get_current_vert_idx(i);
    for (int j=0; j < 3; ++j)
      EXPECT_EQ(vertsI[idx][j], vert->get_point_orig(0)[j]);
  }
}

TEST(ActiveVertData, RemoteElements)
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

  mesh::MeshEntityPtr centerVert = find_closest_vert(mesh1, utils::Point(0.5, 0.5));

  opt::impl::ActiveVertData active(mesh1, centerVert);

  int numLocalVerts = 9;
  int numLocalTris = 12;
  EXPECT_EQ(active.get_num_verts(), numLocalVerts);
  EXPECT_EQ(active.get_num_elements(), numLocalTris);

  active.add_remote_vert({1, 0});
  active.add_remote_vert({1, 1});
  active.add_remote_element({0, 9, 10});
  int numVerts = 11;
  int numTris = 13;

  active.get_unique_verts()[numLocalVerts]   = {1, 0, 0};
  active.get_unique_verts()[numLocalVerts+1] = {0, 1, 0};
  active.finish_initialization();

  EXPECT_EQ(active.get_num_verts(), numVerts);
  EXPECT_EQ(active.get_num_elements(), numTris);
  EXPECT_EQ(active.get_current_vert_idx(numTris-1), 0);  //TODO: is current vert indices getting updated?
  std::array<int, 3> expectedIds = {0, 9, 10}; 
  EXPECT_EQ(active.get_element_vert_ids(numTris-1), expectedIds);
  std::array<utils::Point, 3> expectedPts = { utils::Point{0.5, 0.5}, utils::Point{1, 0}, utils::Point{0, 1}};
  EXPECT_EQ(active.get_element_verts(numTris-1), expectedPts);
  expectedPts = std::array<utils::Point, 3>{ utils::Point{0.5, 0.5}, utils::Point{1, 0}, utils::Point{0, 1}};
  EXPECT_EQ(active.get_element_verts_orig(numTris-1), expectedPts);

  // modify local coordinates
  // check get_element_verts/get_element_verts/orig
  auto& localVerts = active.get_local_verts();
  for (auto& vert : localVerts)
    vert->set_point_orig(0, vert->get_point_orig(0) + utils::Point(0.5, 0.5));

  for (int tri=0; tri < numLocalTris; ++tri)
  {
    std::array<int, 3> vertIndices  = active.get_element_vert_ids(tri);
    std::array<utils::Point, 3> pts = active.get_element_verts(tri);

    for (int i=0; i < 3; ++i)
    {
      expect_near(pts[i], localVerts[vertIndices[i]]->get_point_orig(0), 1e-13);
    }
  }

  active.get_unique_verts()[numLocalVerts] = utils::Point(2, 0);
  active.get_unique_verts()[numLocalVerts+1] = utils::Point(0, 2);

  std::array<utils::Point, 3> pts      = active.get_element_verts(numTris-1);
  std::array<utils::Point, 3> pts_orig = active.get_element_verts_orig(numTris-1);

  expect_near(pts[1], {2, 0}, 1e-13);
  expect_near(pts[2], {0, 2}, 1e-13);
  expect_near(pts_orig[1], {1, 0}, 1e-13);
  expect_near(pts_orig[2], {0, 1}, 1e-13);


  for (int i=0; i < numLocalVerts; ++i)
    EXPECT_EQ(active.get_vert_owner(i), get_owner_remote(mesh1, localVerts[i]));

  mesh::RemoteSharedEntity remote1{1, 0};
  mesh::RemoteSharedEntity remote2{1, 1};
  EXPECT_EQ(active.get_vert_owner(numLocalVerts), remote1);
  EXPECT_EQ(active.get_vert_owner(numLocalVerts + 1), remote2);
}


TEST(ActiveVertData, ExtraLocalVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);

  mesh::MeshEntityPtr centerVert = find_closest_vert(mesh1, utils::Point(1.0/3, 1.0/3));

  opt::impl::ActiveVertData active(mesh1, centerVert);

  mesh::MeshEntityPtr extraLocalVert = find_closest_vert(mesh1, utils::Point(1.0, 1.0/3));
  active.add_local_vert(extraLocalVert);

  int numLocalVerts = 10;
  int numLocalTris = 12;
  EXPECT_EQ(active.get_num_verts(), numLocalVerts);
  EXPECT_EQ(active.get_num_local_verts(), numLocalVerts);
  EXPECT_EQ(active.get_num_elements(), numLocalTris);
  EXPECT_EQ(active.get_unique_verts().size(), size_t(numLocalVerts));
  EXPECT_EQ(active.get_local_verts().size(), size_t(numLocalVerts));
  EXPECT_EQ(active.get_local_verts().back(), extraLocalVert);
  EXPECT_EQ(active.get_points_orig().size(), size_t(numLocalVerts));
  EXPECT_EQ(active.get_points_orig().back().y, extraLocalVert->get_point_orig(0).y);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
