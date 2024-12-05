#include "gtest/gtest.h"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_middle_mesh/active_vert_container.hpp"
#include "stk_middle_mesh/boundary_fixture.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/point.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "stk_middle_mesh/create_sharing_from_verts.hpp"

namespace {

using stk::middle_mesh::opt::impl::ActiveVertData;
using stk::middle_mesh::utils::Point;
using namespace stk::middle_mesh;

class AlwaysTrue
{
  public:
    bool operator()(stk::middle_mesh::mesh::MeshEntityPtr) { return true; }
};

ActiveVertData& get_closest_patch(stk::middle_mesh::mesh::impl::ActiveVertContainer& container, const Point& pt)
{
  double minDist = std::numeric_limits<double>::max();
  ActiveVertData* minPatch = nullptr;
  for (auto& patch : container.get_active_verts())
  {
    auto disp = patch.get_current_vert()->get_point_orig(0) - pt;
    double dist = std::sqrt(dot(disp, disp));
    if (dist < minDist)
    {
      minDist = dist;
      minPatch = &patch;
    }
  }

  STK_ThrowRequireMsg(minPatch != nullptr, "get_closest_patch failed to find patch");
  return *minPatch;
}

bool has_vert(ActiveVertData& patch, const Point& pt, double tol)
{
  for (auto& patchPoint : patch.get_unique_verts())
  {
    auto disp = patchPoint - pt;
    double dist = std::sqrt(dot(disp, disp));
    if (dist < tol)
      return true;
  }

  return false;
}

stk::middle_mesh::mesh::MeshEntityPtr find_closest_vert(std::shared_ptr<stk::middle_mesh::mesh::Mesh> mesh, const Point& pt, double tol)
{
  stk::middle_mesh::mesh::MeshEntityPtr entity = nullptr;
  double min_dist = tol;
  for (auto& e : mesh->get_vertices())
  {
    auto disp = e->get_point_orig(0) - pt;
    double dist = std::sqrt(dot(disp, disp));
    if (dist < min_dist)
    {
      min_dist = dist;
      entity = e;
    }
  }

  return entity;
}

void increment_vert(stk::middle_mesh::mesh::MeshEntityPtr vert, const Point& pt)
{
  auto vertPt = vert->get_point_orig(0);
  vert->set_point_orig(0, vertPt + pt);
}

void expect_near(const Point& lhs, const Point& rhs, double tol)
{
  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(lhs[i], rhs[i], tol);
}

}

TEST(ActiveVertContainer, PointsOrig)
{
  if (stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  stk::middle_mesh::mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const Point& pt) { return pt; };
  auto mesh = create_mesh(spec, func);
  stk::middle_mesh::mesh::impl::BoundaryFixture filter(mesh);

  stk::middle_mesh::mesh::impl::ActiveVertContainer container(mesh, filter, AlwaysTrue());

  for (auto& patch : container.get_active_verts())
  {
    auto& pts = patch.get_unique_verts();
    auto& ptsOrig = patch.get_points_orig();

    EXPECT_EQ(pts.size(), ptsOrig.size());
    for (size_t i=0; i < pts.size(); ++i)
      expect_near(pts[i], ptsOrig[i], 1e-13);
  }
}

TEST(ActiveVertContainer, 2Procs)
{
  if (stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  stk::middle_mesh::mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const Point& pt) { return pt; };
  auto mesh = create_mesh(spec, func);
  stk::middle_mesh::mesh::impl::BoundaryFixture filter(mesh);

  stk::middle_mesh::mesh::impl::ActiveVertContainer container(mesh, filter, AlwaysTrue());

  int myrank = stk::middle_mesh::utils::impl::comm_rank(MPI_COMM_WORLD);

  for (auto& patch : container.get_active_verts())
  {
    EXPECT_EQ(patch.get_num_verts(), 9);
  }

  double delta = 0.25;
  std::vector<Point> centerPts;
  if (myrank == 0)
  {
    centerPts = { {0.25, 0.25}, {0.5, 0.25}, {0.25, 0.5}, {0.5, 0.5}, {0.25, 0.75}, {0.5, 0.75}};
  } else
  {
    centerPts = { {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75} };
  }

  EXPECT_EQ(container.get_active_verts().size(), centerPts.size());
  for (auto& centerPt : centerPts)
  {
    std::cout << "\nchecking patch centered at " << centerPt << std::endl;
    ActiveVertData& patch1 = get_closest_patch(container, centerPt);
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(-delta, -delta), 1e-12));
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(0,      -delta), 1e-12));
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(delta,  -delta), 1e-12));

    EXPECT_TRUE(has_vert(patch1, centerPt + Point(-delta, 0), 1e-12));
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(0,      0), 1e-12));
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(delta,  0), 1e-12));   

    EXPECT_TRUE(has_vert(patch1, centerPt + Point(-delta, delta), 1e-12));
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(0,      delta), 1e-12));
    EXPECT_TRUE(has_vert(patch1, centerPt + Point(delta,  delta), 1e-12));         
  }

  if (myrank == 0)
  {
    increment_vert(find_closest_vert(mesh, {0.5, 0},    1e-12), {0.01, 0}); 
    increment_vert(find_closest_vert(mesh, {0.5, 0.25}, 1e-12), {0.02, 0});
    increment_vert(find_closest_vert(mesh, {0.5, 0.50}, 1e-12), {0.03, 0});
    increment_vert(find_closest_vert(mesh, {0.5, 0.75}, 1e-12), {0.04, 0});
    increment_vert(find_closest_vert(mesh, {0.5, 1.00}, 1e-12), {0.05, 0});
  } else
  {
    increment_vert(find_closest_vert(mesh, {0.75, 0.00}, 1e-12), {0.01, 0}); 
    increment_vert(find_closest_vert(mesh, {0.75, 0.25}, 1e-12), {0.02, 0}); 
    increment_vert(find_closest_vert(mesh, {0.75, 0.50}, 1e-12), {0.03, 0}); 
    increment_vert(find_closest_vert(mesh, {0.75, 0.75}, 1e-12), {0.04, 0}); 
    increment_vert(find_closest_vert(mesh, {0.75, 1.00}, 1e-12), {0.05, 0});
  }

  container.update_remote_coords();

  if (myrank == 0)
  {
    ActiveVertData& patch1 = get_closest_patch(container, {0.52, 0.25});
    EXPECT_TRUE(has_vert(patch1, {0.76, 0.00}, 1e-12));
    EXPECT_TRUE(has_vert(patch1, {0.77, 0.25}, 1e-12));
    EXPECT_TRUE(has_vert(patch1, {0.78, 0.50}, 1e-12));

    ActiveVertData& patch2 = get_closest_patch(container, {0.53, 0.50});
    EXPECT_TRUE(has_vert(patch2, {0.77, 0.25}, 1e-12));
    EXPECT_TRUE(has_vert(patch2, {0.78, 0.50}, 1e-12));  
    EXPECT_TRUE(has_vert(patch2, {0.79, 0.75}, 1e-12));

    ActiveVertData& patch3 = get_closest_patch(container, {0.54, 0.75});  
    EXPECT_TRUE(has_vert(patch3, {0.78, 0.50}, 1e-12));  
    EXPECT_TRUE(has_vert(patch3, {0.79, 0.75}, 1e-12));
    EXPECT_TRUE(has_vert(patch3, {0.80, 1.00}, 1e-12));
  } else
  { 
    ActiveVertData& patch1 = get_closest_patch(container, {0.75, 0.25});    
    EXPECT_TRUE(has_vert(patch1, {0.51, 0},    1e-12));
    EXPECT_TRUE(has_vert(patch1, {0.52, 0.25}, 1e-12));
    EXPECT_TRUE(has_vert(patch1, {0.53, 0.50}, 1e-12));

    ActiveVertData& patch2 = get_closest_patch(container, {0.75, 0.50});
    EXPECT_TRUE(has_vert(patch2, {0.52, 0.25}, 1e-12));
    EXPECT_TRUE(has_vert(patch2, {0.53, 0.50}, 1e-12));
    EXPECT_TRUE(has_vert(patch2, {0.54, 0.75}, 1e-12));

    ActiveVertData& patch3 = get_closest_patch(container, {0.75, 0.75});
    EXPECT_TRUE(has_vert(patch3, {0.53, 0.50}, 1e-12));
    EXPECT_TRUE(has_vert(patch3, {0.54, 0.75}, 1e-12));
    EXPECT_TRUE(has_vert(patch3, {0.55, 1.00}, 1e-12));    
  }
}

TEST(ActiveVertContainer, VertNotConnectedByLocalElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  auto mesh = mesh::make_empty_mesh(MPI_COMM_WORLD);

  int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
  if (myrank == 0)
  {
    auto v1 = mesh->create_vertex({0,    0.5, 0});
    auto v2 = mesh->create_vertex({0.5,  0.5, 0});
    auto v3 = mesh->create_vertex({0.75, 0.5, 0});
    auto v4 = mesh->create_vertex({0,    1.0, 0});
    auto v5 = mesh->create_vertex({0.5,  1.0, 0});
    auto v6 = mesh->create_vertex({1.0,  1.0, 0});
    auto v7 = mesh->create_vertex({1.25, 0.5, 0});
    auto v8 = mesh->create_vertex({1.0,  0.0, 0});

    mesh->create_quad_from_verts(v1, v2, v5, v4);
    mesh->create_quad_from_verts(v2, v3, v6, v5);
    mesh->create_quad_from_verts(v8, v7, v6, v3);

    v1->add_remote_shared_entity({1, 3});
    v2->add_remote_shared_entity({1, 4});
    v3->add_remote_shared_entity({1, 5});
    v8->add_remote_shared_entity({1, 2});    

  } else
  {
    auto v1 = mesh->create_vertex({0,    0.0, 0});
    auto v2 = mesh->create_vertex({0.5,  0.0, 0});
    auto v3 = mesh->create_vertex({1.0,  0.0, 0});

    auto v4 = mesh->create_vertex({0,    0.5, 0});
    auto v5 = mesh->create_vertex({0.5,  0.5, 0});
    auto v6 = mesh->create_vertex({0.75, 0.5, 0});

    mesh->create_quad_from_verts(v1, v2, v5, v4);
    mesh->create_quad_from_verts(v2, v3, v6, v5);

    v4->add_remote_shared_entity({0, 0});
    v5->add_remote_shared_entity({0, 1});
    v6->add_remote_shared_entity({0, 2});
    v3->add_remote_shared_entity({0, 7});    
  }

  mesh::impl::CreateSharingFromVert edgeSharing(mesh);
  edgeSharing.create_sharing_from_verts();

  auto filter = [&](mesh::MeshEntityPtr vert) 
  { 
    if (myrank == 0)
    {
      return vert->get_id() == 1;
    } else
    {
      return vert->get_id() == 4;
    }
  };

  auto isValid = [](mesh::MeshEntityPtr vert) { return true; };
  mesh::impl::ActiveVertContainer activeVerts(mesh, filter, isValid);

  if (myrank == 0)
  {
    EXPECT_EQ(activeVerts.get_active_verts().size(), 1u);

    auto& activeVert = activeVerts.get_active_verts()[0];

    EXPECT_EQ(activeVert.get_num_local_verts(), 7);
    EXPECT_EQ(activeVert.get_num_verts(), 9);
  }
}
