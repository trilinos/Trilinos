#include "gtest/gtest.h"
#include "stk_middle_mesh/destination_field_synchronizer.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/utils.hpp"

using namespace stk::middle_mesh;

namespace {

mesh::MeshEntityPtr find_closest_entity(std::shared_ptr<mesh::Mesh> mesh, int dim, const utils::Point& centroid)
{
  double minDist = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEntity = nullptr;
  for (auto& entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      utils::Point pt = mesh::compute_centroid(entity);
      utils::Point disp = centroid - pt;
      double dist = std::sqrt(dot(disp, disp));
      if (dist < minDist)
      {
        minDist = dist;
        minEntity = entity;
      }
    }

  return minEntity;
}

std::vector<int> get_values(mesh::VariableSizeFieldPtr<int> fieldPtr, mesh::MeshEntityPtr entity, int node)
{
  auto& field = *fieldPtr;
  std::vector<int> vals(field(entity, node).begin(), field(entity, node).end());
  return vals;
}

void expect_eq(const std::vector<int>& lhs, const std::vector<int>& rhs)
{
  EXPECT_EQ(lhs, rhs);
}

}

TEST(DestinationFieldSynchronizer, 2Procs)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.xmin = 0;   spec.ymin = 0;
  spec.xmax = 1;   spec.ymax = 1;
  spec.numelX = 2; spec.numelY = 2;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = mesh::impl::create_mesh(spec, f);

  auto scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
  if (myrank == 0)
  {
    auto el1 = find_closest_entity(mesh, 2, {0.25, 0.25});
    auto el2 = find_closest_entity(mesh, 2, {0.25, 0.75});
    scatterSpec->add_destination(el1, 0);
    scatterSpec->add_destination(el2, 1);
  } else
  {
    auto el1 = find_closest_entity(mesh, 2, {0.75, 0.25});
    auto el2 = find_closest_entity(mesh, 2, {0.75, 0.75});
    scatterSpec->add_destination(el1, 2);
    scatterSpec->add_destination(el2, 3);    
  }

  mesh::impl::DestinationFieldSynchronizer gatherer(mesh, scatterSpec);
  auto gatheredDestField = gatherer.synchronize();

  if (myrank == 0)
  {
    auto v1 = find_closest_entity(mesh, 0, {0,   0,   0});
    auto v2 = find_closest_entity(mesh, 0, {0.5, 0,   0});
    auto v3 = find_closest_entity(mesh, 0, {0,   0.5, 0});
    auto v4 = find_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto v5 = find_closest_entity(mesh, 0, {0,   1,   0});
    auto v6 = find_closest_entity(mesh, 0, {0.5, 1,   0});
    auto el1 = find_closest_entity(mesh, 2, {0.25, 0.25, 0});
    auto el2 = find_closest_entity(mesh, 2, {0.25, 0.75, 0});

    expect_eq(get_values(gatheredDestField, v1, 0), {0});
    expect_eq(get_values(gatheredDestField, v2, 0), {0, 2});
    expect_eq(get_values(gatheredDestField, v3, 0), {0, 1});
    expect_eq(get_values(gatheredDestField, v4, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v5, 0), {1});
    expect_eq(get_values(gatheredDestField, v6, 0), {1, 3});

    expect_eq(get_values(gatheredDestField, el1->get_down(0), 0), {0});
    expect_eq(get_values(gatheredDestField, el1->get_down(1), 0), {0, 2});
    expect_eq(get_values(gatheredDestField, el1->get_down(2), 0), {0, 1});
    expect_eq(get_values(gatheredDestField, el1->get_down(3), 0), {0});

    expect_eq(get_values(gatheredDestField, el2->get_down(1), 0), {1, 3});
    expect_eq(get_values(gatheredDestField, el2->get_down(2), 0), {1});
    expect_eq(get_values(gatheredDestField, el2->get_down(3), 0), {1});    

  } else
  {
    auto v1 = find_closest_entity(mesh, 0, {0.5, 0,   0});
    auto v2 = find_closest_entity(mesh, 0, {1.0, 0,   0});
    auto v3 = find_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto v4 = find_closest_entity(mesh, 0, {1.0, 0.5, 0});
    auto v5 = find_closest_entity(mesh, 0, {0.5, 1,   0});
    auto v6 = find_closest_entity(mesh, 0, {1.0, 1,   0});
    auto el1 = find_closest_entity(mesh, 2, {0.75, 0.25, 0});
    auto el2 = find_closest_entity(mesh, 2, {0.75, 0.75, 0});


    expect_eq(get_values(gatheredDestField, v1, 0), {0, 2});
    expect_eq(get_values(gatheredDestField, v2, 0), {2});
    expect_eq(get_values(gatheredDestField, v3, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v4, 0), {2, 3});
    expect_eq(get_values(gatheredDestField, v5, 0), {1, 3});
    expect_eq(get_values(gatheredDestField, v6, 0), {3});

    expect_eq(get_values(gatheredDestField, el1->get_down(0), 0), {2});
    expect_eq(get_values(gatheredDestField, el1->get_down(1), 0), {2});
    expect_eq(get_values(gatheredDestField, el1->get_down(2), 0), {2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(3), 0), {0, 2});

    expect_eq(get_values(gatheredDestField, el2->get_down(1), 0), {3});
    expect_eq(get_values(gatheredDestField, el2->get_down(2), 0), {3});
    expect_eq(get_values(gatheredDestField, el2->get_down(3), 0), {1, 3});
  }
}

TEST(DestinationFieldSynchronizer, 4ProcsMultipleDestinations)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.xmin = 0;   spec.ymin = 0;
  spec.xmax = 1;   spec.ymax = 1;
  spec.numelX = 2; spec.numelY = 2;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = mesh::impl::create_mesh(spec, f);

  auto scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
  if (myrank == 0)
  {
    auto el1 = find_closest_entity(mesh, 2, {0.25, 0.25});
    scatterSpec->add_destination(el1, 0);
    scatterSpec->add_destination(el1, 1);
  } else if (myrank == 1)
  {
    auto el1 = find_closest_entity(mesh, 2, {0.75, 0.25});
    scatterSpec->add_destination(el1, 1);
    scatterSpec->add_destination(el1, 2);    
  } else if (myrank == 2)
  {
    auto el1 = find_closest_entity(mesh, 2, {0.25, 0.75});
    scatterSpec->add_destination(el1, 2);
    scatterSpec->add_destination(el1, 3);      
  } else if (myrank == 3)
  {
    auto el1 = find_closest_entity(mesh, 2, {0.75, 0.75});
    scatterSpec->add_destination(el1, 3);
    scatterSpec->add_destination(el1, 1);      
  }

  mesh::impl::DestinationFieldSynchronizer gatherer(mesh, scatterSpec);
  auto gatheredDestField = gatherer.synchronize();

  if (myrank == 0)
  {
    auto v1 = find_closest_entity(mesh, 0, {0,   0,   0});
    auto v2 = find_closest_entity(mesh, 0, {0.5, 0,   0});
    auto v3 = find_closest_entity(mesh, 0, {0,   0.5, 0});
    auto v4 = find_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto el1 = find_closest_entity(mesh, 2, {0.25, 0.25, 0});

    expect_eq(get_values(gatheredDestField, v1, 0), {0, 1});
    expect_eq(get_values(gatheredDestField, v2, 0), {0, 1, 2});
    expect_eq(get_values(gatheredDestField, v3, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v4, 0), {0, 1, 2, 3});

    expect_eq(get_values(gatheredDestField, el1->get_down(0), 0), {0, 1});
    expect_eq(get_values(gatheredDestField, el1->get_down(1), 0), {0, 1, 2});
    expect_eq(get_values(gatheredDestField, el1->get_down(2), 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(3), 0), {0, 1});
  } else if (myrank == 1)
  {
    auto v1 = find_closest_entity(mesh, 0, {0.5, 0,   0});
    auto v2 = find_closest_entity(mesh, 0, {1.0, 0,   0});
    auto v3 = find_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto v4 = find_closest_entity(mesh, 0, {1.0, 0.5, 0});
    auto el1 = find_closest_entity(mesh, 2, {0.75, 0.25, 0});

    expect_eq(get_values(gatheredDestField, v1, 0), {0, 1, 2});
    expect_eq(get_values(gatheredDestField, v2, 0), {1, 2});
    expect_eq(get_values(gatheredDestField, v3, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v4, 0), {1, 2, 3});

    expect_eq(get_values(gatheredDestField, el1->get_down(0), 0), {1, 2});
    expect_eq(get_values(gatheredDestField, el1->get_down(1), 0), {1, 2});
    expect_eq(get_values(gatheredDestField, el1->get_down(2), 0), {1, 2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(3), 0), {0, 1, 2});
  } else if (myrank == 2)
  {
    auto v1 = find_closest_entity(mesh, 0, {0.0, 0.5, 0});
    auto v2 = find_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto v3 = find_closest_entity(mesh, 0, {0.0, 1.0, 0});
    auto v4 = find_closest_entity(mesh, 0, {0.5, 1.0, 0});
    auto el1 = find_closest_entity(mesh, 2, {0.25, 0.75, 0});

    expect_eq(get_values(gatheredDestField, v1, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v2, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v3, 0), {2, 3});
    expect_eq(get_values(gatheredDestField, v4, 0), {1, 2, 3});

    expect_eq(get_values(gatheredDestField, el1->get_down(0), 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(1), 0), {1, 2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(2), 0), {2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(3), 0), {2, 3});    
  } else if (myrank == 3)
  {
    auto v1 = find_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto v2 = find_closest_entity(mesh, 0, {1.0, 0.5, 0});
    auto v3 = find_closest_entity(mesh, 0, {0.5, 1.0, 0});
    auto v4 = find_closest_entity(mesh, 0, {1.0, 1.0, 0});
    auto el1 = find_closest_entity(mesh, 2, {0.75, 0.75, 0});

    expect_eq(get_values(gatheredDestField, v1, 0), {0, 1, 2, 3});
    expect_eq(get_values(gatheredDestField, v2, 0), {1, 2, 3});
    expect_eq(get_values(gatheredDestField, v3, 0), {1, 2, 3});
    expect_eq(get_values(gatheredDestField, v4, 0), {1, 3});

    expect_eq(get_values(gatheredDestField, el1->get_down(0), 0), {1, 2, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(1), 0), {1, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(2), 0), {1, 3});
    expect_eq(get_values(gatheredDestField, el1->get_down(3), 0), {1, 2, 3});     
  }
}