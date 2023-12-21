#include "gtest/gtest.h"
#include "stk_middle_mesh/entity_sorted_by_owner.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/utils.hpp"

using namespace stk::middle_mesh;

TEST(EntitySortedByOwner, Access)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };
  std::shared_ptr<mesh::Mesh> mesh = create_mesh(spec, func);

  int ownerCommSize = 2;
  mesh::impl::EntitySortedByOwner container(ownerCommSize);

  mesh::MeshEntityPtr v0 = mesh->get_vertices()[0];
  mesh::MeshEntityPtr v1 = mesh->get_vertices()[1];
  mesh::MeshEntityPtr v2 = mesh->get_vertices()[2];

  mesh::MeshEntityPtr v3 = mesh->get_vertices()[3];
  mesh::MeshEntityPtr v4 = mesh->get_vertices()[4];
  mesh::MeshEntityPtr v5 = mesh->get_vertices()[5];  

  container.insert({0, 0}, v1);
  container.insert({0, 1}, v0);
  container.insert({0, 2}, v2);

  EXPECT_EQ(container.get_value({0, 2}), v2);
  EXPECT_EQ(container.get_value({0, 1}), v0);
  EXPECT_EQ(container.get_value({0, 0}), v1);

  container.insert({1, 1}, v4);
  container.insert({1, 3}, v3);
  container.insert({1, 5}, v5);

  EXPECT_EQ(container.get_value({1, 5}), v5);
  EXPECT_EQ(container.get_value({1, 3}), v3);
  EXPECT_EQ(container.get_value({1, 1}), v4);  
  EXPECT_EQ(container.get_value({1, 0}), nullptr);
  EXPECT_EQ(container.get_value({1, 2}), nullptr);
  EXPECT_EQ(container.get_value({1, 4}), nullptr);
  EXPECT_EQ(container.get_value({1, 6}), nullptr);
}