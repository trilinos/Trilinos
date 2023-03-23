#include "mesh.hpp"
#include "variable_size_field.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

TEST(VariableSizeField, ValuesOneNode)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);
  auto v3 = mesh->create_vertex(2, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v2, 0, 1);
  field.insert(v2, 0, 2);
  field.insert(v1, 0, 3);
  field.insert(v3, 0, 4);
  field.insert(v3, 0, 5);
  field.insert(v3, 0, 6);
  field.insert(v1, 0, 7);
  field.insert(v1, 0, 8);
  field.insert(v1, 0, 9);

  EXPECT_EQ(field.get_num_nodes(0), 1);
  EXPECT_EQ(field.get_num_nodes(1), 0);
  EXPECT_EQ(field.get_num_nodes(2), 0);
  EXPECT_EQ(field.get_num_comp(v1, 0), 4);
  EXPECT_EQ(field.get_num_comp(v2, 0), 2);
  EXPECT_EQ(field.get_num_comp(v3, 0), 3);

  EXPECT_EQ(field(v1, 0, 0), 3);
  EXPECT_EQ(field(v1, 0, 1), 7);
  EXPECT_EQ(field(v1, 0, 2), 8);
  EXPECT_EQ(field(v1, 0, 3), 9);
  EXPECT_EQ(field(v2, 0, 0), 1);
  EXPECT_EQ(field(v2, 0, 1), 2);
  EXPECT_EQ(field(v3, 0, 0), 4);
  EXPECT_EQ(field(v3, 0, 1), 5);
  EXPECT_EQ(field(v3, 0, 2), 6);

  field.insert(v2, 0, 10);

  EXPECT_EQ(field(v1, 0, 0), 3);
  EXPECT_EQ(field(v1, 0, 1), 7);
  EXPECT_EQ(field(v1, 0, 2), 8);
  EXPECT_EQ(field(v1, 0, 3), 9);
  EXPECT_EQ(field(v2, 0, 0), 1);
  EXPECT_EQ(field(v2, 0, 1), 2);
  EXPECT_EQ(field(v2, 0, 2), 10);
  EXPECT_EQ(field(v3, 0, 0), 4);
  EXPECT_EQ(field(v3, 0, 1), 5);
  EXPECT_EQ(field(v3, 0, 2), 6);
}

TEST(VariableSizeField, ValuesTwoNode2)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(2, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v1, 0, 1);
  field.insert(v2, 0, 2);
  field.insert(v1, 1, 3);
  field.insert(v2, 1, 4);
  field.insert(v1, 0, 5);

  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 5);
  EXPECT_EQ(field(v1, 1, 0), 3);
  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 1, 0), 4);
}

TEST(VariableSizeField, AddEntity)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v1, 0, 1);
  field.insert(v1, 0, 2);
  field.insert(v2, 0, 3);

  auto v3 = mesh->create_vertex(2, 0);
  EXPECT_EQ(field.get_num_comp(v3, 0), 0);

  field.insert(v3, 0, 4);

  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 2);
  EXPECT_EQ(field(v2, 0, 0), 3);
  EXPECT_EQ(field(v3, 0, 0), 4);
}

TEST(VariableSizeField, DeleteEntity)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1   = mesh->create_vertex(0, 0, 0);
  auto v2   = mesh->create_vertex(1, 0, 0);
  auto v3   = mesh->create_vertex(1, 1, 0);
  auto v4   = mesh->create_vertex(1, 0, 0);
  auto tri1 = mesh->create_triangle_from_verts(v1, v3, v4);
  mesh->create_triangle_from_verts(v2, v3, v4);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v2, 0, 2);
  field.insert(v1, 0, 1);
  field.insert(v3, 0, 3);
  field.insert(v4, 0, 4);

  mesh->delete_face(tri1);
  mesh->condense_arrays();

  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v3, 0, 0), 3);
  EXPECT_EQ(field(v4, 0, 0), 4);

  field.insert(v2, 0, 5);
  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 5);
  EXPECT_EQ(field(v3, 0, 0), 3);
  EXPECT_EQ(field(v4, 0, 0), 4);
}

TEST(VariableSizeField, MoveToEndRegression)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v1, 0, 1);
  field.insert(v2, 0, 2);
  field.insert(v1, 0, 3);
  field.insert(v2, 0, 4);

  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);
  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 4);
}

TEST(VariableSizeField, Clear)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1                            = mesh->create_vertex(0, 0, 0);
  auto v2                            = mesh->create_vertex(1, 0, 0);
  auto v3                            = mesh->create_vertex(2, 0, 0);
  auto v4                            = mesh->create_vertex(3, 0, 0);
  std::array<MeshEntityPtr, 4> verts = {v1, v2, v3, v4};

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(2, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v2, 0, 0);
  field.insert(v2, 0, 1);
  field.insert(v2, 1, 2);
  field.insert(v2, 1, 3);
  field.insert(v2, 1, 4);

  field.clear(0);
  for (auto& vert : verts)
  {
    EXPECT_EQ(field.get_num_comp(vert, 0), 0);
    EXPECT_EQ(field.get_num_comp(vert, 1), 0);
  }
}

TEST(VariableSizeField, ClearThenAddValues)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  auto& field   = *fieldPtr;

  field.insert(v1, 0, 1);
  field.insert(v2, 0, 2);
  field.insert(v1, 0, 3);
  field.insert(v2, 0, 4);

  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);
  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 4);

  field.clear(0);

  field.insert(v1, 0, 1);
  field.insert(v2, 0, 2);
  field.insert(v1, 0, 3);
  field.insert(v2, 0, 4);

  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);
  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 4);
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
