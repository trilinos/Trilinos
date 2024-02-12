#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/variable_size_field.hpp"
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

TEST(VariableSizeField, InsertMoveToFront)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);
  auto v3 = mesh->create_vertex(2, 0, 0);
  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));

  fieldPtr->insert(v1, 0, 1);
  fieldPtr->insert(v2, 0, 2);
  fieldPtr->insert(v3, 0, 3);
  fieldPtr->insert(v1, 0, 4);
  fieldPtr->insert(v2, 0, 5);

  EXPECT_EQ(fieldPtr->get_num_comp(v1, 0), 2);
  EXPECT_EQ(fieldPtr->get_num_comp(v2, 0), 2);
  EXPECT_EQ(fieldPtr->get_num_comp(v3, 0), 1);

  auto& field = *fieldPtr;
  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 4);
  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 5);
  EXPECT_EQ(field(v3, 0, 0), 3);
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

TEST(VariableSizeField, Iteration)
{

  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));

  fieldPtr->insert(v1, 0, 1);
  fieldPtr->insert(v2, 0, 2);
  fieldPtr->insert(v1, 0, 3);
  fieldPtr->insert(v2, 0, 4);
  fieldPtr->insert(v1, 0, 5);
  fieldPtr->insert(v2, 0, 6);

  std::vector<int> v1Vals, v2Vals;
  auto& field = *fieldPtr;
  for (auto& v : field(v1, 0))
    v1Vals.push_back(v);

  for (auto& v : field(v2, 0))
    v2Vals.push_back(v);

  std::vector<int> v1ValsExpected = {1, 3, 5};
  std::vector<int> v2ValsExpected = {2, 4, 6};
  EXPECT_EQ(v1Vals, v1ValsExpected);
  EXPECT_EQ(v2Vals, v2ValsExpected);

  v1Vals.clear();
  v2Vals.clear();
  const auto& fieldConst = *fieldPtr;
  for (auto& v : fieldConst(v1, 0))
    v1Vals.push_back(v);

  for (auto& v : fieldConst(v2, 0))
    v2Vals.push_back(v);

  EXPECT_EQ(v1Vals, v1ValsExpected);
  EXPECT_EQ(v2Vals, v2ValsExpected);
}

TEST(VariableSizeField, ResizeLarger)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  fieldPtr->insert(v1, 0, 1);
  fieldPtr->insert(v2, 0, 2);
  fieldPtr->insert(v1, 0, 3);
  fieldPtr->insert(v2, 0, 4);
  fieldPtr->insert(v1, 0, 5);
  fieldPtr->insert(v2, 0, 6);

  fieldPtr->resize(v1, 0, 5, 666);
  auto& field = *fieldPtr;
  EXPECT_EQ(field.get_num_comp(v1, 0), 5);
  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);
  EXPECT_EQ(field(v1, 0, 2), 5);
  EXPECT_EQ(field(v1, 0, 3), 666);
  EXPECT_EQ(field(v1, 0, 4), 666);

  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 4);
  EXPECT_EQ(field(v2, 0, 2), 6);
}

TEST(VariableSizeField, ResizeSmaller)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0);

  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));
  fieldPtr->insert(v1, 0, 1);
  fieldPtr->insert(v1, 0, 3);
  fieldPtr->insert(v1, 0, 5);
  fieldPtr->insert(v1, 0, 7);


  fieldPtr->insert(v2, 0, 2);
  fieldPtr->insert(v2, 0, 4);
  fieldPtr->insert(v2, 0, 6);

  fieldPtr->resize(v1, 0, 2);
  auto& field = *fieldPtr;
  EXPECT_EQ(field.get_num_comp(v1, 0), 2);
  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);

  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 4);
  EXPECT_EQ(field(v2, 0, 2), 6);

  fieldPtr->insert(v1, 0, 9);
  fieldPtr->insert(v1, 0, 11);
  EXPECT_EQ(field.get_num_comp(v1, 0), 4);
  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);  
  EXPECT_EQ(field(v1, 0, 2), 9);
  EXPECT_EQ(field(v1, 0, 3), 11);

  fieldPtr->insert(v1, 0, 13);
  EXPECT_EQ(field.get_num_comp(v1, 0), 5);
  EXPECT_EQ(field(v1, 0, 0), 1);
  EXPECT_EQ(field(v1, 0, 1), 3);  
  EXPECT_EQ(field(v1, 0, 2), 9);
  EXPECT_EQ(field(v1, 0, 3), 11);  
  EXPECT_EQ(field(v1, 0, 4), 13);  

  EXPECT_EQ(field(v2, 0, 0), 2);
  EXPECT_EQ(field(v2, 0, 1), 4);
  EXPECT_EQ(field(v2, 0, 2), 6);  

}


TEST(VariableSizeField, OneThenTwoStorage)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0); 
  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));

  for (int i=0; i < 4; ++i)
  {
    fieldPtr->insert(v1, 0, 0);
  } 

  for (int i=0; i < 4; ++i)
  {
    fieldPtr->insert(v2, 0, 1);
  }   

  EXPECT_EQ(mesh::compute_storage_efficiency(fieldPtr), 1.0);
}

TEST(VariableSizeField, TwoThenOneStorage)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0); 
  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));

  for (int i=0; i < 100; ++i)
  {
    fieldPtr->insert(v2, 0, 0);
  } 

  for (int i=0; i < 100; ++i)
  {
    fieldPtr->insert(v1, 0, 1);
  }   

  EXPECT_EQ(mesh::compute_storage_efficiency(fieldPtr), 1.0);
}

TEST(VariableSizeField, AlternatingPatternStorage)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0); 
  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));

  for (int i=0; i < 100; ++i)
  {
    fieldPtr->insert(v1, 0, 0);
    fieldPtr->insert(v2, 0, 1);
  } 

  std::cout << "storage usage: " << fieldPtr->get_storage_usage() << std::endl;
  EXPECT_GE(mesh::compute_storage_efficiency(fieldPtr), 0.5);
}

TEST(VariableSizeField, ThreeElementCyclicPatternStorage)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);

  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(1, 0, 0); 
  auto v3 = mesh->create_vertex(2, 0, 0); 
  auto fieldPtr = create_variable_size_field<int>(mesh, FieldShape(1, 0, 0));

  for (int i=0; i < 100; ++i)
  {
    fieldPtr->insert(v1, 0, 0);
    fieldPtr->insert(v2, 0, 1);
    fieldPtr->insert(v3, 0, 2);
  } 

  std::cout << "storage usage: " << fieldPtr->get_storage_usage() << std::endl;
  EXPECT_GE(mesh::compute_storage_efficiency(fieldPtr), 0.5);
}


//TODO: test iteration after resize

} // namespace impl
} // namespace middle_mesh
} // namespace stk
