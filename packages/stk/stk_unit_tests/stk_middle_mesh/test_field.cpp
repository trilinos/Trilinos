#include "gtest/gtest.h"

#include <iostream>

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

int get_value(const Field<int>& field, MeshEntityPtr e, const int node, const int comp)
{
  int dim = get_type_dimension(e->get_type());
  int id  = e->get_id();
  return id * (field.get_num_nodes(dim) * field.get_num_comp()) + node * field.get_num_comp() + comp;
}

void test_field(std::shared_ptr<Mesh> mesh, const FieldShape& fshape, Field<int>& field)
{
  // test getting and setting values
  for (int dim = 0; dim < 3; ++dim)
    for (auto& e : mesh->get_mesh_entities(dim))
      for (int node = 0; node < fshape.count[0]; ++node)
        for (int comp = 0; comp < field.get_num_comp(); ++comp)
          field(e, node, comp) = get_value(field, e, node, comp);

  for (int dim = 0; dim < 3; ++dim)
    for (auto& e : mesh->get_mesh_entities(dim))
      for (int node = 0; node < fshape.count[0]; ++node)
        for (int comp = 0; comp < field.get_num_comp(); ++comp)
          EXPECT_EQ(field(e, node, comp), get_value(field, e, node, comp));
}

TEST(FieldShape, operators)
{
  EXPECT_TRUE(FieldShape(1, 1, 1) == FieldShape(1, 1, 1));
  EXPECT_FALSE(FieldShape(1, 1, 0) ==  FieldShape(1, 1, 1));
}

TEST(Field, Values)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  FieldShape fshape(2, 3, 4);
  auto fieldPtr     = create_field<int>(mesh, fshape, 5, -1);
  Field<int>& field = *fieldPtr;

  // test initial value
  for (int dim = 0; dim < 3; ++dim)
    for (auto& e : mesh->get_mesh_entities(dim))
      for (int node = 0; node < fshape.count[0]; ++node)
        for (int comp = 0; comp < field.get_num_comp(); ++comp)
          EXPECT_EQ(field(e, node, comp), -1);

  test_field(mesh, fshape, field);
}

TEST(Field, Grow)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  int nelemX = 4, nelemY = 0;
  double xmin = 0, xmax = 1, ymin = 0, ymax = 1;
  double dx = (xmax - xmin) / nelemX;
  double dy = (ymax - ymin) / nelemY;
  std::vector<MeshEntityPtr> verts((nelemX + 1) * (nelemY + 1));

  auto getIdx = [nelemY](const int i, const int j) { return j + i * (nelemY + 1); };

  auto mesh = make_empty_mesh();
  FieldShape fshape(2, 3, 4);
  auto fieldPtr     = create_field<int>(mesh, fshape, 5, -1);
  Field<int>& field = *fieldPtr;

  for (int i = 0; i < (nelemX + 1); ++i)
    for (int j = 0; j < (nelemY + 1); ++j)
      verts[getIdx(i, j)] = mesh->create_vertex(xmin + i * dx, ymin + j * dy);

  for (int i = 0; i < nelemX; ++i)
    for (int j = 0; j < nelemY; ++j)
    {
      auto v1 = verts[getIdx(i, j)];
      auto v2 = verts[getIdx(i + 1, j)];
      auto v3 = verts[getIdx(i + 1, j + 1)];
      auto v4 = verts[getIdx(i, j + 1)];
      mesh->create_quad_from_verts(v1, v2, v3, v4);
    }

  test_field(mesh, fshape, field);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
