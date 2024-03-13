#include "gtest/gtest.h"
#include "stk_middle_mesh/field_shared_reduction.hpp"
#include "stk_middle_mesh/create_mesh.hpp"

using namespace stk::middle_mesh;

TEST(FieldSharedReduction, Sum)
{
  mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };
  std::shared_ptr<mesh::Mesh> mesh = create_mesh(spec, func);

  auto fieldPtr = mesh::create_field<int>(mesh, mesh::FieldShape(2, 2, 0), 3);
  auto& field = *fieldPtr;

  for (int dim=0; dim < 2; ++dim)
    for (auto entity : mesh->get_mesh_entities(dim))
    {
      int idx = 0;
      for (int i=0; i < 2; ++i)
        for (int j=0; j < 3; ++j)
          field(entity, i, j) = entity->get_id() + idx++;
    }

  mesh::ReductionOpSum<int> op;
  mesh::FieldSharedReduction<int> reducer(fieldPtr, op);
  reducer.reduce();

  for (int dim=0; dim < 2; ++dim)
  {
    for (auto entity : mesh->get_mesh_entities(dim))
    {
      int baseVal = entity->get_id();
      for (int i=0; i < entity->count_remote_shared_entities(); ++i)
      {
        baseVal += entity->get_remote_shared_entity(i).remoteId;
      }

      int idx = 0;
      for (int i=0; i < 2; ++i)
        for (int j=0; j < 3; ++j)
        {
          EXPECT_EQ(field(entity, i, j), baseVal + (entity->count_remote_shared_entities() + 1)*idx);
          ++idx;
        }
    }
  }


}