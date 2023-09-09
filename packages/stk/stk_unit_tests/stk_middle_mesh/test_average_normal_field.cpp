#include "gtest/gtest.h"
#include "stk_middle_mesh/predicates/average_normal_field.hpp"
#include "stk_middle_mesh/create_mesh.hpp"

using namespace stk::middle_mesh;

TEST(AveragedNormalField, Plane)
{
  mesh::impl::MeshSpec spec;
  spec.xmin = 0; spec.ymin = 0;
  spec.xmax = 1; spec.ymax = 1;
  spec.numelX = 4, spec.numelY = 4;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = mesh::impl::create_mesh(spec, f);

  predicates::impl::AveragedNormalField averagedNormalField(mesh);
  auto& normalField = *(averagedNormalField.get_field());

  double edgeLength = (spec.xmax - spec.xmin)/spec.numelX;
  for (auto& vert : mesh->get_vertices())
  {
    if (vert)
    {
      utils::Point normal = normalField(vert, 0, 0);
      
      EXPECT_NEAR(normal.x, 0.0, 1e-13);
      EXPECT_NEAR(normal.y, 0.0, 1e-13);
      EXPECT_NEAR(normal.z, edgeLength, 1e-13);
    }
  }

}