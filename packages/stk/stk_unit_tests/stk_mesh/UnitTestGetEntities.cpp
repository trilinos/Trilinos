#include "gtest/gtest.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

class GetEntitiesTest : public stk::unit_test_util::MeshFixture { };

TEST_F(GetEntitiesTest, get_num_entities)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);
    EXPECT_EQ(9u, stk::mesh::get_num_entities(get_bulk()));
  }
}

}
