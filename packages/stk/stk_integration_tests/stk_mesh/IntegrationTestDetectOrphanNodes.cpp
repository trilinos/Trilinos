#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include "stk_unit_test_utils/getOption.h"
#include "stk_mesh/base/GetEntities.hpp"

namespace
{

class OrphanedNodeMesh : public stk::unit_test_util::MeshFixture {};

bool doOrphanededNodesExist(const stk::mesh::BulkData& bulk)
{
  std::ostringstream os;
  os << "For processor " << bulk.parallel_rank() << std::endl;
  stk::mesh::EntityVector nodes;
  stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::NODE_RANK), nodes);
  bool isOrphanedMesh = false;
  for(stk::mesh::Entity node : nodes)
  {
    if(bulk.num_elements(node) == 0 )
    {
      isOrphanedMesh = true;
      os << "\tNode " << bulk.identifier(node) << " is an orphan\n";
    }
  }

  if(isOrphanedMesh)
    std::cerr << os.str();

  return isOrphanedMesh;
}

TEST_F(OrphanedNodeMesh, detectOrphanedNodes)
{
  std::string filename = stk::unit_test_util::get_option("-i", "generated:1x1x100");
  setup_mesh(filename, stk::mesh::BulkData::NO_AUTO_AURA);

  EXPECT_TRUE(!doOrphanededNodesExist(get_bulk()));
}

}
