#include "stk_mesh/base/GetEntities.hpp"  // for get_selected_entities
#include "stk_mesh/base/SideSetUtil.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_io/StkIoUtils.hpp>

namespace
{

class StkToolsB : public stk::unit_test_util::MeshFixture {};

TEST_F(StkToolsB, GetBlockIdsForSpecifiedSideset)
{
  const std::string unNamed = "mesh not specified";
  const std::string meshName = stk::unit_test_util::get_option("-i", unNamed);
  STK_ThrowRequireMsg(meshName!=unNamed, "Please specify mesh with -i option.");
  setup_mesh(meshName, stk::mesh::BulkData::NO_AUTO_AURA);

  int invalidSideset = -1;
  int sideset = stk::unit_test_util::get_command_line_option("-s", invalidSideset);
  STK_ThrowRequireMsg(sideset!=invalidSideset, "Please specify sideset with -s.");

  std::string sidesetName = "surface_" + std::to_string(sideset);
  stk::mesh::Part* ss = get_meta().get_part(sidesetName);
  STK_ThrowRequire(ss!=nullptr);
  stk::mesh::EntityVector sides;
  stk::mesh::get_selected_entities(*ss, get_bulk().buckets(get_meta().side_rank()),sides);

  std::ofstream outfile("test.out.txt", std::ofstream::app);

  std::set<int> blockIds;
  for (stk::mesh::Entity side : sides)
  {
    unsigned numElements = get_bulk().num_elements(side);
    if (numElements > 1)
    {
      outfile << "This is an internal side" << std::endl;
    }
    const stk::mesh::Entity *elements = get_bulk().begin_elements(side);
    for (unsigned i = 0; i < numElements; ++i)
    {
      const stk::mesh::Part *block = stk::mesh::get_element_block_part(get_bulk(), elements[i]);
      blockIds.insert(block->id());
    }
  }
  for (int id : blockIds)
  {
    outfile << "blockId: " << id << std::endl;
  }

}
}
