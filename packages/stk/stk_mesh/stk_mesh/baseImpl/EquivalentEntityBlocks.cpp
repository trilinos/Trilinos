#include "EquivalentEntityBlocks.hpp"
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, stk::mesh::Entity entity2, std::vector<stk::mesh::PartOrdinal>& ordinalScratchSpace1, std::vector<stk::mesh::PartOrdinal>& ordinalScratchSpace2)
{
    get_element_block_part_ordinals(entity1, bulkData, ordinalScratchSpace1);
    get_element_block_part_ordinals(entity2, bulkData, ordinalScratchSpace2);
    return ordinalScratchSpace1==ordinalScratchSpace2;
}

bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, const std::vector<stk::mesh::PartOrdinal>& other_element_part_ordinals, std::vector<stk::mesh::PartOrdinal>& scratch_space)
{
    get_element_block_part_ordinals(entity1, bulkData, scratch_space);
    return scratch_space == other_element_part_ordinals;
}

void get_element_block_part_ordinals(stk::mesh::Entity element, const stk::mesh::BulkData& bulkData, std::vector<PartOrdinal>& partOrdinalsElementBlock)
{
  partOrdinalsElementBlock.clear();
  const stk::mesh::PartVector& parts = bulkData.bucket(element).supersets();
  for(const stk::mesh::Part* part : parts) {
    if (stk::mesh::is_element_block(*part)) {
      partOrdinalsElementBlock.push_back(part->mesh_meta_data_ordinal());
    }
  }
}

}
}
}
