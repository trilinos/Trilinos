#include "EquivalentEntityBlocks.hpp"
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, stk::mesh::Entity entity2)
{
    std::vector<stk::mesh::PartOrdinal> blockOrdinalsForEntity1 = get_element_block_part_ordinals(entity1, bulkData);
    std::vector<stk::mesh::PartOrdinal> blockOrdinalsForEntity2 = get_element_block_part_ordinals(entity2, bulkData);
    return blockOrdinalsForEntity1==blockOrdinalsForEntity2;
}

bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, const std::vector<stk::mesh::PartOrdinal>& other_element_part_ordinals)
{
    std::vector<stk::mesh::PartOrdinal> blockOrdinalsForEntity1 = get_element_block_part_ordinals(entity1, bulkData);
    return blockOrdinalsForEntity1 == other_element_part_ordinals;
}

std::vector<PartOrdinal> get_element_block_part_ordinals(stk::mesh::Entity element, const stk::mesh::BulkData& bulkData)
{
  std::vector<PartOrdinal> partOrdinals;
    bulkData.bucket(element).supersets(partOrdinals);
    std::vector<PartOrdinal> partOrdinalsElementBlock;
    for(PartOrdinal part_ordinal : partOrdinals)
        if (stk::mesh::is_element_block(bulkData.mesh_meta_data().get_part(part_ordinal)))
            partOrdinalsElementBlock.push_back(part_ordinal);
    return partOrdinalsElementBlock;
}


}
}
}
