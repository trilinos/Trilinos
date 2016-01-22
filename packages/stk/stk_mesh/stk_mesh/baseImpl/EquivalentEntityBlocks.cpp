#include "EquivalentEntityBlocks.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool is_element_block(const stk::mesh::Part &part)
{
    return (part.primary_entity_rank() == stk::topology::ELEMENT_RANK && part.id() > 0 );
}

void get_element_blocks_from_parts(stk::mesh::PartVector &element_blocks, const stk::mesh::PartVector &partVector)
{
    for(stk::mesh::Part *part : partVector) {
        if(is_element_block(*part))
            element_blocks.push_back(part);
    }
}

void get_element_blocks_for_entity(const stk::mesh::BulkData &bulkData, stk::mesh::Entity entity, stk::mesh::PartVector &element_blocks)
{
    const stk::mesh::PartVector &partVector = bulkData.bucket(entity).supersets();
    get_element_blocks_from_parts(element_blocks, partVector);
}

bool are_blocks_equivalent(stk::mesh::PartVector &blocksForEntity1, stk::mesh::PartVector &blocksForEntity2)
{
    stk::util::sort_and_unique(blocksForEntity1);
    stk::util::sort_and_unique(blocksForEntity2);
    return (blocksForEntity1 == blocksForEntity2);
}

void get_element_blocks(const stk::mesh::MetaData &metaData, stk::mesh::PartVector &element_blocks)
{
    const stk::mesh::PartVector partVector = metaData.get_mesh_parts();
    get_element_blocks_from_parts(element_blocks, partVector);
}

bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, stk::mesh::Entity entity2)
{
    stk::mesh::PartVector blocksForEntity1, blocksForEntity2;
    get_element_blocks_for_entity(bulkData, entity1, blocksForEntity1);
    get_element_blocks_for_entity(bulkData, entity2, blocksForEntity2);
    return are_blocks_equivalent(blocksForEntity1, blocksForEntity2);
}

}
}
}
