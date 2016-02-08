
#ifndef _EQUIVALENT_ENTITY_BLOCKS_H_
#define _EQUIVALENT_ENTITY_BLOCKS_H_

#include <vector>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace impl {

bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, stk::mesh::Entity entity2);
bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, const std::vector<stk::mesh::PartOrdinal>& other_element_part_ordinals);
bool is_element_block(const stk::mesh::Part &part);
std::vector<PartOrdinal> get_element_block_part_ordinals(stk::mesh::Entity element, const stk::mesh::BulkData& bulkData);

}
}
}

#endif
