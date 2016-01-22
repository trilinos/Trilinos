
#ifndef _EQUIVALENT_ENTITY_BLOCKS_H_
#define _EQUIVALENT_ENTITY_BLOCKS_H_

#include <vector>

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace impl {

void get_element_blocks(const stk::mesh::MetaData &metaData, std::vector<stk::mesh::Part*> &element_blocks);
bool are_entity_element_blocks_equivalent(const stk::mesh::BulkData& bulkData, stk::mesh::Entity entity1, stk::mesh::Entity entity2);
bool is_element_block(const stk::mesh::Part &part);

}
}
}

#endif
