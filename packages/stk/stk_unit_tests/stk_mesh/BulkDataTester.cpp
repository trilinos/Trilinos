#include "BulkDataTester.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc


namespace stk { namespace mesh { namespace unit_test {

bool BulkDataTester::is_entity_in_ghosting_comm_map(stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = this->entity_key(entity);
    bool is_entity_in_aura_comm_map = !this->entity_comm_map(entityKey, this->aura_ghosting()).empty();
    return is_entity_in_aura_comm_map;
}

} } } // namespace stk mesh unit_test
