#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------


std::ostream &
print_entity_key( std::ostream & os , const MetaData & meta_data ,
                  unsigned type , EntityId id )
{
  const std::string & name = meta_data.entity_type_name( type );
  return os << name << "[" << id << "]" ;
}

std::ostream &
print_entity_key( std::ostream & os , const MetaData & meta_data ,
                  const EntityKey & key )
{
  const unsigned type   = entity_type(key);
  const EntityId id = entity_id(key);
  return print_entity_key( os , meta_data , type , id );
}

//----------------------------------------------------------------------

Entity::Entity( const EntityKey & arg_key )
  : m_key( arg_key ), m_relation(), m_bucket(), m_bucket_ord(0),
    m_owner_rank(0), m_sync_count(0),
    m_sharing(std::vector<EntityProc>().begin(),std::vector<EntityProc>().end())
{}

Entity::~Entity()
{}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

