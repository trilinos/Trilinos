/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>
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
  : m_key( arg_key ),
    m_relation(),
    m_comm(),
    m_bucket( NULL ),
    m_bucket_ord(0),
    m_owner_rank(0),
    m_sync_count(0)
{}

Entity::~Entity()
{}

//----------------------------------------------------------------------

PairIterRelation Entity::relations( unsigned rank ) const
{
  std::vector<Relation>::const_iterator i = m_relation.begin();
  std::vector<Relation>::const_iterator e = m_relation.end();

  if ( rank ) {
    const Relation::raw_attr_type lo_attr = Relation::attribute( rank , 0 );
    i = std::lower_bound( i , e , lo_attr , LessRelation() );
  }

  const Relation::raw_attr_type hi_attr = Relation::attribute( rank + 1 , 0 );
  e = std::lower_bound( i , e , hi_attr , LessRelation() );

  return PairIterRelation( i , e );
}

//----------------------------------------------------------------------

PairIterEntityComm Entity::sharing() const
{
  typedef std::vector< EntityCommInfo > EntityComm ;

  EntityComm::const_iterator i = m_comm.begin();
  EntityComm::const_iterator e = m_comm.end();

  e = std::lower_bound( i , e , EntityCommInfo(1,0) );

  return PairIterEntityComm( i , e );
}

PairIterEntityComm Entity::comm( const Ghosting & sub ) const
{
  typedef std::vector< EntityCommInfo > EntityComm ;

  const EntityCommInfo s_begin( sub.ordinal() ,     0 );
  const EntityCommInfo s_end(   sub.ordinal() + 1 , 0 );

  EntityComm::const_iterator i = m_comm.begin();
  EntityComm::const_iterator e = m_comm.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  return PairIterEntityComm( i , e );
}

bool Entity::insert( const EntityCommInfo & val )
{
  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( m_comm.begin() , m_comm.end() , val );

  const bool result = i == m_comm.end() || val != *i ;

  if ( result ) { m_comm.insert( i , val ); }

  return result ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

