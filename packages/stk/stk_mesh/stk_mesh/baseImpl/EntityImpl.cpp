/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//#include <stddef.h>
//#include <stdexcept>
//#include <iostream>
//#include <sstream>
#include <algorithm>

#include <stk_mesh/base/Ghosting.hpp>

//#include <stk_mesh/base/BulkData.hpp>
//#include <stk_mesh/base/MetaData.hpp>

#include <stk_mesh/baseImpl/EntityImpl.hpp>

namespace stk {
namespace mesh {
namespace impl {

EntityImpl::~EntityImpl()
{
}


EntityImpl::EntityImpl( const EntityKey & arg_key )
  : m_key(arg_key),
    m_relation(),
    m_comm(),
    m_bucket( NULL ),
    m_bucket_ord(0),
    m_owner_rank(0),
    m_sync_count(0),
    m_mod_log( LogNoChange )
{
}


PairIterRelation EntityImpl::relations( unsigned rank ) const
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


PairIterEntityComm EntityImpl::sharing() const
{
  EntityCommInfoVector::const_iterator i = m_comm.begin();
  EntityCommInfoVector::const_iterator e = m_comm.end();

  e = std::lower_bound( i , e , EntityCommInfo(1,0) );

  return PairIterEntityComm( i , e );
}


PairIterEntityComm EntityImpl::comm( const Ghosting & sub ) const
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


bool EntityImpl::insert( const EntityCommInfo & val )
{
  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( m_comm.begin() , m_comm.end() , val );

  const bool result = i == m_comm.end() || val != *i ;

  if ( result ) { m_comm.insert( i , val ); }

  return result ;
}


bool EntityImpl::erase( const EntityCommInfo & val )
{
  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( m_comm.begin() , m_comm.end() , val );

  const bool result = i != m_comm.end() && val == *i ;

  if ( result ) { m_comm.erase( i ); }

  return result ;
}


bool EntityImpl::erase( const Ghosting & ghost )
{
  typedef std::vector< EntityCommInfo > EntityComm ;

  const EntityCommInfo s_begin( ghost.ordinal() ,     0 );
  const EntityCommInfo s_end(   ghost.ordinal() + 1 , 0 );

  EntityComm::iterator i = m_comm.begin();
  EntityComm::iterator e = m_comm.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  const bool result = i != e ;

  if ( result ) { m_comm.erase( i , e ); }

  return result ;
}


void EntityImpl::comm_clear_ghosting()
{
  std::vector< EntityCommInfo >::iterator j = m_comm.begin();
  while ( j != m_comm.end() && j->ghost_id == 0 ) { ++j ; }
  m_comm.erase( j , m_comm.end() );
}


void EntityImpl::comm_clear()
{ m_comm.clear(); }

void EntityImpl::log_clear()
{ m_mod_log = LogNoChange ; }


void EntityImpl::log_created()
{
  if ( LogNoChange == m_mod_log ) {
    m_mod_log = LogCreated ;
  }
}


void EntityImpl::log_modified()
{
  if ( LogNoChange == m_mod_log ) {
    m_mod_log = LogModified ;
  }
}


} // namespace impl 
} // namespace mesh
} // namespace stk

