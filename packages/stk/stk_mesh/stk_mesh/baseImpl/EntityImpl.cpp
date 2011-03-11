/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <stk_mesh/base/Ghosting.hpp>

#include <stk_mesh/baseImpl/EntityImpl.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

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
    m_mod_log( EntityLogCreated )
{
}

PairIterRelation EntityImpl::relations( unsigned rank ) const
{
  std::vector<Relation>::const_iterator i = m_relation.begin();
  std::vector<Relation>::const_iterator e = m_relation.end();

  //Nodes
  if ( rank != 0 ) {
    const Relation::raw_relation_id_type lo_attr = Relation::raw_relation_id( rank , 0 );
    i = std::lower_bound( i , e , lo_attr , LessRelation() );
  }

  const Relation::raw_relation_id_type hi_attr = Relation::raw_relation_id( rank + 1 , 0 );
  e = std::lower_bound( i , e , hi_attr , LessRelation() );

  return PairIterRelation( i , e );
}

PairIterEntityComm EntityImpl::sharing() const
{
  EntityCommInfoVector::const_iterator i = m_comm.begin();
  EntityCommInfoVector::const_iterator e = m_comm.end();

  e = std::lower_bound( i , e , EntityCommInfo(1,     // ghost id, 1->aura
                                               0 ) ); // proc

  // Contains everything up the first aura comm (IE, only contains shared comms)
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

  if ( result ) {
    m_comm.insert( i , val );
  }

  return result ;
}

bool EntityImpl::erase( const EntityCommInfo & val )
{
  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( m_comm.begin() , m_comm.end() , val );

  const bool result = i != m_comm.end() && val == *i ;

  if ( result ) {
    m_comm.erase( i );
  }

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

  if ( result ) {
    m_comm.erase( i , e );
  }

  return result ;
}

void EntityImpl::comm_clear_ghosting()
{
  std::vector< EntityCommInfo >::iterator j = m_comm.begin();
  while ( j != m_comm.end() && j->ghost_id == 0 ) { ++j ; }
  m_comm.erase( j , m_comm.end() );
}


void EntityImpl::comm_clear()
{
  m_comm.clear();
}

bool EntityImpl::marked_for_destruction() const
{
  // The original implementation of this method checked bucket capacity. In
  // order to ensure that the addition of EntityLogDeleted does not change
  // behavior, we put error check here.
  ThrowErrorMsgIf((bucket().capacity() == 0) != (m_mod_log == EntityLogDeleted),
      "Inconsistent destruction state; " <<
      "destroyed entities should be in the nil bucket and vice versa.\n" <<
      "Problem is with entity: " <<
      print_entity_key( MetaData::get( bucket() ), key() ) <<
      "\nWas in nil bucket: " << (bucket().capacity() == 0) << ", " <<
      "was in destroyed state: " << (m_mod_log == EntityLogDeleted) );

  return m_mod_log == EntityLogDeleted;
}

namespace {

inline bool is_degenerate_relation ( const Relation &r1 , const Relation &r2 )
{
  return r1.raw_relation_id() == r2.raw_relation_id() && r1.entity() != r2.entity() ;
}

}

void EntityImpl::log_resurrect()
{
  ThrowErrorMsgIf( EntityLogDeleted != m_mod_log,
      "Trying to resurrect non-deleted entity: " <<
      print_entity_key( MetaData::get( bucket() ), key() ) );

  m_mod_log = EntityLogModified;
  m_bucket = NULL;
}

void EntityImpl::log_modified_and_propagate()
{
  // If already in modified state, return
  if (m_mod_log != EntityLogNoChange) {
    return;
  }

  // mark this entity as modified
  m_mod_log = EntityLogModified;

  // recurse on related entities w/ higher rank
  EntityRank rank_of_original_entity = entity_rank();
  for ( PairIterRelation irel = relations() ; irel.first != irel.second ; ) {
    --irel.second;
    Entity & entity = *(irel.second->entity());
    if ( rank_of_original_entity >= entity.entity_rank() ) {
      break; //we're done
    }
    else if ( entity.log_query() == EntityLogNoChange ) {
      entity.m_entityImpl.log_modified_and_propagate();
    }
  }

}

void EntityImpl::log_created_parallel_copy()
{
  if ( EntityLogCreated == m_mod_log ) {
    m_mod_log = EntityLogModified ;
  }
}

bool EntityImpl::destroy_relation( Entity& e_to, const RelationIdentifier local_id )
{
  bool destroyed_relations = false;
  for ( std::vector<Relation>::iterator
        i = m_relation.begin() ; i != m_relation.end() ; ++i ) {
    if ( i->entity() == & e_to && i->identifier() == local_id ) {
      i = m_relation.erase( i ); // invalidates iterators, but we're breaking so it's OK
      destroyed_relations = true;
      break;
    }
  }
  return destroyed_relations;
}

bool EntityImpl::declare_relation( Entity & e_to,
                                   const RelationIdentifier local_id,
                                   unsigned sync_count,
                                   bool is_back_relation )
{
  const MetaData & meta_data = MetaData::get( bucket() );

  const Relation new_relation( e_to , local_id );

  const std::vector<Relation>::iterator rel_end   = m_relation.end();
        std::vector<Relation>::iterator rel_begin = m_relation.begin();
        std::vector<Relation>::iterator lower;

  lower = std::lower_bound( rel_begin , rel_end , new_relation , LessRelation() );

  // The ordering of the Relations allows for two situations that do
  // not arise often in meshes.  The first situation is 2 relations between
  // e_from and e_to with the same kind but different local_ids.  This
  // can happen if, for example, a triangle should be used as a quad.  In
  // this case, one node of the triangle must be two different local nodes of
  // the quad.  This situation is a valid state of mesh entities.

  // The second situation involves malformed stencils.  Given e_from, e_to1,
  // and e_to2, e_to1 and eto2 can share a relation with e_from with the same
  // kind and local_id.  This can arise, for instance, if an edge has three
  // nodes.  The local_id 1 of the edge may point to two different nodes.
  // This situation is disallowed in the mesh.  We now check for it.

  // "Degenerate" -> case where we have two relations whose attributes
  // (rel id + rel rank) match but point to different entities. It's
  // OK for back-relations to be degenerate because there's nothing
  // wrong with a node having several back-relations (with similar id)
  // to different elements.

  // Check for bad degenerate relations (degenerate forward relations)
  // Cannot be degenerate relation if there are no prior relations
  if ( !m_relation.empty() && !is_back_relation ) {
    // Since LessRelation takes the related entity into account, we must check
    // the result of lower_bound AND the iter before to be sure this isn't a
    // bad degenerate relation.
    std::vector<Relation>::iterator start, end;
    start = (lower == rel_begin) ? rel_begin : lower - 1;
    end   = (lower == rel_end)   ? rel_end   : lower + 1;

    for (std::vector<Relation>::iterator itr = start; itr != end; ++itr) {
      ThrowErrorMsgIf( is_degenerate_relation ( new_relation , *itr ),
                       "Could not declare relation from " <<
                       print_entity_key( meta_data, key() ) << " to " <<
                       print_entity_key( meta_data, e_to.key() ) << ", with id " <<
                       local_id << ". Relation already exists to " <<
                       print_entity_key( meta_data, itr->entity()->key() ));
    }
  }

  bool not_already_exists = (rel_end == lower) ||
    ( !is_back_relation && new_relation.raw_relation_id() != lower->raw_relation_id() ) ||
    ( is_back_relation && new_relation != *lower );

  // If the relation does not already exist, we add it
  if (not_already_exists) {
    lower = m_relation.insert( lower , new_relation );

    set_sync_count( sync_count );

    return true;
  }
  else {
    return false;
  }
}

} // namespace impl
} // namespace mesh
} // namespace stk

