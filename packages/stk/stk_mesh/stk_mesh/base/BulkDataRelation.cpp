/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Trace.hpp>

namespace stk {
namespace mesh {

void set_field_relations( Entity e_from ,
                          Entity e_to ,
                          const unsigned ident )
{
  const std::vector<FieldRelation> & field_rels =
    MetaData::get(e_from).get_field_relations();

  for ( std::vector<FieldRelation>::const_iterator
        j = field_rels.begin() ; j != field_rels.end() ; ++j ) {

    const FieldRelation & fr = *j ;

    void ** const ptr = (void**) field_data( * fr.m_root , e_from );

    if ( ptr ) {

      void * const src = field_data( * fr.m_target , e_to );

      const size_t number =
        field_data_size(*fr.m_root,e_from) / sizeof(void*);

      const size_t offset =
         (*fr.m_function)( e_from.entity_rank() ,
                           e_to.entity_rank() , ident );

      if ( offset < number ) {
        ptr[ offset ] = src ;
      }
    }
  }
}

namespace {

void clear_field_relations( Entity e_from ,
                            const unsigned type ,
                            const unsigned ident )
{
  const std::vector<FieldRelation> & field_rels =
    MetaData::get(e_from).get_field_relations();

  for ( std::vector<FieldRelation>::const_iterator
        j = field_rels.begin() ; j != field_rels.end() ; ++j ) {

    const FieldRelation & fr = *j ;

    void ** const ptr = (void**) field_data( * fr.m_root , e_from );

    if ( ptr ) {

      const size_t number =
        field_data_size(*fr.m_root,e_from) / sizeof(void*);

      const size_t offset =
        (*fr.m_function)( e_from.entity_rank() , type , ident );

      if ( offset < number ) {
        ptr[ offset ] = NULL ;
      }
    }
  }
}

} // empty namespace

//----------------------------------------------------------------------

void BulkData::require_valid_relation( const char action[] ,
                                       const BulkData & mesh ,
                                       const Entity e_from ,
                                       const Entity e_to )
{
  const bool error_mesh_from = & mesh != & BulkData::get(e_from);
  const bool error_mesh_to   = & mesh != & BulkData::get(e_to);
  const bool error_type      = e_from.entity_rank() <= e_to.entity_rank();
  const bool error_nil_from  = !e_from.is_valid();
  const bool error_nil_to    = !e_to.is_valid();

  if ( error_mesh_from || error_mesh_to || error_type ||
       error_nil_from || error_nil_to ) {
    std::ostringstream msg ;

    msg << "Could not " << action << " relation from entity "
        << print_entity_key(e_from) << " to entity "
        << print_entity_key(e_to) << "\n";

    ThrowErrorMsgIf( error_mesh_from || error_mesh_to,
                     msg.str() << (error_mesh_from ? "e_from" : "e_to" ) <<
                     " not member of this mesh");
    ThrowErrorMsgIf( error_nil_from  || error_nil_to,
                     msg.str() << (error_mesh_from ? "e_from" : "e_to" ) <<
                     " was destroyed");
    ThrowErrorMsgIf( error_type, msg.str() <<
                     "A relation must be from higher to lower ranking entity");
  }
}

//----------------------------------------------------------------------

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id )
{
  TraceIfWatching("stk::mesh::BulkData::declare_relation", LOG_ENTITY, e_from.key());
  TraceIfWatchingDec("stk::mesh::BulkData::declare_relation", LOG_ENTITY, e_to.key(), 1);
  DiagIfWatching(LOG_ENTITY, e_from.key(),
                 "from: " << e_from << ";  " <<
                 "to: " << e_to << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, e_to.key(),
                 "from: " << e_from << ";  " <<
                 "to: " << e_to << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "declare" , *this , e_from , e_to );

  // TODO: Don't throw if exact relation already exists, that should be a no-op.
  // Should be an exact match if relation of local_id already exists (e_to should be the same).
  m_entity_repo.declare_relation( e_from, e_to, local_id, m_sync_count);

  OrdinalVector add , empty ;

  // Deduce and set new part memberships:

  induced_part_membership( e_from, empty, e_to.entity_rank(), local_id, add );

  internal_change_entity_parts( e_to , add , empty );

  set_field_relations( e_from , e_to , local_id );
}

//----------------------------------------------------------------------

void BulkData::declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel )
{
  require_ok_to_modify();

  const unsigned erank = entity.entity_rank();

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity e = i->entity();
    const unsigned n = i->identifier();
    if ( e.entity_rank() < erank ) {
      declare_relation( entity , e , n );
    }
    else if ( erank < e.entity_rank() ) {
      declare_relation( e , entity , n );
    }
    else {
      ThrowErrorMsg("Given entities of the same entity rank. entity is " <<
                    print_entity_key(entity));
    }
  }
}

//----------------------------------------------------------------------

bool BulkData::destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
  TraceIfWatching("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, e_from.key());
  TraceIfWatchingDec("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, e_to.key(), 1);
  DiagIfWatching(LOG_ENTITY, e_from.key(),
                 "from: " << e_from << ";  " <<
                 "to: " << e_to << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, e_to.key(),
                 "from: " << e_from << ";  " <<
                 "to: " << e_to << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "destroy" , *this , e_from , e_to );

  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  if ( parallel_size() < 2 || entity_comm_sharing(e_to.key()).empty() ) {

    //------------------------------
    // 'keep' contains the parts deduced from kept relations
    // 'del'  contains the parts deduced from deleted relations
    //        that are not in 'keep'
    // Only remove these part memberships the entity is not shared.
    // If the entity is shared then wait until modificaton_end_synchronize.
    //------------------------------

    OrdinalVector del, keep, empty;

    // For all relations that are *not* being deleted, add induced parts for
    // these relations to the 'keep' vector
    for ( PairIterRelation i = e_to.relations(); !i.empty(); ++i ) {
      if (e_to.entity_rank() < i->entity_rank()) { // Need to look at back rels only
        if ( !( i->entity() == e_from && i->identifier() == local_id ) ) {
          induced_part_membership( i->entity(), empty, e_to.entity_rank(),
                                   i->identifier(), keep,
                                   false /*Do not look at supersets*/);
        }
      }
    }

    // Find the relation this is being deleted and add the parts that are
    // induced from that relation (and that are not in 'keep') to 'del'
    for ( PairIterRelation i = e_from.relations() ; !i.empty() ; ++i ) {
      if ( i->entity() == e_to && i->identifier() == local_id ) {
        induced_part_membership( e_from, keep, e_to.entity_rank(),
                                 i->identifier(), del,
                                 false /*Do not look at supersets*/);
        clear_field_relations( e_from , e_to.entity_rank() ,
                               i->identifier() );
        break; // at most 1 relation can match our specification
      }
    }

    if ( !del.empty() ) {
      internal_change_entity_parts( e_to , empty , del );
    }
  }
  else {
    // Just clear the field, part membership will be handled by modification end
    for ( PairIterRelation i = e_from.relations() ; !i.empty() ; ++i ) {
      if ( i->entity() == e_to && i->identifier() == local_id ) {
        clear_field_relations( e_from , e_to.entity_rank() ,
                               i->identifier() );
        break; // at most 1 relation can match our specification
      }
    }
  }

  //delete relations from the entities
  return m_entity_repo.destroy_relation( e_from, e_to, local_id);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

