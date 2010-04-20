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
#include <assert.h>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>

namespace stk {
namespace mesh {

namespace {

void assert_valid_relation( const char method[] ,
                            const BulkData & mesh ,
                            const Entity   & e_from ,
                            const Entity   & e_to )
{
  const bool error_mesh_from = & mesh != & e_from.bucket().mesh();
  const bool error_mesh_to   = & mesh != & e_to.bucket().mesh();
  const bool error_type      = e_from.entity_rank() <= e_to.entity_rank();
  const bool error_nil_from  = 0 == e_from.bucket().capacity();
  const bool error_nil_to    = 0 == e_to.bucket().capacity();

  if ( error_mesh_from || error_mesh_to || error_type ||
       error_nil_from || error_nil_to ) {
    std::ostringstream msg ;
    msg << method << "( from " ;
    print_entity_key( msg , mesh.mesh_meta_data(), e_from.key() );
    if ( error_mesh_from ) {
      msg << " NOT MEMBER OF THIS MESH" ;
    }
    if ( error_nil_from ) {
      msg << " WAS DESTROYED" ;
    }
    msg << " , to " ;
    print_entity_key( msg , mesh.mesh_meta_data(), e_to.key() );
    if ( error_mesh_to ) {
      msg << " NOT MEMBER OF THIS MESH" ;
    }
    if ( error_nil_to ) {
      msg << " WAS DESTROYED" ;
    }
    msg << " ) FAILED" ;
    if ( error_type ) {
      msg << " A relation must be from higher to lower ranking entity" ;
    }
    throw std::runtime_error( msg.str() );
  }
}

bool is_degenerate_relation ( const Relation &r1 , const Relation &r2 )
{
  return r1.attribute() == r2.attribute() && r1.entity() != r2.entity() ;
}

void set_field_relations( Entity & e_from ,
                          Entity & e_to ,
                          const unsigned ident )
{
  const std::vector<FieldRelation> & field_rels =
    e_from.bucket().mesh().mesh_meta_data().get_field_relations();

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

void clear_field_relations( Entity & e_from ,
                            const unsigned type ,
                            const unsigned ident )

{
  const std::vector<FieldRelation> & field_rels =
    e_from.bucket().mesh().mesh_meta_data().get_field_relations();

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

}

//----------------------------------------------------------------------

void BulkData::declare_relation( Entity & e_from ,
                                 Entity & e_to ,
                                 const unsigned local_id )
{
  static const char method[] = "stk::mesh::BulkData::declare_relation" ;

  assert_ok_to_modify( method );

  assert_valid_relation( method , *this , e_from , e_to );

  const Relation forward( e_to , local_id );

  const std::vector<Relation>::iterator fe = e_from.m_relation.end();
        std::vector<Relation>::iterator fi = e_from.m_relation.begin();

  fi = std::lower_bound( fi , fe , forward , LessRelation() );

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

  bool found_degenerate_relation = false;
  EntityKey  degenerate_key;
  if ( fi != fe )
     {
     bool  downstream = fi->entity_rank() < e_from.entity_rank();
     if ( is_degenerate_relation ( forward , *fi ) && downstream )
        {
        found_degenerate_relation = true;
        degenerate_key = fi->entity()->key();
        }
     }
  if ( fi != e_from.m_relation.begin() )
     {
     --fi;
     bool  downstream = fi->entity_rank() < e_from.entity_rank();
     if ( is_degenerate_relation ( forward , *fi ) && downstream )
        {
        found_degenerate_relation = true;
        degenerate_key = fi->entity()->key();
        }
     ++fi;
     }
  if ( found_degenerate_relation )
     {
     std::ostringstream msg ;
     msg << method << "( from " ;
     print_entity_key( msg , m_mesh_meta_data, e_from.key() );
     msg << " , to " ;
     print_entity_key( msg , m_mesh_meta_data, e_to.key() );
     msg << " , id " << local_id ;
     msg << " ) FAILED ";
     msg << " Relation already exists to " ;
     print_entity_key( msg , m_mesh_meta_data, degenerate_key );
     throw std::runtime_error( msg.str() );
     }

  // If the relation is not degenerate, we add it and its converse

  if ( fe == fi || forward.attribute() != fi->attribute() ) {

    // A new relation and its converse

    const Relation converse( e_from , local_id );
    const std::vector<Relation>::iterator ce = e_to.m_relation.end();
          std::vector<Relation>::iterator ci = e_to.m_relation.begin();

    ci = std::lower_bound( ci , ce , converse , LessRelation() );

    if ( ce == ci || converse != *ci ) {
      fi = e_from.m_relation.insert( fi , forward );
      ci = e_to  .m_relation.insert( ci , converse );

      e_from.m_sync_count = m_sync_count ;
      e_to  .m_sync_count = m_sync_count ;

      PartVector add , empty ;

      // Deduce and set new part memberships:

      induced_part_membership( e_from, empty,
                               e_to.entity_rank(), local_id, add );

      internal_change_entity_parts( e_to , add , empty );

      set_field_relations( e_from , e_to , local_id );
    }
    else {
     /* this is unreachable unless a friend of bulk data creates a half-edge
        in the relationship graph. */
      std::ostringstream msg ;
      msg << method << "( from "
          << print_entity_key( msg , m_mesh_meta_data, e_from.key() )
          << " , to "
          << print_entity_key( msg , m_mesh_meta_data, e_to.key() )
          << " , id " << local_id
          << " ) FAILED"
          << " Internal error - converse relation already exists" ;
      throw std::runtime_error( msg.str() );
    }
  }

  // m_transaction_log.modify_entity ( e_from );
  // m_transaction_log.modify_sole_entity ( e_to );
}

//----------------------------------------------------------------------

void BulkData::declare_relation( Entity & entity ,
                                 const std::vector<Relation> & rel )
{
  static const char method[] = "stk::mesh::BulkData::declare_relation" ;

  assert_ok_to_modify( method );

  const unsigned etype = entity.entity_rank();

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity & e = * i->entity();
    const unsigned n = i->identifier();
    if ( e.entity_rank() < etype ) {
      declare_relation( entity , e , n );
    }
    else if ( etype < e.entity_rank() ) {
      declare_relation( e , entity , n );
    }
    else {
      throw std::runtime_error( std::string("stk::mesh::BulkData::declare_relation FAILED: given entities of the same entity type") );
    }
  }
}

//----------------------------------------------------------------------

void BulkData::destroy_relation( Entity & e_from , Entity & e_to )
{
  static const char method[]= "stk::mesh::BulkData::destroy_relation" ;

  assert_ok_to_modify( method );

  assert_valid_relation( method , *this , e_from , e_to );


  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  PartVector del , keep ;

  for ( std::vector<Relation>::iterator
        i = e_to.m_relation.end() ; i != e_to.m_relation.begin() ; ) {
    --i ;
    if ( i->entity() == & e_from ) {
      i = e_to.m_relation.erase( i );
    }
    else if ( e_to.entity_rank() < i->entity_rank() ) {
      induced_part_membership( * i->entity(), del, e_to.entity_rank(),
                               i->identifier(), keep );
    }
  }

  for ( std::vector<Relation>::iterator
        i = e_from.m_relation.end() ; i != e_from.m_relation.begin() ; ) {
    --i ;
    if ( i->entity() == & e_to ) {

      induced_part_membership( e_from, keep, e_to.entity_rank(),
                               i->identifier(), del );

      clear_field_relations( e_from , e_to.entity_rank() ,
                             i->identifier() );

      i = e_from.m_relation.erase( i );
    }
  }

  //------------------------------
  // 'keep' contains the parts deduced from kept relations
  // 'del'  contains the parts deduced from deleted relations
  //        that are not in 'keep'
  // Only remove these part memberships the entity is not shared.
  // If the entity is shared then wait until modificaton_end_synchronize.
  //------------------------------

  if ( ! del.empty() && (parallel_size() < 2 || e_to.sharing().empty()) ) {
    PartVector add ;
    internal_change_entity_parts( e_to , add , del );
  }

  // Mark e_from and e_to as modified
  // m_transaction_log.modify_entity ( e_from );
  // m_transaction_log.modify_sole_entity ( e_to );

}

//----------------------------------------------------------------------
// Deduce propagation of part membership changes to a 'from' entity
// to the related 'to' entities.  There can be both additions and
// removals.

void BulkData::internal_propagate_part_changes(
  Entity           & entity ,
  const PartVector & removed )
{
  const unsigned etype = entity.entity_rank();

  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel ) {
    const unsigned rel_type  = rel->entity_rank();
    const unsigned rel_ident = rel->identifier();

    if ( rel_type < etype ) { // a 'to' entity

      Entity & e_to = * rel->entity();

      PartVector to_del , to_add , empty ;

      // Induce part membership from this relationship to
      // pick up any additions.
      induced_part_membership( entity, empty,
                               rel_type, rel_ident, to_add );

      if ( ! removed.empty() ) {
        // Something was removed from the 'from' entity,
        // deduce what may have to be removed from the 'to' entity.

        // Deduce parts for 'e_to' from all upward relations.
        // Any non-parallel part that I removed that is not deduced for
        // 'e_to' must be removed from 'e_to'

        for ( PairIterRelation
              to_rel = e_to.relations(); ! to_rel.empty() ; ++to_rel ) {
          if ( e_to.entity_rank() < to_rel->entity_rank() &&
               & entity != to_rel->entity() /* Already did this entity */ ) {
            // Relation from to_rel->entity() to e_to
            induced_part_membership( * to_rel->entity(), empty,
                                     e_to.entity_rank(),
                                     to_rel->identifier(),
                                     to_add );
          }
        }

        for ( PartVector::const_iterator
              j = removed.begin() ; j != removed.end() ; ++j ) {
          if ( ! contain( to_add , **j ) ) {
            induced_part_membership( **j, etype, rel_type, rel_ident, to_del );
          }
        }
      }

      if ( parallel_size() < 2 || e_to.sharing().empty() ) {
        // Entirely local, ok to remove memberships now
        internal_change_entity_parts( e_to , to_add , to_del );
      }
      else {
        // Shared, do not remove memberships now.
        // Wait until modification_end.
        internal_change_entity_parts( e_to , to_add , empty );
      }

      set_field_relations( entity, e_to, rel_ident );
    }
    else if ( etype < rel_type ) { // a 'from' entity
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel_ident );
    }
  }
}

//----------------------------------------------------------------------

void BulkData::internal_propagate_relocation( Entity & entity )
{
  const unsigned etype = entity.entity_rank();
  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel ) {
    const unsigned rel_type = rel->entity_rank();
    if ( rel_type < etype ) {
      Entity & e_to = * rel->entity();

      set_field_relations( entity, e_to, rel->identifier() );
    }
    else if ( etype < rel_type ) {
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel->identifier() );
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

