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

//----------------------------------------------------------------------

void BulkData::require_valid_relation( const char action[] ,
                                       const BulkData & mesh ,
                                       const Entity e_from ,
                                       const Entity e_to )
{
  const bool error_type      = mesh.entity_rank(e_from) <= mesh.entity_rank(e_to);
  const bool error_nil_from  = !mesh.is_valid(e_from);
  const bool error_nil_to    = !mesh.is_valid(e_to);

  if ( error_type || error_nil_from || error_nil_to ) {
    std::ostringstream msg ;

    msg << "Could not " << action << " relation from entity "
        << print_entity_key(MetaData::get(mesh), mesh.entity_key(e_from)) << " to entity "
        << print_entity_key(MetaData::get(mesh), mesh.entity_key(e_to)) << "\n";

    ThrowErrorMsgIf( error_nil_from  || error_nil_to,
                     msg.str() << ", entity was destroyed");
    ThrowErrorMsgIf( error_type, msg.str() <<
                     "A relation must be from higher to lower ranking entity");
  }
}


//----------------------------------------------------------------------
bool BulkData::internal_declare_relation(Entity e_from, Entity e_to,
                                         RelationIdentifier local_id,
                                         unsigned sync_count, bool is_back_relation,
                                         Permutation permut)
{
  TraceIfWatching("stk::mesh::BuilkData::declare_relation", LOG_ENTITY, entity_key(e_from));

  const MeshIndex& idx = mesh_index(e_from);

  bool modified = idx.bucket->declare_relation(idx.bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id),
                                               permut);

  if (modified) {
    set_synchronized_count( e_from, sync_count );
  }
  return modified;
}

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut)
{
  TraceIfWatching("stk::mesh::BulkData::declare_relation", LOG_ENTITY, entity_key(e_from));
  TraceIfWatchingDec("stk::mesh::BulkData::declare_relation", LOG_ENTITY, entity_key(e_to), 1);
  DiagIfWatching(LOG_ENTITY, entity_key(e_from),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, entity_key(e_to),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "declare" , *this , e_from , e_to );

  // TODO: Don't throw if exact relation already exists, that should be a no-op.
  // Should be an exact match if relation of local_id already exists (e_to should be the same).
  bool is_converse = false;
  bool caused_change_fwd = internal_declare_relation(e_from, e_to, local_id, m_sync_count,
                                                     is_converse, permut);

  //TODO: check connectivity map
  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {

    const bool higher_order_relation = stk::topology::ELEMENT_RANK < entity_rank(e_from);
    if (    higher_order_relation
         || m_bucket_repository.connectivity_map()(entity_rank(e_to), entity_rank(e_from)) != stk::mesh::INVALID_CONNECTIVITY_TYPE
       )
    {
      // the setup for the converse relationship works slightly differently
      is_converse = true;
      internal_declare_relation(e_to, e_from, local_id, m_sync_count, is_converse, permut );
    }
  }

  // It is critical that the modification be done AFTER the relations are
  // added so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    this->modified(e_to);
    this->modified(e_from);
  }

  OrdinalVector add , empty ;

  // Deduce and set new part memberships:

  induced_part_membership(*this, e_from, empty, entity_rank(e_to), add );

  internal_change_entity_parts( e_to , add , empty );
}

//----------------------------------------------------------------------

void BulkData::declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel )
{
  require_ok_to_modify();

  const unsigned erank = entity_rank(entity);

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity e = i->entity();
    const unsigned n = i->relation_ordinal();
    const Permutation permut = static_cast<Permutation>(i->getOrientation());
    if ( entity_rank(e) < erank ) {
      declare_relation( entity , e , n, permut );
    }
    else if ( erank < entity_rank(e) ) {
      declare_relation( e , entity , n, permut );
    }
    else {
      ThrowErrorMsg("Given entities of the same entity rank. entity is " <<
                    identifier(entity));
    }
  }
}

//----------------------------------------------------------------------

bool BulkData::destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
  TraceIfWatching("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, entity_key(e_from));
  TraceIfWatchingDec("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, entity_key(e_to), 1);
  DiagIfWatching(LOG_ENTITY, entity_key(e_from),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, entity_key(e_to),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "destroy" , *this , e_from , e_to );

  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();
  const EntityRank e_to_entity_rank = entity_rank(e_to);
  const EntityRank e_from_entity_rank = entity_rank(e_from);

  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  m_check_invalid_rels = false; // OK to have gaps when deleting

  if ( parallel_size() < 2 || entity_comm_sharing(entity_key(e_to)).empty() ) {

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
    EntityVector temp_entities;
    std::vector<ConnectivityOrdinal> temp_ordinals;
    Entity const* rel_entities = NULL;
    ConnectivityOrdinal const* rel_ordinals = NULL;
    int num_rels = 0;
    for (EntityRank irank = e_to_entity_rank + 1; irank < end_rank; ++irank)
    {
      if (connectivity_map().valid(e_to_entity_rank, irank)) {
        num_rels     = num_connectivity(e_to, irank);
        rel_entities = begin(e_to, irank);
        rel_ordinals = begin_ordinals(e_to, irank);
      }
      else {
        num_rels     = get_connectivity(*this, e_to, irank, temp_entities, temp_ordinals);
        rel_entities = &*temp_entities.begin();
        rel_ordinals = &*temp_ordinals.begin();
      }

      for (int j = 0; j < num_rels; ++j)
      {
        ThrowAssertMsg(is_valid(rel_entities[j]), "Error, entity " << e_to.local_offset() << " has invalid back-relation for ordinal: "
                       << rel_ordinals[j] << " to rank: " << irank << ", target entity is: " << rel_entities[j].local_offset());
        if ( !(rel_entities[j] == e_from && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) ) )
        {
          induced_part_membership(*this, rel_entities[j], empty, e_to_entity_rank, keep,
                                  false /*Do not look at supersets*/);
        }
      }
    }

    // Find the relation this is being deleted and add the parts that are
    // induced from that relation (and that are not in 'keep') to 'del'
    {
      size_t num_rels = num_connectivity(e_from, e_to_entity_rank);
      Entity const *rel_entities = begin(e_from, e_to_entity_rank);
      ConnectivityOrdinal const *rel_ordinals = begin_ordinals(e_from, e_to_entity_rank);
      for (size_t j = 0; j < num_rels; ++j)
      {
        if ( rel_entities[j] == e_to && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) )
        {
          induced_part_membership(*this, e_from, keep, e_to_entity_rank, del,
                                  false /*Do not look at supersets*/);
          break; // at most 1 relation can match our specification
        }
      }
    }

    if ( !del.empty() ) {
      internal_change_entity_parts( e_to , empty , del );
    }
  }

  //delete relations from the entities
  bool caused_change_fwd = bucket(e_from).destroy_relation(e_from, e_to, local_id);

  // Relationships should always be symmetrical
  if ( caused_change_fwd &&
       (e_to_entity_rank > stk::topology::ELEMENT_RANK || e_from_entity_rank > stk::topology::ELEMENT_RANK ||
        connectivity_map().valid(entity_rank(e_to), entity_rank(e_from))) ) {
    bucket(e_to).destroy_relation(e_to, e_from, local_id);
  }

  // It is critical that the modification be done AFTER the relations are
  // changed so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    this->modified(e_to);
    this->modified(e_from);
  }

  m_check_invalid_rels = true;

  return caused_change_fwd;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk
