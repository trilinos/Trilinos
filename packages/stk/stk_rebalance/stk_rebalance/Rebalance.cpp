/*--------------------------------------------------------------------*/
/*    Copyright 2001, 2008, 2009, 2010 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001,2002 Sandia Corporation, Albuquerque, NM.

#include <memory>
#include <stdexcept>
#include <vector>
#include <string>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>

using namespace stk;
using namespace stk::rebalance;

namespace {

bool balance_comm_spec_domain( Partition * partition,
                               mesh::EntityProcVec & rebal_spec )
{
  bool rebalancingHasOccurred = false;
  {
    int num_elems = partition->num_elems();
    int tot_elems;
    all_reduce_sum(partition->parallel(), &num_elems, &tot_elems, 1);

    if (tot_elems) {
      partition->determine_new_partition(rebalancingHasOccurred);
    }
  }
  if (rebalancingHasOccurred) partition->get_new_partition(rebal_spec);

  return rebalancingHasOccurred;
}


/*
 * Traversing the migrating elements in reverse order produces a simplistic
 * attempt at lowest-rank element proc greedy partitioning of dependents
 * which seems to often work in practice.  Some logic could be added here
 * as needed to enforce more deterministic dependent partitions.
 */

void rebalance_dependent_entities( const mesh::BulkData    & bulk_data ,
                                   const Partition         * partition,
                                   const mesh::EntityRank  & dep_rank,
                                   mesh::EntityProcVec     & entity_procs,
                                   const stk::mesh::EntityRank rank)
{

  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(bulk_data);
  const stk::mesh::EntityRank element_rank = (rank != stk::mesh::InvalidEntityRank) ? rank :
                                             fem_meta.element_rank();

  if (dep_rank == element_rank) return;
  // Create a map of ids of migrating elements to their owner proc and a vector of the migrating elements.
  std::map<mesh::EntityId, unsigned> elem_procs;
  mesh::EntityVector owned_moving_elems;
  mesh::EntityProcVec::iterator ep_iter = entity_procs.begin(),
                                 ep_end = entity_procs.end();
  for( ; ep_end != ep_iter; ++ep_iter ) {
    if( element_rank == ep_iter->first->entity_rank() )
    {
      const mesh::EntityId elem_id = ep_iter->first->identifier();
      elem_procs[elem_id] = ep_iter->second;
      owned_moving_elems.push_back(ep_iter->first);
    }
  }
  // TODO: Determine if this "dumb" greedy approach is adequate and the cost/benefit
  //       of doing something more sophisticated

  // This reverse traversal of elements overwrites assignment of procs for
  // dependents resulting in the last assignment winning.

  // For all dep-rank entities related to migrating elements, pack their info in to
  // dep_entity_procs.
  std::map<mesh::EntityId, unsigned> dep_entity_procs;
  mesh::EntityVector::reverse_iterator r_iter = owned_moving_elems.rbegin(),
                                        r_end = owned_moving_elems.rend();
  for( ; r_end != r_iter; ++r_iter )
  {
    const mesh::EntityId elem_id = (*r_iter)->identifier();
    mesh::EntityVector related_entities;
    mesh::EntityVector elems(1);
    elems[0] = *r_iter;
    stk::mesh::get_entities_through_relations(elems, dep_rank, related_entities);
    for( size_t j = 0; j < related_entities.size(); ++j ) {
      dep_entity_procs[related_entities[j]->identifier()] = elem_procs[elem_id];
    }
  }


  std::map<mesh::EntityId, unsigned>::const_iterator c_iter = dep_entity_procs.begin(),
                                                      c_end = dep_entity_procs.end();
  for( ; c_end != c_iter; ++c_iter )
  {
    mesh::Entity * de = bulk_data.get_entity( dep_rank, c_iter->first );
    if( parallel_machine_rank(partition->parallel()) == de->owner_rank() )
    {
      stk::mesh::EntityProc dep_proc(de, c_iter->second);
      entity_procs.push_back(dep_proc);
    }
  }
}


bool full_rebalance(mesh::BulkData  & bulk_data ,
                    Partition       * partition,
                    const stk::mesh::EntityRank rank)
{
  mesh::EntityProcVec cs_elem;
  bool rebalancingHasOccurred =  balance_comm_spec_domain( partition, cs_elem );

  if(rebalancingHasOccurred && partition->partition_dependents_needed() )
  {
    stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(bulk_data);

    const stk::mesh::EntityRank node_rank = fem_meta.node_rank();
    const stk::mesh::EntityRank edge_rank = fem_meta.edge_rank();
    const stk::mesh::EntityRank face_rank = fem_meta.face_rank();
    const stk::mesh::EntityRank elem_rank = fem_meta.element_rank();
    const stk::mesh::EntityRank cons_rank = elem_rank+1;

    // Don't know the rank of the elements rebalanced, assume all are dependent.
    rebalance_dependent_entities( bulk_data, partition, node_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != edge_rank && rank != edge_rank)
      rebalance_dependent_entities( bulk_data, partition, edge_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != face_rank && rank != face_rank)
      rebalance_dependent_entities( bulk_data, partition, face_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != elem_rank && rank != elem_rank)
      rebalance_dependent_entities( bulk_data, partition, elem_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != cons_rank && rank != cons_rank)
      rebalance_dependent_entities( bulk_data, partition, cons_rank, cs_elem, rank );
  }

  if ( rebalancingHasOccurred )
  {
    bulk_data.modification_begin();
    bulk_data.change_entity_owner( cs_elem );
    bulk_data.modification_end();
  }

  //: Finished
  return rebalancingHasOccurred;
}
} // namespace


bool stk::rebalance::rebalance(mesh::BulkData   & bulk_data  ,
                               const mesh::Selector  & selector ,
                               const VectorField     * rebal_coord_ref ,
                               const ScalarField     * rebal_elem_weight_ref ,
                               Partition & partition,
                               const stk::mesh::EntityRank rank)
{
  stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk_data);
  const stk::mesh::EntityRank element_rank = (rank != stk::mesh::InvalidEntityRank) ? rank :
                                             fem_meta.element_rank();

  mesh::EntityVector rebal_elem_ptrs;
  mesh::EntityVector entities;

  mesh::get_selected_entities(selector,
                              bulk_data.buckets(element_rank),
                              entities);

  for (mesh::EntityVector::iterator iA = entities.begin() ; iA != entities.end() ; ++iA ) {
    if(rebal_elem_weight_ref)
    {
      double * const w = mesh::field_data( *rebal_elem_weight_ref, **iA );
      ThrowRequireMsg( NULL != w,
        "Rebalance weight field is not defined on entities but should be defined on all entities.");
      // Should this be a throw instead???
      if ( *w <= 0.0 ) {
        *w = 1.0 ;
      }
    }
    rebal_elem_ptrs.push_back( *iA );
  }

  (&partition)->set_mesh_info(
      rebal_elem_ptrs,
      rebal_coord_ref,
      rebal_elem_weight_ref);

  bool rebalancingHasOccurred = full_rebalance(bulk_data, &partition, rank);

  return rebalancingHasOccurred;
}
