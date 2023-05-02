// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.




#include <memory>
#include <stdexcept>
#include <vector>
#include <string>

#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
//#include <stk_mesh/base/FieldData.hpp>
//#include <stk_mesh/fem/FEMMetaData.hpp>

#include <percept/stk_rebalance/Rebalance.hpp>
#include <percept/stk_rebalance_utils/RebalanceUtils.hpp>
#include <percept/stk_rebalance/Partition.hpp>

namespace stk {
  namespace rebalance {

bool Rebalance::balance_comm_spec_domain( Partition * partition,
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

void Rebalance::rebalance_dependent_entities( const mesh::BulkData    & bulk_data ,
                                   const Partition         * partition,
                                   const mesh::EntityRank  & dep_rank,
                                   mesh::EntityProcVec     & entity_procs,
                                   const stk::mesh::EntityRank rank)
{

  size_t init_size = entity_procs.size();
  //stk::mesh::MetaData & fem_meta = stk::mesh::MetaData::get(bulk_data);
  const stk::mesh::EntityRank element_rank = (rank != stk::mesh::InvalidEntityRank) ? rank :
    stk::topology::ELEMENT_RANK;

  if (bulk_data.parallel_rank() == 0)
    {
      std::cout << "Reblance::rebalance_dependent_entities rank_to_reb= " << rank
                << " dep_rank = " << dep_rank << std::endl;
    }
  // why??? srk 06/06/14
  if (dep_rank == element_rank) return;
  // Create a map of ids of migrating elements to their owner proc and a vector of the migrating elements.
  std::map<mesh::EntityId, unsigned> elem_procs;
  mesh::EntityVector owned_moving_elems;
  mesh::EntityProcVec::iterator ep_iter = entity_procs.begin(),
                                 ep_end = entity_procs.end();
  for( ; ep_end != ep_iter; ++ep_iter ) {
    if( element_rank == bulk_data.entity_rank(ep_iter->first) )
    {
      const mesh::EntityId elem_id = bulk_data.identifier(ep_iter->first);
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
    const mesh::EntityId elem_id = bulk_data.identifier(*r_iter);
    mesh::EntityVector related_entities;
    mesh::EntityVector elems(1);
    elems[0] = *r_iter;
    stk::mesh::get_entities_through_relations(bulk_data, elems, dep_rank, related_entities);
    for( size_t j = 0; j < related_entities.size(); ++j ) {
      //if (bulk_data.owner_rank(related_entities[j]) == bulk_data.parallel_rank())
        {
          dep_entity_procs[bulk_data.identifier(related_entities[j])] = elem_procs[elem_id];
        }
    }
  }


  std::map<mesh::EntityId, unsigned>::const_iterator c_iter = dep_entity_procs.begin(),
                                                      c_end = dep_entity_procs.end();
  for( ; c_end != c_iter; ++c_iter )
  {
    mesh::Entity  de = bulk_data.get_entity( dep_rank, c_iter->first );
    if( parallel_machine_rank(partition->parallel()) == bulk_data.parallel_owner_rank(de))
    {
      stk::mesh::EntityProc dep_proc(de, c_iter->second);
      entity_procs.push_back(dep_proc);
    }
  }

  size_t final_size = entity_procs.size();
  if (0 && bulk_data.parallel_rank() == 0)
    {
      std::cout << "rebalance_dependent_entities: for dep_rank = " << dep_rank << " rank= " << rank << " init_size= " << init_size
                << " final_size= " << final_size << std::endl;
    }
}


bool Rebalance::full_rebalance(mesh::BulkData  & bulk_data ,
                    Partition       * partition,
                    const stk::mesh::EntityRank rank)
{
  mesh::EntityProcVec cs_elem;
  bool rebalancingHasOccurred =  balance_comm_spec_domain( partition, cs_elem );

  if(rebalancingHasOccurred && partition->partition_dependents_needed() )
  {
    //stk::mesh::MetaData & fem_meta = stk::mesh::MetaData::get(bulk_data);

    const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
    const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
    const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
    const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
    // note, this could also be the FAMILY_TREE_RANK...
    const stk::mesh::EntityRank constraint_rank = static_cast<stk::mesh::EntityRank>(elem_rank+1);

    // Don't know the rank of the elements rebalanced, assume all are dependent.
    rebalance_dependent_entities( bulk_data, partition, node_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != edge_rank && rank != edge_rank)
      rebalance_dependent_entities( bulk_data, partition, edge_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != face_rank && rank != face_rank)
      rebalance_dependent_entities( bulk_data, partition, face_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != elem_rank && rank != elem_rank)
      rebalance_dependent_entities( bulk_data, partition, elem_rank, cs_elem, rank );
    if (stk::mesh::InvalidEntityRank != constraint_rank && rank != constraint_rank)
      rebalance_dependent_entities( bulk_data, partition, constraint_rank, cs_elem, rank );
  }

  if ( rebalancingHasOccurred )
  {
    bulk_data.change_entity_owner( cs_elem );
  }

  //: Finished
  return rebalancingHasOccurred;
}


bool Rebalance::rebalance(mesh::BulkData   & bulk_data  ,
                               const mesh::Selector  & selector ,
                               const stk::mesh::FieldBase     * rebal_coord_ref ,
                               const stk::mesh::FieldBase     * rebal_elem_weight_ref ,
                               Partition & partition,
                               const stk::mesh::EntityRank rank,
                               stk::mesh::EntityVector *entities_to_rebalance ,
                               bool do_rebal,
                               bool debug)
{
  //stk::mesh::MetaData &fem_meta = stk::mesh::MetaData::get(bulk_data);
  const stk::mesh::EntityRank element_rank = (rank != stk::mesh::InvalidEntityRank) ? rank :
    stk::topology::ELEMENT_RANK;

  mesh::EntityVector rebal_elem_ptrs;
  mesh::EntityVector entities;

  const mesh::MetaData &meta_data = bulk_data.mesh_meta_data();
  if (entities_to_rebalance)
    {
      entities = *entities_to_rebalance;
    }
  else
    {
      mesh::get_selected_entities(selector & meta_data.locally_owned_part(),
                                  bulk_data.buckets(element_rank),
                                  entities);
    }

  stk::mesh::Selector locally_owned(meta_data.locally_owned_part());

  //std::cout << "element_rank= " << element_rank << " selector= " << selector << std::endl;
  //stk::rebalance::check_ownership(bulk_data, entities, "here 1");
  for (mesh::EntityVector::iterator iA = entities.begin() ; iA != entities.end() ; ++iA ) {
    if (!locally_owned(bulk_data.bucket(*iA)))
      throw std::runtime_error("Rebalance::rebalance: entities must be locally owned");
    if(rebal_elem_weight_ref)
    {
      double * const w = static_cast<double*>(mesh::field_data( *rebal_elem_weight_ref, *iA ));
      STK_ThrowRequireMsg( NULL != w,
        "Rebalance weight field is not defined on entities but should be defined on all entities.");
      // Should this be a throw instead???
      if ( *w <= 0.0 ) {
        *w = 1.0 ;
      }
    }
    rebal_elem_ptrs.push_back( *iA );
  }

  partition.set_mesh_info(
      rebal_elem_ptrs,
      rebal_coord_ref,
      rebal_elem_weight_ref);

  if (debug)
    {
      const int  print_stats = 1;
      int        nentity     = 0;
      double     entity_wgt  = 0;
      int        ncuts       = 0;
      double     cut_wgt     = 0;
      int        nboundary   = 0;
      int        nadj        = 0;
      const int ierr = partition.evaluate (print_stats, &nentity, &entity_wgt, &ncuts, &cut_wgt, &nboundary, &nadj);
      if (0 && bulk_data.parallel_rank() == 0)
        {
          std::cout <<" Information returned from the Zoltan evaluate function:"<<std::endl;
          std::cout <<" Error Code:             :"<<ierr      <<std::endl;
          std::cout <<" Number of entities:     :"<<nentity   <<std::endl;
          std::cout <<" Number of cuts:         :"<<ncuts     <<std::endl;
          std::cout <<" Cut Weight:             :"<<cut_wgt   <<std::endl;
          std::cout <<" Number on Boundary:     :"<<nboundary <<std::endl;
          std::cout <<" Number Adjancent:       :"<<nadj      <<std::endl;
        }
    }


  bool rebalancingHasOccurred = false;
  if (do_rebal)
    full_rebalance(bulk_data, &partition, rank);
  return rebalancingHasOccurred;
}

}
}
