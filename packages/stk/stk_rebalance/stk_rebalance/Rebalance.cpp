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

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>


#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>

using namespace stk;
using namespace stk::rebalance;


namespace {

bool balance_comm_spec_domain( Partition & partition,
                               std::vector<mesh::EntityProc> & rebal_spec )
{
  bool rebalancingHasOccurred = false;
  {
    int num_elems = partition.num_elems();
    int tot_elems;
    all_reduce_sum(partition.parallel(), &num_elems, &tot_elems, 1);

    if (tot_elems) {
      partition.determine_new_partition(rebalancingHasOccurred);
    }
  }
  if (rebalancingHasOccurred) partition.get_new_partition(rebal_spec);

  return rebalancingHasOccurred;
}


bool full_rebalance(mesh::BulkData        & bulk_data ,
                    Partition             & partition)
{
  std::vector<mesh::EntityProc> cs_elem;
  bool rebalancingHasOccurred =  balance_comm_spec_domain( partition,
                                                           cs_elem );

  if ( rebalancingHasOccurred ) {
    bulk_data.modification_begin();
    bulk_data.change_entity_owner( cs_elem );
    bulk_data.modification_end();
  }

  //: Finished
  return rebalancingHasOccurred;
}

}

// ------------------------------------------------------------------------------

bool stk::rebalance::rebalance_needed(mesh::BulkData &    bulk_data,
                                 const mesh::Field<double> & load_measure,
                                 ParallelMachine    comm,
				 double & imbalance_threshold)
{
  // Need to make load_measure optional with weights defaulting to 1.0. ??

  if ( imbalance_threshold < 1.0 ) return true;

  double my_load = 0.0;

  const mesh::MetaData & meta_data = bulk_data.mesh_meta_data();
  const mesh::TopologicalMetaData & topo_data = mesh::TopologicalMetaData::find_TopologicalMetaData(meta_data);

  mesh::EntityVector local_elems;
  mesh::Selector select_owned( meta_data.locally_owned_part() );

  // Determine imbalance based on current element decomposition
  mesh::get_selected_entities(select_owned,
                              bulk_data.buckets(topo_data.element_rank),
                              local_elems);

  for(mesh::EntityVector::iterator elem_it = local_elems.begin(); elem_it != local_elems.end(); ++elem_it)
  {
    double * load_val = mesh::field_data(load_measure, **elem_it);
    my_load += *load_val;
  }

  double max_load = my_load;
  double tot_load = my_load;

  all_reduce(comm, ReduceMax<1>(&max_load) & ReduceSum<1>(&tot_load));

  const int proc_size = parallel_machine_size(comm);
  const double avg_load = tot_load / proc_size;

  bool need_to_rebalance =  max_load / avg_load > imbalance_threshold;
  imbalance_threshold = max_load / avg_load;

  return need_to_rebalance;
}

bool stk::rebalance::rebalance(mesh::BulkData        & bulk_data  ,
                          const mesh::Selector  & selector ,
                          const VectorField * rebal_coord_ref ,
                          const ScalarField * rebal_elem_weight_ref ,
                          Partition & partition)
{
  mesh::EntityVector rebal_elem_ptrs;
  mesh::EntityVector entities;
  const mesh::TopologicalMetaData & topo_data = mesh::TopologicalMetaData::find_TopologicalMetaData(bulk_data.mesh_meta_data());
  mesh::get_selected_entities(selector,
                              bulk_data.buckets(topo_data.element_rank),
                              entities);

  for (mesh::EntityVector::iterator iA = entities.begin() ; iA != entities.end() ; ++iA ) {
    if(rebal_elem_weight_ref)
    {
      double * const w = mesh::field_data( *rebal_elem_weight_ref, **iA );
      ThrowRequire( NULL != w );
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

  bool rebalancingHasOccurred = full_rebalance(bulk_data, partition);

  return rebalancingHasOccurred;
}
