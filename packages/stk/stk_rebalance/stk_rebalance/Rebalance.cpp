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


#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>

using namespace stk;
using namespace rebalance;


namespace {


bool balance_comm_spec_domain (Partition & partition,
                               std::vector<mesh::EntityProc> & rebal_spec)
{
  bool rebalancingHasOccurred = false;
  {
    int num_elems = partition.num_elems();
    int tot_elems;
    //all_reduce_sum(comm, &num_elems, &tot_elems, 1);
    tot_elems = num_elems;
    
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

// ----------------------------------------------------------------------------
}

// ------------------------------------------------------------------------------

bool rebalance::rebalance_needed(mesh::BulkData &    bulk_data,
                                 mesh::MetaData &    meta_data,
                                 const mesh::Field<double> & load_measure,
                                 ParallelMachine    comm,
				 const double  imbalance_threshold)
{
  if ( imbalance_threshold < 1.0 ) return true;

  double my_load = 0.0;

  mesh::EntityVector local_nodes;
  mesh::Selector select_owned(meta_data.locally_owned_part());
  mesh::get_selected_entities(select_owned,
                              bulk_data.buckets(mesh::Node),
                              local_nodes);

  for(mesh::EntityVector::iterator elem_it = local_nodes.begin(); elem_it != local_nodes.end(); ++elem_it)
  {
    double * load_val = mesh::field_data(load_measure, **elem_it);
    my_load += *load_val;
  }

  double max_load = my_load=my_load;
  double tot_load = my_load=0;

  all_reduce(comm, ReduceMax<1>(&my_load));
  all_reduce_sum(comm, &my_load, &tot_load, 1);

  const int   proc_size = parallel_machine_size(comm);
  const double avg_load = tot_load / proc_size;

  return ( max_load / avg_load > imbalance_threshold );
}

bool rebalance::rebalance(mesh::BulkData        & bulk_data  ,
                          const mesh::Selector  & selector ,
                          const mesh::Field<double> * rebal_coord_ref ,
                          const mesh::Field<double> * rebal_elem_weight_ref ,
                          Partition & partition)
{
  mesh::EntityVector rebal_elem_ptrs;
  mesh::EntityVector entities;
  mesh::get_selected_entities(selector,
                              bulk_data.buckets(mesh::Node),
                              entities);

  for (mesh::EntityVector::iterator iA = entities.begin() ; iA != entities.end() ; ++iA ) {
    if(rebal_elem_weight_ref)
    {
      double * const w = mesh::field_data( *rebal_elem_weight_ref, **iA );
      ThrowRequire( NULL != w );
      if ( *w <= 0.0 ) {
        *w = 1.0 ;
      }
    }
    rebal_elem_ptrs.push_back( *iA );
  }

  partition.replace_mesh(
      rebal_elem_ptrs,
      rebal_coord_ref,
      rebal_elem_weight_ref);

  bool rebalancingHasOccurred =  full_rebalance(bulk_data, partition);

  return rebalancingHasOccurred;
}
