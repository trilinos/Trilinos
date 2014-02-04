/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_sum
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/Bucket.hpp"     // for has_superset, Bucket
#include "stk_mesh/base/Types.hpp"      // for EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { class Part; } }



namespace stk {
namespace mesh {

bool comm_mesh_counts( const BulkData & M ,
                       std::vector<size_t> & counts ,
                       bool local_flag )
{
  const size_t zero = 0 ;

  // Count locally owned entities

  const MetaData & S = MetaData::get(M);
  const EntityRank entity_rank_count = static_cast<EntityRank>(S.entity_rank_count());
  const size_t     comm_count        = entity_rank_count + 1 ;

  std::vector<size_t> local(  comm_count , zero );
  std::vector<size_t> global( comm_count , zero );

  ParallelMachine comm = M.parallel();
  Part & owns = S.locally_owned_part();

  for ( EntityRank i = stk::topology::NODE_RANK ; i < entity_rank_count ; ++i ) {
    const std::vector<Bucket*> & ks = M.buckets( i );

    std::vector<Bucket*>::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( has_superset( **ik , owns ) ) {
        local[i] += (*ik)->size();
      }
    }
  }

  local[ entity_rank_count ] = local_flag ;

  all_reduce_sum( comm , & local[0] , & global[0] , comm_count );

  counts.assign( global.begin() , global.begin() + entity_rank_count );

  return 0 < global[ entity_rank_count ] ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

