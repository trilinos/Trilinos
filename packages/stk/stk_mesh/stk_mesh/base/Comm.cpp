/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>


namespace stk {
namespace mesh {

bool comm_mesh_counts( BulkData & M ,
                       std::vector<size_t> & counts ,
                       bool local_flag )
{
  const size_t zero = 0 ;

  // Count locally owned entities

  const MetaData & S = MetaData::get(M);
  const unsigned entity_rank_count = S.entity_rank_count();
  const size_t   comm_count        = entity_rank_count + 1 ;

  std::vector<size_t> local(  comm_count , zero );
  std::vector<size_t> global( comm_count , zero );

  ParallelMachine comm = M.parallel();
  Part & owns = S.locally_owned_part();

  for ( unsigned i = 0 ; i < entity_rank_count ; ++i ) {
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

