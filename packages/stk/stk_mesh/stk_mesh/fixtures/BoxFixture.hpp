/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Mesh_Fixtures_BoxFixture_hpp
#define Stk_Mesh_Fixtures_BoxFixture_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>   // for entity_rank_names, MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityRank
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <stk_mesh/fixtures/FixtureNodeSharing.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

static const size_t spatial_dimension = 3;
/**
 * A fixture that creates a "box" mesh of hexes
 */
class  BoxFixture {
public:
  BoxFixture(stk::ParallelMachine pm = MPI_COMM_WORLD,
             unsigned block_size = 1000,
             const std::vector<std::string>& entity_names = stk::mesh::entity_rank_names());

  ~BoxFixture () {}

  MetaData & fem_meta () { return m_fem_meta; }
  BulkData & bulk_data () { return m_bulk_data; }

  int  comm_size() const { return m_comm_size; }
  int  comm_rank() const { return m_comm_rank; }

  Part &get_elem_part() const { return m_elem_part; }

  stk::topology get_elem_topology() const { return m_elem_topology; }

  typedef int BOX[3][2];

  /**
   *  Generate a box mesh which is globally ( ngx X ngy X ngz )
   *  elements where:
   *    ngx = root_box[0][1] - root_box[0][0] ;
   *    ngy = root_box[1][1] - root_box[1][0] ;
   *    ngz = root_box[2][1] - root_box[2][0] ;
   *
   *  The box is partitioned via recursive coordinate bisection
   *  and the resulting local box are given in 'local_box'.
   */
  void generate_boxes (const BOX root_box,
                             BOX local_box);

  Entity get_new_entity ( EntityRank rank , EntityId parallel_dependent_id );

protected:
  MetaData m_fem_meta;
  BulkData m_bulk_data;

  int m_comm_rank;
  int m_comm_size;

  NodeToProcsMMap m_nodes_to_procs;

  Part &m_elem_part;
  stk::topology m_elem_topology;

  BulkData::BulkDataSyncState m_previous_state;

  /**
   * Recursively split a box into ( up - ip ) sub-boxes
   */
  static void box_partition( int ip , int up , int axis ,
                             const BOX box ,
                             BOX p_box[] );

private:

  BoxFixture();
  BoxFixture( const BoxFixture & );
  BoxFixture & operator = ( const BoxFixture & );

  void fill_node_map(int proc_rank, const BOX root_box);

};

} // namespace fixtures
} // namespace mesh
} // namespace stk

#endif
