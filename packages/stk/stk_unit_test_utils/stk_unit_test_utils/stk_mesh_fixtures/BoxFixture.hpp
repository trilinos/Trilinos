// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef Stk_Mesh_Fixtures_BoxFixture_hpp
#define Stk_Mesh_Fixtures_BoxFixture_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/MetaData.hpp>   // for entity_rank_names, MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityRank
#include <string>                       // for string
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp"

namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {

static const size_t spatial_dimension = 3;
/**
 * A fixture that creates a "box" mesh of hexes
 */
class BoxFixture {
public:
  BoxFixture(stk::ParallelMachine pm = MPI_COMM_WORLD,
             stk::mesh::BulkData::AutomaticAuraOption autoAuraOption = stk::mesh::BulkData::AUTO_AURA,
             unsigned bucketCapacity = stk::mesh::get_default_bucket_capacity(),
             const std::vector<std::string>& entity_names = stk::mesh::entity_rank_names());

  ~BoxFixture () {}

  MetaData & fem_meta () { return m_fem_meta; }
  stk::unit_test_util::BulkDataTester & bulk_data () { return m_bulk_data; }

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

protected:
  MetaData m_fem_meta;
  stk::unit_test_util::BulkDataTester m_bulk_data;

  int m_comm_rank;
  int m_comm_size;

  NodeToProcsMMap m_nodes_to_procs;

  Part &m_elem_part;
  stk::topology m_elem_topology;

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

namespace simple_fields {

static const size_t spatial_dimension = 3;
/**
 * A fixture that creates a "box" mesh of hexes
 */
class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
BoxFixture {
public:
  BoxFixture(stk::ParallelMachine pm = MPI_COMM_WORLD,
             stk::mesh::BulkData::AutomaticAuraOption autoAuraOption = stk::mesh::BulkData::AUTO_AURA,
             unsigned bucketCapacity = stk::mesh::get_default_bucket_capacity(),
             const std::vector<std::string>& entity_names = stk::mesh::entity_rank_names());

  ~BoxFixture () {}

  MetaData & fem_meta () { return m_fem_meta; }
  stk::unit_test_util::BulkDataTester & bulk_data () { return m_bulk_data; }

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

protected:
  MetaData m_fem_meta;
  stk::unit_test_util::BulkDataTester m_bulk_data;

  int m_comm_rank;
  int m_comm_size;

  NodeToProcsMMap m_nodes_to_procs;

  Part &m_elem_part;
  stk::topology m_elem_topology;

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

} // namespace simple_fields

} // namespace fixtures
} // namespace mesh
} // namespace stk

#endif
