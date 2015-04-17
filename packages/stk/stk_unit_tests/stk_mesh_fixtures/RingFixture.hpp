// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STK_MESH_FIXTURES_RING_FIXTURE_HPP
#define STK_MESH_FIXTURES_RING_FIXTURE_HPP

#include <stddef.h>                     // for size_t
#include <unit_tests/BulkDataTester.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId, PartVector
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <vector>                       // for vector
namespace stk { namespace mesh { class Part; } }




namespace stk {
namespace mesh {
namespace fixtures {

/**
 * Creates a ring mesh (circular loop of elements and nodes). Note that we create
 * a part for each locally owned element. This fixture is 1d, so elements are just lines.
 */

class RingFixture {
 public:
  const int             m_spatial_dimension;
  MetaData              m_meta_data;
  stk::mesh::unit_test::BulkDataTester m_bulk_data;
  PartVector            m_element_parts ;
  Part &                m_element_part_extra ;
  const size_t          m_num_element_per_proc ;
  std::vector<EntityId> m_node_ids , m_element_ids ;
  Part &                m_beam_2_part;

  RingFixture( stk::ParallelMachine pm ,
               unsigned num_element_per_proc = 10 ,
               bool use_element_parts = false,
               enum stk::mesh::BulkData::AutomaticAuraOption auto_aura_option = stk::mesh::BulkData::AUTO_AURA);

  ~RingFixture() {}

  /**
   * Generate a simple loop of mesh entities:
   * node[i] : element[i] : node[ ( i + 1 ) % node.size() ]
   */
  void generate_mesh();

  /**
   * Make sure that element->owner_rank() == element->node[1]->owner_rank()
   */
  void fixup_node_ownership(BulkData::modification_optimization mod_optimize = BulkData::MOD_END_SORT);

 protected:

  typedef std::multimap<EntityId, int> NodeToProcsMMap;
  NodeToProcsMMap m_nodes_to_procs;

 private:

   RingFixture();
   RingFixture( const RingFixture & );
   RingFixture & operator = ( const RingFixture & );

   void fill_node_map(int proc_rank);
};

}
}
}

#endif
