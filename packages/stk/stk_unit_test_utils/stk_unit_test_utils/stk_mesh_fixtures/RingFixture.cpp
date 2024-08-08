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

#include "mpi.h"                        // for ompi_communicator_t
#include "stk_topology/topology.hpp"    // for topology, etc
#include <ostream>                      // for ostringstream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityProc, PartVector, etc
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/RingFixture.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine


namespace stk {
namespace mesh {
namespace fixtures {

/**
 * Creates a ring mesh (circular loop of elements and nodes). Note that we create
 * a part for each locally owned element.
 */

RingFixture::RingFixture( stk::ParallelMachine pm ,
                          unsigned num_element_per_proc ,
                          bool use_element_parts,
                          enum stk::mesh::BulkData::AutomaticAuraOption auto_aura_option)
  : m_spatial_dimension(2),
    m_meta_data( m_spatial_dimension ),
    m_bulk_data( m_meta_data, pm, auto_aura_option ),
    m_element_parts(),
    m_element_part_extra( m_meta_data.declare_part("element_extra" , stk::topology::ELEMENT_RANK ) ),
    m_num_element_per_proc( num_element_per_proc ),
    m_node_ids(),
    m_element_ids(),
    m_beam_2_part( m_meta_data.declare_part_with_topology("beam_2_part", stk::topology::BEAM_2 ) )
{

  if ( use_element_parts ) {
    m_element_parts.resize( num_element_per_proc );
    for ( unsigned i = 0 ; i < num_element_per_proc ; ++i ) {
      std::ostringstream name ;
      name << "ElementPart_" << i ;
      m_element_parts[i] = & m_meta_data.declare_part( name.str() , stk::topology::ELEMENT_RANK );
      //m_element_parts[i] = & m_meta_data.declare_part_with_topology(name.str(), stk::topology::BEAM_2);
    }
  }
}

void RingFixture::generate_mesh( )
{
  const int p_rank          = m_bulk_data.parallel_rank();
  const int p_size          = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_element_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );

  m_node_ids.resize( id_total );
  m_element_ids.resize( id_total );
  std::vector<unsigned> local_count ;

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    m_node_ids[i] = i + 1;
    m_element_ids[i] = i + 1;
  }

  m_nodes_to_procs.clear();

  for (int proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  // Create a loop of elements:
  {
    const PartVector no_parts ;
    PartVector add_parts ;
    add_parts.push_back(&m_beam_2_part);

    if ( ! m_element_parts.empty() ) { add_parts.resize(2); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! m_element_parts.empty() ) {
        add_parts[1] = m_element_parts[ i % m_element_parts.size() ];
      }

      Entity e_node_0 = m_bulk_data.declare_node(m_node_ids[n0], no_parts);
      Entity e_node_1 = m_bulk_data.declare_node(m_node_ids[n1], no_parts);
      Entity e_element   = m_bulk_data.declare_element(m_element_ids[i], add_parts);
      m_bulk_data.declare_relation( e_element , e_node_0 , 0 );
      m_bulk_data.declare_relation( e_element , e_node_1 , 1 );
      DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, m_node_ids[n0], e_node_0);
      DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, m_node_ids[n1], e_node_1);
    }
  }
}

void RingFixture::fill_node_map(int p_rank)
{
  //const int p_rank          = m_bulk_data.parallel_rank();
  const int p_size          = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_element_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );

  // Create a loop of elements:
  {
    const PartVector no_parts ;
    PartVector add_parts ;
    add_parts.push_back(&m_beam_2_part);

    if ( ! m_element_parts.empty() ) { add_parts.resize(2); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! m_element_parts.empty() ) {
        add_parts[1] = m_element_parts[ i % m_element_parts.size() ];
      }

      AddToNodeProcsMMap(m_nodes_to_procs, m_node_ids[n0] , p_rank);
      AddToNodeProcsMMap(m_nodes_to_procs, m_node_ids[n1] , p_rank);
    }
  }
}

void RingFixture::fixup_node_ownership()
{
  const int p_rank          = m_bulk_data.parallel_rank();
  const int p_size          = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_element_per_proc ;
  const unsigned id_begin   = nPerProc * p_rank ;

  // Make sure that element->owner_rank() == element->node[1]->owner_rank()

  if ( 1 < p_size ) {
    std::vector<EntityProc> change ;
    Entity const e_node_0 = m_bulk_data.get_entity( stk::topology::NODE_RANK , m_node_ids[id_begin] );
    if ( p_rank == m_bulk_data.parallel_owner_rank(e_node_0) ) {
      EntityProc entry ;
      entry.first = e_node_0 ;
      entry.second = ( p_rank + p_size - 1 ) % p_size ;
      change.push_back( entry );
    }
    m_bulk_data.change_entity_owner(change);
  }
}

namespace simple_fields {

RingFixture::RingFixture( stk::ParallelMachine pm ,
                          unsigned num_element_per_proc ,
                          bool use_element_parts,
                          enum stk::mesh::BulkData::AutomaticAuraOption auto_aura_option)
  : m_spatial_dimension(2),
    m_meta_data( m_spatial_dimension ),
    m_bulk_data( m_meta_data, pm, auto_aura_option ),
    m_element_parts(),
    m_element_part_extra( m_meta_data.declare_part("element_extra" , stk::topology::ELEMENT_RANK ) ),
    m_num_element_per_proc( num_element_per_proc ),
    m_node_ids(),
    m_element_ids(),
    m_beam_2_part( m_meta_data.declare_part_with_topology("beam_2_part", stk::topology::BEAM_2 ) )
{

  if ( use_element_parts ) {
    m_element_parts.resize( num_element_per_proc );
    for ( unsigned i = 0 ; i < num_element_per_proc ; ++i ) {
      std::ostringstream name ;
      name << "ElementPart_" << i ;
      m_element_parts[i] = & m_meta_data.declare_part( name.str() , stk::topology::ELEMENT_RANK );
      //m_element_parts[i] = & m_meta_data.declare_part_with_topology(name.str(), stk::topology::BEAM_2);
    }
  }
}

void RingFixture::generate_mesh( )
{
  const int p_rank          = m_bulk_data.parallel_rank();
  const int p_size          = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_element_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );

  m_node_ids.resize( id_total );
  m_element_ids.resize( id_total );
  std::vector<unsigned> local_count ;

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    m_node_ids[i] = i + 1;
    m_element_ids[i] = i + 1;
  }

  m_nodes_to_procs.clear();

  for (int proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  // Create a loop of elements:
  {
    const PartVector no_parts ;
    PartVector add_parts ;
    add_parts.push_back(&m_beam_2_part);

    if ( ! m_element_parts.empty() ) { add_parts.resize(2); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! m_element_parts.empty() ) {
        add_parts[1] = m_element_parts[ i % m_element_parts.size() ];
      }

      Entity e_node_0 = m_bulk_data.declare_node(m_node_ids[n0], no_parts);
      Entity e_node_1 = m_bulk_data.declare_node(m_node_ids[n1], no_parts);
      Entity e_element   = m_bulk_data.declare_element(m_element_ids[i], add_parts);
      m_bulk_data.declare_relation( e_element , e_node_0 , 0 );
      m_bulk_data.declare_relation( e_element , e_node_1 , 1 );
      stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, m_node_ids[n0], e_node_0);
      stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, m_node_ids[n1], e_node_1);
    }
  }
}

void RingFixture::fill_node_map(int p_rank)
{
  //const int p_rank          = m_bulk_data.parallel_rank();
  const int p_size          = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_element_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );

  // Create a loop of elements:
  {
    const PartVector no_parts ;
    PartVector add_parts ;
    add_parts.push_back(&m_beam_2_part);

    if ( ! m_element_parts.empty() ) { add_parts.resize(2); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! m_element_parts.empty() ) {
        add_parts[1] = m_element_parts[ i % m_element_parts.size() ];
      }

      stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, m_node_ids[n0] , p_rank);
      stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, m_node_ids[n1] , p_rank);
    }
  }
}

void RingFixture::fixup_node_ownership()
{
  const int p_rank          = m_bulk_data.parallel_rank();
  const int p_size          = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_element_per_proc ;
  const unsigned id_begin   = nPerProc * p_rank ;

  // Make sure that element->owner_rank() == element->node[1]->owner_rank()

  if ( 1 < p_size ) {
    std::vector<EntityProc> change ;
    Entity const e_node_0 = m_bulk_data.get_entity( stk::topology::NODE_RANK , m_node_ids[id_begin] );
    if ( p_rank == m_bulk_data.parallel_owner_rank(e_node_0) ) {
      EntityProc entry ;
      entry.first = e_node_0 ;
      entry.second = ( p_rank + p_size - 1 ) % p_size ;
      change.push_back( entry );
    }
    m_bulk_data.change_entity_owner(change);
  }
}

} // namespace simple_fields

}
}
}
