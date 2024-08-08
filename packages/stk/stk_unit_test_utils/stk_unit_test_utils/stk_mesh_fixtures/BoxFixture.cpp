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

#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Types.hpp>      // for EntityId, PartVector, etc
#include <stk_unit_test_utils/stk_mesh_fixtures/BoxFixture.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {

BoxFixture::BoxFixture( stk::ParallelMachine pm ,
                        stk::mesh::BulkData::AutomaticAuraOption autoAuraOption,
                        unsigned bucketCapacity,
                        const std::vector<std::string>& entity_names )
  : m_fem_meta ( spatial_dimension, entity_names ),
    m_bulk_data ( (bucketCapacity != stk::mesh::get_default_bucket_capacity()) ?
                    stk::unit_test_util::BulkDataTester(m_fem_meta, pm, autoAuraOption, false,
                                                        std::unique_ptr<stk::mesh::FieldDataManager>(),
                                                        bucketCapacity, bucketCapacity)
                  : stk::unit_test_util::BulkDataTester(m_fem_meta, pm, autoAuraOption, false) ),
    m_comm_rank( stk::parallel_machine_rank( pm ) ),
    m_comm_size( stk::parallel_machine_size( pm ) ),
    m_elem_part( m_fem_meta.declare_part_with_topology("elem_part", stk::topology::HEX_8) ),
    m_elem_topology( stk::topology::HEX_8 )
{
}

void BoxFixture::generate_boxes( const BOX   root_box,
                                       BOX   local_box )
{
  const int p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();

  m_nodes_to_procs.clear();

  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank, root_box);
  }

  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;

  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  local_box[0][0] = p_box[ p_rank ][0][0] ;
  local_box[0][1] = p_box[ p_rank ][0][1] ;
  local_box[1][0] = p_box[ p_rank ][1][0] ;
  local_box[1][1] = p_box[ p_rank ][1][1] ;
  local_box[2][0] = p_box[ p_rank ][2][0] ;
  local_box[2][1] = p_box[ p_rank ][2][1] ;

  // Create elements:

  std::vector<unsigned> local_count ;

  const stk::mesh::PartVector no_parts ;
  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&m_elem_part);

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    EntityId node_ids[8];
    node_ids[0] = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[1] = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[2] = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[3] = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[4]= 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[5] = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[6]= 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[7]= 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    const EntityId elem_id =  1 + i + j * ngx + k * ngx * ngy;
    Entity elem  = m_bulk_data.declare_element(elem_id , elem_parts);

    Entity nodes[8];
    for (int en_i = 0; en_i < 8; ++en_i) {
      nodes[en_i] = m_bulk_data.declare_node(node_ids[en_i] , no_parts);
      m_bulk_data.declare_relation(elem, nodes[en_i], en_i);
      DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, node_ids[en_i], nodes[en_i]);
    }
  }
  }
  }

  delete[] p_box ;
}


void BoxFixture::fill_node_map(int proc_rank, const BOX root_box)
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;

  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  BOX local_box;
  local_box[0][0] = p_box[ proc_rank ][0][0] ;
  local_box[0][1] = p_box[ proc_rank ][0][1] ;
  local_box[1][0] = p_box[ proc_rank ][1][0] ;
  local_box[1][1] = p_box[ proc_rank ][1][1] ;
  local_box[2][0] = p_box[ proc_rank ][2][0] ;
  local_box[2][1] = p_box[ proc_rank ][2][1] ;

  const stk::mesh::PartVector no_parts ;
  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&m_elem_part);

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    EntityId node_ids[8];
    node_ids[0] = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[1] = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[2] = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[3] = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[4]= 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[5] = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[6]= 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[7]= 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    for (int en_i = 0; en_i < 8; ++en_i) {
      AddToNodeProcsMMap(m_nodes_to_procs, node_ids[en_i], proc_rank);
     }
  }
  }
  }

  delete[] p_box ;
}

void BoxFixture::box_partition( int ip , int up , int axis ,
                                const BOX box ,
                                BOX p_box[] )
{
  const int np = up - ip ;
  if ( 1 == np ) {
    p_box[ip][0][0] = box[0][0] ; p_box[ip][0][1] = box[0][1] ;
    p_box[ip][1][0] = box[1][0] ; p_box[ip][1][1] = box[1][1] ;
    p_box[ip][2][0] = box[2][0] ; p_box[ip][2][1] = box[2][1] ;
  }
  else {
    const int n = box[ axis ][1] - box[ axis ][0] ;
    const int np_low = np / 2 ;  /* Rounded down */
    const int np_upp = np - np_low ;

    const int n_upp =
        static_cast<int>(
            static_cast<double>(n) * ( static_cast<double>(np_upp) / static_cast<double>(np) )
            );
    const int n_low = n - n_upp ;
    const int next_axis = ( axis + 2 ) % 3 ;

    if ( np_low ) { /* P = [ip,ip+np_low) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      dbox[ axis ][1] = dbox[ axis ][0] + n_low ;

      box_partition( ip, ip + np_low, next_axis,
                     static_cast<const int (*)[2]>(dbox), p_box );
    }

    if ( np_upp ) { /* P = [ip+np_low,ip+np_low+np_upp) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      ip += np_low ;
      dbox[ axis ][0] += n_low ;
      dbox[ axis ][1]  = dbox[ axis ][0] + n_upp ;

      box_partition( ip, ip + np_upp, next_axis,
                     static_cast<const int (*)[2]>(dbox), p_box );
    }
  }
}

namespace simple_fields {

BoxFixture::BoxFixture( stk::ParallelMachine pm ,
                        stk::mesh::BulkData::AutomaticAuraOption autoAuraOption,
                        unsigned bucketCapacity,
                        const std::vector<std::string>& entity_names )
  : m_fem_meta ( spatial_dimension, entity_names ),
    m_bulk_data ( (bucketCapacity != stk::mesh::get_default_bucket_capacity()) ?
                    stk::unit_test_util::BulkDataTester(m_fem_meta, pm, autoAuraOption, false,
                                                        std::unique_ptr<stk::mesh::FieldDataManager>(),
                                                        bucketCapacity, bucketCapacity)
                  : stk::unit_test_util::BulkDataTester(m_fem_meta, pm, autoAuraOption, false) ),
    m_comm_rank( stk::parallel_machine_rank( pm ) ),
    m_comm_size( stk::parallel_machine_size( pm ) ),
    m_elem_part( m_fem_meta.declare_part_with_topology("elem_part", stk::topology::HEX_8) ),
    m_elem_topology( stk::topology::HEX_8 )
{
}

void BoxFixture::generate_boxes( const BOX   root_box,
                                       BOX   local_box )
{
  const int p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();

  m_nodes_to_procs.clear();

  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank, root_box);
  }

  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;

  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  local_box[0][0] = p_box[ p_rank ][0][0] ;
  local_box[0][1] = p_box[ p_rank ][0][1] ;
  local_box[1][0] = p_box[ p_rank ][1][0] ;
  local_box[1][1] = p_box[ p_rank ][1][1] ;
  local_box[2][0] = p_box[ p_rank ][2][0] ;
  local_box[2][1] = p_box[ p_rank ][2][1] ;

  // Create elements:

  std::vector<unsigned> local_count ;

  const stk::mesh::PartVector no_parts ;
  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&m_elem_part);

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    EntityId node_ids[8];
    node_ids[0] = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[1] = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[2] = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[3] = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[4]= 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[5] = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[6]= 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[7]= 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    const EntityId elem_id =  1 + i + j * ngx + k * ngx * ngy;
    Entity elem  = m_bulk_data.declare_element(elem_id , elem_parts);

    Entity nodes[8];
    for (int en_i = 0; en_i < 8; ++en_i) {
      nodes[en_i] = m_bulk_data.declare_node(node_ids[en_i] , no_parts);
      m_bulk_data.declare_relation(elem, nodes[en_i], en_i);
      stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, node_ids[en_i], nodes[en_i]);
    }
  }
  }
  }

  delete[] p_box ;
}


void BoxFixture::fill_node_map(int proc_rank, const BOX root_box)
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;

  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  BOX local_box;
  local_box[0][0] = p_box[ proc_rank ][0][0] ;
  local_box[0][1] = p_box[ proc_rank ][0][1] ;
  local_box[1][0] = p_box[ proc_rank ][1][0] ;
  local_box[1][1] = p_box[ proc_rank ][1][1] ;
  local_box[2][0] = p_box[ proc_rank ][2][0] ;
  local_box[2][1] = p_box[ proc_rank ][2][1] ;

  const stk::mesh::PartVector no_parts ;
  stk::mesh::PartVector elem_parts;
  elem_parts.push_back(&m_elem_part);

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    EntityId node_ids[8];
    node_ids[0] = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[1] = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[2] = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[3] = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    node_ids[4]= 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[5] = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[6]= 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    node_ids[7]= 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    for (int en_i = 0; en_i < 8; ++en_i) {
      stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, node_ids[en_i], proc_rank);
     }
  }
  }
  }

  delete[] p_box ;
}

void BoxFixture::box_partition( int ip , int up , int axis ,
                                const BOX box ,
                                BOX p_box[] )
{
  const int np = up - ip ;
  if ( 1 == np ) {
    p_box[ip][0][0] = box[0][0] ; p_box[ip][0][1] = box[0][1] ;
    p_box[ip][1][0] = box[1][0] ; p_box[ip][1][1] = box[1][1] ;
    p_box[ip][2][0] = box[2][0] ; p_box[ip][2][1] = box[2][1] ;
  }
  else {
    const int n = box[ axis ][1] - box[ axis ][0] ;
    const int np_low = np / 2 ;  /* Rounded down */
    const int np_upp = np - np_low ;

    const int n_upp =
        static_cast<int>(
            static_cast<double>(n) * ( static_cast<double>(np_upp) / static_cast<double>(np) )
            );
    const int n_low = n - n_upp ;
    const int next_axis = ( axis + 2 ) % 3 ;

    if ( np_low ) { /* P = [ip,ip+np_low) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      dbox[ axis ][1] = dbox[ axis ][0] + n_low ;

      box_partition( ip, ip + np_low, next_axis,
                     static_cast<const int (*)[2]>(dbox), p_box );
    }

    if ( np_upp ) { /* P = [ip+np_low,ip+np_low+np_upp) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;

      ip += np_low ;
      dbox[ axis ][0] += n_low ;
      dbox[ axis ][1]  = dbox[ axis ][0] + n_upp ;

      box_partition( ip, ip + np_upp, next_axis,
                     static_cast<const int (*)[2]>(dbox), p_box );
    }
  }
}

}  // namespace simple_fields

}
}
}
