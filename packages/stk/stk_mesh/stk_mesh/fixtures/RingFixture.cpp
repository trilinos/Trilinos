/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/fixtures/RingFixture.hpp>
#include <ostream>                      // for ostringstream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityProc, EntityId, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include "stk_topology/topology.hpp"    // for topology, etc





namespace stk {
namespace mesh {
namespace fixtures {

/**
 * Creates a ring mesh (circular loop of elements and nodes). Note that we create
 * a part for each locally owned element.
 */

RingFixture::RingFixture( stk::ParallelMachine pm ,
                          unsigned num_element_per_proc ,
                          bool use_element_parts )
  : m_spatial_dimension(2),
    m_meta_data( m_spatial_dimension ),
    m_bulk_data( m_meta_data, pm, 100 ),
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
      Entity e_node_0 = m_bulk_data.declare_entity( stk::topology::NODE_RANK , m_node_ids[n0] , no_parts );
      Entity e_node_1 = m_bulk_data.declare_entity( stk::topology::NODE_RANK , m_node_ids[n1] , no_parts );
      Entity e_element   = m_bulk_data.declare_entity( stk::topology::ELEMENT_RANK , m_element_ids[i] , add_parts );
      m_bulk_data.declare_relation( e_element , e_node_0 , 0 );
      m_bulk_data.declare_relation( e_element , e_node_1 , 1 );
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
    m_bulk_data.change_entity_owner( change );
  }
}

}
}
}
