/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/fixtures/RingFixture.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/GetEntities.hpp>


#include <Shards_BasicTopologies.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

/**
 * Creates a ring mesh (circular loop of edges and nodes). Note that we create
 * a part for each locally owned edge. Also note that, since we are in 1D, 
 * the "edges" are actually elements.
 */

RingFixture::RingFixture( stk::ParallelMachine pm ,
                          unsigned num_edge_per_proc ,
                          bool use_edge_parts )
  : m_spatial_dimension(1),
    m_meta_data( fem::entity_rank_names(m_spatial_dimension) ),
    m_bulk_data( m_meta_data, pm, 100 ),
    m_fem( m_meta_data, m_spatial_dimension ),
    m_edge_parts(),
    m_edge_part_extra( declare_part(m_meta_data, "edge_extra" , fem::element_rank(m_fem) ) ),
    m_num_edge_per_proc( num_edge_per_proc ),
    m_node_ids(),
    m_edge_ids()
{
  if ( use_edge_parts ) {
    m_edge_parts.resize( num_edge_per_proc );
    for ( unsigned i = 0 ; i < num_edge_per_proc ; ++i ) {
      std::ostringstream name ;
      name << "EdgePart_" << i ;
      m_edge_parts[i] = & declare_part(m_meta_data,  name.str() , fem::element_rank(m_fem) );
    }
  }
}

void RingFixture::generate_mesh( )
{
  const unsigned p_rank     = m_bulk_data.parallel_rank();
  const unsigned p_size     = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_edge_per_proc ;
  const unsigned id_total   = nPerProc * p_size ;
  const unsigned id_begin   = nPerProc * p_rank ;
  const unsigned id_end     = nPerProc * ( p_rank + 1 );

  m_node_ids.resize( id_total );
  m_edge_ids.resize( id_total );
  std::vector<unsigned> local_count ;

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    m_node_ids[i] = i + 1;
    m_edge_ids[i] = i + 1;
  }

  // Create a loop of edges:
  {
    const PartVector no_parts ;
    PartVector add_parts ;

    if ( ! m_edge_parts.empty() ) { add_parts.resize(1); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! m_edge_parts.empty() ) {
        add_parts[0] = m_edge_parts[ i % m_edge_parts.size() ];
      }
      Entity & e_node_0 = m_bulk_data.declare_entity( 0 , m_node_ids[n0] , no_parts );
      Entity & e_node_1 = m_bulk_data.declare_entity( 0 , m_node_ids[n1] , no_parts );
      Entity & e_edge   = m_bulk_data.declare_entity( 1 , m_edge_ids[i] , add_parts );
      m_bulk_data.declare_relation( e_edge , e_node_0 , 0 );
      m_bulk_data.declare_relation( e_edge , e_node_1 , 1 );
    }
  }
}

void RingFixture::fixup_node_ownership()
{
  const unsigned p_rank     = m_bulk_data.parallel_rank();
  const unsigned p_size     = m_bulk_data.parallel_size();
  const unsigned nPerProc   = m_num_edge_per_proc ;
  const unsigned id_begin   = nPerProc * p_rank ;

  // Make sure that edge->owner_rank() == edge->node[1]->owner_rank()

  if ( 1 < p_size ) {
    std::vector<EntityProc> change ;
    Entity * const e_node_0 = m_bulk_data.get_entity( 0 , m_node_ids[id_begin] );
    if ( p_rank == e_node_0->owner_rank() ) {
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
