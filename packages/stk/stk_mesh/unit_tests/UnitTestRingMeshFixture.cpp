#include <unit_tests/UnitTestRingMeshFixture.hpp>

#include <unit_tests/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

using namespace stk;
using namespace stk::mesh;

RingMeshFixture::RingMeshFixture(ParallelMachine pm) : m_node_ids(),
                                                       m_edge_ids(),
                                                       m_meta_data(NULL),
                                                       m_bulk_data(NULL)
{
  enum { nPerProc = 10 };

  m_meta_data = new MetaData( fem_entity_type_names() );
  m_meta_data->commit();

  const size_t num_entities_per_bucket = 100;
  m_bulk_data = new BulkData( *m_meta_data, pm, num_entities_per_bucket );

  PartVector no_parts;

  generate_loop(no_parts , true /* aura */, nPerProc, m_node_ids, m_edge_ids );
}

RingMeshFixture::~RingMeshFixture()
{
  delete m_bulk_data;
  delete m_meta_data;
}

void RingMeshFixture::generate_loop(
  const PartVector      & edge_parts ,
  const bool              generate_aura ,
  const unsigned          nPerProc ,
  std::vector<EntityId> & node_ids ,
  std::vector<EntityId> & edge_ids )
{
  const unsigned p_rank = m_bulk_data->parallel_rank();
  const unsigned p_size = m_bulk_data->parallel_size();
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;
  const unsigned n_extra = generate_aura && 1 < p_size ? 2 : 0 ;

  node_ids.resize( id_total );
  edge_ids.resize( id_total );
  std::vector<unsigned> local_count ;

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    node_ids[i] = i + 1;
    edge_ids[i] = i + 1;
  }

  // Create a loop of edges:
  {
    const PartVector no_parts ;
    PartVector add_parts ;

    if ( ! edge_parts.empty() ) { add_parts.resize(1); }

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      const unsigned n0 = i ;
      const unsigned n1 = ( i + 1 ) % id_total ;
      if ( ! edge_parts.empty() ) {
        add_parts[0] = edge_parts[ i % edge_parts.size() ];
      }
      Entity & e_node_0 = m_bulk_data->declare_entity( 0 , node_ids[n0] , no_parts );
      Entity & e_node_1 = m_bulk_data->declare_entity( 0 , node_ids[n1] , no_parts );
      Entity & e_edge   = m_bulk_data->declare_entity( 1 , edge_ids[i] , add_parts );
      m_bulk_data->declare_relation( e_edge , e_node_0 , 0 );
      m_bulk_data->declare_relation( e_edge , e_node_1 , 1 );
    }
  }

  Selector select_owned( m_bulk_data->mesh_meta_data().locally_owned_part() );
  Selector select_used( m_bulk_data->mesh_meta_data().locally_used_part() );
  Selector select_all(  m_bulk_data->mesh_meta_data().universal_part() );

  count_entities( select_used , *m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[Node] == nLocalNode );
  STKUNIT_ASSERT( local_count[Edge] == nLocalEdge );

  std::vector<Entity*> all_nodes;
  get_entities(*m_bulk_data, Node, all_nodes);

  unsigned num_selected_nodes =
      count_selected_entities( select_used, m_bulk_data->buckets(Node) );
  STKUNIT_ASSERT( num_selected_nodes == local_count[Node] );

  std::vector<Entity*> universal_nodes;
  get_selected_entities( select_all, m_bulk_data->buckets(Node), universal_nodes );
  STKUNIT_ASSERT( universal_nodes.size() == all_nodes.size() );

  m_bulk_data->modification_end();

  // Verify declarations and sharing two end nodes:

  count_entities( select_used , *m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  if ( 1 < p_size ) {
    const std::vector<EntityProc> & shared = m_bulk_data->shared_entities();

    STKUNIT_ASSERT_EQUAL( shared.size() , size_t(2) );

    const unsigned n0 = id_end < id_total ? id_begin : 0 ;
    const unsigned n1 = id_end < id_total ? id_end : id_begin ;

    STKUNIT_ASSERT( shared[0].first->identifier() == node_ids[n0] );
    STKUNIT_ASSERT( shared[1].first->identifier() == node_ids[n1] );
    STKUNIT_ASSERT_EQUAL( shared[0].first->sharing().size() , size_t(1) );
    STKUNIT_ASSERT_EQUAL( shared[1].first->sharing().size() , size_t(1) );
  }

  // Test no-op first:

  std::vector<EntityProc> change ;

  STKUNIT_ASSERT( m_bulk_data->modification_begin() );
  m_bulk_data->change_entity_owner( change );
  m_bulk_data->modification_end();

  count_entities( select_used , *m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  count_entities( select_all , *m_bulk_data , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

  // Make sure that edge->owner_rank() == edge->node[1]->owner_rank() 
  if ( 1 < p_size ) {
    Entity * const e_node_0 = m_bulk_data->get_entity( 0 , node_ids[id_begin] );
    if ( p_rank == e_node_0->owner_rank() ) {
      EntityProc entry ;
      entry.first = e_node_0 ;
      entry.second = ( p_rank + p_size - 1 ) % p_size ;
      change.push_back( entry );
    }
    STKUNIT_ASSERT( m_bulk_data->modification_begin() );
    m_bulk_data->change_entity_owner( change );
    m_bulk_data->modification_end();

    count_entities( select_all , *m_bulk_data , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

    count_entities( select_used , *m_bulk_data , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    count_entities( select_owned , *m_bulk_data , local_count );
    STKUNIT_ASSERT( local_count[0] == nPerProc );
    STKUNIT_ASSERT( local_count[1] == nPerProc );
  }
}
