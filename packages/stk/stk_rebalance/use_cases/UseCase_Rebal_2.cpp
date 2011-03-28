/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Rebal_2.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

#include <stk_rebalance_utils/RebalanceUtils.hpp>

//----------------------------------------------------------------------

using namespace stk::mesh::fixtures;

typedef stk::mesh::Field<double> ScalarField ;

namespace stk {
namespace rebalance {
namespace use_cases {

namespace {

  void sum_element_weights_through_relations( stk::mesh::EntityVector & elements, 
      ScalarField & field, const std::vector<stk::mesh::EntityRank> & ranks )
  {
    for( size_t i = 0; i < ranks.size(); ++i )
    {
      for( size_t ielem = 0; ielem < elements.size(); ++ielem )
      {
        stk::mesh::Entity * elem = elements[ielem];
        double * elem_weight = field_data( field , *elem );
        const stk::mesh::PairIterRelation rel = elem->relations( ranks[i] );
        const unsigned num_entities = rel.size();

        for ( unsigned j = 0 ; j < num_entities ; ++j ) 
        {
          stk::mesh::Entity & entity = * rel[j].entity();
          const double * entity_weight = field_data( field , entity );
          elem_weight[0] += entity_weight[0];
        }
      }
    }
  }

} // namespace

bool test_heavy_nodes( stk::ParallelMachine pm )
{
  const size_t nx = 3;
  const size_t ny = 3;
  const size_t nz = 3;

  const unsigned p_size = stk::parallel_machine_size(pm);
  const unsigned p_rank = stk::parallel_machine_rank(pm);

  stk::mesh::fixtures::HexFixture fixture(pm, nx, ny, nz);

  stk::mesh::fem::FEMMetaData & fem_meta  = fixture.m_fem_meta;
  stk::mesh::BulkData & bulk  = fixture.m_bulk_data;

  // Put weights field on all entities
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();
  const stk::mesh::EntityRank face_rank = fem_meta.face_rank();
  const stk::mesh::EntityRank edge_rank = fem_meta.edge_rank();
  const stk::mesh::EntityRank node_rank = fem_meta.node_rank();
  ScalarField & weight_field( fem_meta.declare_field< ScalarField >( "entity_weights" ) );
  stk::mesh::put_field(weight_field , element_rank , fem_meta.universal_part() );
  stk::mesh::put_field(weight_field , face_rank , fem_meta.universal_part() );
  stk::mesh::put_field(weight_field , edge_rank , fem_meta.universal_part() );
  stk::mesh::put_field(weight_field , node_rank , fem_meta.universal_part() );

  fem_meta.commit();

  // Configure our mesh on proc 0
  std::vector<stk::mesh::EntityId> my_element_ids;
  if( 0 == p_rank )
  {
    for ( unsigned i = 0 ; i < nx*ny*nz; ++i )
      my_element_ids.push_back(i+1);
  }

  fixture.generate_mesh(my_element_ids);

  // create faces and edges for the mesh
  stk::mesh::PartVector add_parts;
  stk::mesh::create_adjacent_entities(bulk, add_parts);

  bulk.modification_begin();

  // Assign entity weights
  if( 0 == p_rank ) 
  {
    // Get the faces on the x=0 plane and give them a characteristic weight
    stk::mesh::EntityVector selected_nodes;
    stk::mesh::EntityVector selected_faces;
    stk::mesh::EntityVector one_face;
    for ( unsigned j = 0 ; j < ny; ++j ) 
      for ( unsigned k = 0 ; k < nz; ++k ) 
      {
        selected_nodes.clear();
        selected_nodes.push_back( fixture.node(0, j,   k  ) );
        selected_nodes.push_back( fixture.node(0, j+1, k  ) );
        selected_nodes.push_back( fixture.node(0, j,   k+1) );
        selected_nodes.push_back( fixture.node(0, j+1, k+1) );
        stk::mesh::get_entities_through_relations(selected_nodes, face_rank, one_face);
        selected_faces.push_back(one_face[0]);
      }

    for( size_t iface = 0; iface < selected_faces.size(); ++iface )
    {
      stk::mesh::Entity * face = selected_faces[iface];
      double * const weight = stk::mesh::field_data( weight_field, *face );
      weight[0] = 10.0;
    }

    // Get the edges on the boundary of the x=0 plane and give them a characteristic weight
    stk::mesh::EntityVector selected_edges;
    stk::mesh::EntityVector one_edge;
    for ( unsigned j = 0 ; j < ny; ++j ) 
    {
      selected_nodes.clear();
      selected_nodes.push_back( fixture.node(0, j,   0) );
      selected_nodes.push_back( fixture.node(0, j+1, 0) );
      stk::mesh::get_entities_through_relations(selected_nodes, edge_rank, one_edge);
      selected_edges.push_back(one_edge[0]);
      selected_nodes.clear();
      selected_nodes.push_back( fixture.node(0, j,   nz) );
      selected_nodes.push_back( fixture.node(0, j+1, nz) );
      stk::mesh::get_entities_through_relations(selected_nodes, edge_rank, one_edge);
      selected_edges.push_back(one_edge[0]);
    }
    for ( unsigned k = 0 ; k < nz; ++k ) 
    {
      selected_nodes.clear();
      selected_nodes.push_back( fixture.node(0, 0, k) );
      selected_nodes.push_back( fixture.node(0, 0, k+1) );
      stk::mesh::get_entities_through_relations(selected_nodes, edge_rank, one_edge);
      selected_edges.push_back(one_edge[0]);
      selected_nodes.clear();
      selected_nodes.push_back( fixture.node(0, ny, k) );
      selected_nodes.push_back( fixture.node(0, ny, k+1) );
      stk::mesh::get_entities_through_relations(selected_nodes, edge_rank, one_edge);
      selected_edges.push_back(one_edge[0]);
    }
    for( size_t iedge = 0; iedge < selected_edges.size(); ++iedge )
    {
      stk::mesh::Entity * edge = selected_edges[iedge];
      double * const weight = stk::mesh::field_data( weight_field, *edge );
      weight[0] = 100.0;
    }

    // Finally, give the corner nodes of the x=0 plane a characteristic weight
    selected_nodes.clear();
    double * weight = stk::mesh::field_data( weight_field, *fixture.node(0, 0, 0) );
    weight[0] = 1000.0;
    weight = stk::mesh::field_data( weight_field, *fixture.node(0, ny, 0) );
    weight[0] = 1000.0;
    weight = stk::mesh::field_data( weight_field, *fixture.node(0, 0, nz) );
    weight[0] = 1000.0;
    weight = stk::mesh::field_data( weight_field, *fixture.node(0, ny, nz) );
    weight[0] = 1000.0;

    // Assign element weights
    for( size_t i = 0; i < my_element_ids.size(); ++i )
    {
      stk::mesh::Entity * elem = bulk.get_entity(element_rank, my_element_ids[i]);
      double * const e_weight = stk::mesh::field_data( weight_field , *elem );
      *e_weight = 1.0;
    }
    //
    // Get the elements on the x=0 plane and sum in weights from relations
    selected_nodes.clear();
    for ( unsigned j = 0 ; j < ny+1; ++j )
      for ( unsigned k = 0 ; k < nz+1; ++k ) 
        selected_nodes.push_back( fixture.node(0, j, k) );

    std::vector<stk::mesh::EntityRank> ranks;
    ranks.push_back(face_rank);
    ranks.push_back(edge_rank);
    ranks.push_back(node_rank);
    stk::mesh::EntityVector selected_elems;
    for ( unsigned j = 0 ; j < ny; ++j ) 
      for ( unsigned k = 0 ; k < nz; ++k ) 
      {
        selected_nodes.clear();
        selected_nodes.push_back( fixture.node(0, j,   k  ) );
        selected_nodes.push_back( fixture.node(0, j+1, k  ) );
        selected_nodes.push_back( fixture.node(0, j,   k+1) );
        selected_nodes.push_back( fixture.node(0, j+1, k+1) );
        stk::mesh::get_entities_through_relations(selected_nodes, element_rank, one_face);
        selected_elems.push_back(one_face[0]);
      }
    sum_element_weights_through_relations(selected_elems, weight_field, ranks);
  }

  bulk.modification_end();

  // Use Zoltan to determine new partition
  Teuchos::ParameterList emptyList;
  stk::rebalance::Zoltan zoltan_partition(pm, fixture.m_spatial_dimension, emptyList);

  stk::rebalance::rebalance(bulk, fem_meta.universal_part(), &fixture.m_coord_field, &weight_field, zoltan_partition);

  const double imbalance_threshold = ( 3 == p_size )? 1.45 // Zoltan does not do so well for 3 procs
                                                    : 1.1; // ... but does pretty well for 2 and 4 procs
  const bool do_rebal = imbalance_threshold < stk::rebalance::check_balance(bulk, &weight_field, element_rank);

  if( 0 == p_rank )
    std::cerr << std::endl 
      << "imbalance_threshold after rebalance = " << imbalance_threshold << ", " << do_rebal << std::endl;


  // Check that we satisfy our threshhold
  bool result = !do_rebal;

  // And verify that all dependent entities are on the same proc as their parent element
  {
    stk::mesh::EntityVector entities;
    stk::mesh::Selector selector = fem_meta.locally_owned_part();

    get_selected_entities(selector, bulk.buckets(node_rank), entities);
    result &= verify_dependent_ownership(element_rank, entities);

    get_selected_entities(selector, bulk.buckets(edge_rank), entities);
    result &= verify_dependent_ownership(element_rank, entities);

    get_selected_entities(selector, bulk.buckets(face_rank), entities);
    result &= verify_dependent_ownership(element_rank, entities);
  }

  return result;
}

} //namespace use_cases
} //namespace rebalance
} //namespace stk


