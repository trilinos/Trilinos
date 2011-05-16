/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <common/gnu_malloc_hooks.hpp>

using namespace stk ;

//----------------------------------------------------------------------
// This file contains the implementation of use-case 13
// The function 'use_case_13_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk_use_cases {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

typedef mesh::Field<double>                    ScalarFieldType ;
typedef mesh::Field<double,mesh::Cartesian>    VectorFieldType ;

// Specification for the aggressive gather pointer-field for elements.

typedef mesh::Field<double*,mesh::ElementNode> ElementNodePointerFieldType ;

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_13_generate_mesh(
  mesh::BulkData & mesh ,
  const unsigned N[] ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  mesh::Part & hex_block );

void use_case_13_generate_sides(
  mesh::BulkData & mesh , const bool skin_only );

void use_case_13_algorithm( mesh::BulkData & ,
                                    const unsigned          side_type ,
                                    const VectorFieldType & side_field ,
                                    const VectorFieldType & elem_field );

//--------------------------------------------------------------------
//
// main driver for use-case 13: heterogeneous element mesh.
//

void use_case_13_driver( MPI_Comm comm )
{
  const unsigned p_rank = parallel_machine_rank( comm );

  reset_malloc_stats();

  const unsigned box_size[3] = { 20 , 20 , 20 };

  //--------------------------------------------------------------------

  if ( ! p_rank ) {
    std::cout << "stk_mesh Use Case #13, begin" << std::endl
              << "  Number Processes = " << parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  {
    //------------------------------------------------------------------
    // Declare the mesh meta data and bulk data.

    mesh::fem::FEMMetaData mesh_meta_data( SpatialDim );
    const stk::mesh::EntityRank element_rank = mesh_meta_data.element_rank();
    const stk::mesh::EntityRank side_rank    = mesh_meta_data.side_rank();
    const stk::mesh::EntityRank edge_rank    = mesh_meta_data.edge_rank();
    mesh::BulkData mesh_bulk_data( mesh_meta_data.get_meta_data(mesh_meta_data) , MPI_COMM_WORLD );

    //--------------------------------
    // Element-block declarations typically occur when reading the
    // mesh-file meta-data, and thus won't usually appear in application code.
    // Declaring the element blocks and associating an element traits
    // with each element block.

    mesh::Part & universal = mesh_meta_data.universal_part();
    mesh::Part & block_hex = mesh_meta_data.declare_part("block_1", element_rank);

    stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
    mesh::fem::set_cell_topology(block_hex, hex_top );

    //--------------------------------
    // Declare coordinates field on all nodes with 3D:

    VectorFieldType & coordinates_field =
      mesh_meta_data.declare_field< VectorFieldType >( "coordinates" );

    stk::mesh::put_field(
      coordinates_field , mesh::fem::FEMMetaData::NODE_RANK , universal , SpatialDim );

    //--------------------------------

    VectorFieldType & face_field =
      mesh_meta_data.declare_field< VectorFieldType >( "face_flux" );

    VectorFieldType & elem_field =
      mesh_meta_data.declare_field< VectorFieldType >( "elem_flux" );

    stk::mesh::put_field(elem_field , element_rank , block_hex , SpatialDim );

    stk::mesh::put_field(face_field , side_rank , universal , SpatialDim );

    //--------------------------------
    // Declare an aggressive "gather" field which is an
    // array of pointers to the element's nodes' coordinate field data.
    // The declaration specifies:
    //
    //     double * elem_node_coord[number_of_nodes]

    ElementNodePointerFieldType & elem_node_coord =
      mesh_meta_data.
        declare_field< ElementNodePointerFieldType >( "elem_node_coord" );

    // Declare that the 'elem_node_coord' pointer field data
    // points to the 'coordinates_field' data on the nodes.

    mesh_meta_data.declare_field_relation(
      elem_node_coord ,
      mesh::fem::get_element_node_stencil(SpatialDim) ,
      coordinates_field );

    // Declare the size of the aggressive "gather" field
    //     double * elem_node_coord[ size = number_of_nodes ]
    // is the number of nodes of the elements.
    // This size is different for each element block.

    stk::mesh::put_field(
        elem_node_coord , element_rank, block_hex , shards::Hexahedron<8> ::node_count );

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    mesh_meta_data.commit();

    //------------------------------------------------------------------
    //------------------------------------------------------------------

    // In a typical app, the mesh would be read from file at this point.
    // But in this use-case, we generate the mesh and initialize
    // field data to use-case defined values.

    use_case_13_generate_mesh(
      mesh_bulk_data ,
      box_size ,
      coordinates_field ,
      elem_node_coord ,
      block_hex );

    use_case_13_generate_sides( mesh_bulk_data , false );


    {
      std::vector<unsigned> count ;
      mesh::Selector selector = mesh_meta_data.locally_owned_part() |
                                mesh_meta_data.globally_shared_part() ;
      count_entities( selector, mesh_bulk_data, count );

      std::cout << "  P" << p_rank << ": Uses {"
                << " Node = " << count[ stk::mesh::fem::FEMMetaData::NODE_RANK ]
                << " Edge = " << count[ edge_rank ]
                << " Face = " << count[ side_rank ]
                << " Elem = " << count[ element_rank ]
                << " }" << std::endl ;
      std::cout.flush();
    }

    //------------------------------------------------------------------

#ifdef USE_GNU_MALLOC_HOOKS
    if (parallel_machine_rank(comm) == 0) {
      double net_alloc = alloc_MB() - freed_MB();
      std::cout << "Mesh creation:" << "\n   Total allocated: "
        << alloc_MB()<<"MB in "<<alloc_blks() << " blocks."
        << "\n   Total freed: " << freed_MB() << "MB in "
        << freed_blks() << " blocks."
        << "\n   Net allocated: "<<net_alloc << "MB."<<std::endl;
    }
#endif

    //------------------------------------------------------------------

    use_case_13_algorithm( mesh_bulk_data , side_rank,
                           face_field , elem_field );

    //------------------------------------------------------------------

  }
}

//--------------------------------------------------------------------

namespace {

/* Determine if the element's side is in the outward orientation. */
bool outward_orientation( const mesh::Entity & elem ,
                          const mesh::Entity & side ,
                          const unsigned side_ord )
{
  const CellTopologyData * const elem_top = mesh::fem::get_cell_topology( elem ).getCellTopologyData();

  const mesh::PairIterRelation elem_nodes = elem.relations( mesh::fem::FEMMetaData::NODE_RANK );
  const mesh::PairIterRelation side_nodes = side.relations( mesh::fem::FEMMetaData::NODE_RANK );

  // PairIterRelation is a typedef for:
  //
  //   PairIter< std::vector<Relation>::const_iterator >
  //
  // Where PairIter is a wrapper for a pair of iterators.

  const CellTopologyData * const side_top = elem_top->side[ side_ord ].topology ;
  const unsigned         * const side_map = elem_top->side[ side_ord ].node ;

  bool good = true ;
  for ( unsigned j = 0 ; good && j < side_top->node_count ; ++j ) {
    good = side_nodes[j].entity() == elem_nodes[ side_map[j] ].entity();
  }
  return good ;
}

}

void use_case_13_algorithm(
  mesh::BulkData & M ,
  const unsigned          side_type ,
  const VectorFieldType & side_field ,
  const VectorFieldType & elem_field )
{
  const mesh::fem::FEMMetaData & meta_data = mesh::fem::FEMMetaData::get(M);

  const stk::mesh::EntityRank element_rank = meta_data.element_rank();

  {
    // Communicate the element field data that we care about
    // from the processor owners of the ghosted elements
    // to   the ghosted copies of the elements.

    std::vector< const mesh::FieldBase * > sync_fields( 1 , & elem_field );

    mesh::communicate_field_data( M.shared_aura() , sync_fields );
  }

  // Get vector of buckets ( entities and field data)
  // for which the sides are all locally owned.

  mesh::Selector select_owned( meta_data.locally_owned_part() );

  const std::vector<mesh::Bucket*> & buckets = M.buckets( side_type );

  for ( std::vector<mesh::Bucket *>::const_iterator
        ik = buckets.begin() ; ik != buckets.end() ; ++ik ) if ( select_owned( **ik ) ) {

    const mesh::Bucket & bucket = **ik ;

    // Number of sides in this bucket of sides and side field data

    const int number = bucket.size();

    double * side_data = field_data( side_field , bucket.begin() );

    for ( int i = 0 ; i < number ; ++i , side_data += 3 ) {

      mesh::Entity & side = bucket[i] ;

      const mesh::PairIterRelation side_elems = side.relations(element_rank);

      if ( side_elems.size() == 2 ) {
        mesh::Entity & elem1 = * side_elems[0].entity();
        mesh::Entity & elem2 = * side_elems[1].entity();

        double * const elem1_data = field_data( elem_field , elem1 );
        double * const elem2_data = field_data( elem_field , elem2 );

        // Which way is the side oriented, natural for #1 or #2 ?

        if ( outward_orientation( elem1, side, side_elems[0].identifier()) ){
          side_data[0] = elem2_data[0] - elem1_data[0] ;
          side_data[1] = elem2_data[1] - elem1_data[1] ;
          side_data[2] = elem2_data[2] - elem1_data[2] ;
        }
        else {
          side_data[0] = elem1_data[0] - elem2_data[0] ;
          side_data[1] = elem1_data[1] - elem2_data[1] ;
          side_data[2] = elem1_data[2] - elem2_data[2] ;
        }
      }
      else {
        // Nothing happening with only one element?
        side_data[0] = 0 ;
        side_data[1] = 0 ;
        side_data[2] = 0 ;
      }
    }
  }
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// Does not work for shells since side neighbor can be ambiguous...

mesh::Entity * get_side_neighbor(
  const CellTopologyData & elem_top , const mesh::Entity & elem , unsigned side_id )
{
  const CellTopologyData & side_top = * elem_top.side[ side_id ].topology ;
  const unsigned * const node_map = elem_top.side[ side_id ].node ;
  const mesh::PairIterRelation elem_nodes = elem.relations( mesh::fem::FEMMetaData::NODE_RANK );

  // Find other element that shares this side...

  mesh::Entity & node = * elem_nodes[ node_map[0] ].entity();

  const mesh::PairIterRelation node_elems = node.relations(elem.entity_rank());

  mesh::Entity * neighbor = NULL ;

  for ( unsigned i = 0 ; neighbor == NULL && i < node_elems.size() ; ++i ) {

    neighbor = node_elems[i].entity();

    const mesh::PairIterRelation neighbor_nodes = neighbor->relations( mesh::fem::FEMMetaData::NODE_RANK );

    if ( & elem == neighbor ) { neighbor = NULL ; }

    for ( unsigned j = 1 ;
          neighbor != NULL && j < side_top.node_count ; ++j ) {

      mesh::Entity * const next_node = elem_nodes[ node_map[j] ].entity();

      // If neighbor does not have node then not this element ...

      bool found = false ;
      for ( unsigned k = 0 ; ! found && k < neighbor_nodes.size() ; ++k ) {
        found = next_node == neighbor_nodes[k].entity();
      }
      if ( ! found ) { neighbor = NULL ; }
    }

#if 0
    if ( NULL != neighbor ) {
      std::cout << "neighbors( " ;
      std::cout << " Element[ " ;
      std::cout << elem.identifier();
      std::cout << " ]{" ;
      for ( int i = 0 ; i < elem_nodes.size() ; ++i ) {
        std::cout << " " << elem_nodes[i].entity()->identifier();
      }
      std::cout << " } , Element[ " ;
      std::cout << neighbor->identifier();
      std::cout << " ]{" ;
      for ( int i = 0 ; i < neighbor_nodes.size() ; ++i ) {
        std::cout << " " << neighbor_nodes[i].entity()->identifier();
      }
      std::cout << " } , Share { " ;
      for ( unsigned j = 0 ; j < side_top.node_count ; ++j ) {
        mesh::Entity * const next_node = elem_nodes[ node_map[j] ].entity();
        std::cout << " " << next_node->identifier();
      }
      std::cout << " } )" ;
      std::cout << std::endl ;
      std::cout.flush();
    }
#endif
  }

  return neighbor ;
}

unsigned determine_local_side_id( const mesh::Entity & elem , mesh::Entity & side )
{
  const CellTopologyData * const elem_top = mesh::fem::get_cell_topology( elem ).getCellTopologyData();

  const mesh::PairIterRelation elem_nodes = elem.relations( mesh::fem::FEMMetaData::NODE_RANK );
  const mesh::PairIterRelation side_nodes = side.relations( mesh::fem::FEMMetaData::NODE_RANK );

  int side_id = -1 ;

  for ( unsigned i = 0 ; side_id == -1 && i < elem_top->side_count ; ++i ) {
    const CellTopologyData & side_top = * elem_top->side[i].topology ;
    const unsigned     * side_map =   elem_top->side[i].node ;

    if ( side_nodes.size() == side_top.node_count ) {

      side_id = i ;

      for ( unsigned j = 0 ;
            side_id == static_cast<int>(i) && j < side_top.node_count ; ++j ) {

        mesh::Entity * const elem_node = elem_nodes[ side_map[j] ].entity();

        bool found = false ;

        for ( unsigned k = 0 ; ! found && k < side_top.node_count ; ++k ) {
          found = elem_node == side_nodes[k].entity();
        }

        if ( ! found ) { side_id = -1 ; }
      }
    }
  }

  if ( side_id < 0 ) {
    std::ostringstream msg ;
    msg << "determine_local_side_id( " ;
    msg << elem_top->name ;
    msg << " , Element[ " ;
    msg << elem.identifier();
    msg << " ]{" ;
    for ( unsigned i = 0 ; i < elem_nodes.size() ; ++i ) {
      msg << " " << elem_nodes[i].entity()->identifier();
    }
    msg << " } , Side[ " ;
    msg << side.identifier();
    msg << " ]{" ;
    for ( unsigned i = 0 ; i < side_nodes.size() ; ++i ) {
      msg << " " << side_nodes[i].entity()->identifier();
    }
    msg << " } ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }

  return static_cast<unsigned>(side_id) ;
}

}

void use_case_13_generate_sides(
  mesh::BulkData & mesh , const bool skin_only )
{
  const mesh::fem::FEMMetaData & meta_data = mesh::fem::FEMMetaData::get(mesh);
  const stk::mesh::EntityRank element_rank = meta_data.element_rank();
  const unsigned p_rank = mesh.parallel_rank();

  mesh.modification_begin();

  // For each element
  //   for each element face
  //     get the neighbor element
  //     if no neighbor and owned then generate the face
  //     else if element has smaller id and
  //             either element or neighbor is local
  //       then generate face and attach to neighbor

  const std::vector<mesh::Bucket*> & es = mesh.buckets(element_rank);

  for ( std::vector<mesh::Bucket*>::const_iterator
        ie = es.begin() ; ie != es.end() ; ++ie ) {

    const mesh::Bucket & bucket = **ie;
    size_t n = bucket.size();

    for ( size_t j = 0; j < n; ++j) {

      // The const_cast indicates that there is an API issue to be addressed:
      //mesh::Entity & element = const_cast<mesh::Entity&>( bucket[j] );
      mesh::Entity & element = bucket[j];

      const CellTopologyData * const elem_top = mesh::fem::get_cell_topology( element ).getCellTopologyData();

      if ( NULL == elem_top ) {
        throw std::runtime_error( std::string("Element has no topology" ) );
      }

      // Generate if either the element or its neighbor is owned...
      // Only generate if the element has a smaller identifier

      for ( unsigned i = 0 ; i < elem_top->side_count ; ++i ) {

        mesh::Entity * const elem_neighbor = get_side_neighbor(*elem_top,element,i);

        const bool element_owned  = p_rank == element.owner_rank();
        const bool neighbor_owned = elem_neighbor &&
                                    elem_neighbor->owner_rank() == p_rank ;

        const bool create_side =
          ( element_owned || neighbor_owned ) &&
          ( ! elem_neighbor ||
            ( ! skin_only &&
              element.identifier() < elem_neighbor->identifier() ) );

        if ( create_side ) {

          const CellTopologyData * const side_top  = elem_top->side[i].topology ;
          const unsigned * const side_node = elem_top->side[i].node ;
          const unsigned         side_type = side_top->dimension ;

          const mesh::EntityId side_id = element.identifier() * 10 + i + 1;

          mesh::PartVector parts ;

          mesh::Entity & side = mesh.declare_entity( side_type, side_id , parts );

          mesh::PairIterRelation rel = element.relations( mesh::fem::FEMMetaData::NODE_RANK );

          for ( unsigned k = 0 ; k < side_top->node_count ; ++k ) {
            mesh::Entity & node = * rel[ side_node[k] ].entity();
            mesh.declare_relation( side , node , k );
          }

          mesh.declare_relation( element , side , i );

          if ( elem_neighbor ) {

            const unsigned other_side_id =
              determine_local_side_id( *elem_neighbor , side );

            mesh.declare_relation( *elem_neighbor , side , other_side_id );
          }
        }
      }
    }
  }

  mesh.modification_end();
}

} // namespace stk_use_cases

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#include <generated/Iogn_GeneratedMesh.h>

namespace stk_use_cases {

void use_case_13_generate_mesh(
  mesh::BulkData & mesh ,
  const unsigned N[] ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & /*elem_node_coord */,
  mesh::Part & hex_block )
{
  mesh.modification_begin();

  const unsigned parallel_size = mesh.parallel_size();
  const unsigned parallel_rank = mesh.parallel_rank();

  double t = 0 ;
  size_t num_hex = 0 ;
  size_t num_shell = 0 ;
  size_t num_nodes = 0 ;
  size_t num_block = 0 ;
  int error_flag = 0 ;

  try {

    Iogn::GeneratedMesh gmesh( N[0], N[1], N[2], parallel_size, parallel_rank );

    num_nodes = gmesh.node_count_proc();
    num_block = gmesh.block_count();

    t = wall_time();

    std::vector<int> node_map( num_nodes , 0 );

    gmesh.node_map( node_map );

    ThrowRequire( num_nodes == node_map.size() );

    {
      for ( size_t i = 1 ; i <= num_block ; ++i ) {
        const size_t                     num_elem = gmesh.element_count_proc(i);
        const std::pair<std::string,int> top_info = gmesh.topology_type(i);

	std::vector<int> elem_map( num_elem , 0 );
	gmesh.element_map( i, elem_map );

        std::vector<int> elem_conn( num_elem * top_info.second );

        gmesh.connectivity( i , elem_conn );

        if ( top_info.second == 8 ) {

          for ( size_t j = 0 ; j < num_elem ; ++j ) {

            const int * const local_node_id = & elem_conn[ j * 8 ] ;

            const stk::mesh::EntityId node_id[8] = {
	      local_node_id[0] ,
	      local_node_id[1] ,
	      local_node_id[2] ,
	      local_node_id[3] ,
	      local_node_id[4] ,
	      local_node_id[5] ,
	      local_node_id[6] ,
	      local_node_id[7]
            };

            const stk::mesh::EntityId elem_id = elem_map[ j ];

            mesh::fem::declare_element( mesh , hex_block , elem_id , node_id );

            ++num_hex ;
          }
        }
      }
    }

    std::vector<double> node_coordinates( 3 * node_map.size() );

    gmesh.coordinates( node_coordinates );

    if ( 3 * node_map.size() != node_coordinates.size() ) {
      std::ostringstream msg ;
      msg << "  P" << mesh.parallel_rank()
          << ": ERROR, node_map.size() = "
          << node_map.size()
          << " , node_coordinates.size() / 3 = "
          << ( node_coordinates.size() / 3 );
      throw std::runtime_error( msg.str() );
    }

    for ( unsigned i = 0 ; i < node_map.size() ; ++i ) {
      const unsigned i3 = i * 3 ;

      mesh::Entity * const node = mesh.get_entity( mesh::fem::FEMMetaData::NODE_RANK , node_map[i] );

      if ( NULL == node ) {
        std::ostringstream msg ;
        msg << "  P:" << mesh.parallel_rank()
            << " ERROR, Node not found: "
            << node_map[i] << " = node_map[" << i << "]" ;
        throw std::runtime_error( msg.str() );
      }

      double * const data = field_data( node_coord , *node );
      data[0] = node_coordinates[ i3 + 0 ];
      data[1] = node_coordinates[ i3 + 1 ];
      data[2] = node_coordinates[ i3 + 2 ];
    }
  }
  catch ( const std::exception & X ) {
    std::cout << "  P:" << mesh.parallel_rank() << ": " << X.what()
              << std::endl ;
    std::cout.flush();
    error_flag = 1 ;
  }
  catch( ... ) {
    std::cout << "  P:" << mesh.parallel_rank()
              << " Caught unknown exception"
              << std::endl ;
    std::cout.flush();
    error_flag = 1 ;
  }

  all_reduce( mesh.parallel() , ReduceMax<1>( & error_flag ) );

  if ( error_flag ) {
    std::string msg( "Failed mesh generation" );
    throw std::runtime_error( msg );
  }

  mesh.modification_end();

  double dt = wall_dtime( t );

  all_reduce( mesh.parallel() , ReduceMax<1>( & dt ) );

  std::cout << "  P" << mesh.parallel_rank()
            << ": Meshed Hex = " << num_hex
            << " , Shell = " << num_shell
            << " , Node = " << num_nodes
            << " in " << dt << " sec"
            << std::endl ;
  std::cout.flush();
}

} // namespace stk_use_cases

//----------------------------------------------------------------------


