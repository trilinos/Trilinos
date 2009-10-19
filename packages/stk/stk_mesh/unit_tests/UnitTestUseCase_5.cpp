
#include <mpi.h>
#include <vector>

#include <Shards_BasicTopologies.hpp>

#include <unit_tests/stk_utest_macros.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>


using namespace stk ;

//----------------------------------------------------------------------
// This file contains the implementation of use-case 5: heterogeneous mesh.
// The function 'use_case_5_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk_unit_tests {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

typedef mesh::Field<double>                 ScalarFieldType ;
typedef mesh::Field<double,mesh::Cartesian> VectorFieldType ;

// Specification for the aggressive gather pointer-field for elements.

typedef mesh::Field<double*,mesh::ElementNode> ElementNodePointerFieldType ;

// Centroid algorithm generic programming functions:

#include <unit_tests/centroid_algorithm.hpp>

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_5_generate_mesh(
  mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  mesh::Part & hex_block ,
  mesh::Part & wedge_block ,
  mesh::Part & tetra_block ,
  mesh::Part & pyramid_block ,
  mesh::Part & quad_shell_block ,
  mesh::Part & tri_shell_block );

//--------------------------------------------------------------------
//
// main driver for use-case 5: heterogeneous element mesh.
//

void use_case_5_driver( MPI_Comm /*comm */)
{
  //--------------------------------------------------------------------
  // Declare the mesh meta data: six element blocks and associated fields

  mesh::MetaData mesh_meta_data( mesh::fem_entity_type_names() );

  //--------------------------------
  // Element-block declarations typically occur when reading the
  // mesh-file meta-data, and thus won't usually appear in application code.
  // Declaring the element blocks and associating an element traits
  // with each element block.

  mesh::Part & universal        = mesh_meta_data.universal_part();
  mesh::Part & block_hex        = mesh_meta_data.declare_part("block_1", mesh::Element );
  mesh::Part & block_wedge      = mesh_meta_data.declare_part("block_2", mesh::Element );
  mesh::Part & block_tet        = mesh_meta_data.declare_part("block_3", mesh::Element );
  mesh::Part & block_pyramid    = mesh_meta_data.declare_part("block_4", mesh::Element );
  mesh::Part & block_quad_shell = mesh_meta_data.declare_part("block_5", mesh::Element );
  mesh::Part & block_tri_shell  = mesh_meta_data.declare_part("block_6", mesh::Element );

  mesh::set_cell_topology< shards::Hexahedron<8>          >( block_hex );
  mesh::set_cell_topology< shards::Wedge<6>               >( block_wedge );
  mesh::set_cell_topology< shards::Tetrahedron<4>         >( block_tet );
  mesh::set_cell_topology< shards::Pyramid<5>             >( block_pyramid );
  mesh::set_cell_topology< shards::ShellQuadrilateral<4>  >( block_quad_shell );
  mesh::set_cell_topology< shards::ShellTriangle<3>       >( block_tri_shell );

  //--------------------------------
  // Declaring fields and put them on the mesh.

  VectorFieldType & coordinates_field =
    mesh_meta_data.declare_field< VectorFieldType >( "coordinates" );

  mesh::put_field( coordinates_field , mesh::Node , universal );


  VectorFieldType & centroid_field =
    mesh_meta_data.declare_field< VectorFieldType >( "centroid" );

  mesh::put_field( centroid_field , mesh::Element , universal );


  ScalarFieldType & temperature_field =
    mesh_meta_data.declare_field< ScalarFieldType >( "temperature" );

  mesh::put_field( temperature_field, mesh::Node, universal );


  // Declare that the 'volume' field exists on elements in the
  // four solid element-blocks (it does not exist for shell elements).
  // If the field existed on all elements, then we would have
  // specified the universal part as we did above for the centroid field.

  ScalarFieldType & volume_field =
    mesh_meta_data.declare_field< ScalarFieldType >( "volume" );

  mesh::put_field( volume_field, mesh::Element, block_hex );
  mesh::put_field( volume_field, mesh::Element, block_wedge );
  mesh::put_field( volume_field, mesh::Element, block_tet );
  mesh::put_field( volume_field, mesh::Element, block_pyramid );

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
    & mesh::element_node_stencil<void> ,
    coordinates_field );

  mesh::put_field(
    elem_node_coord , mesh::Element , block_hex , shards::Hexahedron<> ::node_count );

  mesh::put_field(
    elem_node_coord , mesh::Element , block_wedge , shards::Wedge<> ::node_count );

  mesh::put_field(
    elem_node_coord , mesh::Element , block_tet , shards::Tetrahedron<> ::node_count );

  mesh::put_field(
    elem_node_coord , mesh::Element , block_pyramid , shards::Pyramid<> ::node_count );

  mesh::put_field(
    elem_node_coord, mesh::Element, block_quad_shell,
    shards::ShellQuadrilateral<> ::node_count);

  mesh::put_field(
    elem_node_coord, mesh::Element, block_tri_shell,
    shards::ShellTriangle<> ::node_count );

  //--------------------------------
  // Commit (finalize) the meta data.  Is now ready to be used
  // in the creation and management of mesh bulk data.

  mesh_meta_data.commit();

  //--------------------------------------------------------------------
  // Declare the mesh bulk data, with 1000 entities per field data chunk.
  // (The parameter '1000' is the application's preferred "workset" size.)

  mesh::BulkData mesh_bulk_data( mesh_meta_data , MPI_COMM_WORLD , 1000 );

  // In a typical app, the mesh would be read from file at this point.
  // But in this use-case, we generate the mesh.

  use_case_5_generate_mesh(
    mesh_bulk_data ,
    coordinates_field ,
    elem_node_coord ,
    block_hex , block_wedge ,
    block_tet , block_pyramid ,
    block_quad_shell , block_tri_shell );

  // Example generic-programming element algorithms.
  // The following calls compute element-centroids for the input element-block.

  // The mapping of application algorithms onto a mesh is a capability
  // of the "legacy" framework, implemented by the Fmwk::Algorithm,
  // Fmwk::WorksetAlgorithm, Fmwk::Mechanics, and Fmwk::Region classes.
  //
  // This functionality is not part of the toolkit's mesh module.
  // However, there may be a separate apply-algorithm-to-mesh 
  // module within the toolkit to provide this functionality.

  centroid_algorithm< shards::Hexahedron<>  >( mesh_bulk_data ,
                                       centroid_field ,
                                       elem_node_coord ,
                                       block_hex );

  centroid_algorithm< shards::Wedge<>  >( mesh_bulk_data ,
                                  centroid_field ,
                                  elem_node_coord ,
                                  block_wedge );

  centroid_algorithm< shards::Tetrahedron<>  >( mesh_bulk_data ,
                                        centroid_field ,
                                        elem_node_coord ,
                                        block_tet );

  centroid_algorithm< shards::Pyramid<>  >( mesh_bulk_data ,
                                    centroid_field ,
                                    elem_node_coord ,
                                    block_pyramid );

  centroid_algorithm< shards::ShellQuadrilateral<>  >( mesh_bulk_data ,
                                               centroid_field ,
                                               elem_node_coord ,
                                               block_quad_shell );

  centroid_algorithm< shards::ShellTriangle<>  >( mesh_bulk_data ,
                                          centroid_field ,
                                          elem_node_coord ,
                                          block_tri_shell );

  //------------------------------
  // Done with use-case example.
  // Now have some unit testing functions.

  centroid_algorithm_unit_test_dimensions< shards::Hexahedron<>  >(
    mesh_bulk_data , centroid_field , elem_node_coord , block_hex );

  centroid_algorithm_unit_test_dimensions< shards::Wedge<>  >(
    mesh_bulk_data , centroid_field , elem_node_coord , block_wedge );

  centroid_algorithm_unit_test_dimensions< shards::Tetrahedron<>  >(
    mesh_bulk_data , centroid_field , elem_node_coord , block_tet );

  centroid_algorithm_unit_test_dimensions< shards::Pyramid<>  >(
    mesh_bulk_data , centroid_field , elem_node_coord , block_pyramid );

  centroid_algorithm_unit_test_dimensions< shards::ShellQuadrilateral<>  >(
    mesh_bulk_data , centroid_field , elem_node_coord , block_quad_shell );

  centroid_algorithm_unit_test_dimensions< shards::ShellTriangle<>  >(
    mesh_bulk_data , centroid_field , elem_node_coord , block_tri_shell );

  //------------------------------
}

//--------------------------------------------------------------------
/*----------------------------------------------------------------------
 * Internal use-case #5 mesh generation.
 *
 * Three hexes, three wedges, three tets, two pyramids,
 * three quad shells, and three triangle shells.
 *
 *  Z = 0 plane:
 *
 *    Y
 *    ^   9      10
 *    !   *-------*                                               
 *    !  / \     / \ 
 *    ! /   \   /   \ 
 *     /     \ /     \ 
 *    *-------*-------*-------*
 *   5|      6|      7|      8|
 *    |       |       |       |
 *    |       |       |       |
 *    *-------*-------*-------*    ----> X
 *    1       2       3       4
 *   
 *  Z = -1 plane:
 *
 *    Y
 *    ^  19      20
 *    !   *-------*
 *    !  / \     / \ 
 *    ! /   \   /   \ 
 *     /     \ /     \ 
 *    *-------*-------*-------*
 *  15|     16|     17|     18|
 *    |       |       |       |
 *    |       |       |       |
 *    *-------*-------*-------*    ----> X
 *   11      12      13      14
 *   
 *
 *  Last node (#21) at Z = -2, translated from node #16
 *----------------------------------------------------------------------*/

enum { node_count = 21 };

enum { number_hex = 3 };
enum { number_wedge = 3 };
enum { number_tetra = 3 };
enum { number_pyramid = 2 };
enum { number_shell_quad = 3 };
enum { number_shell_tri = 3 };

namespace {

static const double node_coord_data[ node_count ][ SpatialDim ] = {
  { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
  { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
  { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
  { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
  { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
  { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
  { 1 , 1 , -2 } };

static const stk::mesh::EntityId hex_node_ids[3][ shards::Hexahedron<> ::node_count ] = {
  { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
  { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
  { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

static const stk::mesh::EntityId wedge_node_ids[3][ shards::Wedge<> ::node_count ] = {
  { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
  { 10 ,  9 ,  6 , 20 , 19 , 16 } ,
  { 16 , 17 , 20 ,  6 ,  7 , 10 } };

static const stk::mesh::EntityId tetra_node_ids[3][ shards::Tetrahedron<> ::node_count ] = {
  { 15 , 19 , 16 , 21 } ,
  { 19 , 20 , 16 , 21 } ,
  { 16 , 20 , 17 , 21 } };

static const stk::mesh::EntityId pyramid_node_ids[2][ shards::Pyramid<> ::node_count ] = {
  { 11 , 15 , 16 , 12 , 21 } ,
  { 12 , 16 , 17 , 13 , 21 } };

static const stk::mesh::EntityId shell_quad_node_ids[3][ shards::ShellQuadrilateral<> ::node_count ]={
  { 9 , 6 , 16 , 19 } ,
  { 6 , 7 , 17 , 16 } ,
  { 7 , 8 , 18 , 17 } };

static const stk::mesh::EntityId shell_tri_node_ids[3][ shards::ShellTriangle<> ::node_count ] ={
  { 19 , 16 , 21 } ,
  { 16 , 17 , 21 } ,
  { 17 , 13 , 21 } };

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void verify_elem_node_coord(
  mesh::Entity & elem ,
  const ElementNodePointerFieldType & elem_node_coord ,
  const VectorFieldType & node_coord ,
  const unsigned node_count )
{
  mesh::PairIterRelation rel = elem.relations( mesh::Node );

  STKUNIT_ASSERT( (unsigned) rel.size() == node_count );

  mesh::EntityArray< ElementNodePointerFieldType >
    elem_node_array( elem_node_coord , elem );

  {
    const unsigned n1 = elem_node_array.dimension<0>();
    STKUNIT_ASSERT( n1 == node_count );
    STKUNIT_ASSERT( (unsigned) elem_node_array.size() == node_count );
  }

  double * const * const elem_data = elem_node_array.contiguous_data();

  for ( unsigned j = 0 ; j < node_count ; ++j ) {
    mesh::Entity & node = * rel[j].entity();

    mesh::EntityArray< VectorFieldType > node_coord_array( node_coord , node );

    {
      const unsigned n1 = node_coord_array.dimension<0>();
      STKUNIT_ASSERT( n1 == (unsigned) SpatialDim );
      STKUNIT_ASSERT( node_coord_array.size() == SpatialDim );
    }

    double * const node_data = node_coord_array.contiguous_data();
    STKUNIT_ASSERT( elem_data[j] == node_data );
  }
}

}

//----------------------------------------------------------------------

void use_case_5_generate_mesh(
  mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  mesh::Part & hex_block ,
  mesh::Part & wedge_block ,
  mesh::Part & tetra_block ,
  mesh::Part & pyramid_block ,
  mesh::Part & quad_shell_block ,
  mesh::Part & tri_shell_block )
{
  stk::mesh::EntityId elem_id = 1 ;

  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id ) {
    mesh::Entity & elem =
      declare_element( mesh, hex_block, elem_id, hex_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord ,
                            shards::Hexahedron<> ::node_count );
  }

  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id ) {
    mesh::Entity & elem =
      declare_element( mesh, wedge_block, elem_id, wedge_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord ,
                            shards::Wedge<> ::node_count );
  }

  for ( unsigned i = 0 ; i < number_tetra ; ++i , ++elem_id ) {
    mesh::Entity & elem =
      declare_element( mesh, tetra_block, elem_id, tetra_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord ,
                            shards::Tetrahedron<> ::node_count );
  }

  for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++elem_id ) {
    mesh::Entity & elem =
      declare_element( mesh, pyramid_block, elem_id, pyramid_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord ,
                            shards::Pyramid<> ::node_count );
  }

  for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++elem_id ) {
    mesh::Entity & elem =
      declare_element( mesh, quad_shell_block, elem_id, shell_quad_node_ids[i]);

    verify_elem_node_coord( elem , elem_node_coord , node_coord ,
                            shards::ShellQuadrilateral<> ::node_count );
  }

  for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++elem_id ) {
    mesh::Entity & elem =
      declare_element( mesh, tri_shell_block, elem_id, shell_tri_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord ,
                            shards::ShellTriangle<> ::node_count );
  }

  for ( unsigned i = 0 ; i < node_count ; ++i ) {
    mesh::Entity * const node = mesh.get_entity( mesh::Node , i + 1 );
    
    STKUNIT_ASSERT( node != NULL );

    double * const coord = field_data( node_coord , *node );

    coord[0] = node_coord_data[i][0] ;
    coord[1] = node_coord_data[i][1] ;
    coord[2] = node_coord_data[i][2] ;
  }

  // No parallel stuff for now
}

} // namespace stk_unit_tests

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestUseCase_5, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;

  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);

  stk_unit_tests::use_case_5_driver(MPI_COMM_WORLD);
}


