/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <vector>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

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
//----------------------------------------------------------------------

namespace stk_unit_tests {

//----------------------------------------------------------------------
// Type declarations and function prototypes.

enum { SpatialDim   = 3 };

typedef mesh::Field<double>                    ScalarFieldType ;
typedef mesh::Field<double,mesh::Cartesian>    VectorFieldType ;
typedef mesh::Field<double*,mesh::ElementNode> ElementNodePointerFieldType ;

void use_case_6_generate_mesh(
  mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  mesh::Part & hex20_block ,
  mesh::Part & wedge15_block ,
  mesh::Part & side_part );

void verify_boundary_field_data( const mesh::BulkData & mesh ,
                                 mesh::Part & side_part ,
                                 VectorFieldType & boundary_field );

#include <unit_tests/centroid_algorithm.hpp>

//----------------------------------------------------------------------

namespace {

// Example using generic programming techniques:

template< class Traits_Full ,
          class Traits_Linear ,
          class PressureField ,
          class VelocityField >
void verify_pressure_velocity_stencil(
  const mesh::BulkData & M ,
  const mesh::Part     & element_part ,
  const mesh::Part     & linear_node_part ,
  const PressureField  & pressure ,
  const VelocityField  & velocity )
{
  typedef Traits_Full   element_traits ;
  typedef Traits_Linear element_linear_traits ;
  typedef typename mesh::FieldTraits< PressureField >::data_type scalar_type ;
  typedef typename mesh::FieldTraits< VelocityField >::data_type scalar_type_2 ;

  StaticAssert< SameType< scalar_type , scalar_type_2 >::value >::ok();

  STKUNIT_ASSERT(
    (int) element_traits::dimension     ==
    (int) element_linear_traits::dimension &&
    (int) element_traits::vertex_count ==
    (int) element_linear_traits::vertex_count &&
    (int) element_traits::edge_count   ==
    (int) element_linear_traits::edge_count &&
    (int) element_traits::face_count   ==
    (int) element_linear_traits::face_count );

  const std::vector<mesh::Bucket*> & buckets = M.buckets( mesh::Element );

  for ( std::vector<mesh::Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {
    mesh::Bucket & bucket = **k ;

    if ( has_superset( bucket, element_part ) ) {

      for ( mesh::Bucket::iterator
            i = bucket.begin() ; i != bucket.end() ; ++i ) {
        mesh::Entity & elem = *i ;

        mesh::PairIterRelation rel = elem.relations( mesh::Node );

        STKUNIT_ASSERT_EQUAL( (unsigned) rel.size() ,
                              (unsigned) element_traits::node_count );

        for ( unsigned j = 0 ; j < element_traits::node_count ; ++j ) {
          mesh::Entity & node = * rel[j].entity();
          const mesh::Bucket & node_bucket = node.bucket();
          mesh::PartVector node_parts ;

          node_bucket.supersets( node_parts );

          scalar_type * const p = mesh::field_data( pressure , node );
          scalar_type * const v = mesh::field_data( velocity , node );

          if ( j < element_linear_traits::node_count ) {
            STKUNIT_ASSERT( mesh::contain( node_parts , linear_node_part ) );
            STKUNIT_ASSERT( p != NULL );
          }
          else {
            STKUNIT_ASSERT( ! mesh::contain( node_parts , linear_node_part ) );
            STKUNIT_ASSERT( p == NULL );
          }

          STKUNIT_ASSERT( v != NULL );
        }
      }
    }
  }
}

}

//----------------------------------------------------------------------
// The 'main' driver for the use case.

void use_case_6_driver( MPI_Comm /*comm */)
{
  STKUNIT_ASSERT( (int) mesh::Cartesian::Size == (int) SpatialDim );

  //--------------------------------------------------------------------
  // MeshBulkData meta data: part and field declarations

  mesh::MetaData meta_data( mesh::fem_entity_type_names() );

  mesh::Part & universal         = meta_data.universal_part();
  mesh::Part & block_hex20       = meta_data.declare_part("block_1",mesh::Element);
  mesh::Part & block_wedge15     = meta_data.declare_part("block_2",mesh::Element);
  mesh::Part & part_vertex_nodes = meta_data.declare_part("vertex_nodes",mesh::Node);
  mesh::Part & side_part         = meta_data.declare_part("sideset_1",mesh::Face);

  // Set the local topology attribute for the homogeneous element block parts:

  mesh::set_cell_topology< shards::Hexahedron<20> >( block_hex20 );
  mesh::set_cell_topology< shards::Wedge<15>      >( block_wedge15 );

  //--------------------------------
  // The vertex nodes of the hex and wedge elements are members
  // of the vertex part; however, the mid-edge nodes are not.
  //
  // Use an element-node stencil to define this relationship.

  // Declare that the Hexahedron<>  nodes of an element in the
  // hex20 element block are members of the linear part.

  meta_data.declare_part_relation(
    block_hex20 ,
    & mesh::element_node_stencil< shards::Hexahedron<8>  > ,
    part_vertex_nodes );

  // Declare that the Wedge<>  nodes of an element in the
  // wedge15 element block are members of the vertex part.

  meta_data.declare_part_relation(
    block_wedge15 ,
    & mesh::element_node_stencil< shards::Wedge<6>  > ,
    part_vertex_nodes );

  //--------------------------------
  // declare various fields an put them on all nodes:

  VectorFieldType & coordinates_field =
    meta_data.declare_field< VectorFieldType >( "coordinates" );

  mesh::put_field( coordinates_field , mesh::Node , universal );


  VectorFieldType & velocity_field =
    meta_data.declare_field< VectorFieldType >( "velocity" );

  mesh::put_field( velocity_field , mesh::Node , universal );


  VectorFieldType & centroid_field =
    meta_data.declare_field< VectorFieldType >( "centroid" );

  mesh::put_field( centroid_field , mesh::Element , universal );


  ScalarFieldType & temperature_field =
    meta_data.declare_field< ScalarFieldType >( "temperature" );

  mesh::put_field( temperature_field, mesh::Node, universal );

  //--------------------------------
  // The pressure field only exists on the vertex nodes:

  ScalarFieldType & pressure_field =
    meta_data.declare_field< ScalarFieldType >( "pressure" );

  mesh::put_field( pressure_field, mesh::Node, part_vertex_nodes);


  // The boundary field only exists on nodes in the sideset part

  VectorFieldType & boundary_field =
    meta_data.declare_field< VectorFieldType >( "boundary" );

  mesh::put_field( boundary_field , mesh::Node , side_part );

  //--------------------------------
  // Declare an aggressive "gather" field which is an
  // array of pointers to the element's nodes' coordinate field data.
  // The declaration specifies:
  //
  //     double * elem_node_coord[number_of_nodes]

  ElementNodePointerFieldType & elem_node_coord =
    meta_data.declare_field<ElementNodePointerFieldType>("elem_node_coord_ptr");

  meta_data.declare_field_relation(
    elem_node_coord ,
    & mesh::element_node_stencil<void> ,
    coordinates_field );

  // Declare the number of pointers in the elem_node_field to be
  // the number of nodes of the elements in the hex20 and wedge15
  // blocks accordingly.

  mesh::put_field(
    elem_node_coord, mesh::Element, block_hex20, shards::Hexahedron<20>::node_count );

  mesh::put_field(
    elem_node_coord, mesh::Element, block_wedge15, shards::Wedge<15>::node_count );

  //--------------------------------
  // Done with part and field (meta data) declarations.

  meta_data.commit();

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  // Declare the mesh bulk data conforming to the meta data
  // and distributed over the parallel machine:

  mesh::BulkData bulk_data( meta_data , MPI_COMM_WORLD );

  // Internal generation of the use case mesh:

  use_case_6_generate_mesh( bulk_data ,
                            coordinates_field ,
                            elem_node_coord ,
                            block_hex20 , block_wedge15 ,
                            side_part );

  centroid_algorithm< shards::Hexahedron<20> >( bulk_data ,
                                        centroid_field ,
                                        elem_node_coord ,
                                        block_hex20 );

  centroid_algorithm< shards::Wedge<15> >( bulk_data ,
                                   centroid_field ,
                                   elem_node_coord ,
                                   block_wedge15 );

  //--------------------------------------------------------------------
  // Internal unit-testing algorithms

  centroid_algorithm_unit_test_dimensions< shards::Hexahedron<20> >(
    bulk_data , centroid_field , elem_node_coord , block_hex20 );

  centroid_algorithm_unit_test_dimensions< shards::Wedge<15> >(
    bulk_data , centroid_field , elem_node_coord , block_wedge15 );

  verify_boundary_field_data( bulk_data ,
                              side_part ,
                              boundary_field );

  verify_pressure_velocity_stencil
    < shards::Hexahedron<20> , shards::Hexahedron<8>  >
    ( bulk_data , block_hex20 , part_vertex_nodes ,
      pressure_field , velocity_field );

  verify_pressure_velocity_stencil
    < shards::Wedge<15> , shards::Wedge<6>  >
    ( bulk_data , block_wedge15 , part_vertex_nodes ,
      pressure_field , velocity_field );

  //--------------------------------------------------------------------
}

//----------------------------------------------------------------------
/*----------------------------------------------------------------------
 * Two hexes and three wedges
 *
 *  Z = 0 plane:
 *
 *    Y
 *    ^
 *    !    58----61----64
 *    !    / \         / \
 *        /   \       /   \ 
 *      46    49    52    55 
 *      /       \   /       \
 *     /         \ /         \ 
 *   31----34----37----40----43
 *    |           |           |
 *    |           |           |
 *   16    19    22    25    28
 *    |           |           |
 *    |           |           |
 *    1-----4-----7----10----13  -----> X
 * 
 *  Z = -1 plane (midplane):
 *
 *    Y
 *    ^
 *    !    59----62----65
 *    !    / \         / \
 *        /   \       /   \ 
 *      47    50    53    56 
 *      /       \   /       \
 *     /         \ /         \ 
 *   32----35----38----41----44
 *    |           |           |
 *    |           |           |
 *   17    20    23    26    29
 *    |           |           |
 *    |           |           |
 *    2-----5-----8----11----14  -----> X
 * 
 *
 *  Z = -2 plane:
 *
 *    Y
 *    ^
 *    !    60----63----66
 *    !    / \         / \
 *        /   \       /   \ 
 *      48    51    54    57 
 *      /       \   /       \
 *     /         \ /         \ 
 *   33----36----39----42----45
 *    |           |           |
 *    |           |           |
 *   18    21    24    27    30
 *    |           |           |
 *    |           |           |
 *    3-----6-----9----12----15  -----> X
 *   
 *----------------------------------------------------------------------*/

enum { node_count   = 66 };
enum { number_hex   = 2 };
enum { number_wedge = 3 };

namespace {

static const double node_coord_data[ node_count ][ SpatialDim ] = {
  {  0 ,  0 ,  0 } , {  0 ,  0 , -1 } , {  0 ,  0 , -2 } ,
  {  1 ,  0 ,  0 } , {  1 ,  0 , -1 } , {  1 ,  0 , -2 } ,
  {  2 ,  0 ,  0 } , {  2 ,  0 , -1 } , {  2 ,  0 , -2 } ,
  {  3 ,  0 ,  0 } , {  3 ,  0 , -1 } , {  3 ,  0 , -2 } ,
  {  4 ,  0 ,  0 } , {  4 ,  0 , -1 } , {  4 ,  0 , -2 } ,

  {  0 ,  1 ,  0 } , {  0 ,  1 , -1 } , {  0 ,  1 , -2 } ,
  {  1 ,  1 ,  0 } , {  1 ,  1 , -1 } , {  1 ,  1 , -2 } ,
  {  2 ,  1 ,  0 } , {  2 ,  1 , -1 } , {  2 ,  1 , -2 } ,
  {  3 ,  1 ,  0 } , {  3 ,  1 , -1 } , {  3 ,  1 , -2 } ,
  {  4 ,  1 ,  0 } , {  4 ,  1 , -1 } , {  4 ,  1 , -2 } ,

  {  0 ,  2 ,  0 } , {  0 ,  2 , -1 } , {  0 ,  2 , -2 } ,
  {  1 ,  2 ,  0 } , {  1 ,  2 , -1 } , {  1 ,  2 , -2 } ,
  {  2 ,  2 ,  0 } , {  2 ,  2 , -1 } , {  2 ,  2 , -2 } ,
  {  3 ,  2 ,  0 } , {  3 ,  2 , -1 } , {  3 ,  2 , -2 } ,
  {  4 ,  2 ,  0 } , {  4 ,  2 , -1 } , {  4 ,  2 , -2 } ,

  {  0.5 , 3 , 0 } , { 0.5 , 3 , -1 } , { 0.5 , 3 , -2 } ,
  {  1.5 , 3 , 0 } , { 1.5 , 3 , -1 } , { 1.5 , 3 , -2 } ,
  {  2.5 , 3 , 0 } , { 2.5 , 3 , -1 } , { 2.5 , 3 , -2 } ,
  {  3.5 , 3 , 0 } , { 3.5 , 3 , -1 } , { 3.5 , 3 , -2 } ,

  {  1 , 4 , 0 } , { 1 , 4 , -1 } , { 1 , 4 , -2 } ,
  {  2 , 4 , 0 } , { 2 , 4 , -1 } , { 2 , 4 , -2 } ,
  {  3 , 4 , 0 } , { 3 , 4 , -1 } , { 3 , 4 , -2 }
};

static const stk::mesh::EntityId wedge_node_ids[3][ shards::Wedge<18>::node_count ] = {
  { 33 , 39 , 60 , 31 , 37 , 58 ,
    36 , 51 , 48 , 32 , 38 , 59 , 34 , 49 , 46 ,
    35 , 50 , 47 },
  { 39 , 45 , 66 , 37 , 43 , 64 ,
    42 , 57 , 54 , 38 , 44 , 65 , 40 , 55 , 52 ,
    41 , 56 , 53 },
  { 66 , 60 , 39 , 64 , 58 , 37 ,
    63 , 51 , 54 , 65 , 59 , 38 , 61 , 49 , 52 ,
    62 , 50 , 53 }
};

static const stk::mesh::EntityId hex_node_ids[2][ shards::Hexahedron<27>::node_count ] = {
  {  1 ,  7 ,  9 ,  3 , 31 , 37 , 39 , 33 ,
     4 ,  8 ,  6 ,  2 , 16 , 22 , 24 , 18 , 34 , 38 , 36 , 32 ,
    20 ,  5 , 35 , 17 , 23 , 19 , 21 } ,
  {  7 , 13 , 15 ,  9 , 37 , 43 , 45 , 39 ,
    10 , 14 , 12 ,  8 , 22 , 28 , 30 , 24 , 40 , 44 , 42 , 38 ,
    26 , 11 , 41 , 23 , 29 , 25 , 27 }
};

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

  STKUNIT_ASSERT_EQUAL( (unsigned) rel.size() , node_count );

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
      STKUNIT_ASSERT( n1 == SpatialDim );
      STKUNIT_ASSERT( node_coord_array.size() == SpatialDim );
    }

    double * const node_data = node_coord_array.contiguous_data();
    STKUNIT_ASSERT( elem_data[j] == node_data );
  }
}

template< class ElemTraits >
void verify_elem_side_node( const stk::mesh::EntityId * const elem_nodes ,
                            const unsigned local_side ,
                            const mesh::Entity & side )
{
  const CellTopologyData * const elem_top = shards::getCellTopologyData< ElemTraits >();

  const CellTopologyData * const side_top = elem_top->side[ local_side ].topology ;
  const unsigned         * const side_node_map = elem_top->side[ local_side ].node ;

  const mesh::PairIterRelation rel = side.relations( mesh::Node );

  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {

    STKUNIT_ASSERT_EQUAL( elem_nodes[ side_node_map[i] ] ,
                          rel[i].entity()->identifier() );
  }
}

}

void verify_boundary_field_data( const mesh::BulkData & mesh ,
                                 mesh::Part & side_part ,
                                 VectorFieldType & boundary_field )
{
  unsigned num_side_nodes = 0 ;

  const std::vector<mesh::Bucket*> & buckets = mesh.buckets( mesh::Node );

  for ( std::vector<mesh::Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {

    mesh::Bucket & bucket = **k ;

    void * data = mesh::field_data( boundary_field , bucket.begin() );

    if ( has_superset( bucket, side_part ) ) {

      STKUNIT_ASSERT( NULL != data );

      num_side_nodes += bucket.size();
    }
    else {
      STKUNIT_ASSERT( NULL == data );
    }
  }

  STKUNIT_ASSERT_EQUAL( num_side_nodes , 20u );
}

//----------------------------------------------------------------------

void use_case_6_generate_mesh(
  mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  mesh::Part & hex20_block ,
  mesh::Part & wedge15_block ,
  mesh::Part & side_part )
{
  stk::mesh::EntityId elem_id = 1 ;
  stk::mesh::EntityId face_id = 1 ;

  mesh::PartVector side_add ; insert( side_add , side_part );

  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id , ++face_id ) {
    mesh::Entity & elem =
      mesh::declare_element( mesh, hex20_block, elem_id, hex_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord , 20 );

    mesh::Entity & face = mesh::declare_element_side( mesh, face_id, elem, 0 );

    verify_elem_side_node< shards::Hexahedron<20> >( hex_node_ids[i] , 0 , face );

    mesh.change_entity_parts( face , side_add );
  }

  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id , ++face_id ) {
    mesh::Entity & elem =
      mesh::declare_element( mesh, wedge15_block, elem_id, wedge_node_ids[i] );

    verify_elem_node_coord( elem , elem_node_coord , node_coord , 15 );

    mesh::Entity & face = mesh::declare_element_side( mesh, face_id , elem , 4 );

    verify_elem_side_node< shards::Wedge<15> >( wedge_node_ids[i] , 4 , face );

    mesh.change_entity_parts( face , side_add );
  }

  for ( unsigned i = 0 ; i < node_count ; ++i ) {
    mesh::Entity * const node = mesh.get_entity( mesh::Node, i + 1 );

    if ( node != NULL ) {
      double * const coord = mesh::field_data( node_coord , *node );

      coord[0] = node_coord_data[i][0] ;
      coord[1] = node_coord_data[i][1] ;
      coord[2] = node_coord_data[i][2] ;
    }
  }

  // No parallel stuff for now
}

//----------------------------------------------------------------------

} // namespace stk_unit_test

STKUNIT_UNIT_TEST(UnitTestUseCase_6, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  stk_unit_tests::use_case_6_driver(MPI_COMM_WORLD);
}

