/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>

#include <use_cases/centroid_algorithm.hpp>

#include <use_cases/UseCase_4.hpp>

enum { SpatialDim   = 3 };

namespace stk_classic{
namespace mesh {
namespace use_cases {

typedef shards::Hexahedron<8>   Hex8;
typedef shards::Hexahedron<27>  Hex27;
typedef shards::Wedge<6>        Wedge6;
typedef shards::Wedge<18>       Wedge18;

UseCase_4_Mesh::UseCase_4_Mesh( stk_classic::ParallelMachine comm, bool doCommit ) :
  m_fem_metaData( SpatialDim )
  , m_bulkData( fem::FEMMetaData::get_meta_data(m_fem_metaData) , comm )
  , m_elem_rank( m_fem_metaData.element_rank() )
  , m_side_rank( m_fem_metaData.side_rank() )
  , m_node_rank( m_fem_metaData.node_rank() )
  , m_block_hex27( stk_classic::mesh::fem::declare_part<Hex27>( m_fem_metaData, "block_1"))
  , m_block_wedge18( stk_classic::mesh::fem::declare_part<Wedge18>( m_fem_metaData, "block_2"))
  , m_part_vertex_nodes( m_fem_metaData.declare_part("vertex_nodes", m_node_rank ))
  , m_side_part(         m_fem_metaData.declare_part("sideset_1", m_side_rank ))
  , m_coordinates_field(m_fem_metaData.declare_field< VectorFieldType >( "coordinates" ))
  , m_velocity_field(m_fem_metaData.declare_field< VectorFieldType >( "velocity" ))
  , m_centroid_field(m_fem_metaData.declare_field< VectorFieldType >( "centroid" ))
  , m_temperature_field(m_fem_metaData.declare_field< ScalarFieldType >( "temperature" ))
  , m_pressure_field(m_fem_metaData.declare_field< ScalarFieldType >( "pressure" ))
  , m_boundary_field(m_fem_metaData.declare_field< VectorFieldType >( "boundary" ))
  , m_element_node_coordinates_field(m_fem_metaData.declare_field< ElementNodePointerFieldType >( "elem_node_coord" ) )
{
  //--------------------------------
  // The vertex nodes of the hex and wedge elements are members
  // of the vertex part; however, the mid-edge nodes are not.
  //
  // Use an element-node stencil to define this relationship.

  // Declare that the Hexahedron<>  nodes of an element in the
  // hex27 element block are members of the linear part.

  m_fem_metaData.declare_part_relation(
    m_block_hex27 ,
    & fem::element_node_stencil< Hex8, SpatialDim > ,
    m_part_vertex_nodes );

  // Declare that the Wedge<>  nodes of an element in the
  // wedge18 element block are members of the vertex part.

  m_fem_metaData.declare_part_relation(
    m_block_wedge18 ,
    & fem::element_node_stencil< Wedge6, SpatialDim > ,
    m_part_vertex_nodes );

  // Where fields exist on the mesh:
  Part & universal = m_fem_metaData.universal_part();

  put_field( m_coordinates_field , m_node_rank , universal );
  put_field( m_velocity_field , m_node_rank , universal );
  put_field( m_centroid_field , m_elem_rank , universal );
  put_field( m_temperature_field, m_node_rank, universal );

  // The pressure field only exists on the vertex nodes:
  put_field( m_pressure_field, m_node_rank, m_part_vertex_nodes );

  // The boundary field only exists on nodes in the sideset part
  put_field( m_boundary_field, m_node_rank, m_side_part );

  // Set up the relationship between nodal coord field and elem node coord field
  m_fem_metaData.declare_field_relation(
    m_element_node_coordinates_field ,
    fem::get_element_node_stencil(SpatialDim) ,
    m_coordinates_field
    );

  put_field( m_element_node_coordinates_field, m_elem_rank, m_block_hex27, Hex27::node_count );
  put_field( m_element_node_coordinates_field, m_elem_rank, m_block_wedge18, Wedge18::node_count );

  if (doCommit)
    m_fem_metaData.commit();
}

UseCase_4_Mesh::~UseCase_4_Mesh()
{ }

//------------------------------------------------------------------------------
// Use-case specific mesh data:

enum { node_count   = 66 };
enum { number_hex   = 2 };
enum { number_wedge = 3 };

namespace {

// Hard coded node coordinate data for all the nodes in the entire mesh
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

// Hard coded wedge node ids for all the wedge nodes in the entire mesh
static const EntityId wedge_node_ids[3][ Wedge18::node_count ] = {
  { 33 , 39 , 60 , 31 , 37 , 58 ,                // vertices
    36 , 51 , 48 , 32 , 38 , 59 , 34 , 49 , 46 , // mid edge nodes
    35 , 50 , 47 },                              // mid side nodes
  { 39 , 45 , 66 , 37 , 43 , 64 ,                // vertices
    42 , 57 , 54 , 38 , 44 , 65 , 40 , 55 , 52 , // mid edge nodes
    41 , 56 , 53 },                              // mid side nodes
  { 66 , 60 , 39 , 64 , 58 , 37 ,                // vertices
    63 , 51 , 54 , 65 , 59 , 38 , 61 , 49 , 52 , // mid edge nodes
    62 , 50 , 53 }                               // mid side nodes

    };

// Hard coded hex node ids for all the hex nodes in the entire mesh
static const EntityId hex_node_ids[2][ Hex27::node_count ] = {
  {  1 ,  7 ,  9 ,  3 , 31 , 37 , 39 , 33 ,                     // vertices
     4 ,  8 ,  6 ,  2 , 16 , 22 , 24 , 18 , 34 , 38 , 36 , 32 , // mid edge nodes
     20 ,  5 , 35 , 17 , 23 , 19 , 21 } ,                       // mid side nodes
  {  7 , 13 , 15 ,  9 , 37 , 43 , 45 , 39 ,                     // vertices
    10 , 14 , 12 ,  8 , 22 , 28 , 30 , 24 , 40 , 44 , 42 , 38 , // mid edge nodes
     26 , 11 , 41 , 23 , 29 , 25 , 27 }                         // mid side nodes

};

}

//------------------------------------------------------------------------------

void UseCase_4_Mesh::populate()
{
  m_bulkData.modification_begin();

  EntityId elem_id = 1;
  EntityId face_id = 1;

  // Declare element with its nodes.
  // Declare a side on an element.  This utility function creates the
  // side, element-side relations, and side-node relations.
  // It will NOT check if the side sandwiched between two elements.

  // Iterate over the number of desired hexs here and declares them
  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id , ++face_id ) {
    Entity & elem =
      fem::declare_element( m_bulkData, m_block_hex27, elem_id, hex_node_ids[i] );

    fem::declare_element_side( m_bulkData, face_id, elem, 0 /*local side id*/, &m_side_part);
  }

  // Iterate over the number of desired wedges here and declares the
  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id , ++face_id ) {
    Entity & elem =
      fem::declare_element( m_bulkData, m_block_wedge18, elem_id, wedge_node_ids[i] );

    fem::declare_element_side( m_bulkData, face_id , elem , 4 /*local side id*/, &m_side_part);
  }

  // For all nodes assign nodal coordinates
  for ( unsigned i = 0 ; i < node_count ; ++i ) {
    Entity * const node = m_bulkData.get_entity( m_node_rank, i + 1 );
    ThrowRequireMsg( node != NULL, i+1 );

    double * const coord = field_data( m_coordinates_field , *node );
    coord[0] = node_coord_data[i][0] ;
    coord[1] = node_coord_data[i][1] ;
    coord[2] = node_coord_data[i][2] ;
  }

  m_bulkData.modification_end();
}

//------------------------------------------------------------------------------

void runAlgorithms( const UseCase_4_Mesh & mesh )
{
  const BulkData & bulkData = mesh.m_bulkData ;
  VectorFieldType & centroid_field = mesh.m_centroid_field ;
  ElementNodePointerFieldType & elem_node_coord = mesh.m_element_node_coordinates_field ;
  Part & block_hex27 = mesh.m_block_hex27 ;
  Part & block_wedge18 = mesh.m_block_wedge18 ;

  // Run the centroid algorithm on the hexes:
  centroid_algorithm< Hex27 >( bulkData ,
                               centroid_field ,
                               elem_node_coord ,
                               block_hex27,
                               mesh.m_elem_rank );

  // Run the centroid algorithm on the wedges:
  centroid_algorithm< Wedge18 >( bulkData ,
                                 centroid_field ,
                                 elem_node_coord ,
                                 block_wedge18,
                                 mesh.m_elem_rank );

}

//------------------------------------------------------------------------------

namespace {

template< class CellTopology >
bool verify_elem_side_node( const EntityId * const elem_nodes ,
                            const unsigned local_side ,
                            const mesh::Entity & side )
{
  bool result = true;
  const CellTopologyData * const elem_top = shards::getCellTopologyData< CellTopology >();

  const CellTopologyData * const side_top = elem_top->side[ local_side ].topology ;
  const unsigned         * const side_node_map = elem_top->side[ local_side ].node ;

  const mesh::PairIterRelation rel = side.relations( fem::FEMMetaData::NODE_RANK );

  // Verify that the node relations are compatible with the cell topology data
  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {

    if ( elem_nodes[ side_node_map[i] ] != rel[i].entity()->identifier() ) {
      std::cerr << "Error!" << std::endl;
      result = false;
    }
  }
  return result;
}

bool verify_boundary_field_data( const BulkData & mesh ,
                                 Part & side_part ,
                                 const VectorFieldType & boundary_field )
{
  bool result = true;

  unsigned num_nodes_on_side_of_mesh  = 0 ;

  const std::vector<Bucket*> & buckets = mesh.buckets( fem::FEMMetaData::NODE_RANK );

  // Iterate over each bucket and sum the number of nodes
  // into num_nodes_on_side_of_mesh for each bucket which is a subset of side_part.
  // this should give the number of nodes on the front face of the mesh (see .hpp file).
  for ( std::vector<Bucket*>::const_iterator
          k = buckets.begin() ; k != buckets.end() ; ++k ) {
    Bucket & bucket = **k ;

    void * data = field_data( boundary_field , bucket.begin() );

    if ( has_superset( bucket, side_part ) ) {
      // Since this node is in the side part it should have a
      // boundary field.
      if ( data == NULL ) {
        std::cerr << "Error!  data == NULL" << std::endl;
        result = false;
      }
      num_nodes_on_side_of_mesh += bucket.size();
    }
    else {
      // Since this node is not in the side part it should not have a
      // boundary field.
      if ( data != NULL ) {
        std::cerr << "Error!  data != NULL" << std::endl;
        result = false;
      }
    }
  }

  // There should be 22 nodes on the front face of the mesh. This is checking to verify that.
  if ( num_nodes_on_side_of_mesh != 22u ) {
    std::cerr << "Error!  num_nodes_on_side_of_mesh == " << num_nodes_on_side_of_mesh << " !=  22!" << std::endl;
    result = false;
  }
  return result;
}

// xCheck on pressure and velocity on each node
template< class ElementTopology ,
          class StencilTopology ,
          class PressureField ,
          class VelocityField >
bool verify_pressure_velocity_stencil(
  const BulkData & M ,
  const Part     & element_part ,
  const Part     & linear_node_part ,
  const PressureField  & pressure ,
  const VelocityField  & velocity,
  EntityRank             element_rank )
{
  typedef typename FieldTraits< PressureField >::data_type scalar_type ;
  typedef typename FieldTraits< VelocityField >::data_type scalar_type_2 ;

  bool result = true;

  StaticAssert< SameType< scalar_type , scalar_type_2 >::value >::ok();

  // Verify element topoloy and stencil topology are compatible
  if ( (int) ElementTopology::dimension != StencilTopology::dimension ) {
    std::cerr << "Error!  ElementTopology::dimension != StencilTopology::dimension!" << std::endl;
    result = false;
  }
  if ( (int) ElementTopology::vertex_count != StencilTopology::vertex_count ) {
    std::cerr << "Error!  ElementTopology::vertex_nodes != StencilTopology::vertex_count!" << std::endl;
    result = false;
  }
  if ( (int) ElementTopology::edge_count != StencilTopology::edge_count ) {
    std::cerr << "Error!  ElementTopology::edge_count != StencilTopology::edge_count!" << std::endl;
    result = false;
  }
  if ( (int) ElementTopology::face_count != StencilTopology::face_count ) {
    std::cerr << "Error!  ElementTopology::face_count != StencilTopology::face_count!" << std::endl;
    result = false;
  }

  const std::vector<Bucket*> & buckets = M.buckets( element_rank );

  // Iterate over all element buckets
  for ( std::vector<Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {
    Bucket & bucket = **k ;

    // Ignoring buckets which do not contain element part.
    if ( has_superset( bucket, element_part ) ) {

      // Iterate over all entities in bucket
      for ( Bucket::iterator
            i = bucket.begin() ; i != bucket.end() ; ++i ) {
        Entity & elem = *i ;

        PairIterRelation rel = elem.relations( fem::FEMMetaData::NODE_RANK );

        // Check that the nodal relation size is the same as the element topology
        // node count.
        if ( (unsigned) rel.size() != (unsigned) ElementTopology::node_count ) {
          std::cerr << "Error!" << std::endl;
          result = false;
        }

        // Iterate over every node and check pressure and velocity at each node.
        for ( unsigned j = 0 ; j < ElementTopology::node_count ; ++j ) {
          Entity & node = * rel[j].entity();

          // Get the parts that this node belongs to.
          const Bucket & node_bucket = node.bucket();
          PartVector node_parts ;
          node_bucket.supersets( node_parts );

          scalar_type * const p = field_data( pressure , node );
          scalar_type * const v = field_data( velocity , node );

          // Check if this node is within StencilTopology::node_count, (Hex - 8 and Wedge - 6)
          // that the node_parts contains the linear_node_part and that
          // the pressure is not NULL. This code relies on the fact that the vertex nodes come first
          // in the relation PairIter.
          if ( j < StencilTopology::node_count ) {

            if ( !contain( node_parts , linear_node_part ) ) {
              std::cerr << "Error! PartVector node_parts does not contain linear_node_part"
                << std::endl;
              result = false;
            }
            // Pressure should not be NULL on the vertices.
            if ( p == NULL ) {
              std::cerr << "Error! Pressure is NULL" << std::endl;
              result = false;
            }
          }
          // If this node is NOT within StencilTopology::node_count,
          // check that node_parts does NOT contain the linear_node_part and that
          // the pressure IS NULL.
          else {
            if ( contain( node_parts , linear_node_part ) ) {
              std::cerr << "Error! PartVector node_parts contains linear_node_part"
                << std::endl;
              result = false;
            }
            if ( p != NULL ) {
              std::cerr << "Error! Pressure is not NULL" << std::endl;
              result = false;
            }
          }

          // Check velocity is not NULL. This should be present on every node.
          if ( v == NULL ) {
            std::cerr << "Error! Velocity is NULL" << std::endl;
            result = false;
          }
        }
      }
    }
  }
  return result;
}

}

//------------------------------------------------------------------------------

bool verifyMesh( const UseCase_4_Mesh & mesh )
{
  bool result = true;
  const BulkData& bulk_data = mesh.m_bulkData ;

  std::vector<Bucket *> element_buckets = bulk_data.buckets( mesh.m_elem_rank );

  // Verify the element node coordinates and side nodes:
  // block_hex27:
  Part & block_hex27 = mesh.m_block_hex27 ;
  const ElementNodePointerFieldType & elem_node_coord = mesh.m_element_node_coordinates_field ;
  const VectorFieldType & node_coord = mesh.m_coordinates_field ;
  result = result &&
    verify_elem_node_coord_by_part(
        block_hex27,
        element_buckets,
        elem_node_coord,
        node_coord,
        27
        );

  // Verify element side node:
  const std::vector<Bucket *> face_buckets = bulk_data.buckets( mesh.m_side_rank );
  Part & side_part = mesh.m_side_part ;
  {
    Selector selector = block_hex27 & side_part;
    std::vector<Entity *> entities;
    get_selected_entities( selector, face_buckets, entities);
    for (unsigned i=0 ; i < entities.size() ; ++i) {
      result = result &&
        verify_elem_side_node< shards::Hexahedron<27> >(
          hex_node_ids[i], 0, *entities[i]
          );
    }
  }

  // block_wedge18:
  Part & block_wedge18 = mesh.m_block_wedge18 ;
  result = result &&
    verify_elem_node_coord_by_part(
        block_wedge18,
        element_buckets,
        elem_node_coord,
        node_coord,
        18
        );

  // Verify element side node:
  {
    Selector selector = block_wedge18 & side_part;
    std::vector<Entity *> entities;
    get_selected_entities( selector, face_buckets, entities);
    for (unsigned i=0 ; i < entities.size() ; ++i) {
      result = result &&
        verify_elem_side_node< shards::Wedge<18> >(
          wedge_node_ids[i], 4, *entities[i]
          );
    }
  }

  // Verify centroid dimensions
  const VectorFieldType & centroid_field = mesh.m_centroid_field ;
  result = result &&
    centroid_algorithm_unit_test_dimensions< shards::Hexahedron<27> >(
        bulk_data , centroid_field , elem_node_coord , block_hex27, mesh.m_elem_rank );

  result = result &&
    centroid_algorithm_unit_test_dimensions< shards::Wedge<18> >(
        bulk_data , centroid_field , elem_node_coord , block_wedge18, mesh.m_elem_rank );

  // Verify boundary field data
  const VectorFieldType & boundary_field = mesh.m_boundary_field ;
  result = result &&
    verify_boundary_field_data( bulk_data ,
                              side_part ,
                              boundary_field );

  // Verify pressure velocity stencil for block_hex27
  Part & part_vertex_nodes = mesh.m_part_vertex_nodes ;
  const ScalarFieldType & pressure_field = mesh.m_pressure_field ;
  const VectorFieldType & velocity_field = mesh.m_velocity_field ;
  result = result &&
    verify_pressure_velocity_stencil
    < shards::Hexahedron<27> , shards::Hexahedron<8>  >
    ( bulk_data , block_hex27 , part_vertex_nodes ,
      pressure_field , velocity_field, mesh.m_elem_rank );

  // Verify pressure velocity stencil for block_wedge18
  result = result &&
    verify_pressure_velocity_stencil
    < shards::Wedge<18> , shards::Wedge<6>  >
    ( bulk_data , block_wedge18 , part_vertex_nodes ,
      pressure_field , velocity_field, mesh.m_elem_rank );

  return result;
}

} //namespace use_cases
} //namespace mesh
} //namespace stk_classic
