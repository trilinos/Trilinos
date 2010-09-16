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

#include <use_cases/UseCase_3.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FEMTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

//----------------------------------------------------------------------

namespace stk{
namespace mesh {
namespace use_cases {

UseCase_3_Mesh::UseCase_3_Mesh( stk::ParallelMachine comm )
  : m_spatial_dimension(3)
  , m_metaData( TopologicalMetaData::entity_rank_names(m_spatial_dimension) )
  , m_bulkData( m_metaData , comm )
  , m_topData( m_metaData, m_spatial_dimension )
  , m_block_hex(        m_topData.declare_part<shards::Hexahedron<8> >( "block_1" ))
  , m_block_wedge(      m_topData.declare_part<shards::Wedge<6> >( "block_2" ))
  , m_block_tet(        m_topData.declare_part<shards::Tetrahedron<4> >( "block_3" ))
  , m_block_pyramid(    m_topData.declare_part<shards::Pyramid<5> >( "block_4" ))
  , m_block_quad_shell( m_topData.declare_part<shards::ShellQuadrilateral<4> >( "block_5" ))
  , m_block_tri_shell(  m_topData.declare_part<shards::ShellTriangle<3> >( "block_6" ))
  , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
  , m_centroid_field(    m_metaData.declare_field< VectorFieldType >( "centroid" ))
  , m_temperature_field( m_metaData.declare_field< ScalarFieldType >( "temperature" ))
  , m_volume_field( m_metaData.declare_field< ScalarFieldType >( "volume" ))
  , m_element_node_coordinates_field( m_metaData.declare_field< ElementNodePointerFieldType >( "elem_node_coord" ))
{

  // Define where fields exist on the mesh:
  Part & universal = m_metaData.universal_part();

  put_field( m_coordinates_field , m_topData.node_rank , universal );
  put_field( m_centroid_field , m_topData.element_rank , universal );
  put_field( m_temperature_field, m_topData.node_rank, universal );
  put_field( m_volume_field, m_topData.element_rank, m_block_hex );
  put_field( m_volume_field, m_topData.element_rank, m_block_wedge );
  put_field( m_volume_field, m_topData.element_rank, m_block_tet );
  put_field( m_volume_field, m_topData.element_rank, m_block_pyramid );

  // Define the field-relation such that the values of the
  // 'element_node_coordinates_field' are pointers to the
  // element's nodal 'coordinates_field'.
  // I.e., let:
  //   double *const* elem_node_coord =
  //     field_data( m_element_node_coordinates_field , element );
  // then
  //     elem_node_coord[n][0..2] is the coordinates of element node 'n'
  //     that are attached to that node.

  m_metaData.declare_field_relation(
    m_element_node_coordinates_field ,
    & element_node_stencil<void> ,
    m_coordinates_field
    );

  put_field( m_element_node_coordinates_field, m_topData.element_rank, m_block_hex, shards::Hexahedron<> ::node_count );
  put_field( m_element_node_coordinates_field, m_topData.element_rank, m_block_wedge, shards::Wedge<> ::node_count );
  put_field( m_element_node_coordinates_field, m_topData.element_rank, m_block_tet, shards::Tetrahedron<> ::node_count );
  put_field( m_element_node_coordinates_field, m_topData.element_rank, m_block_pyramid, shards::Pyramid<> ::node_count );
  put_field( m_element_node_coordinates_field, m_topData.element_rank, m_block_quad_shell, shards::ShellQuadrilateral<> ::node_count);
  put_field( m_element_node_coordinates_field, m_topData.element_rank, m_block_tri_shell, shards::ShellTriangle<> ::node_count );

  m_metaData.commit();
}

UseCase_3_Mesh::~UseCase_3_Mesh()
{ }

//------------------------------------------------------------------------------
// Use case specific mesh generation data:

enum { SpatialDim = 3 };
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

static const EntityId hex_node_ids[3][ shards::Hexahedron<> ::node_count ] = {
  { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
  { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
  { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

static const EntityId wedge_node_ids[3][ shards::Wedge<> ::node_count ] = {
  { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
  { 10 ,  9 ,  6 , 20 , 19 , 16 } ,
  { 16 , 17 , 20 ,  6 ,  7 , 10 } };

static const EntityId tetra_node_ids[3][ shards::Tetrahedron<> ::node_count ] = {
  { 15 , 19 , 16 , 21 } ,
  { 19 , 20 , 16 , 21 } ,
  { 16 , 20 , 17 , 21 } };

static const EntityId pyramid_node_ids[2][ shards::Pyramid<> ::node_count ] = {
  { 11 , 15 , 16 , 12 , 21 } ,
  { 12 , 16 , 17 , 13 , 21 } };

static const EntityId shell_quad_node_ids[3][ shards::ShellQuadrilateral<> ::node_count ]={
  { 9 , 6 , 16 , 19 } ,
  { 6 , 7 , 17 , 16 } ,
  { 7 , 8 , 18 , 17 } };

static const EntityId shell_tri_node_ids[3][ shards::ShellTriangle<> ::node_count ] ={
  { 19 , 16 , 21 } ,
  { 16 , 17 , 21 } ,
  { 17 , 13 , 21 } };

}

//------------------------------------------------------------------------------

void UseCase_3_Mesh::populate()
{
  m_bulkData.modification_begin();

  EntityId elem_id = 1;

  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id ) {
    declare_element( m_bulkData, m_block_hex, elem_id, hex_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id ) {
    declare_element( m_bulkData, m_block_wedge, elem_id, wedge_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_tetra ; ++i , ++elem_id ) {
    declare_element( m_bulkData, m_block_tet, elem_id, tetra_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++elem_id ) {
    declare_element( m_bulkData, m_block_pyramid, elem_id, pyramid_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++elem_id ) {
    declare_element( m_bulkData, m_block_quad_shell, elem_id, shell_quad_node_ids[i]);
  }

  for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++elem_id ) {
    declare_element( m_bulkData, m_block_tri_shell, elem_id, shell_tri_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < node_count ; ++i ) {
    Entity * const node = m_bulkData.get_entity( m_topData.node_rank , i + 1 );
    double * const coord = field_data( m_coordinates_field , *node );
    coord[0] = node_coord_data[i][0] ;
    coord[1] = node_coord_data[i][1] ;
    coord[2] = node_coord_data[i][2] ;
  }

  m_bulkData.modification_end();
  // No parallel stuff for now
}

//------------------------------------------------------------------------------

bool verify_elem_node_coord_3(
  mesh::Entity & elem ,
  const ElementNodePointerFieldType & elem_node_coord ,
  const VectorFieldType & node_coord ,
  const unsigned node_count )
{
  bool result = true;
  mesh::PairIterRelation rel = elem.relations( NodeRank );

  if( (unsigned) rel.size() != node_count ) {
    std::cerr << "Error!  relation size == " << rel.size() << " != "
      << node_count << " == node count" << std::endl;
    result = false;
  }

  mesh::EntityArray< ElementNodePointerFieldType >
    elem_node_array( elem_node_coord , elem );

  {
    const unsigned n1 = elem_node_array.dimension<0>();
    if( n1 != node_count ) {
      std::cerr << "Error!  element node array dimension<0> == " << n1 << " != "
        << node_count << " == node count" << std::endl;
      result = false;
    }
    if ( (unsigned) elem_node_array.size() != node_count ) {
      std::cerr << "Error!  element node array size == "
        << elem_node_array.size() << " != " << node_count << " == node count"
        << std::endl;
      result = false;
    }
  }

  double * const * const elem_data = elem_node_array.contiguous_data();

  for ( unsigned j = 0 ; j < node_count ; ++j ) {
    mesh::Entity & node = * rel[j].entity();

    mesh::EntityArray< VectorFieldType > node_coord_array( node_coord , node );

    {
      const unsigned n1 = node_coord_array.dimension<0>();
      if( n1 != (unsigned) SpatialDim ) {
        std::cerr << "Error!  node coord array dimension<0> == " << n1 << " != "
          << SpatialDim << " == SpatialDim" << std::endl;
        result = false;
      }
      if( node_coord_array.size() != SpatialDim ) {
        std::cerr << "Error!  node coord array size == "
          << node_coord_array.size() << " != " << SpatialDim
          << " == SpatialDim" << std::endl;
        result = false;
      }
    }

    double * const node_data = node_coord_array.contiguous_data();
    if( elem_data[j] != node_data ) {
      std::cerr << "Error!  elem_data[" << j << "] == " << elem_data[j]
        << " != " << node_data << " node_data" << std::endl;
      result = false;
    }
  }
  return result;
}


bool verify_elem_node_coord_by_part_3(
    Part & part,
    const std::vector<Bucket *> & bucket_vector,
    const ElementNodePointerFieldType & elem_node_coord,
    const VectorFieldType & node_coord,
    const unsigned node_count )
{
  Selector selector(part);
  std::vector<Entity *> entities;
  get_selected_entities( selector, bucket_vector, entities);
  std::vector<Entity *>::iterator entity_it = entities.begin();
  bool result = true;
  for ( ; entity_it != entities.end() ; ++entity_it ) {
    result = result &&
      verify_elem_node_coord_3(
          **entity_it , elem_node_coord , node_coord , node_count
          );
  }
  return result;
}


bool verifyMesh( const UseCase_3_Mesh & mesh )
{
  bool result = true;

  const BulkData & bulkData = mesh.m_bulkData ;
  const VectorFieldType & node_coord = mesh.m_coordinates_field ;
  const ElementNodePointerFieldType & elem_node_coord  =
    mesh.m_element_node_coordinates_field ;

  std::vector<Bucket *> element_buckets = bulkData.buckets( mesh.m_topData.element_rank );

  // Verify entities in each part are set up correcty:

  // hex_block:
  Part & hex_block = mesh.m_block_hex ;
  result = result &&
    verify_elem_node_coord_by_part_3(
      hex_block,
      element_buckets,
      elem_node_coord,
      node_coord,
      shards::Hexahedron<> ::node_count
      );

  // wedge_block:
  Part & wedge_block = mesh.m_block_wedge ;
  result = result &&
    verify_elem_node_coord_by_part_3(
      wedge_block,
      element_buckets,
      elem_node_coord,
      node_coord,
      shards::Wedge<> ::node_count
      );

  // tetra_block
  Part & tetra_block = mesh.m_block_tet ;
  result = result &&
    verify_elem_node_coord_by_part_3(
      tetra_block,
      element_buckets,
      elem_node_coord,
      node_coord,
      shards::Tetrahedron<> ::node_count
      );

  // pyramid_block
  Part & pyramid_block = mesh.m_block_pyramid ;
  result = result &&
    verify_elem_node_coord_by_part_3(
      pyramid_block,
      element_buckets,
      elem_node_coord,
      node_coord,
      shards::Pyramid<> ::node_count
      );

  // quad_shell_block
  Part & quad_shell_block = mesh.m_block_quad_shell ;
  result = result &&
    verify_elem_node_coord_by_part_3(
      quad_shell_block,
      element_buckets,
      elem_node_coord,
      node_coord,
      shards::ShellQuadrilateral<> ::node_count
      );

  // tri_shell_block
  Part & tri_shell_block = mesh.m_block_tri_shell ;
  result = result &&
    verify_elem_node_coord_by_part_3(
      tri_shell_block,
      element_buckets,
      elem_node_coord,
      node_coord,
      shards::ShellTriangle<> ::node_count
      );

  // Check that all the nodes were allocated.
  for ( unsigned i = 0 ; i < node_count ; ++i ) {
    Entity * const node = bulkData.get_entity( mesh.m_topData.node_rank , i + 1 );
    if ( node == NULL ) {
      std::cerr << "Error!  Invalid null pointer for node returned from "
        << "bulkData.get_entity( m_topData.node_rank, " << i+1 << " ) " << std::endl;
      result = false;
    }
  }

  return result;
}

} //namespace use_cases
} //namespace mesh
} //namespace stk
