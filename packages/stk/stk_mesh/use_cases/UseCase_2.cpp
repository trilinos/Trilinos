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

#include <use_cases/UseCase_2.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>

using stk::mesh::fem::NODE_RANK;

//----------------------------------------------------------------------

namespace stk{
namespace mesh {
namespace use_cases {

// Functions to generate the use case mesh information:
namespace {
void usecase_2_elem_node_ids( stk::mesh::EntityId elem_id ,
                              stk::mesh::EntityId node_ids[] );
void usecase_2_node_coordinates( stk::mesh::EntityId node_id ,
                                 double coord[] );

const double TEMPERATURE_VALUE_GOLD = 98.6;
const double VOLUME_VAL = 1.0; // A known value.
}

//----------------------------------------------------------------------

enum { field_data_chunk_size = 10 };
enum { SpatialDim = 3 };

typedef shards::Hexahedron<8>  ElementTraits ;

//----------------------------------------------------------------------
// Populate the mesh meta data.
// We have chosen to do it in the use case constructor;
// however, it could be done incrementally outside of
// the constructor.

UseCase_2_Mesh::UseCase_2_Mesh( stk::ParallelMachine comm ) :
  m_metaData( fem::entity_rank_names(SpatialDim) )
  , m_bulkData(  m_metaData , comm , field_data_chunk_size )
  , m_fem(  m_metaData, SpatialDim )
  , m_partLeft(  declare_part<shards::Hexahedron<8> >( m_metaData, "block_left"))
  , m_partRight( declare_part<shards::Hexahedron<8> >( m_metaData, "block_right"))
  , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
  , m_temperature_field( m_metaData.declare_field< ScalarFieldType >( "temperature" ))
  , m_volume_field( m_metaData.declare_field< ScalarFieldType >( "volume" ))
  , m_elem_rank( stk::mesh::fem::element_rank(m_fem) )
  , m_side_rank( stk::mesh::fem::side_rank(m_fem) )
  , m_edge_rank( stk::mesh::fem::edge_rank(m_fem) )
{
  // Put the coordinates and temperature field on all nodes

  stk::mesh::Part & universal = m_metaData.universal_part();
  stk::mesh::put_field( m_coordinates_field , NODE_RANK , universal , SpatialDim );
  stk::mesh::put_field( m_temperature_field, NODE_RANK, universal );

  // Put the volume field on all elements:
  stk::mesh::put_field( m_volume_field , m_elem_rank , universal );

  // Done populating the mesh meta data.
  // Commit the meta data: this locks out changes,
  // verifies consistency of potentially complex meta data relationships,
  // and allows the internal data structures to be optimized
  // for subsquent use by mesh bulk data.

  m_metaData.commit();
}

UseCase_2_Mesh::~UseCase_2_Mesh()
{ }

//----------------------------------------------------------------------
// Populate the mesh bulk data.
// The mesh meta data must be complete and commited before
// the mesh bulk data can be modified.

void UseCase_2_Mesh::populate( unsigned nleft , unsigned nright )
{
  //------------------------------
  { // Generate the elements and nodes

    m_bulkData.modification_begin(); // Begin modifying the mesh

    stk::mesh::EntityId curr_elem_id = 1 ;
    stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];

    // Note declare_element expects a cell topology
    // to have been attached to m_partLeft.

    // Declare nleft elements
    for ( unsigned j = 0 ; j < nleft ; ++j , ++curr_elem_id ) {
      usecase_2_elem_node_ids( curr_elem_id , node_ids );
      stk::mesh::declare_element( m_bulkData, m_partLeft, curr_elem_id, node_ids );
    }

    // Note declare_element expects a cell topology
    // to have been attached to m_partRight.

    // Declare nright elements
    for ( unsigned j = 0 ; j < nright ; ++j , ++curr_elem_id ) {
      usecase_2_elem_node_ids( curr_elem_id , node_ids );
      stk::mesh::declare_element( m_bulkData, m_partRight, curr_elem_id, node_ids );
    }

    // Done modifying the mesh.
    // Modifications on the local parallel process are communicated
    // among processes, verified for consistency, and changes to
    // parallel shared/ghosted mesh entities are synchronized.
    m_bulkData.modification_end();
  }
  //------------------------------
  { // Assign Nodal Field Data
    // The following operations are entirely parallel local
    // and can be parallel inconsistent.  A user may communicate
    // field data as needed to force consistency.


    // Iterate over all nodes by getting all node buckets and iterating over each bucket.
    const std::vector<stk::mesh::Bucket*> & node_buckets =
      m_bulkData.buckets( NODE_RANK );

    for ( std::vector<stk::mesh::Bucket*>::const_iterator
          node_bucket_it = node_buckets.begin() ;
          node_bucket_it != node_buckets.end() ; ++node_bucket_it ) {

      const stk::mesh::Bucket & bucket = **node_bucket_it;

      // Fill the nodal coordinates.
      // Create a multidimensional array view of the
      // nodal coordinates field data for this bucket of nodes.
      // The array is two dimensional ( Cartesian X NumberNodes )
      // and indexed by ( 0..2 , 0..NumberNodes-1 )

      stk::mesh::BucketArray<VectorFieldType>
        coordinates_array( m_coordinates_field, bucket );

      const int num_nodes_in_bucket = coordinates_array.dimension(1);

      ThrowRequireMsg( coordinates_array.dimension(0) == SpatialDim , "Expected coordinates_array.dimension(0) "
                       << coordinates_array.dimension(0) << " to equal SpatialDim " << SpatialDim );

      // For each node in the bucket populate its nodal coordinates.
      for ( int i=0 ; i < num_nodes_in_bucket ; ++i ) {
        const unsigned node_id = bucket[i].identifier();
        usecase_2_node_coordinates( node_id, & coordinates_array(0,i) );
      }

      // Fill the nodal temperature field.
      // Create a multidimensional array view of the
      // nodal temperature field data for this bucket of nodes.
      // The array is one dimensional ( NumberNodes )
      // and indexed by ( 0..NumberNodes-1 )

      stk::mesh::BucketArray<ScalarFieldType>
        temperature_array( m_temperature_field, bucket );

      const int num_temps = temperature_array.dimension(0);

      ThrowRequireMsg( num_nodes_in_bucket == num_temps , "Expected num_temps "
                       << num_temps << " to equal num_nodes_in_bucket " << num_nodes_in_bucket );

      // For each node in the bucket assign the temperature field to a constant.
      for ( int i=0 ; i < num_nodes_in_bucket ; ++i) {
        temperature_array(i) = TEMPERATURE_VALUE_GOLD ;
      }
    }
  }
  //------------------------------
  { // Assign Element Field Data
    // The following operations are entirely parallel local
    // and can be parallel inconsistent.  A user may communicate
    // field data as needed to force consistency.

    const std::vector<stk::mesh::Bucket*> & elem_buckets =
      m_bulkData.buckets( m_elem_rank );

    // Volume field:

    // Iterate over all element buckets and populate volume field.
    for ( std::vector<stk::mesh::Bucket*>::const_iterator
          element_bucket_it = elem_buckets.begin();
          element_bucket_it != elem_buckets.end() ; ++element_bucket_it ) {

      Bucket & element_bucket = **element_bucket_it;

      // Fill the element volume field.
      // Create a multidimensional array view of the
      // element volume field data for this bucket of elements.
      // The array is one dimensional ( NumberElements )
      // and indexed by ( 0..NumberElements-1 )

      stk::mesh::BucketArray<ScalarFieldType>
        volume_array( m_volume_field, element_bucket );

      const unsigned num_elements = volume_array.dimension(0);

      ThrowRequireMsg( element_bucket.size() == num_elements , "Expected element_bucket.size() "
                       << element_bucket.size() << " to equal num_elements " << num_elements );

      // Populate volume field
      for ( unsigned volume_index=0 ; volume_index < num_elements ; ++volume_index) {
        volume_array(volume_index) = VOLUME_VAL;
      }
    }
  }
}

//------------------------------------------------------------------------------

bool verifyMesh( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright )
{
  bool result = false;
  result = verifyCellTopology(mesh);
  result = result && verifyEntityCounts(mesh,nleft,nright);
  result = result && verifyRelations(mesh,nleft,nright);
  result = result && verifyFields(mesh);
  return result;
}


bool verifyCellTopology( const UseCase_2_Mesh & mesh )
{
  // Verify that cell topology for left and right parts is a hex

  stk::mesh::Part & partLeft = mesh.m_partLeft ;
  stk::mesh::Part & partRight = mesh.m_partRight ;

  bool result = true;
  const CellTopologyData * left_cell_topology = stk::mesh::fem::get_cell_topology( partLeft ).getCellTopologyData();
  if (left_cell_topology != shards::getCellTopologyData< shards::Hexahedron<8> >()) {
    std::cerr << "Error, the left cell topology is not what we asked for!" << std::endl;
    std::cerr << "It is " << left_cell_topology->name<< " rather than a hex" << std::endl;
    result = false;
  }

  const CellTopologyData * right_cell_topology = stk::mesh::fem::get_cell_topology( partRight ).getCellTopologyData();
  if (right_cell_topology != shards::getCellTopologyData< shards::Hexahedron<8> >()) {
    std::cerr << "Error, the right cell topology is not what we asked for!" << std::endl;
    std::cerr << "It is " << right_cell_topology->name<< " rather than a hex" << std::endl;
    result = false;
  }
  return result;
}

bool verifyEntityCounts( const UseCase_2_Mesh & mesh,
                         unsigned nleft, unsigned nright )
{
  // Verify that we have the expected numbers of entities in each part

  const stk::mesh::BulkData & bulkData = mesh.m_bulkData ;
  stk::mesh::Part & partLeft = mesh.m_partLeft ;
  stk::mesh::Part & partRight = mesh.m_partRight ;

  // Compute expected number of entities for each entity rank
  const unsigned expected_num_left_nodes = (nleft+1)*4;
  const unsigned expected_num_right_nodes = (nright+1)*4;
  const unsigned expected_num_edges = 0;
  const unsigned expected_num_faces = 0;

  bool result = true;
  std::vector<unsigned> entity_counts;

  // Create a 'Selector' to select mesh entities
  // (nodes, edges, faces, elements)
  // that are members of the left element block..
  // Use the selector to count those entities.
  stk::mesh::Selector selector_left(partLeft);
  stk::mesh::count_entities( selector_left, bulkData , entity_counts );

  // Verify expected number of entities for left part
  if ( entity_counts[NODE_RANK]         != expected_num_left_nodes ||
       entity_counts[mesh.m_edge_rank] != expected_num_edges ||
       entity_counts[mesh.m_side_rank] != expected_num_faces ||
       entity_counts[mesh.m_elem_rank] != nleft ) {
    std::cerr<< "Error, the left entity counts are incorrect!" << std::endl;
    result = false;
  }

  // Create a 'Selector' to select mesh entities
  // (nodes, edges, faces, elements)
  // that are members of the right element block..
  // Use the selector to count those entities.
  stk::mesh::Selector selector_right(partRight);
  stk::mesh::count_entities( selector_right, bulkData , entity_counts );

  // Verify expected number of entities for right part
  if ( entity_counts[NODE_RANK]         != expected_num_right_nodes ||
       entity_counts[mesh.m_edge_rank] != expected_num_edges ||
       entity_counts[mesh.m_side_rank] != expected_num_faces ||
       entity_counts[mesh.m_elem_rank] != nright ) {
    std::cerr<< "Error, the right counts are incorrect!" << std::endl;
    result = false;
  }
  return result;
}

bool verifyRelations( const UseCase_2_Mesh & mesh,
                      unsigned nleft, unsigned nright )
{
  const stk::mesh::BulkData & bulkData = mesh.m_bulkData ;
  stk::mesh::Part & partLeft  = mesh.m_partLeft ;
  stk::mesh::Part & partRight = mesh.m_partRight ;
  std::vector<stk::mesh::Part*> both_parts;
  both_parts.push_back(&partLeft);
  both_parts.push_back(&partRight);

  bool result = true;

  const std::vector<stk::mesh::Bucket*> & all_elem_buckets = bulkData.buckets( mesh.m_elem_rank );

  // For both left and right parts
  // Verify the element-node relationships for each part:
  for( std::vector<stk::mesh::Part*>::const_iterator iter_part = both_parts.begin();  iter_part < both_parts.end(); iter_part++)
  {
    std::vector<stk::mesh::Bucket*> selected_elem_buckets ;

    // From all element buckets select part buckets
    stk::mesh::Selector selector(**iter_part);
    stk::mesh::get_buckets( selector , all_elem_buckets ,
                                            selected_elem_buckets );

    // Iterate over part buckets
    for ( size_t j = 0 ; j < selected_elem_buckets.size() ; ++j ) {

      const stk::mesh::Bucket & elem_bucket = *selected_elem_buckets[j] ;

      // Iterate over all elements in this bucket
      for ( size_t i = 0 ; i < elem_bucket.size() ; ++i ) {
        stk::mesh::Entity & elem = elem_bucket[i] ;

        // Query the node ids for this element.
        stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];
        usecase_2_elem_node_ids( elem.identifier() , node_ids );

        // Pair of iterators for all of the element's relations.
        // This class has convenience functions for size and indexing.
        //   rel.size() == std::distance( rel.first , rel.second );
        //   rel[i] == *( rel.first + i );
        stk::mesh::PairIterRelation rel = elem.relations();

        // Verify that the number of nodes in this element is correct.
        if( shards::Hexahedron<8> ::node_count != rel.size() ) {
          std::cerr << "Error, number of relations is incorrect! It is " << rel.size()
                    << std::endl;
          result = false;
        }

        // Verify the nodes of this element
        // have the correct relation-identifiers and
        // are members of the block.
        for ( unsigned k = 0 ; k < shards::Hexahedron<8> ::node_count ; ++k ) {
          stk::mesh::Entity & rel_node = * rel[k].entity();
          if ( node_ids[k] != rel_node.identifier() ||
               ! rel_node.bucket().member(**iter_part) ) {
              std::cerr << "Error, an element's node is just plain wrong!"
                        << std::endl;
              result = false;
          }
        }
      }
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// An example of a template function for gathering (copying)
// field data from the nodes of an element.
template< class field_type >
bool gather_field_data( unsigned expected_num_rel, const field_type & field ,
                        const stk::mesh::Entity & entity ,
                        typename stk::mesh::FieldTraits< field_type >::data_type * dst,
                        stk::mesh::EntityRank entity_rank )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type T ;

  stk::mesh::PairIterRelation rel = entity.relations( entity_rank );

  bool result = expected_num_rel == (unsigned) rel.size();

  if ( result ) {
    // Iterate over field data for each related entity and copy data
    // into src for one entity at a time
    T * const dst_end = dst + SpatialDim * expected_num_rel ;
    for ( const T * src ;
          ( dst < dst_end ) &&
          ( src = field_data( field , * rel->entity() ) ) ;
          ++rel , dst += SpatialDim ) {
      stk::Copy<SpatialDim>( dst , src );
    }
    result = dst == dst_end ;
  }
  return result ;
}

namespace {
// This is based on the local node numbering mapping to relative coordinates.
// See header file for local node numbering.
const double ELEMENT_LOCAL_COORDINATES_GOLD[ ElementTraits::node_count ][ SpatialDim ] =
  { { 0, 0, 0 } ,
    { 1, 0, 0 } ,
    { 1, 1, 0 } ,
    { 0, 1, 0 } ,
    { 0, 0, 1 } ,
    { 1, 0, 1 } ,
    { 1, 1, 1 } ,
    { 0, 1, 1 }   };

}

bool verifyFields( const UseCase_2_Mesh & mesh )
{
  bool result = true;

  const VectorFieldType & coordinates_field = mesh.m_coordinates_field ;
  const ScalarFieldType & volume_field = mesh.m_volume_field ;
  const stk::mesh::BulkData & bulkData = mesh.m_bulkData ;

  // All element buckets:
  const std::vector<stk::mesh::Bucket*> & elem_buckets =
    bulkData.buckets( mesh.m_elem_rank );

  // Verify coordinates_field by gathering the nodal coordinates
  // from each element's nodes.

  // Iterate over all element buckets
  for ( std::vector<stk::mesh::Bucket*>::const_iterator
        element_bucket_it = elem_buckets.begin();
        element_bucket_it != elem_buckets.end() ; ++element_bucket_it ) {

    const stk::mesh::Bucket& bucket = **element_bucket_it;
    const size_t num_elements_in_bucket = bucket.size();

    // Iterate over all elements in bucket
    for( size_t bucket_index = 0; bucket_index < num_elements_in_bucket; ++bucket_index) {
      const stk::mesh::Entity & elem = bucket[bucket_index] ;

      // Gather nodal coordinates of this element
      double elem_coord[ ElementTraits::node_count ][ SpatialDim ];
      const bool gather_result =
        gather_field_data
        ( ElementTraits::node_count, coordinates_field , elem , & elem_coord[0][0], NODE_RANK );

      if ( gather_result == false ) {
        std::cerr << "Error!" << std::endl;
        result = false;
      }

      // Take the global nodal coordinates of the nodes in this element and
      // and convert them to local relative coordinates

      // Compute nodal coordinates of first node in element
      double base[3] ;
      base[0] = elem_coord[0][0] ;
      base[1] = elem_coord[0][1] ;
      base[2] = elem_coord[0][2] ;

      // Convert the global coordinates for each node by subtracting the
      // coordinates of the base node
      for ( unsigned j = 0 ; j < ElementTraits::node_count ; ++j ) {
        elem_coord[j][0] -= base[0] ;
        elem_coord[j][1] -= base[1] ;
        elem_coord[j][2] -= base[2] ;
      }

      // Check that computed local nodal coordinates match gold local
      // coordinates.
      for (int node_index=0 ; node_index<ElementTraits::node_count ; ++node_index ) {
        for (int coord_index=0 ; coord_index<SpatialDim ; ++coord_index) {
          double local_elem_coord = elem_coord[node_index][coord_index];
          double gold_local_elem_coord = ELEMENT_LOCAL_COORDINATES_GOLD[node_index][coord_index];
          if ( local_elem_coord != gold_local_elem_coord ) {
            std::cerr << "Error!  elem_coord[" << node_index << "]"
              << "[" << coord_index << "] == " << local_elem_coord
              << " != " << gold_local_elem_coord
              << " == ELEMENT_LOCAL_COORDINATES_GOLD[" << node_index << "]"
              << "[" << coord_index << "]!" << std::endl;
            result = false;
          }
        }
      }
    }

    // Verify volume_field

    stk::mesh::BucketArray<ScalarFieldType> volume_array(
        volume_field,
        **element_bucket_it
        );

    // For all elements volume field should be VOLUME_VAL
    for ( unsigned volume_index=0 ; volume_index < num_elements_in_bucket ; ++volume_index) {
      if ( volume_array(volume_index) != VOLUME_VAL ) {
        std::cerr << "Error!  volume_array(" << volume_index << ") == "
          << volume_array(volume_index) << " != " << VOLUME_VAL
          << " == VOLUME_VAL " << std::endl;
        result = false;
      }
    }
  }

  {
    // Verify temperature_field on the nodes

    const ScalarFieldType & temperature_field = mesh.m_temperature_field ;

    // Get all node buckets
    const std::vector<stk::mesh::Bucket*> & node_buckets =
      bulkData.buckets( NODE_RANK );

    // Iterate over all node buckets
    for ( std::vector<stk::mesh::Bucket*>::const_iterator
          node_bucket_it = node_buckets.begin();
          node_bucket_it != node_buckets.end() ; ++node_bucket_it ) {
      const stk::mesh::Bucket & bucket = **node_bucket_it;
      stk::mesh::BucketArray<ScalarFieldType> temperature_array(
          temperature_field,
          bucket
          );

      // For all nodes in bucket temperature field should be
      // TEMPERATURE_VALUE_GOLD
      int num_nodes_in_bucket = temperature_array.dimension(0);
      for ( int i=0 ; i < num_nodes_in_bucket ; ++i) {
        if (temperature_array(i) != TEMPERATURE_VALUE_GOLD) {
          std::cerr << "Error!  temperature_array("<<i<<") == "
            << temperature_array(i) << " != " << TEMPERATURE_VALUE_GOLD
            << " == TEMPERATURE_VALUE_GOLD" << std::endl;
          result = false;
        }
      }
    }
  }
  return result;
}

//------------------------------------------------------------------------------

namespace {

//Given an element id compute the ids of the associated nodes.
void usecase_2_elem_node_ids( stk::mesh::EntityId elem_id ,
                              stk::mesh::EntityId node_ids[] )
{
  ThrowRequireMsg( elem_id != 0,
                   "usecase_2_elem_node_ids: ERROR, elem_id ("
                   << elem_id << ") must be greater than 0.");

  const unsigned nodes_per_side = 4;
  const unsigned base = ( elem_id - 1 ) * nodes_per_side ;
  node_ids[0] = base + 1 ;
  node_ids[1] = base + 5 ;
  node_ids[2] = base + 6 ;
  node_ids[3] = base + 2 ;
  node_ids[4] = base + 4 ;
  node_ids[5] = base + 8 ;
  node_ids[6] = base + 7 ;
  node_ids[7] = base + 3 ;
}

// Given a node_id compute its spatial coordinates
void usecase_2_node_coordinates( stk::mesh::EntityId node_id , double coord[] )
{
  // i_length is the same as the number of the side it occurs in
  // i_plane is the position of the node in the side
  const unsigned i_length = ( node_id - 1 ) / 4 ;
  const unsigned i_plane  = ( node_id - 1 ) % 4 ;

  coord[0] = i_length ;
  coord[1] = i_plane == 1 || i_plane == 2 ? 1.0 : 0.0 ;
  coord[2] = i_plane == 2 || i_plane == 3 ? 1.0 : 0.0 ;
}

}

} //namespace use_cases
} //namespace mesh
} //namespace stk


