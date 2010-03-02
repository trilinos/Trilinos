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

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

//----------------------------------------------------------------------

namespace stk{
namespace mesh {
namespace use_cases {

enum {
  Node = 0,
  Edge = 1,
  Face = 2,
  Element = 3
};

enum { SpaceDim = 3 };

typedef shards::Hexahedron<8>  ElementTraits ;

UseCase_2_MetaData::UseCase_2_MetaData(const std::vector<std::string> & entity_type_names)
  : m_metaData(entity_type_names),
    m_partLeft(m_metaData.declare_part( "block_left", Element )),
    m_partRight(m_metaData.declare_part( "block_right", Element )),
    m_coordinates_field(m_metaData.declare_field< VectorFieldType >( "coordinates" )),
    m_temperature_field(m_metaData.declare_field< ScalarFieldType >( "temperature" )),
    m_volume_field(m_metaData.declare_field< ScalarFieldType >( "volume" ))
{
  // Declare the intention of the parts to be for hex 8 elements.
  // This attaches this cell topology to these parts:
  stk::mesh::set_cell_topology< shards::Hexahedron<8>  >( m_partLeft );
  stk::mesh::set_cell_topology< shards::Hexahedron<8>  >( m_partRight );

  // Field restrictions:
  stk::mesh::Part & universal = m_metaData.universal_part();
  stk::mesh::put_field(
    m_coordinates_field , stk::mesh::Node , universal , SpaceDim );

  stk::mesh::put_field( m_temperature_field, stk::mesh::Node, universal );
  stk::mesh::put_field( m_volume_field , stk::mesh::Element , universal );
  
  m_metaData.commit();

}

UseCase_2_Mesh::~UseCase_2_Mesh() 
{
}

std::vector<std::string> get_entity_type_names_2()
{
  // Note:  stk::mesh::fem_entity_type_names provides a default set of these too.
  std::vector<std::string> entity_type_names(4);
  entity_type_names[0] = "Node";
  entity_type_names[1] = "Edge";
  entity_type_names[2] = "Face";
  entity_type_names[3] = "Element";
  return entity_type_names;
}

UseCase_2_Mesh::UseCase_2_Mesh( stk::ParallelMachine comm )
  : UseCase_2_MetaData(get_entity_type_names_2()),
    m_bulkData(m_metaData,comm,field_data_chunk_size)
{
}

void get_elem_node_ids_2( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] )
{
  if ( elem_id == 0 ) {
    std::cerr << "use_case_2, elem_node_ids: ERROR, elem_id ("
        << elem_id << ") must be greater than 0." << std::endl;
    return;
  }

  const unsigned base = ( elem_id - 1 ) * 4 ;
  node_ids[0] = base + 1 ;
  node_ids[1] = base + 5 ;
  node_ids[2] = base + 6 ;
  node_ids[3] = base + 2 ;
  node_ids[4] = base + 4 ;
  node_ids[5] = base + 8 ;
  node_ids[6] = base + 7 ;
  node_ids[7] = base + 3 ;
}

void node_coordinates( unsigned node_id , double coord[] )
{
  const unsigned i_length = ( node_id - 1 ) / 4 ;
  const unsigned i_plane  = ( node_id - 1 ) % 4 ;

  coord[0] = i_length ;
  coord[1] = i_plane == 1 || i_plane == 2 ? 1.0 : 0.0 ;
  coord[2] = i_plane == 2 || i_plane == 3 ? 1.0 : 0.0 ;
}


void populate( UseCase_2_Mesh & mesh , unsigned nleft , unsigned nright )
{
  stk::mesh::BulkData & bulkData = mesh.modifiableBulkData();
  stk::mesh::Part & partLeft = mesh.partLeft();
  stk::mesh::Part & partRight = mesh.partRight();

  // Generate Mesh
  stk::mesh::EntityId elem_id = 1 ;
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];

  for ( unsigned j = 0 ; j < nleft ; ++j , ++elem_id ) {
    get_elem_node_ids_2( elem_id , node_ids );
    // Note declare_element expects a cell topology to have been attached to partLeft.
    stk::mesh::declare_element( bulkData , partLeft , elem_id , node_ids );
  }

  for ( unsigned j = 0 ; j < nright ; ++j , ++elem_id ) {
    get_elem_node_ids_2( elem_id , node_ids );
    // Note declare_element expects a cell topology to have been attached to partRight.
    stk::mesh::declare_element( bulkData , partRight , elem_id , node_ids );
  }
  bulkData.modification_end();

  // Assign Field Data:
  VectorFieldType & coordinates_field = mesh.coordinates_field();
  ScalarFieldType & temperature_field = mesh.temperature_field();
  ScalarFieldType & volume_field = mesh.volume_field();

  const std::vector<stk::mesh::Bucket*> & node_buckets = bulkData.buckets( stk::mesh::Node );

  std::vector<stk::mesh::Bucket*>::const_iterator node_bucket_it = node_buckets.begin();

  for ( ; node_bucket_it != node_buckets.end() ; ++node_bucket_it ) {
    const stk::mesh::Bucket & bucket = **node_bucket_it;
    
    // Coordinates field:
    stk::mesh::BucketArray<VectorFieldType> coordinates_array( coordinates_field, bucket );
    int num_coords = coordinates_array.dimension(1);
    for ( int i=0 ; i < num_coords ; ++i ) {
      const unsigned node_id = bucket[i].identifier();
      node_coordinates( node_id, & coordinates_array(0,i) );
    }

    // Temperature field:
    stk::mesh::BucketArray<ScalarFieldType> temperature_array( temperature_field, bucket );
    double temp_val = 98.6;
    int num_temps = temperature_array.dimension(0);
    for ( int i=0 ; i < num_temps ; ++i) {
      temperature_array(i) = temp_val;
    }

  }

  const std::vector<stk::mesh::Bucket*> & elem_buckets = bulkData.buckets( stk::mesh::Element );

  // Volume field:
  std::vector<stk::mesh::Bucket*>::const_iterator element_bucket_it = elem_buckets.begin();
  const double volume_val = 1.0;
  for ( ; element_bucket_it != elem_buckets.end() ; ++element_bucket_it ) {
    stk::mesh::BucketArray<ScalarFieldType> volume_array( volume_field, **element_bucket_it );
    int num_volumes = volume_array.dimension(0);
    for ( int volume_index=0 ; volume_index < num_volumes ; ++volume_index) {
      volume_array(volume_index) = volume_val;
    }
  }

}


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
  stk::mesh::Part & partLeft = mesh.partLeft();
  stk::mesh::Part & partRight = mesh.partRight();

  bool result = true;
  const CellTopologyData * left_cell_topology = stk::mesh::get_cell_topology( partLeft );
  if (left_cell_topology != shards::getCellTopologyData< shards::Hexahedron<8> >()) {
    std::cerr << "Error, the left cell topology is not what we asked for!" << std::endl;
    result = false;
  }

  const CellTopologyData * right_cell_topology = stk::mesh::get_cell_topology( partRight );
  if (right_cell_topology != shards::getCellTopologyData< shards::Hexahedron<8> >()) {
    std::cerr << "Error, the right cell topology is not what we asked for!" << std::endl;
    result = false;
  }
  return result;
}


bool verifyEntityCounts( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright )
{
  const stk::mesh::BulkData & bulkData = mesh.bulkData();
  stk::mesh::Part & partLeft = mesh.partLeft();
  stk::mesh::Part & partRight = mesh.partRight();

  const unsigned expected_num_left_nodes = (nleft+1)*4;
  const unsigned expected_num_right_nodes = (nright+1)*4;
  const unsigned expected_num_edges = 0;
  const unsigned expected_num_faces = 0;

  bool result = true;
  std::vector<unsigned> entity_counts;
  stk::mesh::Selector selector_left(partLeft);
  stk::mesh::count_entities( selector_left, bulkData , entity_counts );
  if (
    (entity_counts[Node] != expected_num_left_nodes) ||
    (entity_counts[Edge] != expected_num_edges) ||
    (entity_counts[Face] != expected_num_faces) ||
    (entity_counts[Element] != nleft)
    ) {
    std::cerr<< "Error, the left entity counts are incorrect!" << std::endl;
    result = false;
  }

  stk::mesh::Selector selector_right(partRight);
  stk::mesh::count_entities( selector_right, bulkData , entity_counts );
  if (
    (entity_counts[Node] != expected_num_right_nodes) ||
    (entity_counts[Edge] != expected_num_edges) ||
    (entity_counts[Face] != expected_num_faces) ||
    (entity_counts[Element] != nright)
    ) {
    std::cerr<< "Error, the right counts are incorrect!" << std::endl;
    result = false;
  }
  return result;
}


bool verifyRelations( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright )
{
  const stk::mesh::BulkData & bulkData = mesh.bulkData();
  stk::mesh::Part & partLeft = mesh.partLeft();
  stk::mesh::Part & partRight = mesh.partRight();

  stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];
  stk::mesh::EntityId elem_id = 1 ;

  bool result = true;
  for ( unsigned j = 0 ; j < nleft ; ++j , ++elem_id ) {
    get_elem_node_ids_2( elem_id , node_ids );

    stk::mesh::Entity * const elem = bulkData.get_entity( Element , elem_id );

    if ( elem == NULL ) {
      std::cerr << "Error, element not found!" << std::endl;
      result = false;
    }
    // Verify this element got into the left block.
    if ( !elem->bucket().member(partLeft) ) {
      std::cerr << "Error, element not a member of left block!" << std::endl;
      result = false;
    }
    stk::mesh::PairIterRelation rel = elem->relations();

    // Verify that the number of nodes in this element is correct.
    if( shards::Hexahedron<8> ::node_count != rel.size() ) {
      std::cerr << "Error, number of relations is incorrect!" << std::endl;
      result = false;
    }

    // Verify the nodes of this element got into the left block.
    for ( unsigned i = 0 ; i < shards::Hexahedron<8> ::node_count ; ++i ) {
      stk::mesh::Entity * const rel_node = rel[i].entity();
      if ( ( node_ids[i] != rel_node->identifier() ) ||
          ( !rel_node->bucket().member(partLeft) )   ) {
          std::cerr << "Error, an element's node is just plain wrong!" << std::endl;
          result = false;
      }
    }
  }

  for ( unsigned j = 0 ; j < nright ; ++j , ++elem_id ) {
    get_elem_node_ids_2( elem_id , node_ids );

    stk::mesh::Entity * const elem = bulkData.get_entity( Element , elem_id );

    if ( elem == NULL ) {
      std::cerr << "Error, element not found!" << std::endl;
      result = false;
    }
    // Verify this element got into the right block.
    if ( !elem->bucket().member(partRight) ) {
      std::cerr << "Error, element not a member of right block!" << std::endl;
      result = false;
    }
    stk::mesh::PairIterRelation rel = elem->relations();

    // Verify that the number of nodes in this element is correct.
    if( shards::Hexahedron<8> ::node_count != rel.size() ) {
      std::cerr << "Error, number of relations is incorrect!" << std::endl;
      result = false;
    }

    // Verify the nodes of this element got into the right block.
    for ( unsigned i = 0 ; i < shards::Hexahedron<8> ::node_count ; ++i ) {
      stk::mesh::Entity * const rel_node = rel[i].entity();
      if ( ( node_ids[i] != rel_node->identifier() ) ||
           ( !rel_node->bucket().member(partRight) )   ) {
        std::cerr << "Error, an element's node is just plain wrong!" << std::endl;
        result = false;
      }
    }
  }
  return result;
}

  

template< unsigned NType , stk::mesh::EntityType EType ,
          unsigned NRel , class field_type >
bool gather_field_data( const field_type & field ,
                        const stk::mesh::Entity     & entity ,
                        typename stk::mesh::FieldTraits< field_type >::data_type * dst )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type T ;

  stk::mesh::PairIterRelation rel = entity.relations( EType );

  bool result = NRel == (unsigned) rel.size();

  if ( result ) {
    T * const dst_end = dst + NType * NRel ;
    for ( const T * src ;
          ( dst < dst_end ) &&
          ( src = field_data( field , * rel->entity() ) ) ;
          ++rel , dst += NType ) {
      stk::Copy<NType>( dst , src );
    }
    result = dst == dst_end ;
  }
  return result ;
}


static const double element_coordinates_gold[ ElementTraits::node_count ][ SpaceDim ] = 
  { { 0, 0, 0 } ,
    { 1, 0, 0 } ,
    { 1, 1, 0 } ,
    { 0, 1, 0 } ,
    { 0, 0, 1 } ,
    { 1, 0, 1 } ,
    { 1, 1, 1 } ,
    { 0, 1, 1 }   };


bool verifyFields( const UseCase_2_Mesh & mesh )
{
  bool result = true;

  const VectorFieldType & coordinates_field = mesh.const_coordinates_field();
  const ScalarFieldType & volume_field = mesh.const_volume_field();

  const stk::mesh::BulkData & bulkData = mesh.bulkData();
  const std::vector<stk::mesh::Bucket*> & elem_buckets = 
    bulkData.buckets( stk::mesh::Element );

  // Verify coordinates_field
  std::vector<stk::mesh::Bucket*>::const_iterator element_bucket_it = 
    elem_buckets.begin();
  for (  ; element_bucket_it != elem_buckets.end() ; ++element_bucket_it ) {

    const stk::mesh::Bucket& bucket = **element_bucket_it;
    const size_t num_buckets = bucket.size();

    for( size_t bucket_index = 0; bucket_index < num_buckets; ++bucket_index) {
      const stk::mesh::Entity & elem = bucket[bucket_index] ;
      double elem_coord[ ElementTraits::node_count ][ SpaceDim ];

      const bool gather_result =
        gather_field_data< SpaceDim , stk::mesh::Node , ElementTraits::node_count >
                         ( coordinates_field , elem , & elem_coord[0][0] );

      if ( gather_result == false ) {
        std::cerr << "Error!" << std::endl;
        result = false;
      }

      double base[3] ; 
      base[0] = elem_coord[0][0] ;
      base[1] = elem_coord[0][1] ;
      base[2] = elem_coord[0][2] ;

      for ( unsigned j = 0 ; j < ElementTraits::node_count ; ++j ) {
        elem_coord[j][0] -= base[0] ;
        elem_coord[j][1] -= base[1] ;
        elem_coord[j][2] -= base[2] ;
      }

      for (int node_index=0 ; node_index<ElementTraits::node_count ; ++node_index ) {
        for (int coord_index=0 ; coord_index<SpaceDim ; ++coord_index) {
          double local_elem_coord = elem_coord[node_index][coord_index];
          double gold_elem_coord = element_coordinates_gold[node_index][coord_index];
          if ( local_elem_coord != gold_elem_coord ) {
            std::cerr << "Error!  elem_coord[" << node_index << "]"
              << "[" << coord_index << "] == " << local_elem_coord 
              << " != " << gold_elem_coord 
              << " == element_coordinates_gold[" << node_index << "]"
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
    const double volume_val_gold = 1.0;
    int num_volumes = volume_array.dimension(0);
    for ( int volume_index=0 ; volume_index < num_volumes ; ++volume_index) {
      if ( volume_array(volume_index) != volume_val_gold ) {
        std::cerr << "Error!  volume_array(" << volume_index << ") == "
          << volume_array(volume_index) << " != " << volume_val_gold 
          << " == volume_val_gold " << std::endl;
        result = false;
      }
    }
  }
  // Verify temperature_field
  const ScalarFieldType & temperature_field = mesh.const_temperature_field();
  const std::vector<stk::mesh::Bucket*> & node_buckets = 
    bulkData.buckets( stk::mesh::Node );
  std::vector<stk::mesh::Bucket*>::const_iterator node_bucket_it = 
    node_buckets.begin();
  const double temperature_value_gold = 98.6;
  for ( ; node_bucket_it != node_buckets.end() ; ++node_bucket_it ) {
    const stk::mesh::Bucket & bucket = **node_bucket_it;
    stk::mesh::BucketArray<ScalarFieldType> temperature_array( 
        temperature_field, 
        bucket 
        );
    int num_temps = temperature_array.dimension(0);
    for ( int i=0 ; i < num_temps ; ++i) {
      if (temperature_array(i) != temperature_value_gold) {
        std::cerr << "Error!  temperature_array("<<i<<") == "
          << temperature_array(i) << " != " << temperature_value_gold 
          << " == temperature_value_gold" << std::endl;
        result = false;
      }
    }
  }
  return result;
}

} //namespace use_cases 
} //namespace mesh 
} //namespace stk


