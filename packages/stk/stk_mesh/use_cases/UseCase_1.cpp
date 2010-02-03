/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <use_cases/UseCase_1.hpp>

//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace use_cases {

enum {
  Node = 0,
  Edge = 1,
  Face = 2,
  Element = 3
};

UseCase_1_MetaData::UseCase_1_MetaData(const std::vector<std::string> & entity_type_names)
  :m_metaData(entity_type_names),
  partLeft(m_metaData.declare_part( "block_left", Element )),
  partRight(m_metaData.declare_part( "block_right", Element ))
{
  // Declare the intention of the parts to be for hex 8 elements.
  // This attaches this cell topology to these parts:
  stk::mesh::set_cell_topology< shards::Hexahedron<8>  >( partLeft );
  stk::mesh::set_cell_topology< shards::Hexahedron<8>  >( partRight );

  m_metaData.commit();

}


UseCase_1_Mesh::~UseCase_1_Mesh() 
{
}


std::vector<std::string> get_entity_type_names_1()
{
  // Note:  stk::mesh::fem_entity_type_names provides a default set of these too.
  std::vector<std::string> entity_type_names(4);
  entity_type_names[0] = "Node";
  entity_type_names[1] = "Edge";
  entity_type_names[2] = "Face";
  entity_type_names[3] = "Element";
  return entity_type_names;
}


UseCase_1_Mesh::UseCase_1_Mesh( stk::ParallelMachine comm )
  : UseCase_1_MetaData(get_entity_type_names_1()),
    m_bulkData(m_metaData,comm,field_data_chunk_size)
{
}


void get_elem_node_ids_1( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] )
{
  if ( elem_id == 0 ) {
    std::cerr << "use_case_1, elem_node_ids: ERROR, elem_id ("
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


void populate( UseCase_1_Mesh & mesh , unsigned nleft , unsigned nright )
{
  stk::mesh::BulkData & bulkData = mesh.modifiableBulkData();
  stk::mesh::Part & partLeft = mesh.partLeft;
  stk::mesh::Part & partRight = mesh.partRight;

  // Generate Mesh
  stk::mesh::EntityId elem_id = 1 ;
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];

  for ( unsigned j = 0 ; j < nleft ; ++j , ++elem_id ) {
    get_elem_node_ids_1( elem_id , node_ids );
    // Note declare_element expects a cell topology to have been attached to partLeft.
    stk::mesh::declare_element( bulkData , partLeft , elem_id , node_ids );
  }

  for ( unsigned j = 0 ; j < nright ; ++j , ++elem_id ) {
    get_elem_node_ids_1( elem_id , node_ids );
    // Note declare_element expects a cell topology to have been attached to partRight.
    stk::mesh::declare_element( bulkData , partRight , elem_id , node_ids );
  }
  bulkData.modification_end();
}


bool verifyMesh( const UseCase_1_Mesh & mesh, unsigned nleft, unsigned nright )
{
  bool result = false;
  result = verifyCellTopology(mesh);
  result = result && verifyEntityCounts(mesh,nleft,nright);
  result = result && verifyRelations(mesh,nleft,nright);
  return result;
}


bool verifyCellTopology( const UseCase_1_Mesh & mesh )
{
  stk::mesh::Part & partLeft = mesh.partLeft;
  stk::mesh::Part & partRight = mesh.partRight;

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


bool verifyEntityCounts( const UseCase_1_Mesh & mesh, unsigned nleft, unsigned nright )
{
  const stk::mesh::BulkData & bulkData = mesh.bulkData();
  stk::mesh::Part & partLeft = mesh.partLeft;
  stk::mesh::Part & partRight = mesh.partRight;

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


bool verifyRelations( const UseCase_1_Mesh & mesh, unsigned nleft, unsigned nright )
{
  const stk::mesh::BulkData & bulkData = mesh.bulkData();
  stk::mesh::Part & partLeft = mesh.partLeft;
  stk::mesh::Part & partRight = mesh.partRight;

  stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];
  stk::mesh::EntityId elem_id = 1 ;

  bool result = true;
  for ( unsigned j = 0 ; j < nleft ; ++j , ++elem_id ) {
    get_elem_node_ids_1( elem_id , node_ids );

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
    get_elem_node_ids_1( elem_id , node_ids );

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

} // namespace use_cases 
} // namespace mesh 
} // namespace stk 


