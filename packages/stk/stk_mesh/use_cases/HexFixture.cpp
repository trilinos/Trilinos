/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/HexFixture.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>


#include <stk_mesh/fem/SkinMesh.hpp>

// The following fixture creates a 3x3x3 mesh

HexFixture::HexFixture(stk::ParallelMachine pm)
  : m_meta_data( stk::mesh::fem_entity_rank_names() )
  , m_bulk_data( m_meta_data , pm )
  , m_hex_part( m_meta_data.declare_part("hex_part", stk::mesh::Element) )
  , m_skin_part( m_meta_data.declare_part("skin_part"))
{
  stk::mesh::set_cell_topology<shards::Hexahedron<8> >(m_hex_part);

  m_meta_data.commit();

  generate_mesh();
}

HexFixture::~HexFixture()
{ }

namespace {
  //I know that there is a simple pattern here for discovering the first node from an element id,
  //but after 5 mins of looking at it I did this and moved on.
  unsigned get_first_node_id_for_element( unsigned element_id) {
    switch(element_id) {
      case  1: return 1;
      case  2: return 2;
      case  3: return 3;
      case  4: return 4+1;
      case  5: return 4+2;
      case  6: return 4+3;
      case  7: return 4+4+1;
      case  8: return 4+4+2;
      case  9: return 4+4+3;
      case 10: return 16+1;
      case 11: return 16+2;
      case 12: return 16+3;
      case 13: return 16+4+1;
      case 14: return 16+4+2;
      case 15: return 16+4+3;
      case 16: return 16+4+4+1;
      case 17: return 16+4+4+2;
      case 18: return 16+4+4+3;
      case 19: return 16+16+1;
      case 20: return 16+16+2;
      case 21: return 16+16+3;
      case 22: return 16+16+4+1;
      case 23: return 16+16+4+2;
      case 24: return 16+16+4+3;
      case 25: return 16+16+4+4+1;
      case 26: return 16+16+4+4+2;
      case 27: return 16+16+4+4+3;
      default: return 0;
    }
}
}

void HexFixture::generate_mesh() {

  //const unsigned num_nodes = 64;
  const unsigned num_elements = 27;
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();


  // Note:  This block of code would normally be replaced with a call to stk_io
  // to generate the mesh.

  {
    const stk::mesh::PartVector no_parts;
    const unsigned first_element = (p_rank * num_elements) / p_size;
    const unsigned end_element = ((p_rank + 1) * num_elements) / p_size;

    stk::mesh::PartVector element_parts;
    element_parts.push_back(&m_hex_part);

    const unsigned num_nodes_per_element = 8;

    const unsigned stencil_for_3x3x3_hex_mesh[num_nodes_per_element] =
    { 0, 1, 5, 4, 16, 17, 21, 20};

    // declare elements
    m_bulk_data.modification_begin();
    for (unsigned i = first_element; i < end_element; ++i) {

      unsigned element_id = i + 1;
      stk::mesh::Entity & element = m_bulk_data.declare_entity( stk::mesh::Element, element_id, element_parts);

      unsigned first_node_id = get_first_node_id_for_element( element_id);

      for ( unsigned ordinal = 0; ordinal < num_nodes_per_element; ++ordinal) {
        unsigned node_id = first_node_id + stencil_for_3x3x3_hex_mesh[ordinal];
        stk::mesh::Entity& node = m_bulk_data.declare_entity(stk::mesh::Node, node_id, no_parts);
        m_bulk_data.declare_relation( element , node , ordinal);
      }
    }
    m_bulk_data.modification_end();
  }

  skin_mesh(m_bulk_data, stk::mesh::Element, &m_skin_part);

}
