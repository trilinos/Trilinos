#ifndef __IBMCPP__
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <sierra/mesh/fixture/hex_fixture.hpp>
#include <sierra/mesh/details/cell_topology.hpp>

using namespace sierra::mesh::details;

namespace sierra {
namespace mesh {
namespace fixture {

hex_fixture::hex_fixture(unsigned nx, unsigned ny, unsigned nz)
  : m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_num_nodes((nx+1)*(ny+1)*(nz+1)),
    m_num_elements(nx*ny*nz),
    m_mesh(),
    m_hex_part(m_mesh.declare_part("hex_top_part",cell_topology(shards::getCellTopologyData<shards::Hexahedron<8> >()))),
    m_node_part(m_mesh.declare_part("node_top_part",cell_topology(shards::getCellTopologyData<shards::Node >()))),
    m_coordinates( )
{}

void hex_fixture::generate_mesh()
{
  std::vector<entity_key> element_keys_on_this_processor;
  entity_rank node_rank(0);
  entity_rank elem_rank(4);
  for(size_t i=0; i<m_num_nodes; ++i) {
    m_mesh.add_entity(entity_property(node_rank));
  }
  for(size_t i=0; i<m_num_elements; ++i) {
    entity_key key = m_mesh.add_entity(entity_property(elem_rank));
    element_keys_on_this_processor.push_back(key);
  }

  generate_mesh(element_keys_on_this_processor);
}

void hex_fixture::node_x_y_z( entity_key entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  size_t key = (size_t)entity_key;

  x = key % (m_nx+1);
  key /= (m_nx+1);

  y = key % (m_ny+1);
  key /= (m_ny+1);

  z = key;
}

void hex_fixture::elem_x_y_z( entity_key entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  size_t key = (size_t)entity_key - m_num_nodes;

  x = key % m_nx;
  key /= m_nx;

  y = key % m_ny;
  key /= m_ny;

  z = key;
}

void hex_fixture::generate_mesh(std::vector<entity_key> & element_keys_on_this_processor)
{
  {
    //sort and unique the input elements
    std::vector<entity_key>::iterator ib = element_keys_on_this_processor.begin();
    std::vector<entity_key>::iterator ie = element_keys_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_keys_on_this_processor.erase(ib, ie);
  }

  {
    // Declare the elements that belong on this process

    std::vector<entity_key>::iterator ib = element_keys_on_this_processor.begin();
    const std::vector<entity_key>::iterator ie = element_keys_on_this_processor.end();
    for (; ib != ie; ++ib) {
      entity_key elem_key = *ib;
      m_mesh.change_entity_parts(elem_key, &m_hex_part, &m_hex_part+1);
      unsigned ix = 0, iy = 0, iz = 0;
      elem_x_y_z(elem_key, ix, iy, iz);

      entity_key elem_node[8] ;

      elem_node[0] = node_index( ix   , iy   , iz   );
      elem_node[1] = node_index( ix+1 , iy   , iz   );
      elem_node[2] = node_index( ix+1 , iy   , iz+1 );
      elem_node[3] = node_index( ix   , iy   , iz+1 );
      elem_node[4] = node_index( ix   , iy+1 , iz   );
      elem_node[5] = node_index( ix+1 , iy+1 , iz   );
      elem_node[6] = node_index( ix+1 , iy+1 , iz+1 );
      elem_node[7] = node_index( ix   , iy+1 , iz+1 );

      for(size_t j=0; j<8u; ++j) {
        m_mesh.add_relation(elem_key, elem_node[j], relation_position(0,j));
        m_mesh.change_entity_parts(elem_node[j], &m_node_part, &m_node_part+1);
        m_mesh.change_entity_parts(elem_node[j], &m_hex_part, &m_hex_part+1);
      }
    }
  }

  populate_field_data();
}

void hex_fixture::populate_field_data()
{
  m_coordinates.update_from_mesh(m_node_part, m_mesh);

  for(size_t i=0; i<m_num_nodes; ++i) {
    entity_key node = static_cast<entity_key>(i);
    unsigned x, y, z;
    node_x_y_z(node, x, y, z);
    Scalar* coords = m_coordinates[m_mesh.get_bucket_location(node)];
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
  }
}

} // fixtures
} // mesh
} // sierra
#endif
