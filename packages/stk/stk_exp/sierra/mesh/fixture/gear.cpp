#ifndef __IBMCPP__
/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cmath>
#include <stdexcept>
#include <limits>
#include <iostream>

#include <sierra/mesh/fixture/gear.hpp>
#include <sierra/mesh/details/cell_topology.hpp>

using namespace sierra::mesh::details;

namespace sierra {
namespace mesh {
namespace fixture {

namespace {

const double PI     = 3.14159265358979;
const double TWO_PI = 2 * PI;

} // namespace

gear::gear(
    double arg_element_size,
    double arg_radius_min,
    double arg_radius_max,
    double arg_height_min,
    double arg_height_max
    ) :
      m_element_size(arg_element_size)
    , m_rad_min(arg_radius_min)
    , m_rad_max(arg_radius_max)
    , m_height_min(arg_height_min)
    , m_height_max(arg_height_max)

    , m_angle_num(static_cast<size_t>(TWO_PI/m_element_size))
    , m_rad_num(static_cast<size_t>(2 +(m_rad_max - m_rad_min)/m_element_size))
    , m_height_num(static_cast<size_t>(2 +(m_height_max - m_height_min)/m_element_size))

    , m_angle_increment( TWO_PI / static_cast<double>(m_angle_num))
    , m_rad_increment( (m_rad_max-m_rad_min) / static_cast<double>(m_rad_num-1))
    , m_height_increment( (m_height_max-m_height_min) / static_cast<double>(m_height_num-1))

    , m_num_elements( m_angle_num * (m_rad_num-1) * (m_height_num-1) )
    , m_num_nodes( m_angle_num * (m_rad_num) * (m_height_num) )

    , m_mesh()
    , m_element_part(m_mesh.declare_part("element_part",entity_rank(3)))
    , m_node_part(m_mesh.declare_part("node_part",entity_rank(0)))
    , m_hex_top_part(m_mesh.declare_part("hex_top_part",cell_topology(shards::getCellTopologyData<shards::Hexahedron<8> >())))
    , m_wedge_top_part(m_mesh.declare_part("wedge_top_part",cell_topology(shards::getCellTopologyData<shards::Wedge<6> >())))
    , m_node_top_part(m_mesh.declare_part("node_top_part",cell_topology(shards::getCellTopologyData<shards::Node >())))
    , m_coordinates()

{}

//
//-----------------------------------------------------------------------------
//

void gear::populate_fields() {

  m_coordinates.update_from_mesh( m_node_part & m_node_top_part, m_mesh);

  //setup the cylindrical_coord_field on the hex nodes
  for ( size_t ir = 0 ; ir < m_rad_num-1; ++ir ) {
    const double rad = m_rad_min + m_rad_increment * ir ;

    for ( size_t ia = 0 ; ia < m_angle_num; ++ia ) {
      const double angle = m_angle_increment * ia ;

      for ( size_t iz = 0 ; iz < m_height_num; ++iz ) {
        const double height = m_height_min + m_height_increment * iz ;

        entity_key node = node_index(iz,ir,ia);
        bucket_location location = m_mesh.get_bucket_location(node);

        double * coords = m_coordinates[location];

        coords[0] = rad * std::cos(angle);
        coords[1] = rad * std::sin(angle);
        coords[2] = height;
      }
    }
  }

  //setup the cylindrical_coord_field on the wedge nodes
  {
    size_t ir = m_rad_num-1;
    //const double rad = m_rad_min + m_rad_increment * ir ;
    const double rad = 1.1*m_rad_max;

    for ( size_t ia = 0 ; ia < m_angle_num; ++ia ) {
      const double angle = m_angle_increment * (ia + ia +1.0)/2.0;

      for ( size_t iz = 0 ; iz < m_height_num; ++iz ) {
        const double height = m_height_min + m_height_increment * iz ;

        entity_key node = node_index(iz,ir,ia);
        bucket_location location = m_mesh.get_bucket_location(node);

        double * coords = m_coordinates[location];

        coords[0] = rad * std::cos(angle);
        coords[1] = rad * std::sin(angle);
        coords[2] = height;
      }
    }
  }
}


//
//-----------------------------------------------------------------------------
//

void gear::generate()
{
  entity_rank node_rank(0);
  entity_rank elem_rank(4);
  for(size_t i=0; i<m_num_nodes; ++i) {
    m_mesh.add_entity(entity_property(node_rank));
  }
  for(size_t i=0; i<m_num_elements; ++i) {
    m_mesh.add_entity(entity_property(elem_rank));
  }
  //setup hex elements
  {
    std::vector<part_key> element_parts;
    element_parts.push_back(m_element_part);
    element_parts.push_back(m_hex_top_part);

    std::vector<part_key> node_parts;
    node_parts.push_back(m_node_part);
    node_parts.push_back(m_node_top_part);
    node_parts.push_back(m_hex_top_part);

    for ( size_t ir = 0 ; ir < m_rad_num -2 ; ++ir ) {
      for ( size_t ia = 0 ; ia < m_angle_num; ++ia ) {
        for ( size_t iz = 0 ; iz < m_height_num -1 ; ++iz ) {

          entity_key  elem = elem_index(iz,ir,ia);
          m_mesh.change_entity_parts(elem, element_parts.begin(), element_parts.end());

          const size_t ia_1 = ( ia + 1 ) % m_angle_num ;
          const size_t ir_1 = ir + 1 ;
          const size_t iz_1 = iz + 1 ;

          entity_key node[8] ;

          node[0] = node_index(iz  , ir  , ia_1 );
          node[1] = node_index(iz  , ir  , ia   );
          node[2] = node_index(iz_1, ir  , ia   );
          node[3] = node_index(iz_1, ir  , ia_1 );
          node[4] = node_index(iz  , ir_1, ia_1 );
          node[5] = node_index(iz  , ir_1, ia   );
          node[6] = node_index(iz_1, ir_1, ia   );
          node[7] = node_index(iz_1, ir_1, ia_1 );

          for ( size_t j = 0 ; j < 8u ; ++j ) {
            m_mesh.add_relation( elem, node[j], relation_position(0,j) );
            m_mesh.change_entity_parts(node[j], node_parts.begin(), node_parts.end());
          }
        }
      }
    }
  }

  //setup wedges elements
  {
    std::vector<part_key> element_parts;
    element_parts.push_back(m_element_part);
    element_parts.push_back(m_wedge_top_part);

    std::vector<part_key> node_parts;
    node_parts.push_back(m_node_part);
    node_parts.push_back(m_node_top_part);
    node_parts.push_back(m_wedge_top_part);

    size_t ir = m_rad_num-2 ;
    for ( size_t ia = 0 ; ia < m_angle_num; ++ia ) {
      for ( size_t iz = 0 ; iz < m_height_num -1 ; ++iz ) {

        entity_key  elem = elem_index(iz,ir,ia);
        m_mesh.change_entity_parts(elem, element_parts.begin(), element_parts.end());

        const size_t ia_1 = ( ia + 1 ) % m_angle_num ;
        const size_t ir_1 = ir + 1 ;
        const size_t iz_1 = iz + 1 ;

        entity_key node[6] ;

        node[0] = node_index(iz  , ir  , ia_1 );
        node[1] = node_index(iz  , ir  , ia   );
        node[2] = node_index(iz  , ir_1, ia   );
        node[3] = node_index(iz_1, ir  , ia_1 );
        node[4] = node_index(iz_1, ir  , ia   );
        node[5] = node_index(iz_1, ir_1, ia   );

        for ( size_t j = 0 ; j < 6u ; ++j ) {
          m_mesh.add_relation( elem, node[j], relation_position(0,j) );
          m_mesh.change_entity_parts(node[j], node_parts.begin(), node_parts.end());
        }
      }
    }
  }

  populate_fields();

}

//
//-----------------------------------------------------------------------------
//


}//namespace fixture
}//namespace mesh
}//namespace sierra
#endif
