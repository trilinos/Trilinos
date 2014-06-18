#ifndef SAMBA_FIXTURES_BOX_FIXTURE_HPP
#define SAMBA_FIXTURES_BOX_FIXTURE_HPP

#include <samba/mesh.hpp>
#include <samba/field.hpp>
#include <samba/rank_field.hpp>

namespace samba { namespace fixtures {

class box_fixture
{
public:

  typedef samba::rank_field< samba::entity_rank::node_type
                            ,double
                            ,samba::spatial_dimension_functor
                           > coordinate_field_type;

  typedef samba::rank_field< samba::entity_rank::node_type
                            ,const double
                            ,samba::spatial_dimension_functor
                           > const_coordinate_field_type;

  box_fixture( uint32_t x = 10u
              ,uint32_t y = 10u
              ,uint32_t z = 10u
              ,connectivity_map map = connectivity_map::default_map()
             )
    : m_mesh(map)
    , m_coordinates(m_mesh)
    , m_x(x)
    , m_y(y)
    , m_z( map.spatial_dimension() == 3u ? z : 0)
    , m_num_nodes( map.spatial_dimension() == 3u
                  ? (m_x+1)*(m_y+1)*(m_z+1)
                  : (m_x+1)*(m_y+1)
                 )
    , m_num_elements( map.spatial_dimension() == 3u
                     ? m_x*m_y*m_z
                     : m_x*m_y
                    )
  {
    if (map.spatial_dimension() == 3u)
      create_hex();
    else
      create_quad();
  }


  samba::mesh mesh() { return m_mesh; }

  const samba::mesh mesh() const  { return m_mesh; }

  coordinate_field_type coordinates()  { return m_coordinates; }

  const_coordinate_field_type coordinates() const  { return m_coordinates; }

  uint32_t num_nodes() const { return m_num_nodes; }
  uint32_t num_elements() const { return m_num_elements; }

private:

  void create_hex()
  {
    samba::entity_key_interval nodes = m_mesh.add_entities( samba::entity_topology::node(), m_num_nodes );
    samba::entity_key_interval hexes = m_mesh.add_entities( samba::entity_topology::hex_8(), m_num_elements );

    const bool create_back_relations = m_mesh.connectivity_map()(entity_rank::node(),entity_rank::element()) != samba::connectivity_kind::invalid();

    for (uint32_t x=0; x<m_x; ++x) {
    for (uint32_t y=0; y<m_y; ++y) {
    for (uint32_t z=0; z<m_z; ++z) {

      samba::connectivity_ordinal ordinal = {0};

      samba::entity_key hex = hexes[hex_offset(x,y,z)];
      samba::entity_key node;

      node = nodes[node_offset( x   , y   , z   )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x+1 , y   , z   )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x+1 , y   , z+1 )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x   , y   , z+1 )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x   , y+1 , z   )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x+1 , y+1 , z   )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x+1 , y+1 , z+1 )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

      node = nodes[node_offset( x   , y+1 , z+1 )];
      m_mesh.add_connectivity(hex,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,hex,ordinal);
      ++ordinal;

    }
    }
    }

    for (size_t i=0; i<nodes.size(); ++i) {
      const samba::entity_key key = nodes[i];
      const double x_coord = i%(m_x+1);
      const double y_coord = (i/(m_x+1))%(m_y+1);
      const double z_coord = (i/((m_x+1)*(m_y+1)));

      m_coordinates[key][0] = x_coord;
      m_coordinates[key][1] = y_coord;
      m_coordinates[key][2] = z_coord;

    }

  }

  void create_quad()
  {
    samba::entity_key_interval nodes = m_mesh.add_entities( samba::entity_topology::node(), m_num_nodes );
    samba::entity_key_interval quads = m_mesh.add_entities( samba::entity_topology::quad_4(), m_num_elements );

    const bool create_back_relations = m_mesh.connectivity_map()(entity_rank::node(),entity_rank::element()) != samba::connectivity_kind::invalid();

    for (uint32_t x=0; x<m_x; ++x) {
    for (uint32_t y=0; y<m_y; ++y) {

      samba::connectivity_ordinal ordinal = {0};

      samba::entity_key quad = quads[quad_offset(x,y)];
      samba::entity_key node;

      node = nodes[node_offset( x   , y  )];
      m_mesh.add_connectivity(quad,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,quad,ordinal);
      ++ordinal;

      node = nodes[node_offset( x+1 , y  )];
      m_mesh.add_connectivity(quad,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,quad,ordinal);
      ++ordinal;

      node = nodes[node_offset( x+1 , y+1 )];
      m_mesh.add_connectivity(quad,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,quad,ordinal);
      ++ordinal;

      node = nodes[node_offset( x   , y+1 )];
      m_mesh.add_connectivity(quad,node,ordinal);
      if (create_back_relations) m_mesh.add_connectivity(node,quad,ordinal);
      ++ordinal;

    }
    }

    for (size_t i=0; i<nodes.size(); ++i) {
      const samba::entity_key key = nodes[i];
      const double x_coord = i%(m_x+1);
      const double y_coord = (i/(m_x+1))%(m_y+1);

      m_coordinates[key][0] = x_coord;
      m_coordinates[key][1] = y_coord;

    }
  }

  size_t node_offset(uint32_t x, uint32_t y, uint32_t z)
  { return (x + ( m_x + 1 ) * ( y + ( m_y + 1 ) * z ));  }

  size_t hex_offset(uint32_t x, uint32_t y, uint32_t z)
  { return (x + m_x * ( y +  m_y * z ));  }

  size_t node_offset(uint32_t x, uint32_t y)
  { return (x + ( m_x + 1 ) *  y );  }

  size_t quad_offset(uint32_t x, uint32_t y)
  { return (x + m_x * y );  }

private:

  samba::mesh m_mesh;
  coordinate_field_type m_coordinates;
  uint32_t m_x;
  uint32_t m_y;
  uint32_t m_z;
  uint32_t m_num_nodes;
  uint32_t m_num_elements;

};


}} //namespace samba::fixtures

#endif //SAMBA_FIXTURES_BOX_FIXTURE_HPP

