#include <algorithm>

//#include <usecases/internal_force/simple_hex_fixture.hpp>
#include "simple_hex_fixture.hpp"

namespace fixture {

simple_hex_fixture::simple_hex_fixture(unsigned nx, unsigned ny, unsigned nz)
  : m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_num_nodes((nx+1)*(ny+1)*(nz+1)),
    m_num_elements(nx*ny*nz),
    m_mesh()
{
}

void simple_hex_fixture::node_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  int key = entity_key;

  x = key % (m_nx+1);
  key /= (m_nx+1);

  y = key % (m_ny+1);
  key /= (m_ny+1);

  z = key;
}

void simple_hex_fixture::elem_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  int key = entity_key;

  x = key % m_nx;
  key /= m_nx;

  y = key % m_ny;
  key /= m_ny;

  z = key;
}

void simple_hex_fixture::generate_mesh()
{
  std::vector<std::pair<int,int> > elem_blocks;
  elem_blocks.push_back(std::make_pair(m_num_elements, (int)8));
  m_mesh.allocate(elem_blocks, m_num_nodes);

  int offset = 0;
  for (int ielem=0; ielem<(int)m_num_elements; ++ielem) {
    unsigned ix = 0, iy = 0, iz = 0;
    elem_x_y_z(ielem, ix, iy, iz);

    m_mesh.elem_node_connectivity[offset++] = node_index( ix   , iy   , iz   );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix+1 , iy   , iz   );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix+1 , iy+1 , iz   );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix   , iy+1 , iz   );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix   , iy   , iz+1 );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix+1 , iy   , iz+1 );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix+1 , iy+1 , iz+1 );
    m_mesh.elem_node_connectivity[offset++] = node_index( ix   , iy+1 , iz+1 );
  }

  //setting up node to element connectivity

  int num_node_elems_conn = 0;
  for ( unsigned inode = 0; inode < m_num_nodes; ++inode) {
    unsigned ix = 0, iy = 0, iz = 0;
    node_x_y_z(inode, ix, iy, iz);

    m_mesh.node_elem_offset[inode] = num_node_elems_conn;

    const bool on_x0_face = ix == 0;
    const bool on_y0_face = iy == 0;
    const bool on_z0_face = iz == 0;

    const bool on_x1_face = m_nx < ix;
    const bool on_y1_face = m_ny < iy;
    const bool on_z1_face = m_nz < iz;


    if ( !(on_x1_face || on_y1_face || on_z1_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix   , iy   , iz   ));
      ++num_node_elems_conn;
    }
    if ( !(on_x0_face || on_y1_face || on_z1_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix-1 , iy   , iz   ));
      ++num_node_elems_conn;
    }
    if ( !(on_x0_face || on_y0_face || on_z1_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix-1 , iy-1 , iz   ));
      ++num_node_elems_conn;
    }
    if ( !(on_x1_face || on_y0_face || on_z1_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix   , iy-1 , iz   ));
      ++num_node_elems_conn;
    }
    if ( !(on_x1_face || on_y1_face || on_z0_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix   , iy   , iz-1 ));
      ++num_node_elems_conn;
    }
    if ( !(on_x0_face || on_y1_face || on_z0_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix-1 , iy   , iz-1 ));
      ++num_node_elems_conn;
    }
    if ( !(on_x0_face || on_y0_face || on_z0_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix-1 , iy-1 , iz-1 ));
      ++num_node_elems_conn;
    }
    if ( !(on_x1_face || on_y0_face || on_z0_face) ) {
      m_mesh.node_elem_ids.push_back( elem_index( ix   , iy-1 , iz-1 ));
      ++num_node_elems_conn;
    }
  }

  m_mesh.node_elem_offset[m_num_nodes] = num_node_elems_conn;



  populate_field_data();
}

void simple_hex_fixture::populate_field_data()
{
  const unsigned spatial_dim = 3;
  m_coord_field.resize(m_num_nodes*spatial_dim);

  for(size_t i=0; i<m_num_nodes; ++i) {
    unsigned x, y, z;
    node_x_y_z(i, x, y, z);
    double* coords = &m_coord_field[i*spatial_dim];
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
  }
}

} // fixtures

