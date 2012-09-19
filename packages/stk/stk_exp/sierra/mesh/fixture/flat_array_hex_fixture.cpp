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

#include <sierra/mesh/fixture/flat_array_hex_fixture.hpp>

namespace sierra {
namespace mesh {
namespace fixture {

flat_array_hex_fixture::flat_array_hex_fixture(unsigned nx, unsigned ny, unsigned nz)
  : m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_num_nodes((nx+1)*(ny+1)*(nz+1)),
    m_num_elements(nx*ny*nz),
    m_mesh()
{
}

void flat_array_hex_fixture::node_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  int key = entity_key;

  x = key % (m_nx+1);
  key /= (m_nx+1);

  y = key % (m_ny+1);
  key /= (m_ny+1);

  z = key;
}

void flat_array_hex_fixture::elem_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  int key = entity_key;

  x = key % m_nx;
  key /= m_nx;

  y = key % m_ny;
  key /= m_ny;

  z = key;
}

void flat_array_hex_fixture::generate_mesh()
{
  std::vector<std::pair<int,int> > elem_blocks;
  elem_blocks.push_back(std::make_pair(m_num_elements, (int)8));
  m_mesh.allocate(elem_blocks);

  int offset = 0;
  for (int ielem=0; ielem<(int)m_num_elements; ++ielem) {
    unsigned ix = 0, iy = 0, iz = 0;
    elem_x_y_z(ielem, ix, iy, iz);

    m_mesh.connectivity_table[offset++] = node_index( ix   , iy   , iz   );
    m_mesh.connectivity_table[offset++] = node_index( ix+1 , iy   , iz   );
    m_mesh.connectivity_table[offset++] = node_index( ix+1 , iy   , iz+1 );
    m_mesh.connectivity_table[offset++] = node_index( ix   , iy   , iz+1 );
    m_mesh.connectivity_table[offset++] = node_index( ix   , iy+1 , iz   );
    m_mesh.connectivity_table[offset++] = node_index( ix+1 , iy+1 , iz   );
    m_mesh.connectivity_table[offset++] = node_index( ix+1 , iy+1 , iz+1 );
    m_mesh.connectivity_table[offset++] = node_index( ix   , iy+1 , iz+1 );
  }

  populate_field_data();
}

void flat_array_hex_fixture::populate_field_data()
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
} // mesh
} // sierra

#endif
