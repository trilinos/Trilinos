/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <sierra/mesh/fixture/array_mesh_hex_fixture.hpp>

namespace sierra {
namespace mesh {
namespace fixture {

array_mesh_hex_fixture::array_mesh_hex_fixture(unsigned nx, unsigned ny, unsigned nz,
                                               bool upward_connectivity)
  : m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_num_nodes((nx+1)*(ny+1)*(nz+1)),
    m_num_elements(nx*ny*nz),
    m_mesh(upward_connectivity),
    m_upward_connectivity(upward_connectivity),
    m_block_index()
{
}

void array_mesh_hex_fixture::node_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  int key = entity_key;

  x = key % (m_nx+1);
  key /= (m_nx+1);

  y = key % (m_ny+1);
  key /= (m_ny+1);

  z = key;
}

void array_mesh_hex_fixture::elem_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const
{
  int key = entity_key;

  x = key % m_nx;
  key /= m_nx;

  y = key % m_ny;
  key /= m_ny;

  z = key;
}

void array_mesh_hex_fixture::generate_mesh()
{
  m_mesh.reserve(m_num_elements, m_num_nodes);

  stk::topology t = stk::topology::HEX_8;
  int block_id = 1;
  m_block_index = m_mesh.add_block(array_mesh::Element, block_id, m_num_elements, t);

  std::vector<int> node_ids(m_num_nodes);
  const int nnodes = m_num_nodes;
  for(int inode=0; inode<nnodes; ++inode) {
    node_ids[inode] = inode;
  }

  m_mesh.add_node_ids(node_ids.begin(), node_ids.end());

  const int nodes_per_elem = t.num_nodes();
  int node_idx[nodes_per_elem];

  std::vector<int> elem_ids(m_num_elements);

  int nelems = m_num_elements;
  for (int ielem=0; ielem<nelems; ++ielem) {
    unsigned ix = 0, iy = 0, iz = 0;

    elem_ids[ielem] = ielem;

    elem_x_y_z(ielem, ix, iy, iz);

    node_idx[0] = node_index( ix   , iy   , iz   );
    node_idx[1] = node_index( ix+1 , iy   , iz   );
    node_idx[2] = node_index( ix+1 , iy+1 , iz   );
    node_idx[3] = node_index( ix   , iy+1 , iz   );
    node_idx[4] = node_index( ix   , iy   , iz+1 );
    node_idx[5] = node_index( ix+1 , iy   , iz+1 );
    node_idx[6] = node_index( ix+1 , iy+1 , iz+1 );
    node_idx[7] = node_index( ix   , iy+1 , iz+1 );

    m_mesh.add_connectivity(m_block_index, ielem, node_idx, node_idx + nodes_per_elem);
  }

  m_mesh.add_element_ids(elem_ids.begin(), elem_ids.end());

  populate_field_data();
}

void array_mesh_hex_fixture::populate_field_data()
{
  const unsigned spatial_dim = 3;
  m_coord_field.resize(m_num_nodes*spatial_dim);

  for(size_t i=0; i<m_num_nodes; ++i) {
    unsigned x=0, y=0, z=0;
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

