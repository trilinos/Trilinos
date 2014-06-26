/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SIERRA_MESH_FIXTURE_ARRAYMESH_HEX_FIXTURE_HPP
#define STK_SIERRA_MESH_FIXTURE_ARRAYMESH_HEX_FIXTURE_HPP

#include <sierra/mesh/array_mesh/array_mesh.hpp>

#include <vector>

namespace sierra {
namespace mesh {
namespace fixture {

/**
 * A 3-dimensional X*Y*Z hex fixture.
 *
 * A coordinate field will be created for all nodes.
 */
class array_mesh_hex_fixture
{
 public:
  /**
    The parameters nx, ny, nz specify the number of elements in each
    spatial dimension. Thus, the number of nodes in each dimension will
    be 1 greater (nx+1, ny+1, nz+1).
   */
  array_mesh_hex_fixture(unsigned nx, unsigned ny, unsigned nz, bool upward_connectivity=false);

  const unsigned                m_nx;
  const unsigned                m_ny;
  const unsigned                m_nz;
  const unsigned                m_num_nodes;
  const unsigned                m_num_elements;
  array_mesh                    m_mesh;
  bool                          m_upward_connectivity;
  array_mesh::BlockIndex   m_block_index;
  std::vector<double>           m_coord_field;

  /**
   * Thinking in terms of a 3D grid of nodes, get the index of the node in
   * the (x, y, z) position.
   */
  int node_index( unsigned x , unsigned y , unsigned z ) const  {
    return x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  int elem_index( unsigned x , unsigned y , unsigned z ) const  {
    return m_num_nodes + x + m_nx * ( y + m_ny * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given its index.
   */
  void node_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void elem_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Create the mesh (into m_mesh).
   */
  void generate_mesh();

  void populate_field_data();

 private:

  array_mesh_hex_fixture();
  array_mesh_hex_fixture( const array_mesh_hex_fixture &);
  array_mesh_hex_fixture & operator = (const array_mesh_hex_fixture &);
};

} // fixtures
} // mesh
} // sierra

#endif

