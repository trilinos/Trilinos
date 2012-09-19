/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef SIERRA_SIERRA_MESH_FIXTURE_HEX_FIXTURE_HPP
#define SIERRA_SIERRA_MESH_FIXTURE_HEX_FIXTURE_HPP

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>

#include <sierra/mesh/details/selected_buckets.hpp>
#include <sierra/mesh/details/constant_size_field.hpp>

namespace sierra {
namespace mesh {
namespace fixture {

/**
 * A 3-dimensional X*Y*Z hex fixture.
 *
 * A coordinate field will be added to all nodes.
 */
class hex_fixture
{
 public:
  typedef double                     Scalar ;
  typedef details::constant_size_field<Scalar,3> CoordinateField;

  /**
    The parameters nx, ny, nz specify the number of elements in each
    spatial dimension. Thus, the number of nodes in each dimension will
    be 1 greater (nx+1, ny+1, nz+1).
   */
  hex_fixture(unsigned nx, unsigned ny, unsigned nz);

  const unsigned                m_nx;
  const unsigned                m_ny;
  const unsigned                m_nz;
  const unsigned                m_num_nodes;
  const unsigned                m_num_elements;
  modifiable_mesh               m_mesh;
  details::part_key             m_hex_part;
  details::part_key             m_node_part;
  CoordinateField               m_coordinates ;

  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  details::entity_key node_index( unsigned x , unsigned y , unsigned z ) const  {
    return static_cast<details::entity_key>(x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z ));
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  details::entity_key elem_index( unsigned x , unsigned y , unsigned z ) const  {
    return static_cast<details::entity_key>(m_num_nodes + x + m_nx * ( y + m_ny * z ));
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y_z( details::entity_key entity_key, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void elem_x_y_z( details::entity_key entity_key, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Create the mesh (into m_mesh).
   */
  void generate_mesh();

  void generate_mesh( std::vector<details::entity_key> & element_keys_on_this_processor );

  void populate_field_data();

 private:

  hex_fixture();
  hex_fixture( const hex_fixture &);
  hex_fixture & operator = (const hex_fixture &);
};

} // fixtures
} // mesh
} // sierra

#endif

