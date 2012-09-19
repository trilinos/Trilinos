/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#ifndef SIERRA_SIERRA_MESH_FIXTURE_GEAR_HPP
#define SIERRA_SIERRA_MESH_FIXTURE_GEAR_HPP

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>

#include <sierra/mesh/details/selected_buckets.hpp>
#include <sierra/mesh/details/constant_size_field.hpp>

namespace sierra {
namespace mesh {
namespace fixture {

class gear {

public:

  typedef details::constant_size_field< double, 3>  CoordinateField;

  gear(
      double arg_element_size  =  0.02,
      double arg_radius_min    =  0.4,
      double arg_radius_max    =  1.5,
      double arg_height_min = -0.4,
      double arg_height_max =  0.4
      );

  const double m_element_size;
  const double m_rad_min, m_rad_max;
  const double m_height_min, m_height_max;

  const size_t m_angle_num;
  const size_t m_rad_num;
  const size_t m_height_num;

  const double m_angle_increment;
  const double m_rad_increment;
  const double m_height_increment;

  const size_t m_num_elements;
  const size_t m_num_nodes;

  modifiable_mesh m_mesh;

  const details::part_key     m_element_part;
  const details::part_key     m_node_part;
  const details::part_key     m_hex_top_part;
  const details::part_key     m_wedge_top_part;
  const details::part_key     m_node_top_part;

  CoordinateField  m_coordinates;

  void generate();

private:

  details::entity_key node_index(
      size_t iz ,       // Thickness index
      size_t ir ,       // Radial index
      size_t ia ) const // Angle index
  {
    return static_cast<details::entity_key>(iz + m_height_num * ( ir + m_rad_num * ia ));
  }

  details::entity_key elem_index(
      size_t iz ,       // Thickness index
      size_t ir ,       // Radial index
      size_t ia ) const // Angle index
  {
    return static_cast<details::entity_key>(m_num_nodes + iz + (m_height_num-1) * ( ir + (m_rad_num-1) * ia ));
  }

  void populate_fields();

  gear(const gear &);
  void operator = (const gear &);
};

}//namespace fixture
}//namespace mesh
}//namespace sierra

#endif // SIERRA_SIERRA_MESH_FIXTURE_GEAR_HPP
