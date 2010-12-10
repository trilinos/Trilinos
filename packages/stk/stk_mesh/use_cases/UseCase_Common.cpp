/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <use_cases/UseCase_Common.hpp>

using stk::mesh::fem::NODE_RANK;
enum { SpatialDim = 3 };

//----------------------------------------------------------------------

namespace stk{
namespace mesh {
namespace use_cases {

bool verify_elem_node_coord(
  mesh::Entity & elem ,
  const ElementNodePointerFieldType & elem_node_coord ,
  const VectorFieldType & node_coord ,
  const unsigned node_count )
{
  bool result = true;
  mesh::PairIterRelation rel = elem.relations( fem::NODE_RANK );

  if( (unsigned) rel.size() != node_count ) {
    std::cerr << "Error!  relation size == " << rel.size() << " != "
      << node_count << " == node count" << std::endl;
    result = false;
  }

  // Field data for the elem_node_coord for elem
  mesh::EntityArray< ElementNodePointerFieldType >
    elem_node_array( elem_node_coord , elem );

  // Checking the size and dimensionality of elem_node_array
  {
    const unsigned n1 = elem_node_array.dimension<0>();
    if( n1 != node_count ) {
      std::cerr << "Error!  element node array dimension<0> == " << n1 << " != "
        << node_count << " == node count" << std::endl;
      result = false;
    }
    if ( (unsigned) elem_node_array.size() != node_count ) {
      std::cerr << "Error!  element node array size == "
        << elem_node_array.size() << " != " << node_count << " == node count"
        << std::endl;
      result = false;
    }
  }

  // Get the raw memory for the entity field array
  double * const * const elem_data = elem_node_array.contiguous_data();

  // Iterating over the nodes in element
  for ( unsigned j = 0 ; j < node_count ; ++j ) {
    mesh::Entity & node = * rel[j].entity();

    // Field data for the nodal coordinate field for node
    mesh::EntityArray< VectorFieldType > node_coord_array( node_coord , node );

    // Checking the size and dimensionality of node_coord_array
    {
      const unsigned n1 = node_coord_array.dimension<0>();
      if( n1 != (unsigned) SpatialDim ) {
        std::cerr << "Error!  node coord array dimension<0> == " << n1 << " != "
          << SpatialDim << " == SpatialDim" << std::endl;
        result = false;
      }
      if( node_coord_array.size() != SpatialDim ) {
        std::cerr << "Error!  node coord array size == "
          << node_coord_array.size() << " != " << SpatialDim
          << " == SpatialDim" << std::endl;
        result = false;
      }
    }

    // Confirming that element data points to nodal coordinate data for
    // this node
    double * const node_data = node_coord_array.contiguous_data();
    if( elem_data[j] != node_data ) {
      std::cerr << "Error!  elem_data[" << j << "] == " << elem_data[j]
        << " != " << node_data << " node_data" << std::endl;
      result = false;
    }
  }
  return result;
}

bool verify_elem_node_coord_by_part(
    stk::mesh::Part & part,
    const std::vector<Bucket *> & bucket_vector,
    const ElementNodePointerFieldType & elem_node_coord,
    const VectorFieldType & node_coord,
    const unsigned node_count )
{
  Selector selector(part);
  std::vector<Entity *> entities;
  get_selected_entities( selector, bucket_vector, entities);
  std::vector<Entity *>::iterator entity_it = entities.begin();
  bool result = true;
  for ( ; entity_it != entities.end() ; ++entity_it ) {
    result = result &&
      verify_elem_node_coord(
          **entity_it , elem_node_coord , node_coord , node_count
          );
  }
  return result;
}

} //namespace use_cases
} //namespace mesh
} //namespace stk
