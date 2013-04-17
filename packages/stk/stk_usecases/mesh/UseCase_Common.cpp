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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>

#include <mesh/UseCase_Common.hpp>

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::mesh::MetaData::NODE_RANK;
enum { SpatialDim = 3 };

}

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace use_cases {

bool verify_elem_node_coord(
  mesh::Entity elem ,
  const VectorFieldType & node_coord ,
  const unsigned node_count )
{
  bool result = true;
  mesh::PairIterRelation rel = elem.relations( NODE_RANK );

  if( (unsigned) rel.size() != node_count ) {
    std::cerr << "Error!  relation size == " << rel.size() << " != "
      << node_count << " == node count" << std::endl;
    result = false;
  }

  // Iterating over the nodes in element
  for ( unsigned j = 0 ; j < node_count ; ++j ) {
    mesh::Entity node = rel[j].entity();

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

  }
  return result;
}

} //namespace use_cases
} //namespace mesh
} //namespace stk
