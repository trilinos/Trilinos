/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Mesh_Use_Cases_UseCase_3_hpp
#define Stk_Mesh_Use_Cases_UseCase_3_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

/** stk_mesh Use Case 3 */

namespace stk {
namespace mesh {
namespace use_cases {

typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
typedef stk::mesh::Field<double*,stk::mesh::ElementNode> ElementNodePointerFieldType ;

/** Use case with mixed element topologies and
 *  field relations to provide fast access to node field data
 *  from an element.
 */

class UseCase_3_Mesh {
public:

  ~UseCase_3_Mesh();

  UseCase_3_Mesh( stk::ParallelMachine comm );

  void populate();

  stk::mesh::MetaData m_metaData;
  stk::mesh::BulkData m_bulkData;

  stk::mesh::Part & m_block_hex;
  stk::mesh::Part & m_block_wedge;
  stk::mesh::Part & m_block_tet;
  stk::mesh::Part & m_block_pyramid;
  stk::mesh::Part & m_block_quad_shell;
  stk::mesh::Part & m_block_tri_shell;

  VectorFieldType & m_coordinates_field;
  VectorFieldType & m_centroid_field;
  ScalarFieldType & m_temperature_field;
  ScalarFieldType & m_volume_field;
  ElementNodePointerFieldType & m_element_node_coordinates_field;
};

bool verifyMesh( const UseCase_3_Mesh & mesh );

} //namespace use_cases
} //namespace mesh
} //namespace stk

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
