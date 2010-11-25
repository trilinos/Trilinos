/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Mesh_Use_Cases_UseCase_4_hpp
#define Stk_Mesh_Use_Cases_UseCase_4_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

/** stk_mesh Use Case 4
 */

namespace stk {
namespace mesh {
namespace use_cases {


typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
typedef stk::mesh::Field<double*,stk::mesh::ElementNode> ElementNodePointerFieldType ;

/** Use case for quadratic element with a mix of linear and quadratic nodes. */

class UseCase_4_Mesh {
public:

  ~UseCase_4_Mesh();

  UseCase_4_Mesh( stk::ParallelMachine comm );

  void populate();

  stk::mesh::MetaData m_metaData;
  stk::mesh::BulkData m_bulkData;
  stk::mesh::DefaultFEM m_fem;

  const EntityRank m_elem_rank;
  const EntityRank m_side_rank;

  stk::mesh::Part & m_block_hex20;
  stk::mesh::Part & m_block_wedge15;
  stk::mesh::Part & m_part_vertex_nodes;
  stk::mesh::Part & m_side_part;

  VectorFieldType & m_coordinates_field;
  VectorFieldType & m_velocity_field;
  VectorFieldType & m_centroid_field;
  ScalarFieldType & m_temperature_field;
  ScalarFieldType & m_pressure_field;
  VectorFieldType & m_boundary_field;
  ElementNodePointerFieldType & m_element_node_coordinates_field;

};

void runAlgorithms( const UseCase_4_Mesh & mesh );

bool verifyMesh( const UseCase_4_Mesh & mesh );

} //namespace use_cases
} //namespace mesh
} //namespace stk

#endif // Stk_Mesh_Use_Cases_UseCase_4_hpp
