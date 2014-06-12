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
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <use_cases/UseCase_Common.hpp>

/** stk_mesh Use Case 4
 *
 * This use case creates the mesh below. We have only labelled the
 * node ids of the vertices visible from the front. Sides are
 * created for the side of the elements that are facing us in this diagram. 
 *
 * This mesh consists of two Hex27 elements side by side with three Wedge18 on top.
 *
 * <PRE>
 *         +60-----66+
 *        /         /\
 *       /         /  \
 *      +58-----64+    \
 *     /  \     /  \    +45
 *    /    \   /    \  /|
 *   /      \ /      \/ |
 *   +31-----+37---43+  |     Z  Y
 *   |       |       |  +15   | /
 *   |       |       | /      |/
 *   |       |       |/       *--X
 *   +-------+-------+
 *   1       7      13
 * </PRE>
 *
 * There are 22 nodes on this front face including mid-edge and
 * mid-side nodes in addition to the labelled vertices.
 *
 */

namespace stk_classic {
namespace mesh {
namespace use_cases {

/** Use case for quadratic element with a mix of linear and quadratic nodes. */

class UseCase_4_Mesh {
public:

  ~UseCase_4_Mesh();

  UseCase_4_Mesh( stk_classic::ParallelMachine comm, bool doCommit=true );

  void populate();

  fem::FEMMetaData m_fem_metaData;
  BulkData m_bulkData;

  const EntityRank m_elem_rank;
  const EntityRank m_side_rank;
  const EntityRank m_node_rank;

  Part & m_block_hex27;
  Part & m_block_wedge18;
  Part & m_part_vertex_nodes;
  Part & m_side_part;

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
} //namespace stk_classic

#endif // Stk_Mesh_Use_Cases_UseCase_4_hpp
