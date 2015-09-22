/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Mesh_Use_Cases_UseCase_2_hpp
#define Stk_Mesh_Use_Cases_UseCase_2_hpp

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

/** stk_mesh Use Case 2
 *
 * This use case creates a mesh containing a chain of elements.It'll
 * split them between left and righte parts with left coming first.The
 * point of this use case is to demonstrate construction of a very simple mesh.
 *
 * Assume the following mesh of 4 hex8 elements.
 *
 *  Global node and element numbering
 * <PRE>
 *      3       7      11      15      19
 *      +-------+-------+-------+-------+
 *     /       /       /       /       /|
 *   4/      8/     12/     16/     20/ |
 *   +-------+-------+-------+-------+  |
 *   |       |       |       |       |  +18        Z  Y
 *   |  e1   |  e2   |  e3   |  e4   | /           | /
 *   |       |       |       |       |/            |/
 *   +-------+-------+-------+-------+             *--X
 *   1       5      9       13      17
 * </PRE>
 *
 *  Local node numbering
 * <PRE>
 *      8       7
 *      +-------+
 *     /       /|
 *   5/      6/ |
 *   +-------+  |
 *   |       |  +3
 *   |  e1   | /
 *   |       |/
 *   +-------+
 *   1       2
 * </PRE>
 */

namespace stk_classic {
namespace mesh {
namespace use_cases {

typedef stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> VectorFieldType ;
typedef stk_classic::mesh::Field<double>                      ScalarFieldType ;

// Two part MetaData with four entity types:
// Node, Edge, Face, Element
// and two parts (partLeft and partRight)
// and three fields (coordinates, temperature, and volume)
class UseCase_2_Mesh
{
public:
  ~UseCase_2_Mesh();

  UseCase_2_Mesh( stk_classic::ParallelMachine comm );

  void populate( unsigned nleft , unsigned nright );

  stk_classic::mesh::fem::FEMMetaData m_fem_metaData;
  stk_classic::mesh::BulkData m_bulkData;
  stk_classic::mesh::Part   & m_partLeft;
  stk_classic::mesh::Part   & m_partRight;
  VectorFieldType   & m_coordinates_field;
  ScalarFieldType   & m_temperature_field;
  ScalarFieldType   & m_volume_field;
  const stk_classic::mesh::EntityRank m_elem_rank;
  const stk_classic::mesh::EntityRank m_side_rank;
  const stk_classic::mesh::EntityRank m_edge_rank;
  const stk_classic::mesh::EntityRank m_node_rank;
};

/**
 * Verify correctness of mesh
 */
bool verifyMesh( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright );

// Helper functions for verifyMesh
bool verifyCellTopology( const UseCase_2_Mesh & mesh );
bool verifyEntityCounts( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright );
bool verifyRelations( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright );
bool verifyFields( const UseCase_2_Mesh & mesh );

} //namespace use_cases
} //namespace mesh
} //namespace stk_classic

#endif // Stk_Mesh_Use_Cases_UseCase_2_hpp

