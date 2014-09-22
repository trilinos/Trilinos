// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef Stk_Mesh_Use_Cases_UseCase_2_hpp
#define Stk_Mesh_Use_Cases_UseCase_2_hpp

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>

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

namespace stk {
namespace mesh {
namespace use_cases {

typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;

// Two part MetaData with four entity types:
// Node, Edge, Face, Element
// and two parts (partLeft and partRight)
// and three fields (coordinates, temperature, and volume)
class UseCase_2_Mesh
{
public:
  ~UseCase_2_Mesh();

  UseCase_2_Mesh( stk::ParallelMachine comm );

  void populate( unsigned nleft , unsigned nright );

  stk::mesh::MetaData m_fem_metaData;
  stk::mesh::BulkData m_bulkData;
  stk::mesh::Part   & m_partLeft;
  stk::mesh::Part   & m_partRight;
  VectorFieldType   & m_coordinates_field;
  ScalarFieldType   & m_temperature_field;
  ScalarFieldType   & m_volume_field;
  const stk::mesh::EntityRank m_elem_rank;
  const stk::mesh::EntityRank m_side_rank;
  const stk::mesh::EntityRank m_edge_rank;
  const stk::mesh::EntityRank m_node_rank;
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
} //namespace stk

#endif // Stk_Mesh_Use_Cases_UseCase_2_hpp

