// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_NE, etc
#include <stddef.h>                     // for size_t
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_topology/topology.hpp>    // for topology, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, FieldBase
#include "stk_mesh/base/Types.hpp"      // for BucketVector
#include <cmath>
#include <string>
namespace stk { namespace mesh { class Part; } }

namespace
{
//-BEGIN
TEST(StkMeshHowTo, iterateElemNodeConnectivity_ForEachEntityWithNodes)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> stkMesh = stk::mesh::MeshBuilder(comm).create();
  // Generate a mesh of unit-cube hexes with a sideset
  const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
  stk::io::fill_mesh(generatedMeshSpecification, *stkMesh);

  typedef stk::mesh::Field<double> CoordinatesField_t;
  CoordinatesField_t const & coord_field =
      *dynamic_cast<CoordinatesField_t const *>(stkMesh->mesh_meta_data().coordinate_field());

  constexpr unsigned nodesPerHex = 8;
  constexpr unsigned spatialDim = 3;
  unsigned count = 0;
  double elementNodeCoords[nodesPerHex][spatialDim] = {
    {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN},
    {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN} };

  stk::mesh::Selector all = stkMesh->mesh_meta_data().universal_part();

  stk::mesh::for_each_entity_run_with_nodes(*stkMesh, stk::topology::ELEM_RANK, all,
    [&](stk::mesh::Entity elem, const stk::mesh::Entity* nodes, size_t numNodesPerEntity) {
      EXPECT_EQ(numNodesPerEntity, nodesPerHex);
      for (unsigned inode = 0; inode < numNodesPerEntity; ++inode) {
        const double *coords = stk::mesh::field_data(coord_field, nodes[inode]);
        elementNodeCoords[inode][0] = coords[0];
        elementNodeCoords[inode][1] = coords[1];
        elementNodeCoords[inode][2] = coords[2];
        ++count;
      }
    });

  const unsigned numElems = 2*2*2;
  const unsigned totalNodesVisited = numElems * nodesPerHex;
  EXPECT_EQ(count, totalNodesVisited);
  EXPECT_FALSE(std::isnan(elementNodeCoords[0][0]));
}

TEST(StkMeshHowTo, iterateConnectivity_General_BulkData)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> stkMesh = stk::mesh::MeshBuilder(comm).create();
  // Generate a mesh of unit-cube hexes with a sideset
  const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
  stk::io::fill_mesh(generatedMeshSpecification, *stkMesh);

  typedef stk::mesh::Field<double> CoordinatesField_t;
  CoordinatesField_t const & coord_field =
      *dynamic_cast<CoordinatesField_t const *>(stkMesh->mesh_meta_data().coordinate_field());

  constexpr unsigned nodesPerHex = 8;
  constexpr unsigned spatialDim = 3;
  unsigned count = 0;
  double elementNodeCoords[nodesPerHex][spatialDim] = {
    {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN},
    {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN} };

  stk::mesh::for_each_entity_run(*stkMesh, stk::topology::ELEM_RANK,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity elem) {
      const stk::mesh::ConnectedEntities nodes = stkMesh->get_connected_entities(elem, stk::topology::NODE_RANK);
      EXPECT_EQ(nodes.size(), nodesPerHex);

      for (unsigned inode = 0; inode < nodes.size(); ++inode) {
        const double *coords = stk::mesh::field_data(coord_field, nodes[inode]);
        elementNodeCoords[inode][0] = coords[0];
        elementNodeCoords[inode][1] = coords[1];
        elementNodeCoords[inode][2] = coords[2];
        ++count;
      }
    });

  const unsigned numElems = 2*2*2;
  const unsigned totalNodesVisited = numElems * nodesPerHex;
  EXPECT_EQ(count, totalNodesVisited);
  EXPECT_FALSE(std::isnan(elementNodeCoords[0][0]));
}

TEST(StkMeshHowTo, iterateConnectivity_Buckets)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { return; }
  std::unique_ptr<stk::mesh::BulkData> stkMesh = stk::mesh::MeshBuilder(comm).create();
  // Generate a mesh of unit-cube hexes with a sideset
  const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
  stk::io::fill_mesh(generatedMeshSpecification, *stkMesh);

  typedef stk::mesh::Field<double> CoordinatesField_t;
  CoordinatesField_t const & coord_field =
      *dynamic_cast<CoordinatesField_t const *>(stkMesh->mesh_meta_data().coordinate_field());

  const stk::mesh::BucketVector &elementBuckets = stkMesh->buckets(stk::topology::ELEMENT_RANK);

  constexpr unsigned nodesPerHex = 8;
  constexpr unsigned spatialDim = 3;
  unsigned count = 0;
  double elementNodeCoords[nodesPerHex][spatialDim] = {
    {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN},
    {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN}, {NAN,NAN,NAN} };

  for (size_t bucketIndex = 0; bucketIndex < elementBuckets.size(); ++bucketIndex) {
    const stk::mesh::Bucket &elemBucket = *elementBuckets[bucketIndex];

    for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex) {
      unsigned numNodes = elemBucket.num_nodes(elemIndex);
      EXPECT_EQ(numNodes, nodesPerHex);
      const stk::mesh::Entity* nodes = elemBucket.begin_nodes(elemIndex);

      for (unsigned inode = 0; inode < numNodes; ++inode) {
        const double *coords = stk::mesh::field_data(coord_field, nodes[inode]);
        elementNodeCoords[inode][0] = coords[0];
        elementNodeCoords[inode][1] = coords[1];
        elementNodeCoords[inode][2] = coords[2];
        ++count;
      }
    }
  }
  const unsigned numElems = 2*2*2;
  const unsigned totalNodesVisited = numElems * nodesPerHex;
  EXPECT_EQ(count, totalNodesVisited);
  EXPECT_FALSE(std::isnan(elementNodeCoords[0][0]));
}
//-END
}
