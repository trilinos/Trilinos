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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp> // for create_faces
#include <stk_mesh/base/CreateEdges.hpp> // for create_edges
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include <stddef.h>

#include <generated/Iogn_DashSurfaceMesh.h>
#include <generated/Iogn_GeneratedMesh.h>
#include <generated/Iogn_DatabaseIO.h>
#include <Ionit_Initializer.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Utils.h>
#include <Ioss_Region.h>

namespace stk { namespace mesh { class BulkData; } }

namespace
{
//Example1 (for documentation)
TEST(StkMeshHowTo, CreateFacesEdgesHex)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { return; }
  stk::io::StkMeshIoBroker stkIo(communicator);

  const std::string generatedFileName = "generated:8x8x8";
  stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  // ============================================================
  //+ EXAMPLE
  //+ Create the faces..
  stk::mesh::create_faces(stkIo.bulk_data());

  //+ Create the edges..
  stk::mesh::create_edges(stkIo.bulk_data());

  // ==================================================
  // VERIFICATION
  stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
  EXPECT_EQ( 512u, entityCounts[stk::topology::ELEMENT_RANK]);
  EXPECT_EQ(1728u, entityCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(1944u, entityCounts[stk::topology::EDGE_RANK]);
  // MAKE SURE FACES ARE HOOKED TO EDGES
  // this should happen if create_faces is called before create_edges
  stk::mesh::BucketVector const & face_buckets = stkIo.bulk_data().buckets(stk::topology::FACE_RANK);
  for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
    stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
    const unsigned num_expected_edges = bucket.topology().num_edges();
    EXPECT_EQ(4u, num_expected_edges);
    for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
      stk::mesh::Entity face = bucket[face_count];
      EXPECT_EQ(num_expected_edges, stkIo.bulk_data().num_edges(face));
    }
  }
}
//End Example1 (for documentation)

//for now, this is just a copy of above with a different order of create_edges and create_faces
TEST(StkMeshHowTo, CreateEdgesFacesHex)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { return; }
  stk::io::StkMeshIoBroker stkIo(communicator);

  const std::string generatedFileName = "generated:8x8x8";
  stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  // ============================================================
  //+ EXAMPLE
  //+ Create the edges..
  stk::mesh::create_edges(stkIo.bulk_data());

  //+ Create the faces..
  stk::mesh::create_faces(stkIo.bulk_data(), true);
  // ==================================================
  // VERIFICATION
  stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
  EXPECT_EQ( 512u, entityCounts[stk::topology::ELEMENT_RANK]);
  EXPECT_EQ(1728u, entityCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(1944u, entityCounts[stk::topology::EDGE_RANK]);
  // MAKE SURE FACES ARE HOOKED TO EDGES
  // this should happen if create_faces is called before create_edges
  stk::mesh::BucketVector const & face_buckets = stkIo.bulk_data().buckets(stk::topology::FACE_RANK);
  for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
    stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
    const unsigned num_expected_edges = bucket.topology().num_edges();
    EXPECT_EQ(4u, num_expected_edges);
    for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
      stk::mesh::Entity face = bucket[face_count];
      EXPECT_EQ(num_expected_edges, stkIo.bulk_data().num_edges(face));
    }
  }
}

//for now, this is just a copy of above with a different order of create_edges and create_faces and false
//passed in to create_faces for connect_faces_to_edges so we don't connect them together
TEST(StkMeshHowTo, CreateEdgesFacesHexNoConnect)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { return; }
  stk::io::StkMeshIoBroker stkIo(communicator);

  const std::string generatedFileName = "generated:8x8x8";
  stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  // ============================================================
  //+ EXAMPLE
  //+ Create the edges..
  stk::mesh::create_edges(stkIo.bulk_data());

  //+ Create the faces..
  stk::mesh::create_faces(stkIo.bulk_data(), false); //false means don't connect faces to edges
  // ==================================================
  // VERIFICATION
  stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
  EXPECT_EQ( 512u, entityCounts[stk::topology::ELEMENT_RANK]);
  EXPECT_EQ(1728u, entityCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(1944u, entityCounts[stk::topology::EDGE_RANK]);
  // MAKE SURE FACES ARE NOT HOOKED TO EDGES
  stk::mesh::BucketVector const & face_buckets = stkIo.bulk_data().buckets(stk::topology::FACE_RANK);
  for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
    stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
    const unsigned num_expected_edges = bucket.topology().num_edges();
    EXPECT_EQ(4u, num_expected_edges);
    for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
      stk::mesh::Entity face = bucket[face_count];
      EXPECT_EQ(0u, stkIo.bulk_data().num_edges(face));  //expect no edges now
    }
  }
}

Iogn::ExodusData createExodusData(int numberOfHexes, stk::mesh::EntityId* globalIds);
void writeExodusFile(Iogn::GeneratedMesh *generatedMesh, const std::string &exodusFileName);

TEST(StkMeshHowTo, UnderstandEdgeAndFaceOrdering)
{
  int oneHex=1;
  stk::mesh::EntityId globalIdsHex[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  Iogn::ExodusData exoData = createExodusData(oneHex, globalIdsHex);
  Iogn::ExodusMesh* exoMesh = new Iogn::ExodusMesh(exoData);
  std::string exodusFileName = "oneElem.exo";
  // No Need To Delete ExoMesh When Using 'writeExodusFile'
  writeExodusFile(exoMesh, exodusFileName);

  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker stkIo(communicator);
  const std::string generatedFileName = exodusFileName;
  stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  stk::mesh::BulkData &bulkData = stkIo.bulk_data();

  // BEGIN EXAMPLE 2
  // ============================================================
  //+ EXAMPLE
  //+ Create the faces..
  stk::mesh::create_faces(bulkData);

  unsigned goldValuesForHexFaceNodesFromStkTopology[6][4] = {
    {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 4, 8, 7}, {1, 5, 8, 4}, {1, 4, 3, 2}, {5, 6, 7, 8} };

  // Lexicographical smallest permutation per face leads from topology ordering (above) for face to ordering below

  unsigned goldValuesForHexFaceNodesFromCreateFaces[6][4] = {
    {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 4, 8, 7}, {1, 4, 8, 5}, {1, 2, 3, 4}, {5, 6, 7, 8} };

  //+ Create the edges..
  stk::mesh::create_edges(bulkData);

  unsigned goldValuesHexEdgeNodesFromStkTopology[12][2] = {
    {1, 2}, {2, 3}, {3, 4}, {4, 1}, {5, 6}, {6, 7}, {7, 8}, {8, 5}, {1, 5}, {2, 6}, {3, 7}, {4, 8} };

  // Lexicographical smallest permutation per edge leads from topology ordering (above) for edge to ordering below

  unsigned goldValuesHexEdgeNodesFromCreateEdges[12][2] = {
    {1, 2}, {2, 3}, {3, 4}, {1, 4}, {5, 6}, {6, 7}, {7, 8}, {5, 8}, {1, 5}, {2, 6}, {3, 7}, {4, 8} };

  // END EXAMPLE 2

  // ==================================================
  // VERIFICATION

  unsigned elementId=1;
  stk::mesh::Entity element=bulkData.get_entity(stk::topology::ELEMENT_RANK, elementId);
  const stk::mesh::Entity *edges = bulkData.begin_edges(element);
  stk::topology hex = stk::topology::HEXAHEDRON_8;
  unsigned edgeNodeIds[2];

  for (unsigned i=0;i<bulkData.num_edges(element);i++)
  {
    const stk::mesh::Entity *edgeNodes = bulkData.begin_nodes(edges[i]);
    EXPECT_EQ(goldValuesHexEdgeNodesFromCreateEdges[i][0], bulkData.identifier(edgeNodes[0]));
    EXPECT_EQ(goldValuesHexEdgeNodesFromCreateEdges[i][1], bulkData.identifier(edgeNodes[1]));

    hex.edge_nodes((stk::mesh::EntityId*)globalIdsHex, i, (unsigned*)edgeNodeIds);
    EXPECT_EQ(goldValuesHexEdgeNodesFromStkTopology[i][0], edgeNodeIds[0]);
    EXPECT_EQ(goldValuesHexEdgeNodesFromStkTopology[i][1], edgeNodeIds[1]);
  }

  const stk::mesh::Entity *faces = bulkData.begin_faces(element);
  unsigned faceNodeIds[4];
  for (unsigned i=0;i<bulkData.num_faces(element);i++)
  {
    const stk::mesh::Entity *faceNodes = bulkData.begin_nodes(faces[i]);
    EXPECT_EQ(goldValuesForHexFaceNodesFromCreateFaces[i][0], bulkData.identifier(faceNodes[0])) << "failed for face " << i << std::endl;
    EXPECT_EQ(goldValuesForHexFaceNodesFromCreateFaces[i][1], bulkData.identifier(faceNodes[1])) << "failed for face " << i << std::endl;
    EXPECT_EQ(goldValuesForHexFaceNodesFromCreateFaces[i][2], bulkData.identifier(faceNodes[2])) << "failed for face " << i << std::endl;
    EXPECT_EQ(goldValuesForHexFaceNodesFromCreateFaces[i][3], bulkData.identifier(faceNodes[3])) << "failed for face " << i << std::endl;

    hex.face_nodes((stk::mesh::EntityId*)globalIdsHex, i, (unsigned*)faceNodeIds);
    EXPECT_EQ(goldValuesForHexFaceNodesFromStkTopology[i][0], faceNodeIds[0]);
    EXPECT_EQ(goldValuesForHexFaceNodesFromStkTopology[i][1], faceNodeIds[1]);
    EXPECT_EQ(goldValuesForHexFaceNodesFromStkTopology[i][2], faceNodeIds[2]);
    EXPECT_EQ(goldValuesForHexFaceNodesFromStkTopology[i][3], faceNodeIds[3]);
  }

  stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
  EXPECT_EQ(1u, entityCounts[stk::topology::ELEMENT_RANK]);
  EXPECT_EQ(6u, entityCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(12u, entityCounts[stk::topology::EDGE_RANK]);

  unlink(exodusFileName.c_str());
}

Iogn::ExodusData createExodusData(int numberOfHexes, stk::mesh::EntityId* globalIds)
{
  int numberOfElementBlocks = numberOfHexes;
  int numberOfNodesInEachElementBlock = 8;
  int globalNumberOfNodes = numberOfElementBlocks * numberOfNodesInEachElementBlock;
  std::vector<double> coordinates(globalNumberOfNodes*3);
  std::vector< std::vector<int> > elementBlockConnectivity(numberOfElementBlocks);
  for(int blockIndex=0; blockIndex < numberOfElementBlocks; blockIndex++)
  {
    elementBlockConnectivity[blockIndex].resize(numberOfNodesInEachElementBlock);
  }
  int numberOfElementsInEachBlock = 1;
  std::vector<int> globalNumberOfElementsInBlock(numberOfElementBlocks, numberOfElementsInEachBlock);
  std::vector<int> localNumberOfElementsInBlock(numberOfElementBlocks, numberOfElementsInEachBlock);
  std::vector<Iogn::Topology> blockTopologicalData(numberOfElementBlocks, Iogn::Hex8);
  int numberOfTotalElements = numberOfElementsInEachBlock * numberOfElementBlocks;
  std::vector<int> globalIdsOfLocalElements(numberOfTotalElements);
  std::vector<int> globalIdsOfLocalNodes(globalNumberOfNodes);
  for(int i=0; i < numberOfTotalElements; i++)
  {
    globalIdsOfLocalElements[i] = i+1;
  }
  for(int i=0; i < globalNumberOfNodes; i++)
  {
    //          globalIdsOfLocalNodes[i] = i+1;
    globalIdsOfLocalNodes[i] = globalIds[i];
  }
  for(int blockIndex=0; blockIndex < numberOfElementBlocks; blockIndex++)
  {
    int elementBlockOffset = numberOfNodesInEachElementBlock * blockIndex;
    for(int i=0; i < numberOfNodesInEachElementBlock; i++)
    {
      elementBlockConnectivity[blockIndex][i] = elementBlockOffset+i+1;
    }
  }

  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1,
    0, 1, 1,
  };

  coordinates.clear();
  coordinates.resize(globalNumberOfNodes*3);
  for (int i=0;i<numberOfElementBlocks;i++)
  {
    int offset = 3*numberOfNodesInEachElementBlock*i;
    for (int j=0;j<numberOfNodesInEachElementBlock;j++)
    {
      coordinates[offset+3*j + 0] = coords[3*j+0]+i;
      coordinates[offset+3*j + 1] = coords[3*j+1];
      coordinates[offset+3*j + 2] = coords[3*j+2];
    }
  }

  return Iogn::ExodusData(coordinates, elementBlockConnectivity, globalNumberOfElementsInBlock, localNumberOfElementsInBlock,
                          blockTopologicalData, globalNumberOfNodes, globalIdsOfLocalElements, globalIdsOfLocalNodes);
}

void writeExodusFile(Iogn::GeneratedMesh *generatedMesh, const std::string &exodusFileName)
{
  Ioss::Init::Initializer io;

  Ioss::DatabaseIO *database = Ioss::IOFactory::create("generated", "2x3x3", Ioss::READ_MODEL, MPI_COMM_WORLD);
  Iogn::DatabaseIO *iognDatabase = dynamic_cast<Iogn::DatabaseIO *>(database);
  iognDatabase->setGeneratedMesh(generatedMesh);

  Ioss::Region* io_region = new Ioss::Region(database);
  stk::io::StkMeshIoBroker meshData;
  std::shared_ptr<Ioss::Region> junk(io_region, [](auto pointerWeWontDelete){});
  meshData.add_mesh_database(junk);
  meshData.create_input_mesh();
  meshData.populate_bulk_data();
  size_t result_file_index = meshData.create_output_mesh(exodusFileName, stk::io::WRITE_RESULTS);
  meshData.write_output_mesh(result_file_index);
  delete io_region;
}

}

