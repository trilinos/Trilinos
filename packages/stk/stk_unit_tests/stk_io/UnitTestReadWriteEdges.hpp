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

#ifndef UNIT_TEST_READ_WRITE_EDGES_HPP
#define UNIT_TEST_READ_WRITE_EDGES_HPP

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/ConnectEdgesImpl.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>

struct ExpectedValues
{
  ExpectedValues()
    : numConnectedEdges(0), globalEdgeCount(0), globalElemCount(0) {}

  ExpectedValues(std::vector<unsigned> numLocalEdgesPerProc_, unsigned numFaces_,
                 unsigned numConnectedEdges_, unsigned globalEdgeCount_, unsigned globalElemCount_)
                 : numLocalEdgesPerProc(numLocalEdgesPerProc_), numConnectedEdges(numConnectedEdges_),
                   globalEdgeCount(globalEdgeCount_), globalElemCount(globalElemCount_) {}

  std::vector<unsigned> numEdgesPerProc;
  std::vector<unsigned> numLocalEdgesPerProc;
  std::vector<unsigned> numFacesPerProc;
  std::vector<unsigned> numLocalFacesPerProc;
  unsigned numConnectedEdges;
  unsigned globalEdgeCount;
  unsigned globalElemCount;
};

class StkEdgeIoTest : public stk::unit_test_util::MeshFixture
{
public:
  StkEdgeIoTest() : stk::unit_test_util::MeshFixture()
  {
  }

  void setup_edge_mesh(unsigned numBlocks);

  void setup_mesh_with_edges(unsigned numBlocks);

  void setup_mesh_with_edges_and_faces(unsigned numBlocks);

  void test_connectivity_to_element(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank entityRank);

  void test_entity_count(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank entityRank,
                         unsigned expectedNumLocalEntities, unsigned expectedNumEntities);

  void test_edges(const stk::mesh::BulkData& bulk);

  void test_faces(const stk::mesh::BulkData& bulk);

  virtual void output_mesh();

  void test_output_mesh();

  void test_output_mesh(stk::mesh::BulkData& bulk);

  virtual void load_output_mesh(stk::mesh::BulkData& bulk);

  void set_expected_values(ExpectedValues& expectedValues_);

  virtual ~StkEdgeIoTest()
  {
    unlink(fileName.c_str());
  }

protected:
  std::string fileName = "output.exo";
  std::string edgePartName = "edgeBlock";
  std::string facePartName = "faceBlock";
  ExpectedValues expectedValues;
};

#endif