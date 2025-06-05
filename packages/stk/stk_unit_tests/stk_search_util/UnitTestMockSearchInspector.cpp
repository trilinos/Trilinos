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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_io/FillMesh.hpp>                             // for fill_mesh
#include "stk_io/WriteMesh.hpp"                            // for write_mesh
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/BulkData.hpp>                      // for BulkData
#include "stk_mesh/base/EntityKey.hpp"                     // for EntityKey
#include "stk_search/FilterCoarseSearch.hpp"               // for ObjectOuts...
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/spmd/GeometricSearch.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"               // for build_mesh
#include "stk_unit_test_utils/TextMesh.hpp"                // for get_full_t...
#include "stk_unit_test_utils/UnitTestSearchUtils.hpp"

#include <gtest/gtest.h>
#include "mpi.h"            // for MPI_COMM_WORLD

#include <algorithm>                                       // for max, copy
#include <cmath>                                           // for sqrt, pow
#include <cstddef>                                         // for size_t
#include <cstdint>                                         // for int64_t
#include <iostream>                                        // for operator<<
#include <memory>                                          // for shared_ptr
#include <stdexcept>                                       // for runtime_error
#include <string>                                          // for string
#include <utility>                                         // for pair
#include <vector>                                          // for vector, swap

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace {

using STKNodeFormatterBase = stk::search::InspectorOutputFormatterBase<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;
using STKNodeFormatter = stk::search::InspectorOutputFormatter<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;
using STKNodeInspectorInfo = stk::search::InspectorInfo<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;

using STKElemFormatterBase = stk::search::InspectorOutputFormatterBase<stk::search::spmd::ElementSendMesh, stk::search::spmd::ElementRecvMesh>;
using STKElemFormatter = stk::search::InspectorOutputFormatter<stk::search::spmd::ElementSendMesh, stk::search::spmd::ElementRecvMesh>;
using STKElemInspectorInfo = stk::search::InspectorInfo<stk::search::spmd::ElementSendMesh, stk::search::spmd::ElementRecvMesh>;

inline std::ostream& operator<<(std::ostream& os, const Kokkos::pair<stk::search::spmd::EntityKeyPair,int> ekey)
{
  os << "{" << ekey.first << ", INDEX: " << ekey.second << "}";
  return os;
}

struct MockOutputNodeFormatter : public STKNodeFormatterBase
{
protected:
  void internal_write_info(std::ostream& os, const STKNodeInspectorInfo& info) const override
  {
    os << "[ ELEMENT(" << info.domainEntityKey << "), PROC(" << info.domainProc << ") ] -> [ NODE("
       <<        info.rangeEntityKey  << "), PROC(" << info.rangeProc  << ") ] ";
    for(const std::string& partName : info.domainParts) {
      os << partName << " ";
    }
    os << std::endl;
  }
};

struct MockOutputElemFormatter : public STKElemFormatterBase
{
protected:
  void internal_write_info(std::ostream& os, const STKElemInspectorInfo& info) const override
  {
    os << "[ ELEMENT(" << info.domainEntityKey << "), PROC(" << info.domainProc << ") ] -> [ NODE("
       <<        info.rangeEntityKey  << "), PROC(" << info.rangeProc  << ") ] ";
    for(const std::string& partName : info.domainParts) {
      os << partName << " ";
    }
    os << std::endl;
  }
};

stk::search::spmd::EntityKeyPair make_spmd_key(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
{
  return stk::search::spmd::make_entity_key_pair(bulk, stk::mesh::EntityKey(rank, id));
}
stk::search::spmd::EntityKeyPair make_spmd_key(stk::mesh::EntityRank rank, stk::mesh::EntityId id)
{
  return stk::search::spmd::EntityKeyPair(stk::mesh::Entity(), stk::mesh::EntityKey(rank, id));
}

template <typename SENDMESH, typename RECVMESH>
void test_inspection(const std::string& fileName,
                     const stk::search::InspectorOutputFormatterBase<SENDMESH, RECVMESH>& formatter,
                     const std::vector<stk::search::InspectorInfo<SENDMESH, RECVMESH>>& expectedInfoVec)
{
  std::ostringstream os;
  formatter.write_info(os, expectedInfoVec);

  std::ifstream inputFile(fileName);
  EXPECT_TRUE(inputFile);
  std::ostringstream ss;
  ss << inputFile.rdbuf();

  EXPECT_EQ(os.str(), ss.str());

  unlink(fileName.c_str());
}

void setup_and_run_serial_mock_inspection(const std::string& fileName)
{
  using STKMockSearch = stk::search::spmd::GeometricSearch<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  STK_ThrowRequire(stk::parallel_machine_size(comm) == 1);

  double parametricTolerance = 0.001;
  double geometricTolerance = 0.1;
  double expansionSum = geometricTolerance;
  double expansionFactor = 0.0;
  double slant = 1.0;
  const stk::search::SearchMethod search_method = stk::search::KDTREE;

  std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_slanted_single_hex_bulk(slant);
  auto sendMesh = stk::unit_test_util::construct_hex_send_mesh(*sendBulk, parametricTolerance);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_single_point_bulk(0.5, 0.5, 0.5);
  auto destMesh = stk::unit_test_util::construct_node_recv_mesh(*recvBulk,
                                                                parametricTolerance, geometricTolerance);

  sendMesh->set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy::EXTRAPOLATE);

  STKMockSearch search(sendMesh, destMesh, "search", MPI_COMM_WORLD, expansionFactor, expansionSum, search_method);

  search.initialize();

  std::vector<stk::search::spmd::NodeRecvMesh::EntityKey> rangeKeys = { make_spmd_key(*recvBulk, stk::topology::NODE_RANK, 1) };

  search.register_inspector(fileName, rangeKeys);
  search.register_output_formatter(std::make_shared<MockOutputNodeFormatter>());

  search.coarse_search();
  search.local_search();
  EXPECT_EQ(1u, search.get_range_to_domain().size());

  search.inspect_user_defined_entities();
}

void test_serial_mock_inspection(const std::string& fileName)
{
  MockOutputNodeFormatter formatter;
  STKNodeInspectorInfo expectedInfo;
  expectedInfo.domainEntityKey = make_spmd_key(stk::topology::ELEM_RANK, 1);
  expectedInfo.rangeEntityKey = make_spmd_key(stk::topology::NODE_RANK, 1);
  expectedInfo.domainProc = 0;
  expectedInfo.rangeProc = 0;
  expectedInfo.domainParts = {"BLOCK_1"};

  test_inspection(fileName, formatter, {expectedInfo});
}

TEST(MockSearchTest, withInspector)
{
  const std::string fileName = "mock_search.txt";

  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    setup_and_run_serial_mock_inspection(fileName);
    test_serial_mock_inspection(fileName);
  }
}

void create_search_and_run_parallel_inspection(const std::string& fileName,
                                               std::shared_ptr<stk::search::spmd::ElementSendMesh> sendMesh,
                                               std::shared_ptr<stk::search::spmd::ElementRecvMesh> destMesh)
{
  using STKMockSearch = stk::search::spmd::GeometricSearch<stk::search::spmd::ElementSendMesh, stk::search::spmd::ElementRecvMesh>;
  using RangeKeyVec = std::vector<stk::search::spmd::ElementRecvMesh::EntityKey>;

  double expansionFactor = 0.0;
  stk::ParallelMachine comm = sendMesh->comm();

  STKMockSearch search(sendMesh, destMesh, "searchtest", comm, expansionFactor, 0.0);
  search.initialize();
  sendMesh->set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy::ABORT);

  stk::mesh::EntityKeyVector inspectedElements = {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1u),
                                                  stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2u)};
  RangeKeyVec rangeKeys;

  destMesh->fill_entity_keys(inspectedElements, rangeKeys);

  search.register_inspector(fileName, rangeKeys);

  search.coarse_search();
  EXPECT_NO_THROW(search.local_search());

  search.inspect_user_defined_entities();
}

STKElemFormatter setup_with_centroid_recv_mesh_and_run_parallel_inspection(const std::string& fileName)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  STK_ThrowRequire(stk::parallel_machine_size(comm) == 2);

  double parametricTolerance = 0.1;
  double geometricTolerance = 0.001;

  std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_mesh(comm);
  stk::io::fill_mesh("generated:1x1x2", *sendBulk);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(comm);
  stk::io::fill_mesh("generated:1x1x2", *recvBulk);

  std::shared_ptr<stk::search::spmd::ElementSendMesh> sendMesh =
      stk::unit_test_util::construct_hex_send_mesh(*sendBulk, parametricTolerance);
  std::shared_ptr<stk::search::spmd::ElementRecvMesh> destMesh =
      stk::unit_test_util::construct_element_centroid_recv_mesh(*recvBulk, parametricTolerance, geometricTolerance);

  create_search_and_run_parallel_inspection(fileName, sendMesh, destMesh);

  return STKElemFormatter(sendMesh, destMesh);
}

STKElemFormatter setup_with_gauss_point_recv_mesh_and_run_parallel_inspection(const std::string& fileName)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  STK_ThrowRequire(stk::parallel_machine_size(comm) == 2);

  double parametricTolerance = 0.1;
  double geometricTolerance = 0.001;

  std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_mesh(comm);
  stk::io::fill_mesh("generated:1x1x2", *sendBulk);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(comm);
  stk::io::fill_mesh("generated:1x1x2", *recvBulk);

  std::shared_ptr<stk::search::spmd::ElementSendMesh> sendMesh =
      stk::unit_test_util::construct_hex_send_mesh(*sendBulk, parametricTolerance);
  std::shared_ptr<stk::search::spmd::ElementRecvMesh> destMesh =
      stk::unit_test_util::construct_hex_gauss_point_recv_mesh(*recvBulk, parametricTolerance, geometricTolerance);

  create_search_and_run_parallel_inspection(fileName, sendMesh, destMesh);

  return STKElemFormatter(sendMesh, destMesh);
}

void test_parallel_inspection(const STKElemFormatter& formatter, const std::string& fileName, unsigned numEvalpoints)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_rank(comm) == 0) {
    STKElemInspectorInfo info;
    std::vector<STKElemInspectorInfo> expectedInfoVec;

    for(unsigned i=0; i<numEvalpoints; ++i) {
      info.domainEntityKey = make_spmd_key(stk::topology::ELEM_RANK, 1u);
      info.rangeEntityKey = Kokkos::make_pair(make_spmd_key(stk::topology::ELEM_RANK, 1u), i);
      info.domainProc = 0;
      info.rangeProc = 0;
      info.domainParts = {"BLOCK_1"};
      expectedInfoVec.push_back(info);
    }

    for(unsigned i=0; i<numEvalpoints; ++i) {
      info.domainEntityKey = make_spmd_key(stk::topology::ELEM_RANK, 2u);
      info.rangeEntityKey = Kokkos::make_pair(make_spmd_key(stk::topology::ELEM_RANK, 2u), i);
      info.domainProc = 1;
      info.rangeProc = 1;
      info.domainParts = {"BLOCK_1"};
      expectedInfoVec.push_back(info);
    }

    test_inspection(fileName, formatter, expectedInfoVec);
  }
}

TEST(MockSearchTest, parallelElementCentroidInspection)
{
  const std::string fileName = "parallel_centroid_search.txt";
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  unsigned numEvalPoints = 1;

  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  STKElemFormatter formatter = setup_with_centroid_recv_mesh_and_run_parallel_inspection(fileName);
  test_parallel_inspection(formatter, fileName, numEvalPoints);
}

TEST(MockSearchTest, parallelElementGaussPointInspection)
{
  const std::string fileName = "parallel_gauss_point_search.txt";
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  unsigned numEvalPoints = 8;

  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  STKElemFormatter formatter = setup_with_gauss_point_recv_mesh_and_run_parallel_inspection(fileName);
  test_parallel_inspection(formatter, fileName, numEvalPoints);
}

}

