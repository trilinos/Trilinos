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

#include "gtest/gtest.h"

#include "util/parallel_search_test_util.hpp"

namespace {

TEST(SplitCommTest, one_to_one)
{
  int nProcsSendMesh = 1, nProcsRecvMesh = 1;
  SplitCommTestUtil sctu(nProcsSendMesh, nProcsRecvMesh);

  if (!sctu.get_status()) { GTEST_SKIP(); }

  EXPECT_EQ(0, sctu.get_local_comm_rank(0));
  EXPECT_EQ(0, sctu.get_local_comm_rank(1));
}

TEST(SplitCommTest, one_to_two)
{
  int nProcsSendMesh = 1, nProcsRecvMesh = 2;
  SplitCommTestUtil sctu(nProcsSendMesh, nProcsRecvMesh);

  if (!sctu.get_status()) { GTEST_SKIP(); }

  EXPECT_EQ(0, sctu.get_local_comm_rank(0));
  EXPECT_EQ(0, sctu.get_local_comm_rank(1));
  EXPECT_EQ(1, sctu.get_local_comm_rank(2));
}

TEST(SplitCommTest, two_to_one)
{
  int nProcsSendMesh = 2, nProcsRecvMesh = 1;
  SplitCommTestUtil sctu(nProcsSendMesh, nProcsRecvMesh);

  if (!sctu.get_status()) { GTEST_SKIP(); }

  EXPECT_EQ(0, sctu.get_local_comm_rank(0));
  EXPECT_EQ(1, sctu.get_local_comm_rank(1));
  EXPECT_EQ(0, sctu.get_local_comm_rank(2));
}

TEST(SearchBoundingBox, CreateBoundingBox)
{
  int commSize;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  if(commSize > 1) { GTEST_SKIP(); }

  stk::middle_mesh::mesh::impl::MeshSpec spec;
  spec.numelX = 1;
  spec.numelY = 1;
  spec.xmin = 0;
  spec.xmax = 1;
  spec.ymin = 0;
  spec.ymax = 1;

  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return stk::middle_mesh::utils::Point(pt.x, pt.y, 0); };
  auto mesh = stk::middle_mesh::mesh::impl::create_mesh(spec, func, MPI_COMM_WORLD, false);

  auto vertices = mesh->get_vertices();
  EXPECT_EQ(4u, vertices.size());

  auto elements = mesh->get_elements();
  EXPECT_EQ(1u, elements.size());

  std::vector<stk::middle_mesh::impl::SearchMesh::BoundingBox> boundingBoxVec;
  stk::middle_mesh::impl::SearchMesh meshA(mesh);
  meshA.fill_bounding_boxes(boundingBoxVec);

  EXPECT_EQ(1u, boundingBoxVec.size());
  auto box = boundingBoxVec[0];
  EXPECT_NEAR(spec.xmin, box.first.get_x_min(), 1e-5);
  EXPECT_NEAR(spec.ymin, box.first.get_y_min(), 1e-5);
  EXPECT_NEAR(spec.xmax, box.first.get_x_max(), 1e-5);
  EXPECT_NEAR(spec.ymax, box.first.get_y_max(), 1e-5);
}


class ParallelSearch : public ::testing::Test {
protected:
  ParallelSearch() {}

  void setup_comm()
  {
    m_splitCommUtil = std::make_shared<SplitCommTestUtil>(m_numProcsSendMesh, m_numProcsRecvMesh);

    m_color = m_splitCommUtil->get_color();
    m_splitComm = m_splitCommUtil->get_comm();
  }

  void setup_mesh(unsigned sendNumElemInX, unsigned sendNumElemInY,
                  unsigned recvNumElemInX, unsigned recvNumElemInY)
  {
    m_sendMeshSpec.numelX = sendNumElemInX;
    m_sendMeshSpec.numelY = sendNumElemInY;
    m_sendMeshSpec.xmin = 0;
    m_sendMeshSpec.xmax = 1;
    m_sendMeshSpec.ymin = 0;
    m_sendMeshSpec.ymax = 1;

    m_recvMeshSpec.numelX = recvNumElemInX;
    m_recvMeshSpec.numelY = recvNumElemInY;
    m_recvMeshSpec.xmin = 0;
    m_recvMeshSpec.xmax = 1;
    m_recvMeshSpec.ymin = 0;
    m_recvMeshSpec.ymax = 1;

    auto func = [&](const stk::middle_mesh::utils::Point& pt) { return stk::middle_mesh::utils::Point(pt.x, pt.y, 0); };

    if (m_color == SplitCommColor::SEND) {
      m_mesh = stk::middle_mesh::mesh::impl::create_mesh(m_sendMeshSpec, func, m_splitComm, false);
    } else if (m_color == SplitCommColor::RECV) {
      m_mesh = stk::middle_mesh::mesh::impl::create_mesh(m_recvMeshSpec, func, m_splitComm, false);
    }
  }

  void setup_search()
  {
    if (m_color == stk::middle_mesh::impl::SplitCommColor::SEND) {
      m_searchSendMesh = std::make_shared<SearchMesh>(m_mesh);
      m_search = std::make_shared<BoundingBoxSearch>(
          m_searchSendMesh, m_searchRecvMesh, "BoundingBoxSearch", MPI_COMM_WORLD);
    } else if (m_color == stk::middle_mesh::impl::SplitCommColor::RECV) {
      m_searchRecvMesh = std::make_shared<SearchMesh>(m_mesh);
      m_search = std::make_shared<BoundingBoxSearch>(
          m_searchSendMesh, m_searchRecvMesh, "BoundingBoxSearch", MPI_COMM_WORLD);
    } else {
      ThrowRequireMsg(false, "Invalid SplitComm color");
    }

    m_search->coarse_search();
  }

  void setup(int numProcsSendMesh, int numProcsRecvMesh,
             unsigned sendNumElemInX, unsigned sendNumElemInY,
             unsigned recvNumElemInX, unsigned recvNumElemInY)
  {
    m_numProcsSendMesh = numProcsSendMesh;
    m_numProcsRecvMesh = numProcsRecvMesh;

    setup_comm();
    setup_mesh(sendNumElemInX, sendNumElemInY, recvNumElemInX, recvNumElemInY);
    setup_search();
  }

  void test_search_result(const unsigned expectedNumSearchPairs)
  {
    int myRankWorld;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankWorld);

    const SearchRelationVec& pairedEntities = m_search->get_range_to_domain();
    const UnpairedRelationVec& unpairedEntities = m_search->get_unpaired_recv_entities();

    if(m_color == SplitCommColor::SEND) {
      EXPECT_EQ(expectedNumSearchPairs, pairedEntities.size());
      EXPECT_EQ(0u, unpairedEntities.size());
    }

    unsigned expectedGlobalRank = 0;

    for(const auto &pairedEntity : pairedEntities) {
      const SearchMesh::EntityProc& recvEntityProc = pairedEntity.first;
      const SearchMesh::EntityProc& sendEntityProc = pairedEntity.second;

      if(m_color == SplitCommColor::RECV) {
        bool isGhost = myRankWorld != (int)recvEntityProc.proc();

        if(!isGhost) {
          EXPECT_EQ(0, recvEntityProc.id());
          EXPECT_EQ(0, sendEntityProc.id());
          EXPECT_EQ(expectedGlobalRank, sendEntityProc.proc());

          expectedGlobalRank++;
        }
      }
    }
  }

  void fill_mesh_scatter(MeshScatterSpec& scatter)
  {
    const SearchRelationVec& pairedEntities = m_search->get_range_to_domain();
    const UnpairedRelationVec& unpairedEntities = m_search->get_unpaired_recv_entities();
    EXPECT_EQ(0u, unpairedEntities.size());

    if(m_color == SplitCommColor::SEND) {
      auto elements = m_search->send_mesh()->get_mesh()->get_elements();

      for(const auto &pairedEntity : pairedEntities) {
        const SearchMesh::EntityProc& recvEntityProc = pairedEntity.first;
        const SearchMesh::EntityProc& sendEntityProc = pairedEntity.second;

        int globalRank = recvEntityProc.proc();
        int splitRank = m_splitCommUtil->get_local_comm_rank(globalRank);

        int entityId = sendEntityProc.id();

        stk::middle_mesh::mesh::MeshEntityPtr entity = elements[entityId];

        scatter.add_destination(entity, splitRank);
      }
    }
  }

  int m_numProcsSendMesh{0};
  int m_numProcsRecvMesh{0};
  SplitCommColor m_color{SplitCommColor::INVALID};

  MPI_Comm m_splitComm{MPI_COMM_NULL};

  stk::middle_mesh::mesh::impl::MeshSpec m_sendMeshSpec;
  stk::middle_mesh::mesh::impl::MeshSpec m_recvMeshSpec;

  std::shared_ptr<BoundingBoxSearch> m_search;

  std::shared_ptr<stk::middle_mesh::mesh::Mesh> m_mesh;

  std::shared_ptr<stk::middle_mesh::impl::SearchMesh> m_searchSendMesh;
  std::shared_ptr<stk::middle_mesh::impl::SearchMesh> m_searchRecvMesh;

  std::shared_ptr<SplitCommTestUtil> m_splitCommUtil;
};

TEST_F(ParallelSearch, one_by_one_to_one_by_one)
{
  unsigned sendNumElemInX = 1,  sendNumElemInY = 1;
  unsigned recvNumElemInX = 1,  recvNumElemInY = 1;

  int numProcsSendMesh = sendNumElemInX*sendNumElemInY;
  int numProcsRecvMesh = recvNumElemInX*recvNumElemInY;

  SplitCommTestUtil sctu(numProcsSendMesh, numProcsRecvMesh);
  if (!sctu.get_status()) { GTEST_SKIP(); }

  setup(numProcsSendMesh, numProcsRecvMesh, sendNumElemInX, sendNumElemInY, recvNumElemInX, recvNumElemInY);

  unsigned expectedNumSearchPairs = 1u;
  test_search_result(expectedNumSearchPairs);

  stk::middle_mesh::impl::MeshScatterSpec sendMeshSpec;
  fill_mesh_scatter(sendMeshSpec);
}

TEST_F(ParallelSearch, one_by_one_to_two_by_two)
{
  unsigned sendNumElemInX = 1,  sendNumElemInY = 1;
  unsigned recvNumElemInX = 2,  recvNumElemInY = 2;

  int numProcsSendMesh = sendNumElemInX*sendNumElemInY;
  int numProcsRecvMesh = recvNumElemInX*recvNumElemInY;

  SplitCommTestUtil sctu(numProcsSendMesh, numProcsRecvMesh);
  if (!sctu.get_status()) { GTEST_SKIP(); }

  setup(numProcsSendMesh, numProcsRecvMesh, sendNumElemInX, sendNumElemInY, recvNumElemInX, recvNumElemInY);

  unsigned expectedNumSearchPairs = 4u;
  test_search_result(expectedNumSearchPairs);

  stk::middle_mesh::impl::MeshScatterSpec sendMeshSpec;
  fill_mesh_scatter(sendMeshSpec);
}

TEST_F(ParallelSearch, two_by_two_to_two_by_one)
{
  unsigned sendNumElemInX = 2,  sendNumElemInY = 2;
  unsigned recvNumElemInX = 2,  recvNumElemInY = 1;

  int numProcsSendMesh = sendNumElemInX*sendNumElemInY;
  int numProcsRecvMesh = recvNumElemInX*recvNumElemInY;

  SplitCommTestUtil sctu(numProcsSendMesh, numProcsRecvMesh);
  if (!sctu.get_status()) { GTEST_SKIP(); }

  setup(numProcsSendMesh, numProcsRecvMesh, sendNumElemInX, sendNumElemInY, recvNumElemInX, recvNumElemInY);

  unsigned expectedNumSearchPairs = 2u;
  test_search_result(expectedNumSearchPairs);

  stk::middle_mesh::impl::MeshScatterSpec sendMeshSpec;
  fill_mesh_scatter(sendMeshSpec);
}

TEST_F(ParallelSearch, two_by_two_to_two_by_two)
{
  unsigned sendNumElemInX = 2,  sendNumElemInY = 2;
  unsigned recvNumElemInX = 2,  recvNumElemInY = 2;

  int numProcsSendMesh = sendNumElemInX*sendNumElemInY;
  int numProcsRecvMesh = recvNumElemInX*recvNumElemInY;

  SplitCommTestUtil sctu(numProcsSendMesh, numProcsRecvMesh);
  if (!sctu.get_status()) { GTEST_SKIP(); }

  setup(numProcsSendMesh, numProcsRecvMesh, sendNumElemInX, sendNumElemInY, recvNumElemInX, recvNumElemInY);

  unsigned expectedNumSearchPairs = 4u;
  test_search_result(expectedNumSearchPairs);

  stk::middle_mesh::impl::MeshScatterSpec sendMeshSpec;
  fill_mesh_scatter(sendMeshSpec);
}

}
