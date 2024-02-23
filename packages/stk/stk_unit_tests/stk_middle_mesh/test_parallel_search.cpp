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
#include "stk_middle_mesh/mesh_scatter_spec.hpp"
#include "stk_middle_mesh/bounding_box_search.hpp"
#include "stk_middle_mesh/create_mesh.hpp"


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

  std::vector<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox::BoundingBox> boundingBoxVec;
  stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox meshA(mesh, MPI_COMM_WORLD);
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
  using BoundingBoxSearch = stk::middle_mesh::search::ElementToElementBoundingBoxSearch;
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

    if (m_color == stk::middle_mesh::search::SplitCommColor::SEND) {
      m_mesh = stk::middle_mesh::mesh::impl::create_mesh(m_sendMeshSpec, func, m_splitComm, false);
    } else if (m_color == stk::middle_mesh::search::SplitCommColor::RECV) {
      m_mesh = stk::middle_mesh::mesh::impl::create_mesh(m_recvMeshSpec, func, m_splitComm, false);
    }
  }

  void setup_search()
  {
    if (m_color == stk::middle_mesh::search::SplitCommColor::SEND) {
      m_searchSendMesh = std::make_shared<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox>(m_mesh, MPI_COMM_WORLD);
      m_search = std::make_shared<stk::middle_mesh::search::ElementToElementBoundingBoxSearch>(
          m_searchSendMesh, m_searchRecvMesh, "BoundingBoxSearch", MPI_COMM_WORLD);
    } else if (m_color == stk::middle_mesh::search::SplitCommColor::RECV) {
      m_searchRecvMesh = std::make_shared<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox>(m_mesh, MPI_COMM_WORLD);
      m_search = std::make_shared<stk::middle_mesh::search::ElementToElementBoundingBoxSearch>(
          m_searchSendMesh, m_searchRecvMesh, "BoundingBoxSearch", MPI_COMM_WORLD);
    } else {
      STK_ThrowRequireMsg(false, "Invalid SplitComm color");
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

    const BoundingBoxSearch::EntityProcRelationVec& pairedEntities = m_search->get_range_to_domain();
    const std::vector<BoundingBoxSearch::BoundingBoxB>& unpairedEntities = m_search->get_unpaired_recv_entities();

    if(m_color == stk::middle_mesh::search::SplitCommColor::SEND) {
      EXPECT_EQ(expectedNumSearchPairs, pairedEntities.size());
      EXPECT_EQ(0u, unpairedEntities.size());
    }

    unsigned expectedGlobalRank = 0;

    for(const auto &pairedEntity : pairedEntities) {
      const stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox::EntityProc& recvEntityProc = pairedEntity.first;
      const stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox::EntityProc& sendEntityProc = pairedEntity.second;

      if(m_color == stk::middle_mesh::search::SplitCommColor::RECV) {
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

  void fill_mesh_scatter(stk::middle_mesh::mesh::impl::MeshScatterSpec& scatter)
  {
    const BoundingBoxSearch::EntityProcRelationVec& pairedEntities = m_search->get_range_to_domain();
    const std::vector<BoundingBoxSearch::BoundingBoxB>& unpairedEntities = m_search->get_unpaired_recv_entities();
    EXPECT_EQ(0u, unpairedEntities.size());

    if(m_color == stk::middle_mesh::search::SplitCommColor::SEND) {
      auto elements = m_search->send_mesh()->get_mesh()->get_elements();

      for(const auto &pairedEntity : pairedEntities) {
        const stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox::EntityProc& recvEntityProc = pairedEntity.first;
        const stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox::EntityProc& sendEntityProc = pairedEntity.second;

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
  stk::middle_mesh::search::SplitCommColor m_color{stk::middle_mesh::search::SplitCommColor::INVALID};

  MPI_Comm m_splitComm{MPI_COMM_NULL};

  stk::middle_mesh::mesh::impl::MeshSpec m_sendMeshSpec;
  stk::middle_mesh::mesh::impl::MeshSpec m_recvMeshSpec;

  std::shared_ptr<stk::middle_mesh::search::ElementToElementBoundingBoxSearch> m_search;

  std::shared_ptr<stk::middle_mesh::mesh::Mesh> m_mesh;

  std::shared_ptr<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox> m_searchSendMesh;
  std::shared_ptr<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox> m_searchRecvMesh;

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

  stk::middle_mesh::mesh::impl::MeshScatterSpec sendMeshSpec(MPI_COMM_WORLD, m_mesh);
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

  stk::middle_mesh::mesh::impl::MeshScatterSpec sendMeshSpec(MPI_COMM_WORLD, m_mesh);
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

  stk::middle_mesh::mesh::impl::MeshScatterSpec sendMeshSpec(MPI_COMM_WORLD, m_mesh);
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

  stk::middle_mesh::mesh::impl::MeshScatterSpec sendMeshSpec(MPI_COMM_WORLD, m_mesh);
  fill_mesh_scatter(sendMeshSpec);
}

TEST(SearchBoundingBox, Regression)
{
  if (stk::middle_mesh::utils::impl::comm_size(MPI_COMM_WORLD) != 1)
  {
    GTEST_SKIP();
  }

  auto mesh1 = stk::middle_mesh::mesh::make_empty_mesh();
  auto v1 = mesh1->create_vertex(1.051722092687432, 0.7641208279802152, 0.7483314773547881);
  auto v2 = mesh1->create_vertex(0.8816778784387097, 1.213525491562421, 0);
  auto v3 = mesh1->create_vertex(0.764120827980215, 1.051722092687432, 0.748331477354788);
  mesh1->create_triangle_from_verts(v1, v2, v3);

  auto mesh2 = stk::middle_mesh::mesh::make_empty_mesh();
  auto v4 = mesh2->create_vertex(0.8032742580699905, 1.105558767447774, 0.6184579612837334);
  auto v5 = mesh2->create_vertex(0.6833453909022464, 1.183422925761223, 0.6185056386059162);
  auto v6 = mesh2->create_vertex(0.7166834109365785, 1.241075309649816, 0.4428266523506837);
  mesh2->create_triangle_from_verts(v4, v5, v6);

  std::vector<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox::BoundingBox> boundingBoxVec1, boundingBoxVec2;
  stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox searchMesh1(mesh1, MPI_COMM_WORLD);
  searchMesh1.fill_bounding_boxes(boundingBoxVec1);  

  stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox searchMesh2(mesh2, MPI_COMM_WORLD);
  searchMesh2.fill_bounding_boxes(boundingBoxVec2);

  EXPECT_TRUE(stk::search::intersects(boundingBoxVec1[0].first, boundingBoxVec2[0].first)); 
}

}
