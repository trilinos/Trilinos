#include "gtest/gtest.h"

#include "stk_middle_mesh/communication_api.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"

using namespace stk::middle_mesh;

namespace {
 

class CommunicationAPIMPMDTest
{
  public:
    CommunicationAPIMPMDTest(MPI_Comm unionComm, MPI_Comm splitComm)
    {
      setup(unionComm, splitComm);
    }

    void runtest()
    {
      send_and_check_data_both_directions(true);
      send_and_check_data_both_directions(false);
    }

  private:

    void setup(MPI_Comm unionComm, MPI_Comm splitComm)
    {
      m_unionComm = unionComm;
      m_splitComm = splitComm;

      if (splitComm == MPI_COMM_NULL)
        return; 

      int commSize = utils::impl::comm_size(splitComm);
      mesh::impl::MeshSpec spec;
      spec.xmin = 0;          spec.ymin = 0;
      spec.xmax = 1;          spec.ymax = 1;
      spec.numelX = commSize, spec.numelY = commSize;

      auto f = [](const utils::Point& pt) { return pt; };
      m_mesh1 = create_mesh(spec, f, splitComm);
    }

    mesh::FieldPtr<mesh::RemoteSharedEntity> create_dest_field()
    {
      if (!m_mesh1)
        return nullptr;

      auto fieldPtr = mesh::create_field<mesh::RemoteSharedEntity>(m_mesh1, mesh::FieldShape(0, 0, 1), 1);

      //int destRank = utils::impl::comm_rank(m_splitComm) + utils::impl::comm_size(m_splitComm);
      int usedCommSize = 2*utils::impl::comm_size(m_splitComm);
      int destRank = (utils::impl::comm_rank(m_unionComm) + utils::impl::comm_size(m_splitComm)) % usedCommSize;
      int numEls = m_mesh1->get_elements().size();
      auto& field = *fieldPtr;
      for (auto el : m_mesh1->get_elements())
      {
        int destLocalId = (el->get_id() + 1) % numEls;
        field(el, 0, 0) = mesh::RemoteSharedEntity{destRank, destLocalId};
      }

      return fieldPtr;
    }

    mesh::FieldPtr<int> create_send_field(mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfoPtr)
    {
      if (!m_mesh1)
        return nullptr;

      auto fieldPtr = mesh::create_field<int>(m_mesh1, mesh::FieldShape(0, 0, m_numNodesPerElement), m_numComponentsPerNode);

      auto& remoteInfo = *remoteInfoPtr;
      auto& field = *fieldPtr;
      for (auto& el : m_mesh1->get_elements())
        if (el)
        {
          mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
          for (int i=0; i < m_numNodesPerElement; ++i)
            for (int j=0; j < m_numComponentsPerNode; ++j)
              field(el, i, j) = get_unique_value(remote.remoteRank, remote.remoteId, i, j);
        }

      return fieldPtr;
    }

    void send_and_check_data_both_directions(bool useFirstArg)
    {
      auto remoteInfo = create_dest_field();
      std::shared_ptr<mesh::Mesh> mesh1 = useFirstArg ? m_mesh1 : nullptr;
      std::shared_ptr<mesh::Mesh> mesh2 = useFirstArg ? nullptr : m_mesh1;
      mesh::FieldPtr<mesh::RemoteSharedEntity> remote1 = useFirstArg ? remoteInfo : nullptr;
      mesh::FieldPtr<mesh::RemoteSharedEntity> remote2 = useFirstArg ? nullptr    : remoteInfo;

      stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(m_unionComm, mesh1, mesh2, remote1, remote2);

      bool amISenderFirstRound = false;
      if (m_splitComm != MPI_COMM_NULL)
        amISenderFirstRound = utils::impl::comm_rank(m_unionComm) < utils::impl::comm_size(m_splitComm);
        
      send_and_check_data(exchanger, remoteInfo,  amISenderFirstRound);
      send_and_check_data(exchanger, remoteInfo, !amISenderFirstRound);
    }

    void send_and_check_data(stk::middle_mesh::MiddleMeshFieldCommunication<int>& exchanger, 
                     mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo, bool amISender)
    {
      mesh::FieldPtr<int> fieldSend, fieldRecv;
      if (amISender) {
        fieldSend = create_send_field(remoteInfo);
        fieldRecv = nullptr;
      } else if (m_splitComm != MPI_COMM_NULL) {
        fieldSend = nullptr;
        fieldRecv = mesh::create_field<int>(m_mesh1, mesh::FieldShape(0, 0, m_numNodesPerElement), m_numComponentsPerNode, -1);
      } else {
        fieldSend = nullptr;
        fieldRecv = nullptr;
      }

      exchanger.start_exchange(fieldSend, fieldRecv);
      exchanger.finish_exchange(fieldRecv);

      if (!amISender)
        check_recv_field(fieldRecv);
    }

    void check_recv_field(mesh::FieldPtr<int> fieldPtr)
    {
      if (!m_mesh1)
        return;

      //TODO: check for extrannious processes
      int myRank = utils::impl::comm_rank(m_unionComm);
      auto& field = *fieldPtr;
      for (auto el : m_mesh1->get_elements())
        for (int i=0; i < m_numNodesPerElement; ++i)
          for (int j=0; j < m_numComponentsPerNode; ++j)
            EXPECT_EQ(field(el, i, j), get_unique_value(myRank, el->get_id(), i, j));
    }

    int get_unique_value(int rank, int elId, int node, int component)
    {
      return rank + elId * m_numNodesPerElement * m_numComponentsPerNode + node * m_numComponentsPerNode + component;
    }

  private:
    MPI_Comm m_unionComm;
    MPI_Comm m_splitComm;
    std::shared_ptr<mesh::Mesh> m_mesh1;

    const int m_numNodesPerElement   = 2;
    const int m_numComponentsPerNode = 3;
};

}

TEST(CommunicationAPIMPMD, Values)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) == 1)
    GTEST_SKIP();

  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  int commRank = utils::impl::comm_rank(MPI_COMM_WORLD);
  
  for (int usedCommSize=2; usedCommSize <= commSize; usedCommSize += 2)
  {
    int color;
    if (commRank < usedCommSize/2)
      color = 0;
    else if (commRank >= usedCommSize/2 && commRank < usedCommSize)
      color = 1;
    else
      color = MPI_UNDEFINED;
    /*
    int color = commRank < i/2 ? 0 : 1;
    if (commSize % 2 != 0 && commRank == (commSize - 1))
      color = MPI_UNDEFINED;
    */

    MPI_Comm splitComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &splitComm);

    CommunicationAPIMPMDTest tester(MPI_COMM_WORLD, splitComm);
    tester.runtest();

    if (splitComm != MPI_COMM_NULL)
      MPI_Comm_free(&splitComm);
  }
}


TEST(CommunicationAPIMPMD, RemoteInfoWrongShape)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.xmin = 0;          spec.ymin = 0;
  spec.xmax = 1;          spec.ymax = 1;
  spec.numelX = 2, spec.numelY = 2;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = create_mesh(spec, f, MPI_COMM_SELF);

  auto remoteInfo = mesh::create_field<mesh::RemoteSharedEntity>(mesh, mesh::FieldShape(1, 0, 0), 1);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(MPI_COMM_WORLD, mesh, nullptr, remoteInfo, nullptr));

  remoteInfo = mesh::create_field<mesh::RemoteSharedEntity>(mesh, mesh::FieldShape(0, 0, 1), 2);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(MPI_COMM_WORLD, mesh, nullptr, remoteInfo, nullptr));
}

TEST(CommunicationAPIMPMD, RemoteInfoNotOnMiddleMesh)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.xmin = 0;          spec.ymin = 0;
  spec.xmax = 1;          spec.ymax = 1;
  spec.numelX = 2, spec.numelY = 2;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh = create_mesh(spec, f, MPI_COMM_SELF);
  auto mesh2 = create_mesh(spec, f, MPI_COMM_SELF);

  auto remoteInfo = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(0, 0, 1), 1);

  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(MPI_COMM_WORLD, mesh, nullptr, remoteInfo, nullptr));
}