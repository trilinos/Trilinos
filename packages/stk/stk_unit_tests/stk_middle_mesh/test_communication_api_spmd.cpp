#include "gtest/gtest.h"

#include "stk_middle_mesh/communication_api.hpp"
#include "stk_middle_mesh/create_mesh.hpp"

using namespace stk::middle_mesh;

namespace {


class CommunicationAPISPMDTest
{
  public:
    CommunicationAPISPMDTest(MPI_Comm comm)
    {
      setup(comm);
    }

    void runtest()
    {
      send_and_check_data_both_directions();
    }

  private:

    void setup(MPI_Comm comm)
    {
      m_comm = comm;

      int commSize = utils::impl::comm_size(comm);
      mesh::impl::MeshSpec spec;
      spec.xmin = 0;          spec.ymin = 0;
      spec.xmax = 1;          spec.ymax = 1;
      spec.numelX = commSize, spec.numelY = commSize;

      auto f = [](const utils::Point& pt) { return pt; };
      m_mesh1 = create_mesh(spec, f, comm);
      m_mesh2 = create_mesh(spec, f, comm);
    }

    mesh::FieldPtr<mesh::RemoteSharedEntity> create_dest_field(std::shared_ptr<mesh::Mesh> mesh)
    {
      auto fieldPtr = mesh::create_field<mesh::RemoteSharedEntity>(mesh, mesh::FieldShape(0, 0, 1), 1);

      int destRank = utils::impl::comm_rank(m_comm);
      int numEls = mesh->get_elements().size();
      auto& field = *fieldPtr;
      for (auto el : mesh->get_elements())
      {
        int destLocalId = (el->get_id() + 1) % numEls;
        field(el, 0, 0) = mesh::RemoteSharedEntity{destRank, destLocalId};
      }

      return fieldPtr;
    }

    mesh::FieldPtr<int> create_send_field(std::shared_ptr<mesh::Mesh> mesh, mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfoPtr)
    {
      auto fieldPtr = mesh::create_field<int>(mesh, mesh::FieldShape(0, 0, m_numNodesPerElement), m_numComponentsPerNode);

      auto& remoteInfo = *remoteInfoPtr;
      auto& field = *fieldPtr;
      for (auto& el : mesh->get_elements())
        if (el)
        {
          mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
          for (int i=0; i < m_numNodesPerElement; ++i)
            for (int j=0; j < m_numComponentsPerNode; ++j)
              field(el, i, j) = get_unique_value(remote.remoteRank, remote.remoteId, i, j);
        }

      return fieldPtr;
    }

    void send_and_check_data_both_directions()
    {
      auto remoteInfo1 = create_dest_field(m_mesh1);
      auto remoteInfo2 = create_dest_field(m_mesh2);
      stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(m_comm, m_mesh1, m_mesh2, remoteInfo1, remoteInfo2);
        
      send_and_check_data(exchanger, remoteInfo1, remoteInfo2, true);
      send_and_check_data(exchanger, remoteInfo1, remoteInfo2, false);
    }

    void send_and_check_data(stk::middle_mesh::MiddleMeshFieldCommunication<int>& exchanger, 
                             mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo1, 
                             mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo2, 
                             bool sendOneToTwo)
    {
      auto sendMesh = sendOneToTwo ? m_mesh1 : m_mesh2;
      auto recvMesh = sendOneToTwo ? m_mesh2 : m_mesh1;
      auto remoteInfoSend = sendOneToTwo ? remoteInfo1 : remoteInfo2;
      auto remoteInfoRecv = sendOneToTwo ? remoteInfo2 : remoteInfo1;



      mesh::FieldPtr<int> fieldSend = create_send_field(sendMesh, remoteInfoSend);
      mesh::FieldPtr<int> fieldRecv = mesh::create_field<int>(recvMesh, mesh::FieldShape(0, 0, m_numNodesPerElement),
                                                              m_numComponentsPerNode, -1);

      exchanger.start_exchange(fieldSend, fieldRecv);
      exchanger.finish_exchange(fieldRecv);
      check_recv_field(fieldRecv);
    }

    void check_recv_field(mesh::FieldPtr<int> fieldPtr)
    {
      int myRank = utils::impl::comm_rank(m_comm);
      auto& field = *fieldPtr;
      for (auto el : fieldPtr->get_mesh()->get_elements())
        for (int i=0; i < m_numNodesPerElement; ++i)
          for (int j=0; j < m_numComponentsPerNode; ++j)
            EXPECT_EQ(field(el, i, j), get_unique_value(myRank, el->get_id(), i, j));
    }

    int get_unique_value(int rank, int elId, int node, int component)
    {
      return rank + elId * m_numNodesPerElement * m_numComponentsPerNode + node * m_numComponentsPerNode + component;
    }

  private:
    MPI_Comm m_comm;
    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;

    const int m_numNodesPerElement   = 2;
    const int m_numComponentsPerNode = 3;
};

}

TEST(CommunicationAPISPMD, Values)
{
  CommunicationAPISPMDTest tester(MPI_COMM_WORLD);
  tester.runtest();
}

TEST(CommunicationAPISPMD, RemoteInfoWrongShape)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MPI_Comm comm = MPI_COMM_WORLD;
  int commSize = utils::impl::comm_size(comm);
  mesh::impl::MeshSpec spec;
  spec.xmin = 0;          spec.ymin = 0;
  spec.xmax = 1;          spec.ymax = 1;
  spec.numelX = commSize, spec.numelY = commSize;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh1 = create_mesh(spec, f, comm);
  auto mesh2 = create_mesh(spec, f, comm);

  auto remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(1, 0, 0), 1);
  auto remoteInfo2 = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(0, 0, 1), 1);

  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2));

  remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(0, 1, 0), 1);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2));

  remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(0, 0, 1), 1);
  remoteInfo2 = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(1, 0, 0), 1);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2));

  remoteInfo2 = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(0, 1, 0), 1);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2));

  remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(0, 0, 1), 2);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2));

  remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(0, 0, 1), 1);
  remoteInfo2 = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(0, 0, 1), 2);
  EXPECT_ANY_THROW(stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2));
}


TEST(CommunicationAPISPMD, FieldsWrongShape)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MPI_Comm comm = MPI_COMM_WORLD;
  int commSize = utils::impl::comm_size(comm);
  mesh::impl::MeshSpec spec;
  spec.xmin = 0;          spec.ymin = 0;
  spec.xmax = 1;          spec.ymax = 1;
  spec.numelX = commSize, spec.numelY = commSize;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh1 = create_mesh(spec, f, comm);
  auto mesh2 = create_mesh(spec, f, comm);

  auto remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(0, 0, 1), 1);
  auto remoteInfo2 = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(0, 0, 1), 1);

  stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2);

  auto field1 = mesh::create_field<int>(mesh1, mesh::FieldShape(1, 0, 0), 1);
  auto field2 = mesh::create_field<int>(mesh2, mesh::FieldShape(0, 0, 1), 1);

  EXPECT_ANY_THROW(exchanger.start_exchange(field1, field2));

  field1 = mesh::create_field<int>(mesh1, mesh::FieldShape(0, 1, 0), 1);
  EXPECT_ANY_THROW(exchanger.start_exchange(field1, field2));

  field1 = mesh::create_field<int>(mesh1, mesh::FieldShape(0, 0, 1), 1);
  field2 = mesh::create_field<int>(mesh2, mesh::FieldShape(1, 0, 0), 1);
  EXPECT_ANY_THROW(exchanger.start_exchange(field1, field2));

  field2 = mesh::create_field<int>(mesh2, mesh::FieldShape(0, 1, 0), 1);
  EXPECT_ANY_THROW(exchanger.start_exchange(field1, field2));
}


TEST(CommunicationAPISPMD, FieldsDifferentShapes)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MPI_Comm comm = MPI_COMM_WORLD;
  int commSize = utils::impl::comm_size(comm);
  mesh::impl::MeshSpec spec;
  spec.xmin = 0;          spec.ymin = 0;
  spec.xmax = 1;          spec.ymax = 1;
  spec.numelX = commSize, spec.numelY = commSize;

  auto f = [](const utils::Point& pt) { return pt; };
  auto mesh1 = create_mesh(spec, f, comm);
  auto mesh2 = create_mesh(spec, f, comm);

  auto remoteInfo1 = mesh::create_field<mesh::RemoteSharedEntity>(mesh1, mesh::FieldShape(0, 0, 1), 1);
  auto remoteInfo2 = mesh::create_field<mesh::RemoteSharedEntity>(mesh2, mesh::FieldShape(0, 0, 1), 1);

  stk::middle_mesh::MiddleMeshFieldCommunication<int> exchanger(comm, mesh1, mesh2, remoteInfo1, remoteInfo2);

  auto field1 = mesh::create_field<int>(mesh1, mesh::FieldShape(0, 0, 1), 2);
  auto field2 = mesh::create_field<int>(mesh2, mesh::FieldShape(0, 0, 1), 1);
  EXPECT_ANY_THROW(exchanger.start_exchange(field1, field2));

  field1 = mesh::create_field<int>(mesh1, mesh::FieldShape(0, 0, 2), 1);
  field2 = mesh::create_field<int>(mesh2, mesh::FieldShape(0, 0, 1), 1);
  EXPECT_ANY_THROW(exchanger.start_exchange(field1, field2));  
}