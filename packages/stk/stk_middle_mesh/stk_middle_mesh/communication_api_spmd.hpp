#ifndef COMMUNICATION_API_SPMD_H
#define COMMUNICATION_API_SPMD_H

#include "mesh.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

#include <stdexcept>
#include <vector>

namespace stk {
namespace middle_mesh {

template <typename T>
class MiddleMeshFieldCommunicationSPMD
{
  public:
    MiddleMeshFieldCommunicationSPMD(std::shared_ptr<mesh::Mesh> middleMesh1, std::shared_ptr<mesh::Mesh> middleMesh2,
                                     mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo1, mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo2);

    void start_exchange(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv);

    void finish_exchange(mesh::FieldPtr<T> fieldRecv);

  private:

    struct DataPlusId
    {
      DataPlusId(const T& data_=T(), int localId_=-1) :
        data(data_),
        localId(localId_)
      {}

      T data;
      int localId;
    };  

    void check_fieldshape(mesh::FieldPtr<T> field);

    void check_field_shapes_same_locally(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv);

    void check_field_shapes_same_globally_debug_only(mesh::FieldPtr<T> field);

    void set_recv_buffer_sizes(mesh::FieldPtr<T> fieldRecv);

    void pack_send_buffers(mesh::FieldPtr<T> fieldSendPtr);

    void complete_receives(mesh::FieldPtr<T> fieldRecvPtr);

    void unpack_buffer(int rank, const std::vector<DataPlusId>& buf, mesh::FieldPtr<T> fieldPtr);


    MPI_Comm m_unionComm;
    std::shared_ptr<mesh::Mesh> m_middleMesh1;
    std::shared_ptr<mesh::Mesh> m_middleMesh2;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo1;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo2;    
    stk::DataExchangeKnownPatternNonBlockingBuffer<DataPlusId> m_exchanger;
};


template <typename T>
MiddleMeshFieldCommunicationSPMD<T>::MiddleMeshFieldCommunicationSPMD(
                                    std::shared_ptr<mesh::Mesh> middleMesh1, std::shared_ptr<mesh::Mesh> middleMesh2,
                                    mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo1, mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo2) :
  m_middleMesh1(middleMesh1),
  m_middleMesh2(middleMesh2),
  m_remoteInfo1(remoteInfo1),
  m_remoteInfo2(remoteInfo2),
  m_exchanger(middleMesh1->get_comm())
{
  for (int i=0; i < 3; ++i)
  {
    ThrowRequire(remoteInfo1->get_field_shape().count[i] == mesh::FieldShape(0, 0, 1).count[i]);
    ThrowRequire(remoteInfo2->get_field_shape().count[i] == mesh::FieldShape(0, 0, 1).count[i]);
  }
  ThrowRequire(remoteInfo1->get_num_comp() == 1);
  ThrowRequire(remoteInfo2->get_num_comp() == 1);

  ThrowRequireMsg(middleMesh1->get_comm() == middleMesh2->get_comm(), "Both middle meshes must be on same communicator");
}


template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::start_exchange(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv)
{
  check_fieldshape(fieldSend);
  check_fieldshape(fieldRecv);
  check_field_shapes_same_locally(fieldSend, fieldRecv);
  check_field_shapes_same_globally_debug_only(fieldSend);
  check_field_shapes_same_globally_debug_only(fieldRecv);
 

  m_exchanger.clear_send_bufs();
  m_exchanger.clear_recv_bufs();

  pack_send_buffers(fieldSend);
  set_recv_buffer_sizes(fieldRecv);

  m_exchanger.start_nonblocking();
}


template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::finish_exchange(mesh::FieldPtr<T> fieldRecv)
{
  complete_receives(fieldRecv);
}
 

template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::check_fieldshape(mesh::FieldPtr<T> field)
{
  if (!field)
    return;

  mesh::FieldShape fshape = field->get_field_shape();

  if (fshape.count[0] != 0)
    throw std::runtime_error("fields with data on vertices are not supported");

  if (fshape.count[1] != 0)
    throw std::runtime_error("fields with data on vertices are not supported");

  if (fshape.count[2] == 0)
    throw std::runtime_error("field must have at least 1 value on faces");
}

template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::check_field_shapes_same_locally(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv)
{
  for (size_t i=0; i < fieldSend->get_field_shape().count.size(); ++i)
  {
    ThrowRequireMsg(fieldSend->get_field_shape().count[i] == fieldRecv->get_field_shape().count[i], "fieldSend and fieldRecv must have same FieldShape");;
  }

  ThrowRequireMsg(fieldSend->get_num_comp() == fieldRecv->get_num_comp(), "fieldSend and FieldRecv must have same number of components");
}


template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::check_field_shapes_same_globally_debug_only(mesh::FieldPtr<T> field)
{
#ifndef NDEBUG
  MPI_Comm comm = field->get_mesh()->get_comm();
  const int root = 0;
  const int nprocs = utils::impl::comm_size(comm);
  mesh::FieldShape fshape;

  if (field) {
    fshape = field->get_field_shape();
  } else
  {
    fshape = mesh::FieldShape(-1, -1, -1);
  }

  int nvals = fshape.count.size();
  int recvBufSize = utils::impl::comm_rank(comm) == root ? nprocs : 0;
  std::vector<mesh::FieldShape> recvBuf(recvBufSize);
  MPI_Gather(fshape.count.data(), nvals, MPI_INT, recvBuf.data(), nvals, MPI_INT, root, comm);

  if (utils::impl::comm_rank(comm) == root)
  {
    mesh::FieldShape fshapeValid;
    for (auto& fshapeRecved : recvBuf)
      if (fshapeRecved.count[0] != -1)
      {
        fshapeValid = fshapeRecved;
        break;
      }

    for (int i=0; i < nprocs; ++i)
      for (int j=0; j < nvals; ++j)
        if (recvBuf[i].count[j] != -1 && recvBuf[i].count[j] != fshapeValid.count[j])
          throw std::runtime_error("fieldshape must be the same on all ranks");
  }
#endif
}

template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::set_recv_buffer_sizes(mesh::FieldPtr<T> fieldRecv)
{
  std::vector<int> recvCounts(utils::impl::comm_size(fieldRecv->get_mesh()->get_comm()), 0);
  auto& remoteInfo              = fieldRecv->get_mesh() == m_remoteInfo1->get_mesh() ? *m_remoteInfo1 : *m_remoteInfo2;
  mesh::FieldShape fshape = fieldRecv->get_field_shape();
  int numCompPerNode            = fieldRecv->get_num_comp();
  int numNodesPerElement        = fshape.count[2];
  for (auto el : fieldRecv->get_mesh()->get_elements())
    if (el)
    {
      mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
      recvCounts[remote.remoteRank] += numNodesPerElement * numCompPerNode;
    }

  for (size_t rank=0; rank < recvCounts.size(); ++rank)
  {
    m_exchanger.get_recv_buf(rank).resize(recvCounts[rank]);
  }
}

template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::pack_send_buffers(mesh::FieldPtr<T> fieldSendPtr)
{
  auto& field            = *fieldSendPtr;
  auto& remoteInfo       = field.get_mesh() == m_remoteInfo1->get_mesh() ? *m_remoteInfo1 : *m_remoteInfo2;
  int numNodesPerElement = field.get_field_shape().count[2];
  int numCompPerNode     = field.get_num_comp();
  for (auto el : field.get_mesh()->get_elements())
    if (el)
    {
      mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
      for (int i=0; i < numNodesPerElement; ++i)
        for (int j=0; j < numCompPerNode; ++j)
          m_exchanger.get_send_buf(remote.remoteRank).emplace_back(field(el, i, j), remote.remoteId);
    }
}


template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::complete_receives(mesh::FieldPtr<T> fieldRecvPtr)
{
  auto f = [&](int rank, const std::vector<DataPlusId>& buf)
  {
    unpack_buffer(rank, buf, fieldRecvPtr);
  };

  m_exchanger.complete_receives(f);
  m_exchanger.complete_sends();
}


template <typename T>
void MiddleMeshFieldCommunicationSPMD<T>::unpack_buffer(int rank, const std::vector<DataPlusId>& buf, mesh::FieldPtr<T> fieldRecvPtr)
{
  std::shared_ptr<mesh::Mesh> middleMeshRecv = fieldRecvPtr->get_mesh();
  int numNodesPerElement = fieldRecvPtr->get_field_shape().count[2];
  int numCompPerNode = fieldRecvPtr->get_num_comp();
  assert(buf.size() % (numNodesPerElement * numCompPerNode) == 0);

  auto& field = *fieldRecvPtr;
  size_t idx = 0;
  while (idx < buf.size())
  {
    auto el = middleMeshRecv->get_elements()[buf[idx].localId];
    for (int i=0; i < numNodesPerElement; ++i)
      for (int j=0; j < numCompPerNode; ++j)
      {      
        field(el, i, j) = buf[idx].data;
        idx++;
      }
  }
}

}  // namespace
}  // namespace

#endif