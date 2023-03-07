#ifndef COMMUNICATION_API_MPMD_H
#define COMMUNICATION_API_MPMD_H

#include "mesh.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

#include <stdexcept>
#include <vector>

namespace stk {
namespace middle_mesh {


template <typename T>
class MiddleMeshFieldCommunicationMPMD
{
  public:
    MiddleMeshFieldCommunicationMPMD(MPI_Comm unionComm, std::shared_ptr<mesh::Mesh> middleMesh,
                                     mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo);

    void start_exchange(mesh::FieldPtr<T> field, bool amISender);

    void finish_exchange(mesh::FieldPtr<T> field, bool amISender);

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

    void check_field_shapes_same_debug_only(mesh::FieldPtr<T> field);

    void check_send_direction_debug_only(bool amISender);

    void set_recv_buffer_sizes(mesh::FieldShape fshape, int numCompPerNode);

    void pack_send_buffers(mesh::FieldPtr<T> fieldPtr);

    void complete_receives(mesh::FieldPtr<T> fieldPtr);

    void unpack_buffer(int rank, const std::vector<DataPlusId>& buf, mesh::FieldPtr<T> fieldPtr);


    MPI_Comm m_unionComm;
    std::shared_ptr<mesh::Mesh> m_middleMesh;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo;
    stk::DataExchangeKnownPatternNonBlockingBuffer<DataPlusId> m_exchanger;
};


template <typename T>
MiddleMeshFieldCommunicationMPMD<T>::MiddleMeshFieldCommunicationMPMD(MPI_Comm unionComm,
                                       std::shared_ptr<mesh::Mesh> middleMesh,
                                       mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo) :
  m_unionComm(unionComm),
  m_middleMesh(middleMesh),
  m_remoteInfo(remoteInfo),
  m_exchanger(unionComm)
{
  if (remoteInfo)
  {
    for (int i=0; i < 3; ++i)
      ThrowRequire(remoteInfo->get_field_shape().count[i] == mesh::FieldShape(0, 0, 1).count[i]);
    ThrowRequire(remoteInfo->get_num_comp() == 1);
    ThrowRequire(remoteInfo->get_mesh() == middleMesh);
  }
}


template <typename T>
void MiddleMeshFieldCommunicationMPMD<T>::start_exchange(mesh::FieldPtr<T> field, bool amISender)
{
  check_fieldshape(field);
  if (m_middleMesh) {
    ThrowRequireMsg(field->get_mesh() == m_middleMesh, "field must be defined on the middle mesh");
  }
  check_field_shapes_same_debug_only(field);
  check_send_direction_debug_only(amISender);


  m_exchanger.clear_send_bufs();
  m_exchanger.clear_recv_bufs();
  if (amISender) {
    pack_send_buffers(field);
  } else if (field) {
    set_recv_buffer_sizes(field->get_field_shape(), field->get_num_comp());
  }

  m_exchanger.start_nonblocking();
}


template <typename T>
void MiddleMeshFieldCommunicationMPMD<T>::finish_exchange(mesh::FieldPtr<T> field, bool amISender)
{
  complete_receives(field);
}
 

template <typename T>
void MiddleMeshFieldCommunicationMPMD<T>::check_fieldshape(mesh::FieldPtr<T> field)
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
void MiddleMeshFieldCommunicationMPMD<T>::check_field_shapes_same_debug_only(mesh::FieldPtr<T> field)
{
#ifndef NDEBUG  
  const int root = 0;
  const int nprocs = utils::impl::comm_size(m_unionComm);
  mesh::FieldShape fshape;

  if (field) {
    fshape = field->get_field_shape();
  } else
  {
    fshape = mesh::FieldShape(-1, -1, -1);
  }

  int nvals = fshape.count.size();
  int recvBufSize = utils::impl::comm_rank(m_unionComm) == root ? nprocs : 0;
  std::vector<mesh::FieldShape> recvBuf(recvBufSize);
  MPI_Gather(fshape.count.data(), nvals, MPI_INT, recvBuf.data(), nvals, MPI_INT, root, m_unionComm);

  if (utils::impl::comm_rank(m_unionComm) == root)
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
void MiddleMeshFieldCommunicationMPMD<T>::check_send_direction_debug_only(bool amISender)
{
#ifndef NDEBUG
  enum class SenderStatus : int
  {
    Sender,
    Receiver,
    Neither
  };

  SenderStatus myStatus;
  if (amISender && m_middleMesh)
    myStatus = SenderStatus::Sender;
  else if (!amISender && m_middleMesh)
    myStatus = SenderStatus::Receiver;
  else
    myStatus = SenderStatus::Neither;

  std::vector<SenderStatus> statuses(utils::impl::comm_size(m_unionComm));

  MPI_Allgather(&myStatus, 1, MPI_INT, statuses.data(), 1, MPI_INT, m_unionComm);

  if (myStatus != SenderStatus::Neither)
  {
    MPI_Group unionCommGroup, meshCommGroup;
    MPI_Comm_group(m_unionComm, &unionCommGroup);
    MPI_Comm_group(m_middleMesh->get_comm(), &meshCommGroup);
    for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
      if (statuses[rank] != SenderStatus::Neither)
      {
        int rankOnMeshComm;
        MPI_Group_translate_ranks(unionCommGroup, 1, &rank, meshCommGroup, &rankOnMeshComm);
        if (rankOnMeshComm == MPI_UNDEFINED)
        {
          ThrowRequireMsg(statuses[rank] != myStatus, "Both application cannot be the sender");
        } else
        {
          ThrowRequireMsg(statuses[rank] == myStatus, "All processes on a given appplication must give the same value for amISender");
        }
      }
  }

#endif
}

template <typename T>
void MiddleMeshFieldCommunicationMPMD<T>::set_recv_buffer_sizes(mesh::FieldShape fshape, int numCompPerNode)
{
  if (!m_middleMesh)
  {
    return;
  }

  std::vector<int> recvCounts(utils::impl::comm_size(m_unionComm), 0);
  auto& remoteInfo = *m_remoteInfo;
  int numNodesPerElement = fshape.count[2];
  for (auto el : m_middleMesh->get_elements())
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
void MiddleMeshFieldCommunicationMPMD<T>::pack_send_buffers(mesh::FieldPtr<T> fieldPtr)
{
  if (!m_middleMesh)
    return;

  auto& field = *fieldPtr;
  auto& remoteInfo = *m_remoteInfo;
  int numNodesPerElement = field.get_field_shape().count[2];
  int numCompPerNode = field.get_num_comp();
  for (auto el : m_middleMesh->get_elements())
    if (el)
    {
      mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
      for (int i=0; i < numNodesPerElement; ++i)
        for (int j=0; j < numCompPerNode; ++j)
        {
          m_exchanger.get_send_buf(remote.remoteRank).emplace_back(field(el, i, j), remote.remoteId);
        }
    }
}


template <typename T>
void MiddleMeshFieldCommunicationMPMD<T>::complete_receives(mesh::FieldPtr<T> fieldPtr)
{
  auto f = [&](int rank, const std::vector<DataPlusId>& buf)
  {
    unpack_buffer(rank, buf, fieldPtr);
  };

  m_exchanger.complete_receives(f);
  m_exchanger.complete_sends();
}


template <typename T>
void MiddleMeshFieldCommunicationMPMD<T>::unpack_buffer(int rank, const std::vector<DataPlusId>& buf, mesh::FieldPtr<T> fieldPtr)
{
  int numNodesPerElement = fieldPtr->get_field_shape().count[2];
  int numCompPerNode = fieldPtr->get_num_comp();
  assert(buf.size() % (numNodesPerElement * numCompPerNode) == 0);

  auto& field = *fieldPtr;
  size_t idx = 0;
  while (idx < buf.size())
  {
    auto el = m_middleMesh->get_elements()[buf[idx].localId];
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