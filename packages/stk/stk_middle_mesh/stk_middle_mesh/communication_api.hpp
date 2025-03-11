#ifndef STK_MIDDLE_MESH_COMMUNICATION_API_H
#define STK_MIDDLE_MESH_COMMUNICATION_API_H

#include "mesh.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

#include <stdexcept>
#include <vector>

namespace stk {
namespace middle_mesh {

template <typename T>
class MiddleMeshFieldCommunication
{
  public:
    MiddleMeshFieldCommunication(MPI_Comm unionComm, std::shared_ptr<mesh::Mesh> middleMesh1, std::shared_ptr<mesh::Mesh> middleMesh2,
                                 mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo1, mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo2);

    void start_exchange(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv);

    void finish_exchange(mesh::FieldPtr<T> fieldRecv);

  private: 

    void check_fieldshape(mesh::FieldPtr<T> field);

    void check_field_shapes_same_locally(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv);

    void check_all_procs_on_field_provided_argument_debug_only(mesh::FieldPtr<T> field);

    void check_field_shapes_same_globally_debug_only(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv);

    mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info(mesh::FieldPtr<T> field);

    void set_recv_buffer_sizes(mesh::FieldPtr<T> fieldRecv);

    void pack_send_buffers(mesh::FieldPtr<T> fieldSendPtr);

    void complete_receives(mesh::FieldPtr<T> fieldRecvPtr);

    void unpack_buffer(stk::CommBuffer& buf, mesh::FieldPtr<T> fieldPtr);


    MPI_Comm m_unionComm;
    std::shared_ptr<mesh::Mesh> m_middleMesh1;
    std::shared_ptr<mesh::Mesh> m_middleMesh2;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo1;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo2;    
    stk::DataExchangeKnownPatternNonBlockingCommBuffer m_exchanger;
};


template <typename T>
MiddleMeshFieldCommunication<T>::MiddleMeshFieldCommunication(MPI_Comm unionComm,
                                    std::shared_ptr<mesh::Mesh> middleMesh1, std::shared_ptr<mesh::Mesh> middleMesh2,
                                    mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo1, mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfo2) :
  m_unionComm(unionComm),
  m_middleMesh1(middleMesh1),
  m_middleMesh2(middleMesh2),
  m_remoteInfo1(remoteInfo1),
  m_remoteInfo2(remoteInfo2),
  m_exchanger(unionComm)
{
  if (middleMesh1)
  {
    STK_ThrowRequireMsg(remoteInfo1, "When passing in middleMesh1, must also pass in remoteInfo1");
    STK_ThrowRequireMsg(remoteInfo1->get_mesh() == middleMesh1, "remoteInfo1 must be defined on middleMesh1");
  }


  if (middleMesh2)
  {
    STK_ThrowRequireMsg(remoteInfo2, "When passing in middleMesh2, must also pass in remoteInfo2");
    STK_ThrowRequireMsg(remoteInfo2->get_mesh() == middleMesh2, "remoteInfo1 must be defined on middleMesh1");
  }
      
  for (int i=0; i < 3; ++i)
  {
    if (remoteInfo1)
      STK_ThrowRequire(remoteInfo1->get_field_shape() == mesh::FieldShape(0, 0, 1));
    
    if (remoteInfo2)
      STK_ThrowRequire(remoteInfo2->get_field_shape() == mesh::FieldShape(0, 0, 1));
  }

  if (remoteInfo1)
    STK_ThrowRequire(remoteInfo1->get_num_comp() == 1);

  if (remoteInfo2)
    STK_ThrowRequire(remoteInfo2->get_num_comp() == 1);
}


template <typename T>
void MiddleMeshFieldCommunication<T>::start_exchange(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv)
{
  check_all_procs_on_field_provided_argument_debug_only(fieldSend);
  check_all_procs_on_field_provided_argument_debug_only(fieldRecv);
  check_fieldshape(fieldSend);
  check_fieldshape(fieldRecv);
  check_field_shapes_same_locally(fieldSend, fieldRecv);
  check_field_shapes_same_globally_debug_only(fieldSend, fieldRecv); 

  m_exchanger.clear_send_bufs();
  m_exchanger.clear_recv_bufs();

  pack_send_buffers(fieldSend);
  set_recv_buffer_sizes(fieldRecv);

  m_exchanger.start_nonblocking();
}


template <typename T>
void MiddleMeshFieldCommunication<T>::finish_exchange(mesh::FieldPtr<T> fieldRecv)
{
  complete_receives(fieldRecv);
}
 

template <typename T>
void MiddleMeshFieldCommunication<T>::check_fieldshape(mesh::FieldPtr<T> field)
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
void MiddleMeshFieldCommunication<T>::check_field_shapes_same_locally(mesh::FieldPtr<T> fieldSend, mesh::FieldPtr<T> fieldRecv)
{
  if (fieldSend && fieldRecv)
  {
    STK_ThrowRequireMsg(fieldSend->get_field_shape() == fieldRecv->get_field_shape(), "fieldSend and fieldRecv must have same FieldShape");
    STK_ThrowRequireMsg(fieldSend->get_num_comp() == fieldRecv->get_num_comp(), "fieldSend and FieldRecv must have same number of components");
  }
}

template <typename T>
void MiddleMeshFieldCommunication<T>::check_all_procs_on_field_provided_argument_debug_only([[maybe_unused]] mesh::FieldPtr<T> field)
{
#ifndef NDEBUG
  if (field)
  {
    MPI_Comm comm = field->get_mesh()->get_comm();
    // this will hang if some processes on comm provided the field but others didn't
    MPI_Barrier(comm);
  }
#endif
}



template <typename T>
void MiddleMeshFieldCommunication<T>::check_field_shapes_same_globally_debug_only([[maybe_unused]] mesh::FieldPtr<T> fieldSend, [[maybe_unused]] mesh::FieldPtr<T> fieldRecv)
{
#ifndef NDEBUG
  MPI_Comm comm = m_unionComm;
  const int root = 0;
  const int nprocs = utils::impl::comm_size(comm);
  mesh::FieldShape fshape;

  if (fieldSend || fieldRecv) {
    if (fieldSend && fieldRecv)
      STK_ThrowRequireMsg(fieldSend->get_field_shape() == fieldRecv->get_field_shape(), "send and receive fields must have same FieldShape");

    
    fshape = fieldSend ? fieldSend->get_field_shape() : fieldRecv->get_field_shape();
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
      if (recvBuf[i].count[0] != -1 && recvBuf[i] != fshapeValid)
        throw std::runtime_error("fieldshape must be the same on all ranks");
  }
#endif
}

template <typename T>
mesh::FieldPtr<mesh::RemoteSharedEntity> MiddleMeshFieldCommunication<T>::get_remote_info(mesh::FieldPtr<T> field)
{
  if (m_remoteInfo1 && field->get_mesh() == m_remoteInfo1->get_mesh())
  {
    return m_remoteInfo1;
  } else if (m_remoteInfo2 && field->get_mesh() == m_remoteInfo2->get_mesh())
  {
    return m_remoteInfo2;
  } else
  {
    throw std::runtime_error("cannot determine remoteInfo");
  }
}

template <typename T>
void MiddleMeshFieldCommunication<T>::set_recv_buffer_sizes(mesh::FieldPtr<T> fieldRecv)
{
  int commSize = utils::impl::comm_size(m_unionComm);
  if (!fieldRecv)
  {
    for (int rank=0; rank < commSize; ++rank)
      m_exchanger.set_recv_buffer_size(rank, 0);
    
    m_exchanger.allocate_recv_buffers();
    return;
  }

  std::vector<int> recvCounts(commSize, 0);
  auto& remoteInfo = *(get_remote_info(fieldRecv));
  for (auto el : fieldRecv->get_mesh()->get_elements())
    if (el)
    {
      mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
      recvCounts[remote.remoteRank]++;
    }

  mesh::FieldShape fshape = fieldRecv->get_field_shape();
  int numCompPerNode            = fieldRecv->get_num_comp();
  int numNodesPerElement        = fshape.count[2];
  for (size_t rank=0; rank < recvCounts.size(); ++rank)
  {
    auto& buf = m_exchanger.get_recv_buf(rank);
    for (int i=0; i < recvCounts[rank]; ++i)
    {
      buf.template pack<int>(0);
      for (int j=0; j < numNodesPerElement; ++j)
        for (int k=0; k < numCompPerNode; ++k)
          buf.pack(T());
    }

    m_exchanger.set_recv_buffer_size(rank, buf.size());
  }

  m_exchanger.allocate_recv_buffers();
}

template <typename T>
void MiddleMeshFieldCommunication<T>::pack_send_buffers(mesh::FieldPtr<T> fieldSendPtr)
{
  if (!fieldSendPtr)
    return;

  auto& field            = *fieldSendPtr;
  auto& remoteInfo       = *(get_remote_info(fieldSendPtr));
  int numNodesPerElement = field.get_field_shape().count[2];
  int numCompPerNode     = field.get_num_comp();

  for (int phase=0; phase < 2; ++phase)
  {
    for (auto el : field.get_mesh()->get_elements())
      if (el)
      {
        mesh::RemoteSharedEntity remote = remoteInfo(el, 0, 0);
        auto& buf = m_exchanger.get_send_buf(remote.remoteRank);

        buf.pack(remote.remoteId);
        for (int i=0; i < numNodesPerElement; ++i)
          for (int j=0; j < numCompPerNode; ++j)
            buf.pack(field(el, i, j));
      }

    if (phase == 0)
    {
      m_exchanger.allocate_send_buffers();
    }
  }
}


template <typename T>
void MiddleMeshFieldCommunication<T>::complete_receives(mesh::FieldPtr<T> fieldRecvPtr)
{
  auto f = [&](int /*rank*/, stk::CommBuffer& buf)
  {
    unpack_buffer(buf, fieldRecvPtr);
  };

  m_exchanger.complete_receives(f);
  m_exchanger.complete_sends();
}


template <typename T>
void MiddleMeshFieldCommunication<T>::unpack_buffer(stk::CommBuffer& buf, mesh::FieldPtr<T> fieldRecvPtr)
{
  std::shared_ptr<mesh::Mesh> middleMeshRecv = fieldRecvPtr->get_mesh();
  int numNodesPerElement = fieldRecvPtr->get_field_shape().count[2];
  int numCompPerNode = fieldRecvPtr->get_num_comp();

  auto& field = *fieldRecvPtr;
  while (buf.remaining() > 0)
  {
    int localId;
    buf.unpack(localId);
    auto el = middleMeshRecv->get_elements()[localId];
    for (int i=0; i < numNodesPerElement; ++i)
      for (int j=0; j < numCompPerNode; ++j)
      {
        T val;  
        buf.unpack(val);
        field(el, i, j) = val;
      }
  }
}

}  // namespace
}  // namespace

#endif