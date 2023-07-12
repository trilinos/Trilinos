#include "stk_util/parallel/ManagedBufferBase.hpp"

namespace stk {

void ManagedCommBufferBase::allocate_send_buffers()
{
  for (size_t i=0; i < m_sendBufs.size(); ++i) {
    auto& storage_i = m_sendBufStorage[i];
    storage_i.resize(m_sendBufs[i].size());
    m_sendBufs[i].set_buffer_ptrs(storage_i.data(), storage_i.data(), storage_i.data() + storage_i.size());
  }

  m_sendBuffersAllocated = true;
}

void ManagedCommBufferBase::clear_send_bufs()
{
  for (size_t i=0; i < m_sendBufs.size(); ++i) {
    m_sendBufs[i] = CommBuffer();
    m_sendBufStorage[i].clear();
  }

  m_sendBuffersAllocated = false;
}

void ManagedCommBufferBase::clear_recv_bufs()
{
  for (size_t i=0; i < m_recvBufs.size(); ++i) {
    m_recvBufs[i] = CommBuffer();
    m_recvBufStorage[i].clear();
  }
}

void ManagedCommBufferBase::deallocate_send_bufs()
{
  for (size_t i=0; i < m_sendBufs.size(); ++i) {
    get_send_buf(i).reset();
    std::vector<unsigned char> tmp;
    m_sendBufStorage[i].swap(tmp);
  }
  m_sendBuffersAllocated = false;
  m_sendBuffersDeallocated = true;
}

void ManagedCommBufferBase::set_recv_buffer_storage()
{
  for (size_t i=0; i < m_recvBufs.size(); ++i) {
    if (m_recvBufStorage[i].size() > 0) {
      set_recv_buffer_storage(i);
    }
  }
}

CommBuffer& ManagedCommBufferBase::set_recv_buffer_storage(int rank)
{
  auto& storage = m_recvBufStorage[rank];
  m_recvBufs[rank].set_buffer_ptrs(storage.data(), storage.data(), storage.data() + storage.size());
  return m_recvBufs[rank];
}

void ManagedCommBufferBase::reset_send_commbufs()
{
  for (size_t i=0; i < m_sendBufs.size(); ++i) {
    get_recv_buf(i).reset();
  }
}

bool ManagedCommBufferBase::get_send_buffers_allocated() const 
{ 
  if (m_sendBuffersAllocated)
    return m_sendBuffersAllocated;
  else
  {
    bool anyDataToSend = false;
    for (auto& buf : m_sendBufs) {
      anyDataToSend = anyDataToSend || buf.size() > 0;
    }

    return !anyDataToSend;
  }
}

}  // namespace
