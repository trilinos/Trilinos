#include "DataExchangeKnownPatternNonBlockingBuffer.hpp"

#ifdef STK_HAS_MPI

namespace stk {

DataExchangeKnownPatternNonBlockingCommBuffer::DataExchangeKnownPatternNonBlockingCommBuffer(
  MPI_Comm comm, int tag_hint) :
  ManagedCommBufferBase(comm),
  m_exchanger(comm, tag_hint)
{}

void DataExchangeKnownPatternNonBlockingCommBuffer::set_recv_buffer_size(int rank, size_t bytes)
{
  m_recvBufStorage[rank].resize(bytes);
}

void DataExchangeKnownPatternNonBlockingCommBuffer::allocate_recv_buffers()
{
  set_recv_buffer_storage();
  m_recvBufsAllocated = true;
}

void DataExchangeKnownPatternNonBlockingCommBuffer::start_nonblocking()
{
  STK_ThrowRequireMsg(get_send_buffers_allocated(), "send buffers must have been allocated before starting the send");
  STK_ThrowRequireMsg(m_recvBufsAllocated, "must allocate recv buffers before starting recvs");

  m_exchanger.start_nonblocking(m_sendBufStorage, m_recvBufStorage);
  set_sends_in_progress(true);
  set_recvs_in_progress(true);
}

void DataExchangeKnownPatternNonBlockingCommBuffer::complete_sends()
{
  m_exchanger.complete_sends();
  set_sends_in_progress(false);
}


}  // namespace

#endif