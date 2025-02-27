#include "stk_util/parallel/CommBuffer.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {

void CommBuffer::pack_overflow() const
{
  STK_ThrowErrorMsg("stk::CommBuffer::pack<T>(...){ overflow by " << remaining() << " bytes. }");
}

void CommBuffer::unpack_overflow() const
{
  STK_ThrowErrorMsg("stk::CommBuffer::unpack<T>(...){ overflow by " << remaining() << " bytes. }");
}

void CommBuffer::set_buffer_ptrs(unsigned char* begin, unsigned char* ptr, unsigned char* end)
{
  m_beg = begin;
  m_ptr = ptr;
  m_end = end;
  m_offset = static_cast<size_t>(m_ptr-m_beg);
}

}
