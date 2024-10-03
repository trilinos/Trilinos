#include "stk_util/parallel/CommBuffer.hpp"

namespace stk {

void CommBuffer::pack_overflow() const
{
#ifndef NDEBUG
  std::ostringstream os ;
  os << "stk::CommBuffer::pack<T>(...){ overflow by " ;
  os << remaining() ;
  os << " bytes. }" ;
  throw std::overflow_error( os.str() );
#endif
}

void CommBuffer::unpack_overflow() const
{
#ifndef NDEBUG
  std::ostringstream os ;
  os << "stk::CommBuffer::unpack<T>(...){ overflow by " ;
  os << remaining();
  os << " bytes. }" ;
  throw std::overflow_error( os.str() );
#endif
}

void CommBuffer::set_buffer_ptrs(unsigned char* begin, unsigned char* ptr, unsigned char* end)
{
  m_beg = begin;
  m_ptr = ptr;
  m_end = end;
  m_offset = static_cast<unsigned>(m_ptr-m_beg);
}

}
