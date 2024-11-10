#ifndef stk_util_parallel_CommBuffer_hpp
#define stk_util_parallel_CommBuffer_hpp

#include <stddef.h>
#include <string>
#include <map>
#include <vector>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, ThrowRequire


namespace stk {

/** Perform collective all-to-all communication with individually
 *  varying message sizes.  The collective operation uses an
 *  all-to-all if the maximum number of sends or receives from
 *  any one processor is greater than the given bounds.
 *
 *  This is a work-horse communication for mesh data structures
 *  with parallel domain decomposition.
 */
class CommSparse;
class CommNeighbors;
class CommBroadcast;

template<typename T>
struct is_pair { static constexpr bool value = false; };

template<template<typename...> class C, typename U, typename V>
struct is_pair<C<U,V>> { static constexpr bool value = std::is_same<C<U,V>, std::pair<U,V>>::value; };

template <typename T>
using IsPair = std::enable_if_t<is_pair<T>::value>;

template <typename T>
using NotPair = std::enable_if_t<!is_pair<T>::value>;

class CommBuffer {
public:

  /** Pack a value to be sent:  buf.pack<type>( value ) */
  template<typename T,
           class = NotPair<T>>
  CommBuffer &pack( const T & value );

  CommBuffer &pack( const std::string & value );

  template<typename P,
           class = IsPair<P>, class = void>
  CommBuffer &pack(const P & value);

  template<typename K, typename V>
  CommBuffer &pack( const std::map<K,V> & value );

  template<typename K>
  CommBuffer &pack( const std::vector<K> & value );

private:
  /** Do not try to pack a pointer for global communication */
  template<typename T> CommBuffer &pack( const T* value ) {
    STK_ThrowAssertMsg(false,"CommBuffer::pack(const T* value) not allowed. Don't pack a pointer for communication!");
    return *this;
  }

public:

  /** Pack an array of values to be sent:  buf.pack<type>( ptr , num ) */
  template<typename T> CommBuffer &pack( const T * value , size_t number );

  /** Unpack a received value:  buf.unpack<type>( value ) */
  template<typename T,
           class = NotPair<T>>
  CommBuffer &unpack( T & value );

  CommBuffer &unpack( std::string& value );

  template<typename P,
           class = IsPair<P>, class = void>
  CommBuffer &unpack( P & value);

  template<typename K, typename V>
  CommBuffer &unpack( std::map<K,V> & value );

  template<typename K>
  CommBuffer &unpack( std::vector<K> & value );

  /** Unpack an array of received values:  buf.unpack<type>( ptr , num ) */
  template<typename T> CommBuffer &unpack( T * value , size_t number );

  /** Peek at a received value (don't advance buffer): buf.peek<type>(value) */
  template<typename T> CommBuffer &peek( T & value );

  /** Peek at an array of received values: buf.peek<type>( ptr , num ) */
  template<typename T> CommBuffer &peek( T * value , size_t number );

  CommBuffer &peek( std::string& value );

  template<typename K, typename V>
  CommBuffer &peek( std::map<K,V> & value );

  /** Skip buffer ahead by a number of values. */
  template<typename T,
           class = NotPair<T>>
  CommBuffer &skip( size_t number );

  /** Skip buffer ahead by a number of values. */
  template<typename T,
           class = IsPair<T>, class = void>
  CommBuffer &skip( size_t number );

  /** Reset the buffer to the beginning so that size() == 0 */
  void reset();

  /** Size, in bytes, of the buffer.
   *  If the buffer is not yet allocated this is zero.
   */
  size_t capacity() const ;

  // TODO - terribly misinforming when used on recv buffer, returns 0!
  /** Size, in bytes, of the buffer that has been processed.
   *  If the buffer is not yet allocated then this is the
   *  number of bytes that has been attempted to pack.
   */
  size_t size() const ;
  void set_size(size_t newsize_bytes);

  /** Size, in bytes, of the buffer remaining to be processed.
   *  Equal to 'capacity() - size()'.  A negative result
   *  indicates either the buffer is not allocated or an
   *  overflow has occurred.  An overflow will have thrown
   *  an exception.
   */
  ptrdiff_t remaining() const ;

  /** Pointer to base of buffer. */
  void * buffer() const ;

  CommBuffer() : m_beg(nullptr), m_ptr(nullptr), m_end(nullptr), m_offset(0) { }

  void set_buffer_ptrs(unsigned char* begin, unsigned char* ptr, unsigned char* end);

private:
  friend class CommSparse ;
  friend class CommNeighbors ;
  friend class CommBroadcast ;

  void pack_overflow() const ;
  void unpack_overflow() const ;

  typedef unsigned char * ucharp ;

  ucharp m_beg ;
  ucharp m_ptr ;
  ucharp m_end ;
  unsigned m_offset;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template implementations for the CommBuffer

template <unsigned long N>
struct CommBufferAlign {
  static size_t align( size_t i ) { i %= N ; return i ? ( N - i ) : 0 ; }
};

template<>
struct CommBufferAlign<1> {
  static size_t align( size_t ) { return 0 ; }
};

template<typename T, class>
inline
CommBuffer &CommBuffer::pack( const T & value )
{
  if constexpr (std::is_same_v<T, std::string>) {
    return pack(value);
  }
  static constexpr auto Size = sizeof(T);
  if ( m_beg ) {
    size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
//std::cout<<"m_beg: "<<(void*)m_beg<<", m_ptr: "<<(void*)m_ptr<<"m_ptr-m_beg: "<<std::distance(m_beg,m_ptr)<<", nalign: "<<nalign<<", Size: "<<Size<<", m_end: "<<(void*)m_end<<std::endl;
    if ( m_end < m_ptr + nalign + Size ) { pack_overflow(); }
    while ( nalign ) { --nalign ; *m_ptr = 0 ; ++m_ptr ; }
    T *tmp = reinterpret_cast<T *>(m_ptr);
    *tmp = value ;
    m_ptr = reinterpret_cast<ucharp>( ++tmp );
  }
  else {
    size_t nalign = CommBufferAlign<Size>::align( m_offset );
    m_offset += nalign + Size ;
  }
  return *this;
}

inline
CommBuffer &CommBuffer::pack( const std::string & value )
{
  size_t length = value.length();
  pack(length);
  pack(value.c_str(), length);
  return *this;
}

template<typename P, class, class>
inline
CommBuffer &CommBuffer::pack(const P & value)
{
  pack(value.first);
  pack(value.second);
  return *this;
}

template<typename K, typename V>
inline
CommBuffer &CommBuffer::pack( const std::map<K,V> & value )
{
  size_t ns = value.size();
  pack(ns);

  for (auto && s : value)
  {
    pack(s.first);
    pack(s.second);
  }

  return *this;
}

template<typename K>
inline
CommBuffer &CommBuffer::pack( const std::vector<K> & value )
{
  pack<unsigned>(value.size());
  for (size_t i=0; i<value.size(); ++i) {
    pack(value[i]);
  }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::pack( const T * value , size_t number )
{
  static constexpr auto Size = sizeof(T);
  if ( m_beg ) {
    size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
    if ( m_end < m_ptr + nalign + number * Size ) { pack_overflow(); }
    while ( nalign ) { --nalign ; *m_ptr = 0 ; ++m_ptr ; }
    T * tmp = reinterpret_cast<T*>(m_ptr);
    while ( number ) { --number ; *tmp = *value ; ++tmp ; ++value ; }
    m_ptr = reinterpret_cast<ucharp>( tmp );
  }
  else {
    size_t nalign = CommBufferAlign<Size>::align( m_offset );
    m_offset += nalign + number * Size ;
  }
  return *this;
}

template<typename T, class>
inline
CommBuffer &CommBuffer::skip( size_t number )
{
  static constexpr auto Size = sizeof(T);
  if ( m_beg ) {
    m_ptr += CommBufferAlign<Size>::align( m_ptr - m_beg ) + Size * number ;
  }
  else {
    m_offset += CommBufferAlign<Size>::align( m_offset ) + Size * number ;
  }
  if ( m_beg && m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

template<typename T, class, class>
inline
CommBuffer &CommBuffer::skip( size_t number )
{
  skip<typename T::first_type>(number);
  skip<typename T::second_type>(number);
  return *this;
}

template<typename T, class>
inline
CommBuffer &CommBuffer::unpack( T & value )
{
  if constexpr (std::is_same_v<T,std::string>) {
    return unpack(value);
  }
  static constexpr auto Size = sizeof(T);
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  value = *tmp ;
  m_ptr = reinterpret_cast<ucharp>( ++tmp );
  if ( m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

inline
CommBuffer &CommBuffer::unpack( std::string & value )
{
  size_t length;
  unpack(length);
  std::vector<char> chars(length);
  unpack(chars.data(), length);
  value.assign(chars.data(), length);
  return *this;
}

template<typename P,
         class, class>
inline
CommBuffer &CommBuffer::unpack( P & value)
{
  unpack(value.first);
  unpack(value.second);
  return *this;
}

template<typename K, typename V>
inline
CommBuffer &CommBuffer::unpack( std::map<K,V> & value )
{
  value.clear();

  size_t ns;
  unpack(ns);

  for (size_t i = 0; i < ns; ++i) {
    K key;
    unpack(key);

    V val;
    unpack(val);

    value[key] = val;
  }
  return *this;
}

template<typename K>
inline
CommBuffer &CommBuffer::unpack( std::vector<K> & value )
{
  unsigned num_items = 0;
  unpack<unsigned>(num_items);
  value.resize(num_items);
  for (unsigned i=0;i<num_items;++i) {
    K val;
    unpack(val);
    value[i] = val;
  }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::unpack( T * value , size_t number )
{
  static constexpr auto Size = sizeof(T);
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  while ( number ) { --number ; *value = *tmp ; ++tmp ; ++value ; }
  m_ptr = reinterpret_cast<ucharp>( tmp );
  if ( m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}
template<typename item>
inline
item unpack(stk::CommBuffer& buf)
{
    item object;
    buf.unpack<item>(object);
    return object;
}

template<typename T>
inline
CommBuffer &CommBuffer::peek( T & value )
{
  ucharp oldPtr = m_ptr;
  unpack<T>(value);
  m_ptr = oldPtr;
  return *this;
}

inline
CommBuffer &CommBuffer::peek( std::string& value )
{
  size_t length;
  peek(length);

  size_t offset = sizeof(size_t);
  std::vector<char> chars(offset+length);
  peek(chars.data(), chars.size());

  value.assign(chars.data() + offset, length);

  return *this;
}

template<typename K, typename V>
inline
CommBuffer &CommBuffer::peek( std::map<K,V> & value )
{
  throw std::runtime_error("Peek not implemented for std::map");
}

template<typename T>
inline
CommBuffer &CommBuffer::peek( T * value , size_t number )
{
  static constexpr auto Size = sizeof(T);
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  while ( number ) { --number ; *value = *tmp ; ++tmp ; ++value ; }
  if ( m_end < reinterpret_cast<ucharp>(tmp) ) { unpack_overflow(); }
  return *this;
}

inline
void CommBuffer::reset()
{ m_ptr = m_beg ; }

inline
size_t CommBuffer::capacity() const
{ return m_end - m_beg ; }

inline
size_t CommBuffer::size() const
{ return m_beg ? static_cast<size_t>(m_ptr - m_beg) : static_cast<size_t>(m_offset) ; }

inline
void CommBuffer::set_size(size_t newsize_bytes)
{ m_beg = nullptr;  m_ptr = nullptr; m_offset = newsize_bytes ; m_end = nullptr; }

inline
ptrdiff_t CommBuffer::remaining() const
{ return m_end - m_ptr ; }

inline
void * CommBuffer::buffer() const
{ return static_cast<void*>( m_beg ); }

} //namespace

#endif
