/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_parallel_ParallelComm_hpp
#define stk_util_parallel_ParallelComm_hpp

#include <cstddef>
#include <iosfwd>
#include <stk_util/parallel/Parallel.hpp>

//------------------------------------------------------------------------

namespace stk_classic {

/** Perform collective all-to-all communication with individually
 *  varying message sizes.  The collective operation uses an
 *  all-to-all if the maximum number of sends or receives from
 *  any one processor is greater than the given bounds.
 *
 *  This is a work-horse communication for mesh data structures
 *  with parallel domain decomposition.
 */
class CommAll ;

/** Pack and unpack buffers for the sparse all-to-all communication.
 */
class CommBuffer ;

/** Given the send sizes determine the receive sizes.
 *  Send and receive size arrays are dimensioned to
 *  the size of the parallel machine.
 *  Return global parallel logical OR of the input local flag.
 *  This parallel reduction is aggregated into the required
 *  communication for determining the sparse sizes.
 *  Output the receive sizes and maximum number of send or
 *  receive messages for a single processor.
 *  A dense all-to-all communication is used if:
 *     num_msg_bound < num_msg_maximum
 *  otherwise a set of point-to-point messages are used.
 */
bool comm_sizes( ParallelMachine ,
                 const unsigned   num_msg_bound ,
                       unsigned & num_msg_maximum ,
                 const unsigned * const send_size ,
                       unsigned * const recv_size ,
                 bool local_flag = false );

/** If the communication is known to be dense.
 */
bool comm_dense_sizes( ParallelMachine ,
                       const unsigned * const send_size ,
                             unsigned * const recv_size ,
                       bool local_flag = false );

//------------------------------------------------------------------------

class CommBuffer {
public:

  /** Pack a value to be sent:  buf.pack<type>( value ) */
  template<typename T> CommBuffer &pack( const T & value );

  /** Pack an array of values to be sent:  buf.pack<type>( ptr , num ) */
  template<typename T> CommBuffer &pack( const T * value , size_t number );

  /** Unpack a received value:  buf.unpack<type>( value ) */
  template<typename T> CommBuffer &unpack( T & value );

  /** Unpack an array of received values:  buf.unpack<type>( ptr , num ) */
  template<typename T> CommBuffer &unpack( T * value , size_t number );

  /** Peek at a received value (don't advance buffer): buf.peek<type>(value) */
  template<typename T> CommBuffer &peek( T & value );

  /** Peek at an array of received values: buf.peek<type>( ptr , num ) */
  template<typename T> CommBuffer &peek( T * value , size_t number );

  /** Skip buffer ahead by a number of values. */
  template<typename T> CommBuffer &skip( size_t number );

  /** Reset the buffer to the beginning so that size() == 0 */
  void reset();

  /** Size, in bytes, of the buffer.
   *  If the buffer is not yet allocated this is zero.
   */
  size_t capacity() const ;

  /** Size, in bytes, of the buffer that has been processed.
   *  If the buffer is not yet allocated then this is the
   *  number of bytes that has been attempted to pack.
   */
  size_t size() const ;

  /** Size, in bytes, of the buffer remaining to be processed.
   *  Equal to 'capacity() - size()'.  A negative result
   *  indicates either the buffer is not allocated or an
   *  overflow has occured.  An overflow will have thrown
   *  an exception.
   */
  ptrdiff_t remaining() const ;

  /** Pointer to base of buffer. */
  void * buffer() const ;

  ~CommBuffer();
  CommBuffer();

private:
  friend class CommAll ;
  friend class CommGather ;
  friend class CommBroadcast ;

  static CommBuffer * allocate( const unsigned, const unsigned * const );
  static void deallocate( const unsigned , CommBuffer * );

  void pack_overflow() const ;
  void unpack_overflow() const ;

  CommBuffer( const CommBuffer & );
  CommBuffer & operator = ( const CommBuffer & );

  typedef unsigned char * ucharp ;

  ucharp m_beg ;
  ucharp m_ptr ;
  ucharp m_end ;
};

//------------------------------------------------------------------------

class CommAll {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  unsigned        parallel_size() const { return m_size ; }
  unsigned        parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for a given processor */
  CommBuffer & send_buffer( unsigned ) const ;

  /** Obtain the message buffer for a given processor */
  CommBuffer & recv_buffer( unsigned ) const ;

  //----------------------------------------
  /** Construct for undefined communication.
   *  No buffers are allocated.
   */
  CommAll();

  /** Allocate send and receive buffers based upon input sizes.
   *  If recv_size == NULL then the receive size
   *  is determined by communicating the send sizes.
   *  Symmetry is given by passing the same data for both
   *  send and receive sizes.
   *  Return global parallel OR of local flags.
   */
  bool allocate_buffers( ParallelMachine ,
                         const unsigned num_msg_bounds ,
                         const unsigned * const send_size ,
                         const unsigned * const recv_size ,
                         const bool local_flag = false );

  //----------------------------------------
  /** Construct for a to-be-sized communication.
   *  Allocate surrogate send buffers to enable
   *  no-op packing for the purpose of send sizing.
   *  Surrogate send scenario:
   *  1) Surrogate send buffers are "packed" for sizing where
   *     packing sizes are recorded but no data is copied.
   *  2) 'allocate_buffers(symmetric,flag)' is called to allocate
   *     buffers.  The symmetric flag guarantees that the send
   *     sizes matches the receive sizes.
   *  3) Send buffers are identically packed; however, this
   *     packing copies data into the send buffers.
   */
  explicit CommAll( ParallelMachine );

  /** Allocate asymmetric communication based upon
   *  sizing from the surrogate send buffer packing.
   *  If symmetric then the receive sizes are guaranteed
   *  to be identical to the send sizes.
   *  Return global parallel OR of local flags.
   */
  bool allocate_buffers( const unsigned num_msg_bounds ,
                         const bool symmetric = false ,
                         const bool local_flag = false );

  //----------------------------------------
  /** Communicate send buffers to receive buffers.  */
  void communicate();

  //----------------------------------------
  /** Swap send and receive buffers leading to reversed communication. */
  void swap_send_recv();

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

  ~CommAll();

private:

  CommAll( const CommAll & );
  CommAll & operator = ( const CommAll & );

  void rank_error( const char * , unsigned ) const ;

  bool allocate_buffers( const unsigned * const send_size ,
                         const unsigned * const recv_size ,
                         bool local_flag );

  ParallelMachine m_comm ;
  unsigned        m_size ;
  unsigned        m_rank ;
  unsigned        m_bound ;
  unsigned        m_max ;
  CommBuffer    * m_send ;
  CommBuffer    * m_recv ;
};

//------------------------------------------------------------------------

class CommBroadcast {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  unsigned        parallel_size() const { return m_size ; }
  unsigned        parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for the root_rank processor */
  CommBuffer & send_buffer();

  /** Obtain the message buffer for the local processor */
  CommBuffer & recv_buffer();

  //----------------------------------------

  CommBroadcast( ParallelMachine , unsigned root_rank );

  void communicate();

  bool allocate_buffer( const bool local_flag = false );

  ~CommBroadcast();

private:

  CommBroadcast();
  CommBroadcast( const CommBroadcast & );
  CommBroadcast & operator = ( const CommBroadcast & );

  ParallelMachine m_comm ;
  unsigned        m_size ;
  unsigned        m_rank ;
  unsigned        m_root_rank ;
  CommBuffer      m_buffer ;
};

//----------------------------------------------------------------------

class CommGather {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  unsigned        parallel_size() const { return m_size ; }
  unsigned        parallel_rank() const { return m_rank ; }

  ~CommGather();

  CommGather( ParallelMachine , unsigned root_rank , unsigned send_size );

  CommBuffer & send_buffer() { return m_send ; }

  void communicate();

  CommBuffer & recv_buffer( unsigned );

  void reset(); 

private:

  CommGather();
  CommGather( const CommBroadcast & );
  CommGather & operator = ( const CommBroadcast & );

  ParallelMachine m_comm ;
  unsigned        m_size ;
  unsigned        m_rank ;
  unsigned        m_root_rank ;
  CommBuffer      m_send ;
  CommBuffer    * m_recv ;
  int           * m_recv_count ;
  int           * m_recv_displ ;
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template implementations for the CommBuffer

namespace stk_classic {

template<unsigned N> struct CommBufferAlign ;

template<>
struct CommBufferAlign<1> {
  static size_t align( size_t ) { return 0 ; }
};

template<unsigned N>
struct CommBufferAlign {
  static size_t align( size_t i ) { i %= N ; return i ? ( N - i ) : 0 ; }
};

template<typename T>
inline
CommBuffer &CommBuffer::pack( const T & value )
{
  enum { Size = sizeof(T) };
  size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  if ( m_beg ) {
    if ( m_end < m_ptr + nalign + Size ) { pack_overflow(); }
    while ( nalign ) { --nalign ; *m_ptr = 0 ; ++m_ptr ; }
    T * tmp = reinterpret_cast<T*>(m_ptr);
    *tmp = value ;
    m_ptr = reinterpret_cast<ucharp>( ++tmp );
  }
  else {
    m_ptr += nalign + Size ;
  }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::pack( const T * value , size_t number )
{
  enum { Size = sizeof(T) };
  size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  if ( m_beg ) {
    if ( m_end < m_ptr + nalign + number * Size ) { pack_overflow(); }
    while ( nalign ) { --nalign ; *m_ptr = 0 ; ++m_ptr ; }
    T * tmp = reinterpret_cast<T*>(m_ptr);
    while ( number ) { --number ; *tmp = *value ; ++tmp ; ++value ; }
    m_ptr = reinterpret_cast<ucharp>( tmp );
  }
  else {
    m_ptr += nalign + number * Size ;
  }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::skip( size_t number )
{
  enum { Size = sizeof(T) };
  m_ptr += CommBufferAlign<Size>::align( m_ptr - m_beg ) + Size * number ;
  if ( m_beg && m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::unpack( T & value )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  value = *tmp ;
  m_ptr = reinterpret_cast<ucharp>( ++tmp );
  if ( m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::unpack( T * value , size_t number )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  while ( number ) { --number ; *value = *tmp ; ++tmp ; ++value ; }
  m_ptr = reinterpret_cast<ucharp>( tmp );
  if ( m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::peek( T & value )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  value = *tmp ;
  if ( m_end < reinterpret_cast<ucharp>(++tmp) ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::peek( T * value , size_t number )
{
  enum { Size = sizeof(T) };
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
{ return m_ptr - m_beg ; }

inline
ptrdiff_t CommBuffer::remaining() const
{ return m_end - m_ptr ; }

inline
void * CommBuffer::buffer() const
{ return static_cast<void*>( m_beg ); }

//----------------------------------------------------------------------
// Inline implementations for the CommAll

inline
CommBuffer & CommAll::send_buffer( unsigned p ) const
{
  if ( m_size <= p ) { rank_error("send_buffer",p); }
  return m_send[p] ;
}

inline
CommBuffer & CommAll::recv_buffer( unsigned p ) const
{
  if ( m_size <= p ) { rank_error("recv_buffer",p); }
  return m_recv[p] ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

