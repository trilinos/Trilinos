/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <Host/KokkosArray_Host_Thread.hpp>

#include <limits>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count );

void host_decrement_not_thread_safe( const void * ptr );

//----------------------------------------------------------------------------

class HostThreadSentinel {
public:
  HostThreadSentinel();
  ~HostThreadSentinel();

  static HostThreadSentinel & singleton();

  unsigned get_thread_count();
  void set_thread( const unsigned , HostThread *);
  HostThread * clear_thread( const unsigned );
  void set_relationships();

  static void resize_reduce( HostThread & , unsigned );

  void resize_reduce( unsigned );

private:

  void clear_relationships();

  unsigned relationships ;
};

//----------------------------------------------------------------------------

HostThread::HostThread()
{
  // Initialize as the only thread:
  m_state        = ThreadInactive ;
  m_fan_count    = 0 ;
  m_thread_rank  = 0 ;
  m_thread_count = 1 ;
  m_gang_rank    = 0 ;
  m_gang_count   = 1 ;
  m_worker_rank  = 0 ;
  m_worker_count = 1 ;
  m_reduce       = 0 ;
  m_gang_tag     = 0 ;
  m_worker_tag   = 0 ;

  for ( unsigned i = 0 ; i < max_fan_count ; ++i ) { m_fan[i] = 0 ; }
}

HostThread::~HostThread()
{
  if ( m_reduce ) {
    std::cerr << "KokkosArray::Impl::HostThread WARNING : destroyed with allocated reduction memory"
              << std::endl ;
  }

  m_state        = ThreadInactive ;
  m_fan_count    = 0 ;
  m_thread_rank  = std::numeric_limits<unsigned>::max();
  m_thread_count = 0 ;
  m_gang_rank    = std::numeric_limits<unsigned>::max();
  m_gang_count   = 0 ;
  m_worker_rank  = std::numeric_limits<unsigned>::max();
  m_worker_count = 0 ;
  m_reduce       = 0 ;
  m_gang_tag     = 0 ;
  m_worker_tag   = 0 ;

  for ( unsigned i = 0 ; i < max_fan_count ; ++i ) { m_fan[i] = 0 ; }
}

void HostThread::resize_reduce( unsigned size )
{
  if ( m_reduce ) {
    host_decrement_not_thread_safe( m_reduce );
    m_reduce = 0 ;
  }

  if ( size ) {

    m_reduce = host_allocate_not_thread_safe( "reduce_scratch_space" , typeid(unsigned char) , 1 , size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( m_reduce );
    unsigned * const end = ptr + size / sizeof(unsigned);

    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

unsigned HostThread::get_thread_count()
{ return HostThreadSentinel::singleton().get_thread_count(); }

void HostThread::set_thread_relationships()
{ HostThreadSentinel::singleton().set_relationships(); }

void HostThread::set_thread( const unsigned entry , HostThread * t )
{ HostThreadSentinel::singleton().set_thread( entry , t ); }

HostThread * HostThread::clear_thread( const unsigned entry )
{ return HostThreadSentinel::singleton().clear_thread( entry ); }

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

HostThread * HostThread::m_thread[ HostThread::max_thread_count ];

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

HostThreadSentinel &
HostThreadSentinel::singleton()
{
  static HostThreadSentinel self ;
  return self ;
}

HostThreadSentinel::HostThreadSentinel()
: relationships(0)
{
  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    HostThread::m_thread[i] = 0 ;
  }
}

HostThreadSentinel::~HostThreadSentinel()
{
  unsigned nonzero_count = 0 ;
  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    if ( 0 != HostThread::m_thread[i] ) ++nonzero_count ;
  }
  if ( nonzero_count ) {
    std::cerr << "KokkosArray::Impl::HostThread WARNING Terminating with "
              << nonzero_count
              << " non-null threads."
              << std::endl ;
  }
}

//----------------------------------------------------------------------------

void HostThreadSentinel::set_thread( const unsigned entry , HostThread * t )
{
  clear_relationships();

  const bool ok_rank = entry < HostThread::max_thread_count ;
  const bool ok_zero = ok_rank && ( 0 == HostThread::m_thread[ entry ] );

  if ( ok_rank && ok_zero ) {
    HostThread::m_thread[ entry ] = t ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::set_thread( "
        << entry << " , ... ) ERROR: " ;
    if ( ! ok_rank ) { msg << " OUT OF BOUNDS" ; }
    else if ( ! ok_zero ) { msg << " ALREADY SET" ; }
    throw std::runtime_error( msg.str() );
  }
}

HostThread * HostThreadSentinel::clear_thread( const unsigned entry )
{
  clear_relationships();

  HostThread * th = 0 ;

  if ( entry < HostThread::max_thread_count ) {
    th = HostThread::m_thread[ entry ] ;
    HostThread::m_thread[ entry ] = 0 ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::clear_thread( "
        << entry << " , ... ) ERROR:  OUT OF BOUNDS" ;
    throw std::runtime_error( msg.str() );
  }

  return th ;
}

//----------------------------------------------------------------------------

unsigned HostThreadSentinel::get_thread_count()
{
  unsigned count = 0 ;

  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    if ( HostThread::m_thread[i] ) ++count ;
  }
  return count ;
}

void HostThreadSentinel::clear_relationships()
{
  if ( relationships ) {
    for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
      if ( HostThread::m_thread[i] ) {

        HostThread & thread = * HostThread::m_thread[i] ;

        // Reset to view self as the only thread:

        thread.m_fan_count    = 0 ;
        thread.m_thread_rank  = 0 ;
        thread.m_thread_count = 1 ;
        thread.m_gang_rank    = 0 ;
        thread.m_gang_count   = 1 ;
        thread.m_worker_rank  = 0 ;
        thread.m_worker_count = 1 ;

        for ( unsigned j = 0 ; j < HostThread::max_fan_count ; ++j ) {
          thread.m_fan[j] = 0 ;
        }
      }
    }
    relationships = 0 ;
  }
}

//----------------------------------------------------------------------------

struct HostThreadCompare {
  bool operator()( const HostThread * const lhs ,
                   const HostThread * const rhs ) const
  {
    return lhs->m_gang_tag   != rhs->m_gang_tag ?
           lhs->m_gang_tag   <  rhs->m_gang_tag : (
           lhs->m_worker_tag != rhs->m_worker_tag ?
           lhs->m_worker_tag <  rhs->m_worker_tag :
           lhs               < rhs );
  }
};

// Set the fan based upon gang and worker rank:

void HostThreadSentinel::set_relationships()
{
  if ( relationships ) return ;

  HostThread * threads[ HostThread::max_thread_count ];

  unsigned thread_rank = 0 ;
  unsigned gang_rank   = 0 ;

  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    if ( HostThread::m_thread[i] != 0 )
      threads[ thread_rank++ ] = HostThread::m_thread[i] ;
  }

  if ( 0 == thread_rank ) return ;

  const unsigned thread_count = thread_rank ;

  HostThread ** const threads_end = threads + thread_count ;

  std::sort( threads , threads_end , HostThreadCompare() );

  { // Assign thread, gang, and worker counts and ranks:

    HostThread ** gang_begin = threads ;
    HostThread ** gang_end   = threads ;

    thread_rank = 0 ;

    while ( gang_end < threads_end ) {

      // Span of gang:
      while ( gang_end < threads_end &&
              (*gang_begin)->m_gang_tag == (*gang_end)->m_gang_tag ) ++gang_end ;

      // Workers within this gang:
      const unsigned worker_count = gang_end - gang_begin ;

      // Assign thread_rank, gang_rank, worker_rank
      // Assign thread_count and worker_count
      for ( unsigned worker_rank = 0 ;
                     worker_rank < worker_count ;
                     ++worker_rank , ++thread_rank ) {

        HostThread & thread = * gang_begin[ worker_rank ];

        // Reassign counts and ranks to be sequential:

        thread.m_thread_count = thread_count ;
        thread.m_thread_rank  = thread_rank ;
        thread.m_gang_count   = 0 ;
        thread.m_gang_rank    = gang_rank ;
        thread.m_worker_count = worker_count ;
        thread.m_worker_rank  = worker_rank ;
      }

      if ( worker_count ) ++gang_rank ;

      gang_begin = gang_end ;
    }

    // Assign gang_count
    for ( unsigned i = 0 ; i < thread_count ; ++i ) {
      threads[i]->m_gang_count = gang_rank ;
    }
  }

  // Assign fan in/out for each thread
  for ( unsigned i = 0 ; i < thread_count ; ++i ) {

    HostThread & thread = * threads[i] ; // i == thread.m_thread_rank

    unsigned fan_count = 0 ;

    // Intra-gang reduction:
    for ( unsigned n = 1 ; thread.m_worker_rank + n < thread.m_worker_count ; n <<= 1 ) {

      if ( n & thread.m_worker_rank ) break ;

      thread.m_fan[ fan_count++ ] = threads[ thread.m_thread_rank + n ];
    }

    // Inter-gang reduction:
    if ( thread.m_worker_rank == 0 ) {

      for ( unsigned n = 1 ; thread.m_gang_rank + n < thread.m_gang_count ; n <<= 1 ) {

        if ( n & thread.m_gang_rank ) break ;

        // Find thread with:
        //   gank_rank   == thread.m_gang_rank + n
        //   worker_rank == 0
        HostThread ** th = threads + i ;
        while ( (*th)->m_gang_rank < thread.m_gang_rank + n ) {
          th += (*th)->m_worker_count ;
        }

        thread.m_fan[ fan_count++ ] = *th ;
      }
    }

    thread.m_fan_count = fan_count ;
  }

  relationships = 1 ;

#if 0

  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    if ( 0 != HostThread::m_thread[i] ) {
      HostThread & thread = * HostThread::m_thread[i] ;

      std::cout << "HostThread[" << i << "] :"
                << " rank[ " << thread.m_thread_rank
                << " / " << thread.m_thread_count << " ]"
                << " gang[ " << thread.m_gang_tag
                << " : " << thread.m_gang_rank
                << " / " << thread.m_gang_count << " ]"
                << " worker[ " << thread.m_worker_tag
                << " : " << thread.m_worker_rank
                << " / " << thread.m_worker_count << " ]"
                ;

      if ( thread.m_fan_count ) {
        std::cout << " fan_ranks[" ;
        for ( unsigned j = 0 ; j < thread.m_fan_count ; ++j ) {
          std::cout << " " << thread.m_fan[j]->m_thread_rank ;
        }
        std::cout << " ]" ;
      }

      std::cout << std::endl ;
    }
  }

#endif

}

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray


