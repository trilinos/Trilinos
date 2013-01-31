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

#include <stddef.h>
#include <limits>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <impl/KokkosArray_Error.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>

namespace KokkosArray {
namespace Impl {
namespace {

//----------------------------------------------------------------------------

const void * upper_bound()
{
  return (const void *) std::numeric_limits<ptrdiff_t>::max();
}

// Fast search for result[-1] <= val < result[0].
// Requires result[max] == upper_bound.
// Start with a binary search until the search range is
// less than LINEAR_LIMIT, then switch to linear search.

int upper_bound( const void * const * const begin , unsigned length ,
                 const void * const val )
{
  typedef const void * const * ptr_type ;

  enum { LINEAR_LIMIT = 32 };

  // precondition: begin[length-1] == upper_bound()

  ptr_type first = begin ;

  while ( LINEAR_LIMIT < length ) {
    unsigned half   = length >> 1 ;
    ptr_type middle = first + half ;

    if ( val < *middle ) {
      length = half ;
    }
    else {
      first   = ++middle ;
      length -= ++half ;
    }
  }

  for ( ; ! ( val < *first ) ; ++first );

  return first - begin ;
}

} // namespace

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void MemoryTracking::Info::print( std::ostream & s ) const
{
  const void * const end =
    ((const unsigned char *) begin ) + size * length ;
  s << "{ "
    << "range[ " << begin << " : " << end << " )"
    << "typeid(" << type->name() << ") "
    << "size(" << size << ") "
    << "length(" << length << ") "
    << "count(" << count << ")  "
    << "label(" << label << ") }" ;
}

//----------------------------------------------------------------------------

bool MemoryTracking::empty() const
{
  return 1 == m_tracking.size();
}

MemoryTracking::MemoryTracking( const std::string & space )
  : m_space( space ), m_tracking(), m_tracking_end()
{
  const void * const max = upper_bound();

  // Sentinal value of end [ max , max )

  m_tracking.reserve(64);
  m_tracking_end.reserve(64);

  Info info ;

  info.begin  = max ;
  info.type   = NULL ;
  info.size   = 0 ;
  info.length = 0 ;
  info.count  = 0 ;

  m_tracking.push_back( info );
  m_tracking_end.push_back( max );
}

MemoryTracking::~MemoryTracking()
{
  try {
    if ( 1 < m_tracking.size() ) {
      std::cerr << m_space << " destroyed with memory leaks:" << std::endl ;
      print( std::cerr , std::string("  ") );
    }
    else if ( 1 != m_tracking.size() ||
              1 != m_tracking_end.size() ||
              m_tracking.back().begin != upper_bound() ||
              m_tracking_end.back() != upper_bound() ) {
      std::cerr << m_space << " corrupted data structure" << std::endl ;
    }
  } catch( ... ) {}
}

void MemoryTracking::track(
  const void           * ptr ,
  const std::type_info * type ,
  const size_t           size ,
  const size_t           length ,
  const std::string      label )
{
  const int offset = upper_bound( & m_tracking_end[0] , m_tracking_end.size() , ptr );
  const std::vector<Info>::iterator i = m_tracking.begin() + offset ;

  // Guaranteed:
  //   i == 0 || m_tracking_end[i-1] <= ptr < m_tracking_end[i]
  //   ptr < m_tracking_end[i]

  if ( ptr < i->begin ) {
    const std::vector<const void *>::iterator j = m_tracking_end.begin() + offset ;

    Info info ;

    info.label  = label ;
    info.begin  = ptr ;
    info.type   = type ;
    info.size   = size ;
    info.length = length ;
    info.count  = 1 ;

    m_tracking.insert( i , info );
    m_tracking_end.insert( j , ((const unsigned char *)ptr) + size * length );
  }
  else {
    std::ostringstream msg ;
    msg << "MemoryTracking::track( "
        << "ptr(" << ptr << ") ,"
        << "typeid(" << type->name() << ") ,"
        << "size(" << size << ") ,"
        << "length(" << length << ") ,"
        << "label(" << label << ") )"
        << " ERROR, already exists as " ;
    i->print( msg );
    throw_runtime_exception( msg.str() );
  }
}

void MemoryTracking::increment( const void * ptr )
{
  const std::vector<Info>::iterator i = m_tracking.begin() +
    upper_bound( & m_tracking_end[0] , m_tracking_end.size() , ptr );

  if ( i->begin <= ptr ) {
    ++( i->count );
  }
  else {
    std::ostringstream msg ;
    msg << "MemoryTracking(" << (void *) this
        << ")::increment( "
        << "ptr(" << ptr << ") ) ERROR, not being tracked" ;
    throw_runtime_exception( msg.str() );
  }
}

void * MemoryTracking::decrement( const void * ptr )
{
  const int offset =
    upper_bound( & m_tracking_end[0] , m_tracking_end.size() , ptr );

  const std::vector<Info>::iterator i = m_tracking.begin() + offset ;

  void * ptr_alloc = 0 ;

  if ( i->begin <= ptr ) {
    --( i->count );

    if ( 0 == i->count ) {
      const std::vector<const void*>::iterator j = m_tracking_end.begin() + offset ;

      ptr_alloc = const_cast<void*>( i->begin );

      m_tracking.erase( i );
      m_tracking_end.erase( j );
    }
  }
  else {
    std::ostringstream msg ;
    msg << "MemoryTracking(" << (void *) this
        << ")::decrement( "
        << "ptr(" << ptr << ") ) ERROR, not being tracked" ;
    throw_runtime_exception( msg.str() );
  }


  return ptr_alloc ;
}

MemoryTracking::Info
MemoryTracking::query( const void * ptr ) const
{
  const std::vector<Info>::const_iterator i = m_tracking.begin() +
    upper_bound( & m_tracking_end[0] , m_tracking_end.size() , ptr );

  return ( i->begin <= ptr ) ? *i : Info();
}

void MemoryTracking::print( std::ostream & s , const std::string & lead ) const
{
  const std::vector<Info>::const_iterator iend = m_tracking.end() - 1 ;

  for ( std::vector<Info>::const_iterator i = m_tracking.begin() ; i != iend ; ++i ) {
    s << lead ;
    i->print( s );
    s << std::endl ;
  }
}

} /* namespace Impl */
} /* namespace KokkosArray */


