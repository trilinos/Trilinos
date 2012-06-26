/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <impl/KokkosArray_MemoryTracking.hpp>

namespace KokkosArray {
namespace Impl {
namespace {

bool contains( const MemoryTracking::Info & block ,
               const void * const ptr )
{
  return block.begin <= ptr && ptr < block.end ;
}

} // namespace

struct LessMemoryTrackingInfo {

  LessMemoryTrackingInfo() {}

  bool operator()( const MemoryTracking::Info & lhs ,
                   const void * const rhs_ptr ) const
  { return lhs.end < rhs_ptr ; }
};

void MemoryTracking::track(
  const void           * ptr ,
  const std::type_info * type ,
  const size_t           size ,
  const size_t           length ,
  const std::string      label )
{
  const LessMemoryTrackingInfo compare ;

  std::vector<Info>::iterator i =
    std::lower_bound( m_tracking.begin() , m_tracking.end() , ptr , compare );

  if ( i != m_tracking.end() && contains( *i , ptr ) ) {
    std::ostringstream msg ;
    msg << "MemoryTracking::track( "
        << "ptr(" << ptr << ") ,"
        << "typeid(" << type->name() << ") ,"
        << "size(" << size << ") ,"
        << "length(" << length << ") ,"
        << "label(" << label << ") )"
        << " ERROR, already exists as { "
        << "begin(" << i->begin << ") ,"
        << "end(" << i->end << ") ,"
        << "typeid(" << i->type->name() << ") ,"
        << "size(" << i->size << ") ,"
        << "length(" << i->length << ") ,"
        << "label(" << i->label << ")"
        << "count(" << i->count << ") }" ;
    throw std::runtime_error( msg.str() );
  }

  Info info ;

  info.label  = label ;
  info.begin  = ptr ;
  info.end    = ((const char *)ptr) + size * length ;
  info.type   = type ;
  info.size   = size ;
  info.length = length ;
  info.count  = 1 ;

  m_tracking.insert( i , info );
}

MemoryTracking::Info
MemoryTracking::increment( const void * ptr )
{
  const LessMemoryTrackingInfo compare ;

  std::vector<Info>::iterator i =
    std::lower_bound( m_tracking.begin() , m_tracking.end() , ptr , compare );

  if ( i == m_tracking.end() || ! contains( *i , ptr ) ) {
    std::ostringstream msg ;
    msg << "MemoryTracking::increment( "
        << "ptr(" << ptr << ") ) ERROR, not being tracked" ;
    throw std::runtime_error(msg.str());
  }

  ++( i->count );

  return *i ;
}

MemoryTracking::Info
MemoryTracking::decrement( const void * ptr )
{
  const LessMemoryTrackingInfo compare ;

  std::vector<Info>::iterator i =
    std::lower_bound( m_tracking.begin() , m_tracking.end() , ptr , compare );

  if ( i == m_tracking.end() || ! contains( *i , ptr ) ) {
    std::ostringstream msg ;
    msg << "MemoryTracking::decrement( "
        << "ptr(" << ptr << ") ) ERROR, not being tracked" ;
    throw std::runtime_error(msg.str());
  }

  const size_t count = --( i->count );

  if ( 0 == count ) {
    Info entry = *i ;
    m_tracking.erase( i );
    return entry ;
  }
  else {
    return *i ;
  }
}

MemoryTracking::Info
MemoryTracking::query( const void * ptr ) const
{
  const LessMemoryTrackingInfo compare ;

  std::vector<Info>::const_iterator i =
    std::lower_bound( m_tracking.begin() , m_tracking.end() , ptr , compare );

  return ( i != m_tracking.end() && contains( *i , ptr ) ) ? *i : Info();
}

void MemoryTracking::print( std::ostream & s , const std::string & lead ) const
{
  for ( std::vector<Info>::const_iterator
        i = m_tracking.begin() ; i != m_tracking.end() ; ++i ) {
    s << lead
      << "{ begin(" << i->begin << "), "
      << "end(" << i->end << "), "
      << "typeid(" << i->type->name() << "), "
      << "size(" << i->size << "), "
      << "length(" << i->length << "), "
      << "label(" << i->label << "), "
      << "count(" << i->count << ") }"
      << std::endl ;
  }
}

} /* namespace Impl */
} /* namespace KokkosArray */


