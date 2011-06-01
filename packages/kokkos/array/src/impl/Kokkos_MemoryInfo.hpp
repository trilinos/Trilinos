/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_MEMORYINFO_HPP
#define KOKKOS_MEMORYINFO_HPP

#include <cstddef>
#include <typeinfo>
#include <set>

namespace Kokkos {
namespace Impl {

class MemoryInfo {
public:
  void                 * m_ptr ;
  const std::type_info * m_type ;
  std::string            m_label ;
  size_t                 m_size ;
  size_t                 m_count ;

  MemoryInfo() : m_ptr(0), m_type(0), m_label(), m_size(0), m_count(0) {}

  MemoryInfo( const MemoryInfo & rhs )
    : m_ptr(   rhs.m_ptr )
    , m_type(  rhs.m_type )
    , m_label( rhs.m_label )
    , m_size(  rhs.m_size )
    , m_count( rhs.m_count )
    {}

  MemoryInfo & operator = ( const MemoryInfo & rhs )
    {
      m_ptr   = rhs.m_ptr ;
      m_type  = rhs.m_type ;
      m_label = rhs.m_label ;
      m_size  = rhs.m_size ;
      m_count = rhs.m_count ;
      return *this ;
    }

  ~MemoryInfo()
    {
      m_ptr   = 0 ;
      m_type  = 0 ;
      m_label.clear();
      m_size = 0 ;
      m_count = 0 ;
    }

  bool operator == ( const MemoryInfo & rhs ) const
  { return m_ptr == rhs.m_ptr ; }

  bool operator < ( const MemoryInfo & rhs ) const
  { return m_ptr < rhs.m_ptr ; }
};

class MemoryInfoSet {
private:
  std::set< MemoryInfo > m_set ;
public:

  bool insert( const MemoryInfo & m )
  {
    std::pair< std::set< MemoryInfo >::iterator , bool > result =
      m_set.insert( m );
    return result.second ;
  }

  bool erase( void * ptr )
  {
    MemoryInfo tmp ;
    tmp.m_ptr = ptr ;
    return 1 == m_set.erase( tmp );
  }

  bool empty() const { return m_set.empty(); }

  void print( std::ostream & ) const ;
};

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_MEMORYINFO_HPP */

