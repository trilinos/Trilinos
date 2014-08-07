/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <Kokkos_Serial.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace {

struct Sentinel {

  void *   m_reduce ;
  void *   m_shared ;
  unsigned m_reduce_size ;
  unsigned m_shared_size ;

  Sentinel() : m_reduce(0), m_shared(0), m_reduce_size(0), m_shared_size(0) {}

  ~Sentinel()
    {
      if ( m_reduce ) { free( m_reduce ); }
      if ( m_shared ) { free( m_shared ); }
    }

  static Sentinel & singleton();
};

Sentinel & Sentinel::singleton()
{
  static Sentinel s ; return s ;
}

}

void * Serial::scratch_memory_space::resize_reduce_scratch( unsigned size )
{
  static Sentinel & s = Sentinel::singleton();

  const unsigned rem = size % Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size ) || ( s.m_reduce_size < size ) ) {

    if ( s.m_reduce ) { free( s.m_reduce ); }
  
    s.m_reduce_size = size ;

    s.m_reduce = malloc( size );
  }

  return s.m_reduce ;
}

void * Serial::scratch_memory_space::resize_shared_scratch( unsigned size )
{
  static Sentinel & s = Sentinel::singleton();

  const unsigned rem = size % Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size ) || ( s.m_shared_size < size ) ) {

    if ( s.m_shared ) { free( s.m_shared ); }
  
    s.m_shared_size = size ;

    s.m_shared = malloc( size );
  }

  return s.m_shared ;
}

void * Serial::scratch_memory_space::get_shmem( const int size ) const
{
  static Sentinel & s = Sentinel::singleton();

  const int offset = m_shmem_iter >> Impl::power_of_two<sizeof(int)>::value ;

  m_shmem_iter += size ;

  if ( int(s.m_shared_size) < m_shmem_iter ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Serial::get_shmem FAILED : exceeded shared memory size" ) );
  }

  return ((int*)s.m_shared) + offset ;
}

} // namespace Kokkos

