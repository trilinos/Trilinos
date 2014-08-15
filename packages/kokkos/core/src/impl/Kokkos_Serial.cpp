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
#include <sstream>
#include <Kokkos_Serial.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace {

struct Sentinel {

  void *   m_scratch ;
  unsigned m_reduce_end ;
  unsigned m_shared_end ;

  Sentinel() : m_scratch(0), m_reduce_end(0), m_shared_end(0) {}

  ~Sentinel() { if ( m_scratch ) { free( m_scratch ); } }

  static Sentinel & singleton();
};

Sentinel & Sentinel::singleton()
{
  static Sentinel s ; return s ;
}

inline
unsigned align( unsigned n )
{
  enum { ALIGN = 0x0100 /* 256 */ , MASK = ALIGN - 1 };
  return ( n + MASK ) & ~MASK ;
}

}

void * Serial::scratch_memory_space::resize( unsigned reduce_size , unsigned shared_size )
{
  static Sentinel & s = Sentinel::singleton();

  reduce_size = align( reduce_size );
  shared_size = align( shared_size );

  if ( ( 0 == reduce_size + shared_size ) ||
       ( s.m_reduce_end < reduce_size ) ||
       ( s.m_shared_end < s.m_reduce_end + shared_size ) ) {

    if ( s.m_scratch ) { free( s.m_scratch ); }
  }

  if ( ( s.m_reduce_end < reduce_size ) ||
       ( s.m_shared_end < s.m_reduce_end + shared_size ) ) {
  
    if ( s.m_reduce_end < reduce_size ) s.m_reduce_end = reduce_size ;
    if ( s.m_shared_end < s.m_reduce_end + shared_size ) s.m_shared_end = s.m_reduce_end + shared_size ;

    s.m_scratch = malloc( s.m_shared_end );
  }

  return s.m_scratch ;
}

Serial::scratch_memory_space::scratch_memory_space()
{
  Sentinel & s = Sentinel::singleton();
  m_shmem_iter = ((char *) s.m_scratch ) + s.m_reduce_end ;
  m_shmem_end  = ((char *) s.m_scratch ) + s.m_shared_end ;
}

void Serial::scratch_memory_space::get_shmem_error()
{
  Kokkos::Impl::throw_runtime_exception( std::string("Serial::get_shmem FAILED : out of memory") );
}

} // namespace Kokkos

