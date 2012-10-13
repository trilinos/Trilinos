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

#ifndef KOKKOSARRAY_CUDAMEMORYSPACE_HPP
#define KOKKOSARRAY_CUDAMEMORYSPACE_HPP

#include <iosfwd>
#include <typeinfo>
#include <string>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

/** \brief  Memory management on the host for devices */

class CudaMemorySpace {
public:

  static int m_memory_view_tracking ;

public:

  typedef size_t size_type ;

  static void * allocate( const std::string    & label ,
                          const std::type_info & scalar_type ,
                          const size_t           scalar_size ,
                          const size_t           scalar_count );

#if ! defined( __CUDA_ARCH__ )
  static void increment( const void * );
  static void decrement( const void * );
#else
  static __device__ void increment( const void * ) {}
  static __device__ void decrement( const void * ) {}
#endif

  static void print_memory_view( std::ostream & );

  /*--------------------------------*/

  static 
  size_t preferred_alignment( size_t scalar_size , size_t scalar_count );

  /*--------------------------------*/

  static void copy_to_device_from_device( void * , const void * , size_t );
  static void copy_to_device_from_host(   void * , const void * , size_t );
  static void copy_to_host_from_device(   void * , const void * , size_t );

  /*--------------------------------*/
};

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
struct DeepCopy ;

template< typename ValueType >
struct DeepCopy< ValueType , Cuda::memory_space , Cuda::memory_space > {
  DeepCopy( ValueType * dst , const ValueType * src , size_t count )
  {
    CudaMemorySpace::copy_to_device_from_device( dst , src , sizeof(ValueType) * count );
  }
};

template< typename ValueType >
struct DeepCopy< ValueType , Cuda::memory_space , Host::memory_space > {
  DeepCopy( ValueType * dst , const ValueType * src , size_t count )
  {
    CudaMemorySpace::copy_to_device_from_host( dst , src , sizeof(ValueType) * count );
  }
};

template< typename ValueType >
struct DeepCopy< ValueType , Host::memory_space , Cuda::memory_space > {
  DeepCopy( ValueType * dst , const ValueType * src , size_t count )
  {
    CudaMemorySpace::copy_to_host_from_device( dst , src , sizeof(ValueType) * count );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_CUDAMEMORYSPACE_HPP */

