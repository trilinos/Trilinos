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

#include <stdlib.h>
#include <iostream>
#include <sstream>

#include <KokkosArray_CudaSpace.hpp>
#include <Cuda/KokkosArray_Cuda_Internal.hpp>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>
#include <impl/KokkosArray_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace {

Impl::MemoryTracking & cuda_space_singleton()
{
  static Impl::MemoryTracking self("KokkosArray::CudaSpace");
  return self ;
}

}

/*--------------------------------------------------------------------------*/

DeepCopy<HostSpace,CudaSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

DeepCopy<CudaSpace,HostSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

DeepCopy<CudaSpace,CudaSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

/*--------------------------------------------------------------------------*/

void * CudaSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  HostSpace::assert_master_thread( "KokkosArray::CudaSpace::allocate" );

  const size_t size = scalar_size * scalar_count ;

  void * ptr = 0 ;

  if ( 0 < scalar_size * scalar_count ) {
    bool ok = true ;

    cudaDeviceSynchronize();

    if ( ok ) ok = cudaSuccess == cudaMalloc( & ptr , size );
    if ( ok ) ok = 0 != ptr ;
    if ( ok ) ok = cudaSuccess == cudaMemset( ptr , 0 , size );
    if ( ok ) ok = cudaSuccess == cudaThreadSynchronize();

    if ( ! ok ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::CudaSpace::allocate( "
          << label
          << " , " << scalar_type.name()
          << " , " << scalar_size
          << " , " << scalar_count
          << " ) FAILED memory allocation" ;
      KokkosArray::Impl::throw_runtime_exception( msg.str() );
    }

    cuda_space_singleton()
      .track( ptr, & scalar_type, scalar_size, scalar_count, label );
  }

  return ptr ;
}

#if ! defined( __CUDA_ARCH__ )

void CudaSpace::increment( const void * ptr )
{
  HostSpace::assert_master_thread( "KokkosArray::CudaSpace::increment" );

  if ( 0 != ptr ) {
    cuda_space_singleton().increment( ptr );
  }
}

void CudaSpace::decrement( const void * ptr )
{
  HostSpace::assert_master_thread( "KokkosArray::CudaSpace::decrement" );

  if ( 0 != ptr ) {

    void * ptr_alloc = cuda_space_singleton().decrement( ptr );

    if ( 0 != ptr_alloc ) {

      cudaDeviceSynchronize();

      const bool failed = cudaSuccess != cudaFree( ptr_alloc );

      if ( failed ) {
        std::string msg("KokkosArray::Impl::CudaSpace::decrement() failed cudaFree");
        KokkosArray::Impl::throw_runtime_exception( msg );
      }
    }
  }
}

#endif

void CudaSpace::print_memory_view( std::ostream & o )
{
  cuda_space_singleton().print( o , std::string("  ") );
}


size_t CudaSpace::preferred_alignment(
  size_t scalar_size , size_t scalar_count )
{
  const size_t align = 0 == MEMORY_ALIGNMENT % scalar_size 
                     ? MEMORY_ALIGNMENT / scalar_size : 0 ;

  if ( align && align < scalar_count && scalar_count % align ) {
    scalar_count += align - scalar_count % align ;
  }

  return scalar_count ;
}

std::string CudaSpace::query_label( const void * p )
{
  const Impl::MemoryTracking::Info info =
    cuda_space_singleton().query( p );

  return info.label ;
}

void CudaSpace::access_error()
{
  const std::string msg("KokkosArray::CudaSpace::access_error attempt to execute Cuda function from non-Cuda space" );

  KokkosArray::Impl::throw_runtime_exception( msg );
}

void CudaSpace::access_error( const void * const ptr )
{
  std::ostringstream msg ;
  msg << "KokkosArray::CudaSpace::access_error:" ;
  msg << " attempt to access Cuda-data labeled(" ;
  msg << query_label( ptr ) ;
  msg << ") from non-Cuda execution" ;
  KokkosArray::Impl::throw_runtime_exception( msg.str() );
}

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray

