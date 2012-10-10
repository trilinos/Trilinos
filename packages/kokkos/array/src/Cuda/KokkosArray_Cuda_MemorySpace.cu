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
#include <stdexcept>
#include <sstream>

#include <KokkosArray_Cuda.hpp>
#include <Cuda/KokkosArray_Cuda_MemorySpace.hpp>
#include <Cuda/KokkosArray_Cuda_Internal.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

namespace {

class CudaMemoryImpl {
public:
  MemoryTracking m_allocations ;

  CudaMemoryImpl();
  ~CudaMemoryImpl();

  static CudaMemoryImpl & singleton();
};

CudaMemoryImpl::CudaMemoryImpl()
  : m_allocations()
{}

CudaMemoryImpl & CudaMemoryImpl::singleton()
{
  static CudaMemoryImpl self ;
  return self ;
}

CudaMemoryImpl::~CudaMemoryImpl()
{
  if ( ! m_allocations.empty() ) {
    std::cerr << "KokkosArray::Cuda memory leaks:" << std::endl ;
    m_allocations.print( std::cerr , std::string("  ") );
  }
}

}

/*--------------------------------------------------------------------------*/

void CudaMemorySpace::copy_to_device_from_device(
  void * dst, const void * src, size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

void CudaMemorySpace::copy_to_device_from_host(
  void * dst, const void * src, size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

void CudaMemorySpace::copy_to_host_from_device(
  void * dst, const void * src, size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}



/*--------------------------------------------------------------------------*/

void * CudaMemorySpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  CudaMemoryImpl & s = CudaMemoryImpl::singleton();

  const size_t size = scalar_size * scalar_count ;

  void * ptr = 0 ;

  if ( 0 < scalar_size * scalar_count ) {
    bool ok = true ;

    if ( ok ) ok = cudaSuccess == cudaMalloc( & ptr , size );
    if ( ok ) ok = 0 != ptr ;
    if ( ok ) ok = cudaSuccess == cudaMemset( ptr , 0 , size );
    if ( ok ) ok = cudaSuccess == cudaThreadSynchronize();

    if ( ! ok ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::CudaMemorySpace::allocate( "
          << label
          << " , " << scalar_type.name()
          << " , " << scalar_size
          << " , " << scalar_count
          << " ) FAILED memory allocation" ;
      throw std::runtime_error( msg.str() );
    }

    s.m_allocations.track( ptr, & scalar_type, scalar_size, scalar_count, label );
  }

  return ptr ;
}

#if ! defined( __CUDA_ARCH__ )

void CudaMemorySpace::increment( const void * ptr )
{
  if ( 0 != ptr ) {
    CudaMemoryImpl & s = CudaMemoryImpl::singleton();

    s.m_allocations.increment( ptr );
  }
}

void CudaMemorySpace::decrement( const void * ptr )
{
  if ( 0 != ptr ) {

    CudaMemoryImpl & s = CudaMemoryImpl::singleton();

    void * ptr_alloc = s.m_allocations.decrement( ptr );

    if ( 0 != ptr_alloc ) {
      const bool failed = cudaSuccess != cudaFree( ptr_alloc );

      if ( failed ) {
        std::string msg("KokkosArray::Impl::CudaMemorySpace::decrement() failed cudaFree");
        throw std::runtime_error( msg );
      }
    }
  }
}

#endif

void CudaMemorySpace::print_memory_view( std::ostream & o )
{
  CudaMemoryImpl & s = CudaMemoryImpl::singleton();

  s.m_allocations.print( o , std::string("  ") );
}


size_t CudaMemorySpace::preferred_alignment(
  size_t scalar_size , size_t scalar_count )
{
  const size_t alignment = Impl::CudaTraits::WarpSize * sizeof(Cuda::size_type);

  // If the array is larger than the warp-alignment
  // then align the count on the warp boundary.

  if ( alignment < scalar_size * scalar_count &&
       0 == alignment % scalar_size ) {
    const size_t align = alignment / scalar_size ;
    const size_t rem   = scalar_count % align ;
    if ( rem ) scalar_count += align - rem ;
  }
  return scalar_count ;
}

/*--------------------------------------------------------------------------*/

} // namespace Impl
} // namespace KokkosArray

