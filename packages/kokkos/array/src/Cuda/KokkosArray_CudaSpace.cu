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

#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_CudaSpace.hpp>

#include <Cuda/KokkosArray_Cuda_Internal.hpp>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>
#include <impl/KokkosArray_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace {

class CudaMemoryTrackingEntry : public Impl::MemoryTrackingEntry
{
public:

  void * const                    ptr_alloc ;
  const unsigned                  size ;
  const unsigned                  count ;
  Impl::cuda_texture_object_type  tex_obj ;

  CudaMemoryTrackingEntry( const std::string & arg_label ,
                           const std::type_info & arg_info ,
                           void * const           arg_ptr ,
                           const unsigned         arg_size ,
                           const unsigned         arg_count )
    : Impl::MemoryTrackingEntry( arg_label , arg_info , arg_ptr , arg_size * arg_count )
    , ptr_alloc( arg_ptr )
    , size( arg_size )
    , count( arg_count )
    , tex_obj( 0 )
    {}

  ~CudaMemoryTrackingEntry();
};

CudaMemoryTrackingEntry::~CudaMemoryTrackingEntry()
{
  cudaDeviceSynchronize();

  if ( tex_obj ) {

  }

  if ( cudaSuccess != cudaFree( ptr_alloc ) ) {
    std::cerr << "cudaFree( " << ptr_alloc << " ) FAILED for " ;
    Impl::MemoryTrackingEntry::print( std::cerr );
  }
}

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

    cuda_space_singleton().insert(
      new CudaMemoryTrackingEntry( label , scalar_type , ptr , scalar_size , scalar_count ) );
  }

  return ptr ;
}

void CudaSpace::increment( const void * ptr )
{
  HostSpace::assert_master_thread( "KokkosArray::CudaSpace::increment" );

  cuda_space_singleton().increment( ptr );
}

void CudaSpace::decrement( const void * ptr )
{
  HostSpace::assert_master_thread( "KokkosArray::CudaSpace::decrement" );

  cuda_space_singleton().decrement( ptr );
}

void CudaSpace::print_memory_view( std::ostream & o )
{
  cuda_space_singleton().print( o , std::string("  ") );
}

//----------------------------------------------------------------------------

size_t CudaSpace::preferred_alignment(
  size_t scalar_size , size_t scalar_count )
{
  const size_t align = 0 == MEMORY_ALIGNMENT % scalar_size 
                     ? MEMORY_ALIGNMENT / scalar_size : 0 ;
  const size_t threshold = align ;

  if ( align && threshold < scalar_count && scalar_count % align ) {
    scalar_count += align - scalar_count % align ;
  }

  return scalar_count ;
}

std::string CudaSpace::query_label( const void * p )
{
  const Impl::MemoryTrackingEntry * entry =
    cuda_space_singleton().query( p );

  return entry ? entry->label : std::string("ERROR NOT FOUND");
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

#if defined( CUDA_VERSION ) && ( 500 <= CUDA_VERSION )

namespace KokkosArray {
namespace Impl {

::cudaTextureObject_t
cuda_texture_object_attach(
  const cudaChannelFormatDesc & desc ,
  const void * const            ptr )
{
  if ( 0 == ptr ) return 0 ;

  const unsigned max_count = 1 << 28 ;

  CudaMemoryTrackingEntry * entry =
    dynamic_cast<CudaMemoryTrackingEntry *>( cuda_space_singleton().query( ptr ) );

  const bool ok_found  = 0 != entry ;
  const bool ok_ptr    = ok_found && ptr == entry->ptr_alloc ;
  const bool ok_count  = ok_found && entry->count < max_count ;

  if ( ok_found && ok_ptr && ok_count ) {

    // Can only create texture object on device architure 3.0 or better

    if ( 0 == entry->tex_obj && 300 <= Cuda::device_arch() ) {

      struct cudaResourceDesc resDesc ;
      struct cudaTextureDesc  texDesc ;

      memset( & resDesc , 0 , sizeof(resDesc) );
      memset( & texDesc , 0 , sizeof(texDesc) );

      resDesc.resType                = cudaResourceTypeLinear ;
      resDesc.res.linear.desc        = desc ;
      resDesc.res.linear.sizeInBytes = entry->size * entry->count ;
      resDesc.res.linear.devPtr      = entry->ptr_alloc ;

      cudaCreateTextureObject( & entry->tex_obj, & resDesc, & texDesc, NULL);
    }
  }
  else {
    std::ostringstream msg ;
    msg << "CudaSpace::texture_object_attach( " << ptr << " ) FAILED: " ;

    if ( ! ok_found ) {
      msg << "Not View allocated" ;
    }
    else if ( ! ok_ptr ) {
      msg << "Not the originally allocated View \"" << entry->label << "\"" ;
    }
    else if ( ! ok_count ) {
      msg << "Cuda texture object limit exceeded "
          << max_count << " <= " << entry->count ;
    }
    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }

  return entry->tex_obj ;
}

} // namespace Impl
} // namespace KokkosArray

#endif


