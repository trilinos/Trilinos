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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <Kokkos_DeviceCuda.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_DeepCopy.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace {

void cuda_safe_call( cudaError e , const char * name )
{
  if ( cudaSuccess != e ) {
    std::ostringstream out ;
    out << name << " error: " << cudaGetErrorString(e);
    throw std::runtime_error( out.str() );
  }
}

#define CUDA_SAFE_CALL( call )  cuda_safe_call( call , # call )

/*--------------------------------------------------------------------------*/

class DeviceCuda_Impl {
public:
  Impl::MemoryInfoSet            m_allocations ;
  struct cudaDeviceProp          m_cudaProp ;
  int                            m_cudaDev ;
  unsigned                       m_maxWarp ;
  DeviceCuda::Traits::WordType * m_reduceScratchSpace ;
  DeviceCuda::Traits::WordType * m_reduceScratchFlag ;

  explicit DeviceCuda_Impl( int cuda_device_id );
  ~DeviceCuda_Impl();

  static DeviceCuda_Impl & singleton( int cuda_device_id = 0 );

  void * allocate_memory(
    const std::string    & label ,
    const std::type_info & type ,
    const size_t member_size ,
    const size_t member_count );

  void deallocate_memory( void * ptr );

private:
  DeviceCuda_Impl();
  DeviceCuda_Impl( const DeviceCuda_Impl & );
  DeviceCuda_Impl & operator = ( const DeviceCuda_Impl & );
};

DeviceCuda_Impl & DeviceCuda_Impl::singleton( int cuda_device_id )
{
  static DeviceCuda_Impl self( cuda_device_id );
  return self ;
}

DeviceCuda_Impl::DeviceCuda_Impl( int cuda_device_id )
  : m_allocations()
  , m_cudaProp()
  , m_cudaDev( cuda_device_id )
  , m_maxWarp( 0 )
  , m_reduceScratchSpace( 0 )
  , m_reduceScratchFlag( 0 )
{
  // Some significant cuda device properties:
  //
  // cudaDeviceProp::major               : Device major number
  // cudaDeviceProp::minor               : Device minor number
  // cudaDeviceProp::multiProcessorCount : number of multiprocessors
  // cudaDeviceProp::sharedMemPerBlock   : capacity of shared memory per block
  // cudaDeviceProp::totalGlobalMem      : capacity of global memory

  enum { n = sizeof(DeviceCuda::Traits::WordType) };

  const DeviceCuda::Traits::WordType zero = 0 ;

  // Device query

  CUDA_SAFE_CALL( cudaGetDevice( & m_cudaDev ) );
  CUDA_SAFE_CALL( cudaGetDeviceProperties( & m_cudaProp , m_cudaDev ) );

  // Maximum number of warps,
  // at most one warp per thread in a warp for reduction.

  m_maxWarp = DeviceCuda::Traits::WarpSize ;
  while ( m_cudaProp.maxThreadsPerBlock <
          DeviceCuda::Traits::WarpSize * m_maxWarp ) {
    m_maxWarp >>= 1 ;
  }

  // Allocate shared memory image for multiblock reduction scratch space

  const size_t sharedWord =
   ( m_cudaProp.sharedMemPerBlock + n - 1 ) / n ;

  m_reduceScratchSpace = (DeviceCuda::Traits::WordType *)
    allocate_memory( std::string("MultiblockReduceScratchSpace") ,
                     typeid( DeviceCuda::Traits::WordType ),
                     sizeof( DeviceCuda::Traits::WordType ),
                     sharedWord + 1 );

  m_reduceScratchFlag = m_reduceScratchSpace + sharedWord ;

  CUDA_SAFE_CALL(
    cudaMemcpy( m_reduceScratchFlag , & zero , n, cudaMemcpyHostToDevice ) );
}

DeviceCuda_Impl::~DeviceCuda_Impl()
{
  deallocate_memory( m_reduceScratchSpace );

  m_reduceScratchSpace = 0 ;
  m_reduceScratchFlag  = 0 ;

  if ( ! m_allocations.empty() ) {
    std::cerr << "Kokkos::DeviceCuda memory leaks:" << std::endl ;
    m_allocations.print( std::cerr );
  }
}

void * DeviceCuda_Impl::allocate_memory(
  const std::string    & label ,
  const std::type_info & type ,
  const size_t member_size ,
  const size_t member_count )
{
  Impl::MemoryInfo tmp ;

  tmp.m_type  = & type ;
  tmp.m_label = label ;
  tmp.m_size  = member_size ;
  tmp.m_count = member_count ;

  CUDA_SAFE_CALL( cudaMalloc( & tmp.m_ptr , member_size * member_count ) );

  const bool ok_alloc  = 0 != tmp.m_ptr ;
  const bool ok_insert = ok_alloc && m_allocations.insert( tmp );

  if ( ! ok_alloc || ! ok_insert ) {
    std::ostringstream msg ;
    msg << "Kokkos::DeviceCuda::allocate_memory( " << label
        << " , " << type.name()
        << " , " << member_size
        << " , " << member_count
        << " ) FAILED " ;
    if ( ok_alloc ) { msg << "memory allocation" ; }
    else            { msg << "with internal error" ; }
    throw std::runtime_error( msg.str() );
  }

  return tmp.m_ptr ;
}

void DeviceCuda_Impl::deallocate_memory( void * ptr )
{
  if ( ! m_allocations.erase( ptr ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::DeviceCuda::deallocate_memory( " << ptr
        << " ) FAILED memory allocated by this device" ;
    throw std::runtime_error( msg.str() );
  }

  CUDA_SAFE_CALL( cudaFree( ptr ) );
}

}

/*--------------------------------------------------------------------------*/

void DeviceCuda::initialize( int cuda_device_id )
{
  DeviceCuda_Impl::singleton( cuda_device_id );
}

DeviceCuda::Traits::WordType *
DeviceCuda::reduce_multiblock_scratch_space()
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();
  return s.m_reduceScratchSpace ;
}

DeviceCuda::Traits::WordType *
DeviceCuda::reduce_multiblock_scratch_flag()
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();
  return s.m_reduceScratchFlag ;
}

DeviceCuda::size_type
DeviceCuda::maximum_warp_count()
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();
  return s.m_maxWarp ;
}

DeviceCuda::size_type
DeviceCuda::maximum_grid_count()
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();
  return s.m_cudaProp.maxGridSize[0];
}

/*--------------------------------------------------------------------------*/

void * DeviceCuda::allocate_memory(
  const std::string    & label ,
  const std::type_info & type ,
  const size_t member_size ,
  const size_t member_count )
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();

  return s.allocate_memory( label ,type , member_size , member_count );
}

void DeviceCuda::deallocate_memory( void * ptr )
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();

  s.deallocate_memory( ptr );
}

void DeviceCuda::print_memory_view( std::ostream & o )
{
  DeviceCuda_Impl & s = DeviceCuda_Impl::singleton();

  s.m_allocations.print( o );
}

/*--------------------------------------------------------------------------*/

unsigned int DeviceCuda::m_launching_kernel = false ;

void DeviceCuda::set_dispatch_functor()
{
  if ( m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceCuda::set_dispatch_functor FAILED: " );
    msg.append( "kernel dispatch is already in progress, " );
    msg.append( "a recursive call or forgotten 'clear_dispatch_kernel" );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = true ;
}

void DeviceCuda::clear_dispatch_functor()
{
  if ( ! m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceCuda::clear_dispatch_functor FAILED: " );
    msg.append( "no kernel dispatch in progress." );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = false ;
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void copy_to_cuda_from_host( void * dst , const void * src ,
                             size_t member_size , size_t member_count )
{
  CUDA_SAFE_CALL(
    cudaMemcpy( dst , src , member_size * member_count , cudaMemcpyHostToDevice ) );

}

void copy_to_host_from_cuda( void * dst , const void * src ,
                             size_t member_size , size_t member_count )
{
  CUDA_SAFE_CALL(
    cudaMemcpy( dst , src , member_size * member_count , cudaMemcpyDeviceToHost ) );
}

}
}

