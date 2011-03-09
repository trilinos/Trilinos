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

#include <stdexcept>
#include <map>
#include <ostream>
#include <sstream>
#include <Kokkos_CudaDevice.hpp>

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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// For the cuda device singleton.

class CudaDeviceImpl {
public:
  std::map<void*,std::string> m_allocations ;
  struct cudaDeviceProp       m_cudaProp ;
  int                         m_cudaDev ;

  // Appropriate cached device information

  CudaDeviceImpl();

  static CudaDeviceImpl & singleton();
};

CudaDeviceImpl::CudaDeviceImpl()
  : m_allocations()
{
  // Appropriate device queries
  m_cudaDev = 0 ;

  CUDA_SAFE_CALL( cudaGetDevice( & m_cudaDev ) );
  CUDA_SAFE_CALL( cudaGetDeviceProperties( & m_cudaProp , m_cudaDev ) );
}

CudaDeviceImpl & CudaDeviceImpl::singleton()
{
  static CudaDeviceImpl * impl = NULL ;
  if ( impl == NULL ) { impl = new CudaDeviceImpl(); }
  return *impl ;
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

__global__
void clear_memory( int * ptr , int count )
{
  const int n = count ;
  const int s = blockDim.x * gridDim.x ;

  int i = threadIdx.x + blockDim.x * blockIdx.x ;
  
  for ( ; i < n ; i += s ) ptr[i] = 0 ;
}

void * CudaDevice::allocate_memory( size_type member_size ,
                                    size_type member_count ,
                                    const std::string & label )
{
  int * ptr_on_device = NULL ;

  dim3 dimGrid(  256 , 1 , 1 );
  dim3 dimBlock( 256 , 1 , 1 );

  // Require member_size be a multiple of word size?

  CUDA_SAFE_CALL( cudaMalloc( & ptr_on_device , member_size * member_count ) );

  const size_type word_count = ( member_size * member_count + sizeof(unsigned) - 1 ) / sizeof(unsigned);

  while ( word_count <= ( dimBlock.x * ( dimGrid.x >> 1 ) ) ) {
    dimGrid.x >>= 1 ;
  }

  clear_memory<<< dimGrid , dimBlock >>>( ptr_on_device , word_count );

  CudaDeviceImpl::singleton().m_allocations[ ptr_on_device ] = label ;

  return ptr_on_device ;
}

void CudaDevice::deallocate_memory( void * ptr_on_device )
{
  if ( 1 != CudaDeviceImpl::singleton().m_allocations.count( ptr_on_device ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::CudaDevice::deallocate_memory("
        << ptr_on_device
        << ") FAILED invalid address" ;
    throw std::runtime_error( msg.str() );
  }

  CudaDeviceImpl::singleton().m_allocations.erase( ptr_on_device );

  CUDA_SAFE_CALL( cudaFree( ptr_on_device ) );
}

void CudaDevice::print_allocations( std::ostream & s ) const
{
  CudaDeviceImpl & impl = CudaDeviceImpl::singleton();

  std::map<void*,std::string>::const_iterator i = impl.m_allocations.begin();
  std::map<void*,std::string>::const_iterator end = impl.m_allocations.end();

  for ( ; i != end ; ++i ) {
    s << i->second << std::endl ;
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/* Define at most two thread blocks per physical multiprocessor */

CudaDevice::size_type
CudaDevice::block_count_max()
{
  CudaDeviceImpl & impl = CudaDeviceImpl::singleton();

  size_type max = impl.m_cudaProp.maxGridSize[0] ;

  // Recomended for good performance 2 * impl.m_cudaProp.multiProcessorCount

  if ( max > 4 * impl.m_cudaProp.multiProcessorCount ) { 
    max = 4 * impl.m_cudaProp.multiProcessorCount ;
  }

  size_type nblock = 1 ;

  while ( ( nblock << 1 ) < max ) { nblock <<= 1 ; }

  return nblock ;
}

/* Final reduction limits: threads per block or shared memory per block ? */

CudaDevice::size_type
CudaDevice::reduction_thread_max( CudaDevice::size_type shmemPerThread )
{
  CudaDeviceImpl & impl = CudaDeviceImpl::singleton();

  size_type nthread = impl.m_cudaProp.maxThreadsPerBlock ;
  size_type maxThreadReduction =
    impl.m_cudaProp.sharedMemPerBlock / shmemPerThread ;

  if ( maxThreadReduction > impl.m_cudaProp.maxThreadsDim[0] ) {
    maxThreadReduction = impl.m_cudaProp.maxThreadsDim[0] ;
  }

  if ( maxThreadReduction > impl.m_cudaProp.maxGridSize[0] ) {
    maxThreadReduction = impl.m_cudaProp.maxGridSize[0] ;
  }

  while ( maxThreadReduction < nthread ) { nthread >>= 1 ; }

  return nthread ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos



