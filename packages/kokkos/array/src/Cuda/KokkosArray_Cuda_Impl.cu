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

/*--------------------------------------------------------------------------*/
/* KokkosArray interfaces */

#include <KokkosArray_Cuda.hpp>
#include <Cuda/KokkosArray_Cuda_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

void cuda_internal_error_throw( cudaError e , const char * name )
{
  std::ostringstream out ;
  out << name << " error: " << cudaGetErrorString(e);
  throw std::runtime_error( out.str() );
}

//----------------------------------------------------------------------------
// Some significant cuda device properties:
//
// cudaDeviceProp::name                : Text label for device
// cudaDeviceProp::major               : Device major number
// cudaDeviceProp::minor               : Device minor number
// cudaDeviceProp::warpSize            : number of threads per warp
// cudaDeviceProp::multiProcessorCount : number of multiprocessors
// cudaDeviceProp::sharedMemPerBlock   : capacity of shared memory per block
// cudaDeviceProp::totalConstMem       : capacity of constant memory
// cudaDeviceProp::totalGlobalMem      : capacity of global memory
// cudaDeviceProp::maxGridSize[3]      : maximum grid size

//
//  Section 4.4.2.4 of the CUDA Toolkit Reference Manual
//
// struct cudaDeviceProp {
//   char name[256];
//   size_t totalGlobalMem;
//   size_t sharedMemPerBlock;
//   int regsPerBlock;
//   int warpSize;
//   size_t memPitch;
//   int maxThreadsPerBlock;
//   int maxThreadsDim[3];
//   int maxGridSize[3];
//   size_t totalConstMem;
//   int major;
//   int minor;
//   int clockRate;
//   size_t textureAlignment;
//   int deviceOverlap;
//   int multiProcessorCount;
//   int kernelExecTimeoutEnabled;
//   int integrated;
//   int canMapHostMemory;
//   int computeMode;
//   int concurrentKernels;
//   int ECCEnabled;
//   int pciBusID;
//   int pciDeviceID;
//   int tccDriver;
// };


namespace {

class CudaInternalDevices {
public:
  enum { MAXIMUM_DEVICE_COUNT = 8 };
  struct cudaDeviceProp  m_cudaProp[ MAXIMUM_DEVICE_COUNT ] ;
  int                    m_cudaDevCount ;

  CudaInternalDevices();

  static const CudaInternalDevices & singleton();
};

CudaInternalDevices::CudaInternalDevices()
{
  // See 'cudaSetDeviceFlags' for host-device thread interaction
  // Section 4.4.2.6 of the CUDA Toolkit Reference Manual

  cudaGetDeviceCount( & m_cudaDevCount );

  for ( int i = 0 ; i < m_cudaDevCount ; ++i ) {
    cudaGetDeviceProperties( m_cudaProp + i , i );
  }
}

const CudaInternalDevices & CudaInternalDevices::singleton()
{
  static CudaInternalDevices self ; return self ;
}

}

//----------------------------------------------------------------------------

class CudaInternal {
private:

  CudaInternal( const CudaInternal & );
  CudaInternal & operator = ( const CudaInternal & );

public:

  typedef Cuda::size_type size_type ;

  std::vector<cudaStream_t> m_streams ;
  int                       m_cudaDev ;
  unsigned                  m_maxWarpCount ;
  unsigned                  m_maxBlock ;
  unsigned                  m_maxSharedWords ;
  size_type               * m_reduceScratchSpace ;
  size_type               * m_reduceScratchFlag ;

  static CudaInternal & raw_singleton();
  static const CudaInternal & singleton();

  void initialize( int cuda_device_id );
  void finalize();

  ~CudaInternal();

  CudaInternal()
    : m_streams() 
    , m_cudaDev( -1 )
    , m_maxWarpCount( 0 )
    , m_maxBlock( 0 ) 
    , m_maxSharedWords( 0 )
    , m_reduceScratchSpace( 0 )
    , m_reduceScratchFlag( 0 )
    {}
};

CudaInternal::~CudaInternal()
{
  if ( m_reduceScratchSpace ) {
    std::cerr << "KokkosArray::Cuda ERROR: Failed to call KokkosArray::Cuda::finalize()"
              << std::endl ;
    std::cerr.flush();
  }
}

CudaInternal & CudaInternal::raw_singleton()
{ static CudaInternal self ; return self ; }

const CudaInternal & CudaInternal::singleton()
{
  const CudaInternal & self = raw_singleton();

  if ( self.m_cudaDev == -1 ) {
    throw std::runtime_error(std::string("CATASTROPHIC FAILURE: Using KokkosArray::Cuda before calling KokkosArray::Cuda::initialize(...)"));
  }

  return self ;
}

void CudaInternal::initialize( int cuda_device_id )
{
  enum { WordSize = sizeof(size_type) };

  const CudaInternalDevices & dev_info = CudaInternalDevices::singleton();

  const bool ok_init = 0 == m_reduceScratchSpace ;
  const bool ok_id   = 0 <= cuda_device_id &&
                            cuda_device_id < dev_info.m_cudaDevCount ;

  // Need device capability 1.3 or better

  const bool ok_dev = ok_id &&
    ( 1 < dev_info.m_cudaProp[ cuda_device_id ].major ||
      2 < dev_info.m_cudaProp[ cuda_device_id ].minor );

  if ( ok_init && ok_dev ) {

    const struct cudaDeviceProp & cudaProp =
      dev_info.m_cudaProp[ cuda_device_id ];

    m_cudaDev = cuda_device_id ;

    CUDA_SAFE_CALL( cudaSetDevice( m_cudaDev ) );
  
    //----------------------------------
    // Maximum number of warps,
    // at most one warp per thread in a warp for reduction.

    // HCE 02/02/2012 :
    // Found bug in CUDA 4.1 that sometimes a kernel launch would fail
    // if the thread count == 1024 and a functor is passed to the kernel.
    // Copying the kernel to constant memory and then launching with
    // thread count == 1024 would work fine.
    // All compute capabilities support at least 16 warps (512 threads).

    m_maxWarpCount = 16 ;

    // m_maxWarpCount = cudaProp.maxThreadsPerBlock / Impl::CudaTraits::WarpSize ;

    if ( Impl::CudaTraits::WarpSize < m_maxWarpCount ) {
      m_maxWarpCount = Impl::CudaTraits::WarpSize ;
    }

    //----------------------------------

    m_maxSharedWords = cudaProp.sharedMemPerBlock / WordSize ;

    //----------------------------------
    {
      // Set the maximum number of blocks to the maximum number of
      // resident blocks per multiprocessor times the number of
      // multiprocessors on the device.  Once a block is active on
      // a multiprocessor then let it do all the work that it can.

      enum { MaxResidentBlocksPerMultiprocessor = 8 };

      const unsigned minBlock = cudaProp.multiProcessorCount *
                                MaxResidentBlocksPerMultiprocessor ;

      m_maxBlock = 1 ;
      while ( m_maxBlock < minBlock &&
              ( m_maxBlock << 1 ) <= (unsigned) cudaProp.maxGridSize[0] ) {
        m_maxBlock <<= 1 ;
      }
    }

    //----------------------------------
    // Allocate a parallel stream for each multiprocessor
    // to support concurrent heterogeneous multi-functor execution.

    m_streams.resize( cudaProp.multiProcessorCount );

    for ( int i = 0 ; i < cudaProp.multiProcessorCount ; ++i ) {
      CUDA_SAFE_CALL( cudaStreamCreate(    & m_streams[i] ) );
      CUDA_SAFE_CALL( cudaStreamSynchronize( m_streams[i] ) );
    }

    //----------------------------------
    // Allocate shared memory image for multiblock reduction scratch space

    const size_t sharedWord =
     ( cudaProp.sharedMemPerBlock + WordSize - 1 ) / WordSize ;

    m_reduceScratchSpace = (size_type *)
      CudaMemorySpace::allocate(
        std::string("MultiblockReduceScratchSpace") ,
        typeid( size_type ),
        sizeof( size_type ),
        sharedWord + 1 );

    m_reduceScratchFlag = m_reduceScratchSpace + sharedWord ;

    CUDA_SAFE_CALL( cudaMemset( m_reduceScratchFlag , 0 , WordSize ) );
  }
  else {

    std::ostringstream msg ;
    msg << "KokkosArray::Cuda::initialize(" << cuda_device_id << ") FAILED" ;

    if ( ! ok_init ) {
      msg << " : Already initialized" ;
    }
    if ( ! ok_id ) {
      msg << " : Device identifier out of range "
          << "[0.." << dev_info.m_cudaDevCount << "]" ;
    }
    else if ( ! ok_dev ) {
      msg << " : Device " ;
      msg << dev_info.m_cudaProp[ cuda_device_id ].major ;
      msg << "." ;
      msg << dev_info.m_cudaProp[ cuda_device_id ].minor ;
      msg << " has insufficient capability, required 1.3 or better" ;
    }
    throw std::runtime_error( msg.str() );
  } 
}

void CudaInternal::finalize()
{
  for ( size_t i = m_streams.size() ; i ; ) {
    --i ;
    CUDA_SAFE_CALL( cudaStreamSynchronize( m_streams[i] ) );
    CUDA_SAFE_CALL( cudaStreamDestroy( m_streams[i] ) );
  }

  CudaMemorySpace::decrement( m_reduceScratchSpace );

  m_streams.clear();
  m_cudaDev            = -1 ;
  m_maxWarpCount       = 0 ;
  m_maxBlock           = 0 ; 
  m_maxSharedWords     = 0 ;
  m_reduceScratchSpace = 0 ;
  m_reduceScratchFlag  = 0 ;
}

//----------------------------------------------------------------------------

Cuda::size_type cuda_internal_maximum_warp_count()
{ return CudaInternal::singleton().m_maxWarpCount ; }

Cuda::size_type cuda_internal_maximum_grid_count()
{ return CudaInternal::singleton().m_maxBlock ; }

Cuda::size_type cuda_internal_maximum_shared_words()
{ return CudaInternal::singleton().m_maxSharedWords ; }

Cuda::size_type cuda_internal_stream_count()
{ return CudaInternal::singleton().m_streams.size(); }

cudaStream_t & cuda_internal_stream( Cuda::size_type i )
{
  CudaInternal & s = CudaInternal::raw_singleton();

  if ( s.m_streams.size() <= i ) {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::cuda_internal_stream( " << i << " ) ERROR "
        << "stream_count = " << s.m_streams.size() ;
    throw std::logic_error( msg.str() );
  }

  return s.m_streams[i] ;
}

Cuda::size_type * cuda_internal_reduce_multiblock_scratch_space()
{ return CudaInternal::singleton().m_reduceScratchSpace ; }

Cuda::size_type * cuda_internal_reduce_multiblock_scratch_flag()
{ return CudaInternal::singleton().m_reduceScratchFlag ; }

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {

Cuda::size_type Cuda::detect_device_count()
{ return Impl::CudaInternalDevices::singleton().m_cudaDevCount ; }

void Cuda::initialize( const Cuda::SelectDevice config )
{ Impl::CudaInternal::raw_singleton().initialize( config.cuda_device_id ); }

void Cuda::finalize()
{ Impl::CudaInternal::raw_singleton().finalize(); }

bool Cuda::sleep() { return false ; }

bool Cuda::wake() { return true ; }

void Cuda::fence()
{ cudaThreadSynchronize(); }

} // namespace KokkosArray

//----------------------------------------------------------------------------

