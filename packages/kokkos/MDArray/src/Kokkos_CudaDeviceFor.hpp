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

#ifndef KOKKOS_CUDADEVICEFOR_HPP
#define KOKKOS_CUDADEVICEFOR_HPP

namespace Kokkos {

class CudaDevice ;

template< class FunctorType , class DeviceType > struct ParallelFor ;

template< class FunctorType >
__global__
static void run_functor_on_cuda(
  const CudaDevice::size_type work_count ,
  const FunctorType functor )
{
  const CudaDevice::size_type work_stride = blockDim.x * gridDim.x ;

  for ( CudaDevice::size_type iwork = threadIdx.x + blockDim.x * blockIdx.x ;
        iwork < work_count ; iwork += work_stride ) {
    functor( iwork );
  }
}

// Partial specialization requires a class:

template< class FunctorType >
struct ParallelFor< FunctorType , CudaDevice > {

  static void run( const CudaDevice::size_type work_count ,
                   const FunctorType & functor )
  {
    unsigned int threadCount = 256 ;
    unsigned int blockCount  = ( work_count + threadCount - 1 ) / threadCount ;

    run_functor_on_cuda<<< blockCount , threadCount >>>( work_count , functor );

    // Reconsider - not necessary, but convenient
    cudaThreadSynchronize();
  }
};

template< typename iType , class FunctorType >
inline
void parallel_for( const iType & work_count , const FunctorType & functor )
{
  ParallelFor< FunctorType , typename FunctorType::device_type >::run( work_count , functor );
}

}

#endif /* KOKKOS_CUDADEVICEFOR_HPP */

