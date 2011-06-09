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

#ifndef KOKKOS_DEVICECUDA_PARALLELFOR_HPP
#define KOKKOS_DEVICECUDA_PARALLELFOR_HPP

#include <DeviceCuda/Kokkos_DeviceCuda_ParallelDriver.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------


template< class DriverType >
__global__

#if KOKKOS_DEVICE_CUDA_USE_CONSTANT_MEMORY 

static void cuda_parallel_for()
{
  // The driver functor has been copied to constant memory
  const DriverType & driver = 
    *((const DriverType *) kokkos_device_cuda_constant_memory_buffer );

#else

static void cuda_parallel_for( const DriverType driver )
{

#endif

  typedef DeviceCuda::size_type size_type ;

  const size_type work_count  = driver.m_work_count ;
  const size_type work_stride = driver.m_work_stride ;

  size_type iwork = threadIdx.x + blockDim.x * blockIdx.x ;

  for ( ; iwork < work_count ; iwork += work_stride ) {
    driver.m_work_functor( iwork );
  }
}

//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelFor< FunctorType , DeviceCuda > {
public:

  const FunctorType           m_work_functor ;
  const DeviceCuda::size_type m_work_count ;  
  const DeviceCuda::size_type m_work_stride ;  

private:

  ParallelFor( const FunctorType & functor ,
               const size_t work_count ,
               const size_t work_stride )
    : m_work_functor( functor )
    , m_work_count(  work_count )
    , m_work_stride( work_stride )
    {}

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

  ParallelFor( const ParallelFor & rhs )
    : m_work_functor( rhs.m_work_functor )
    , m_work_count(   rhs.m_work_count )
    , m_work_stride(  rhs.m_work_stride )
    {}

public:

  static void execute( const size_t work_count ,
                       const FunctorType & functor )
  {
    typedef ParallelFor< FunctorType , DeviceCuda > self_type ;

    enum { WarpSize = DeviceCudaTraits::WarpSize };

    const size_t grid_max = DeviceCuda::maximum_grid_count();

    const dim3 block( WarpSize * DeviceCuda::maximum_warp_count() , 1 , 1 );

    dim3 grid( ( ( work_count + block.x - 1 ) / block.x ) , 1 , 1 );

    if ( grid_max < grid.x ) grid.x = grid_max ;

    DeviceCuda::set_dispatch_functor();

    const self_type driver( functor , work_count , block.x * grid.x );

    DeviceCuda::clear_dispatch_functor();

#if KOKKOS_DEVICE_CUDA_USE_CONSTANT_MEMORY 

  // Copy functor to constant memory on the device
  cudaMemcpyToSymbol( kokkos_device_cuda_constant_memory_buffer , & driver , sizeof(self_type) );

  // Invoke the driver function on the device
  cuda_parallel_for< self_type > <<< grid , block >>>();

#else

  cuda_parallel_for< self_type > <<< grid , block >>>( driver );

#endif
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

#include <Kokkos_DeviceClear_macros.hpp>

#endif /* KOKKOS_DEVICECUDA_PARALLELFOR_HPP */

