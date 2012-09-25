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

#ifndef KOKKOSARRAY_CUDA_PARALLELFOR_HPP
#define KOKKOSARRAY_CUDA_PARALLELFOR_HPP

#if defined( __CUDACC__ )

#include <KokkosArray_ParallelFor.hpp>

#include <Cuda/KokkosArray_Cuda_Parallel.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelFor< FunctorType , Cuda , Cuda::size_type > {
public:

  const FunctorType     m_work_functor ;
  const Cuda::size_type m_work_count ;  

private:

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

public:

  inline
  __device__
  void operator()(void) const
  {
    const Cuda::size_type work_stride = blockDim.x * gridDim.x ;

    Cuda::size_type iwork = threadIdx.x + blockDim.x * blockIdx.x ;

    for ( ; iwork < m_work_count ; iwork += work_stride ) {
      m_work_functor( iwork );
    }
  }


  ParallelFor( const Cuda::size_type  work_count ,
               const FunctorType    & functor )
    : m_work_functor( functor )
    , m_work_count(  work_count )
    {
      const Cuda::size_type grid_max = cuda_internal_maximum_grid_count();

      const dim3 block( CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1, 1);

      dim3 grid( ( ( work_count + block.x - 1 ) / block.x ) , 1 , 1 );

      if ( grid_max < grid.x ) grid.x = grid_max ;

      CudaParallelLaunch< ParallelFor >( *this , grid , block , 0 );
    }
};

//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelFor< FunctorType , Cuda , CudaWorkConfig > {
public:

  ParallelFor( const CudaWorkConfig & work_config ,
               const FunctorType    & functor )
  {
    const dim3 grid( work_config.grid[0] ,
                     work_config.grid[1] ,
                     work_config.grid[2] );

    const dim3 block( work_config.block[0] ,
                      work_config.block[1] ,
                      work_config.block[2] );

    CudaParallelLaunch< ParallelFor >( functor , grid , block , work_config.shared );
  }
};

} // namespace Impl
} // namespace KokkosArray

#endif /* defined( __CUDACC__ ) */

#endif /* KOKKOSARRAY_CUDA_PARALLELFOR_HPP */

