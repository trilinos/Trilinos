/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOS_CUDA_PARALLEL_HPP
#define KOKKOS_CUDA_PARALLEL_HPP

#if defined( __CUDACC__ )

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_CudaExec.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelFor< FunctorType , WorkSpec /* size_t */ , Cuda > {
private:

  const FunctorType     m_functor ;
  const Cuda::size_type m_work ;  

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

public:

  inline
  __device__
  void operator()(void) const
  {
    const Cuda::size_type work_stride = blockDim.x * gridDim.x ;

    for ( Cuda::size_type
            iwork = threadIdx.x + blockDim.x * blockIdx.x ;
            iwork < m_work ;
            iwork += work_stride ) {
      m_functor( iwork );
    }
  }

  ParallelFor( const FunctorType  & functor ,
               const size_t         work )
    : m_functor( functor )
    , m_work(    work )
    {
      const dim3 block( CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1, 1);
      const dim3 grid( std::min( ( m_work + block.x - 1 ) / block.x , cuda_internal_maximum_grid_count() ) , 1 , 1 );

      CudaParallelLaunch< ParallelFor >( *this , grid , block , 0 );
    }
};

template< class FunctorType >
class ParallelFor< FunctorType , ParallelWorkRequest , Cuda > {
private:

  const FunctorType          m_functor ;
  const ParallelWorkRequest  m_work ;
  const int                  m_shmem ;

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

public:

  inline
  __device__
  void operator()(void) const
  {
    CudaExec exec( 0 , m_shmem );
    m_functor( Cuda( exec ) );
  }

  ParallelFor( const FunctorType         & functor ,
               const ParallelWorkRequest &  work )
    : m_functor( functor )
    , m_work( std::min( work.league_size , size_t(cuda_internal_maximum_grid_count()) ) ,
              std::min( work.team_size ,   size_t(CudaTraits::WarpSize * cuda_internal_maximum_warp_count()) ) )
    , m_shmem( FunctorShmemSize< FunctorType >::value( functor ) )
    {
      const dim3 grid(  m_work.league_size , 1 , 1 );
      const dim3 block( m_work.team_size , 1, 1 );

      CudaParallelLaunch< ParallelFor >( *this , grid , block , m_shmem );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , CudaWorkConfig , Cuda > {
public:

  const FunctorType m_work_functor ;

  inline
  __device__
  void operator()(void) const
  {
    Cuda::size_type iwork = threadIdx.x + blockDim.x * (
                            threadIdx.y + blockDim.y * (
                            threadIdx.z + blockDim.z * (
                            blockIdx.x + gridDim.x * (
                            blockIdx.y + gridDim.y * (
                            blockIdx.z )))));

    m_work_functor( iwork );
  }

  ParallelFor( const FunctorType    & functor ,
               const CudaWorkConfig & work_config )
  : m_work_functor( functor )
  {
    const dim3 grid( work_config.grid[0] ,
                     work_config.grid[1] ,
                     work_config.grid[2] );

    const dim3 block( work_config.block[0] ,
                      work_config.block[1] ,
                      work_config.block[2] );

    CudaParallelLaunch< ParallelFor >( *this , grid , block , work_config.shared );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */

