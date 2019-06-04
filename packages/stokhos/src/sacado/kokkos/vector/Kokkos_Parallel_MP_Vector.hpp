// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_PARALLEL_MP_VECTOR_HPP
#define KOKKOS_PARALLEL_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_Core.hpp"

//----------------------------------------------------------------------------
// Kokkos execution policies useful for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

/*!
 * \brief Team-based parallel work configuration for Sacado::MP::Vector
 */
template< class ExecSpace >
struct MPVectorWorkConfig {

  typedef MPVectorWorkConfig execution_policy ;
  typedef ExecSpace          execution_space ;
  typedef void               work_tag ;

  size_t range;
  size_t team;
  size_t shared;

  MPVectorWorkConfig( const size_t range_,
                      const size_t team_,
                      const size_t shared_ = 0 ) :
    range(range_), team(team_), shared(shared_) {}
};

namespace Impl {

#if defined( KOKKOS_ENABLE_THREADS )
// Specialization of ParallelFor<> for MPVectorWorkConfig and Threads
// The default implementation ignores the team size and uses the standard
// work-range implementation.  In the future maybe we should try and use
// hyperthreads in a useful way.  That would require:
//   -- interpreting the team-size differently, rather as the sacado size
//   -- determining the vector size of the architecture
//   -- laying out the threads differently to use hyperthreads across the
//      the sacado dimension
template< class FunctorType >
class ParallelFor< FunctorType , MPVectorWorkConfig< Threads > > :
  public ParallelFor< FunctorType , Kokkos::RangePolicy< Threads > > {
  typedef Kokkos::RangePolicy< Threads > Policy ;
public:
  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< Threads > & work_config ) :
    ParallelFor< FunctorType , Policy >( functor ,
                                         Policy( 0, work_config.range ) ) {}
};
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
// Specialization of ParallelFor<> for MPVectorWorkConfig and OpenMP
// The default implementation ignores the team size and uses the standard
// work-range implementation.  In the future maybe we should try and use
// hyperthreads in a useful way.  That would require:
//   -- interpreting the team-size differently, rather as the sacado size
//   -- determining the vector size of the architecture
//   -- laying out the threads differently to use hyperthreads across the
//      the sacado dimension
template< class FunctorType >
class ParallelFor< FunctorType , MPVectorWorkConfig< OpenMP > > :
  public ParallelFor< FunctorType , Kokkos::RangePolicy< OpenMP > > {
  typedef Kokkos::RangePolicy< OpenMP > Policy ;
public:
  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< OpenMP > & work_config ) :
    ParallelFor< FunctorType , Policy >( functor ,
                                         Policy( 0, work_config.range ) ) {}
};
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
// Specialization of ParallelFor<> for MPVectorWorkConfig and Serial
// The default implementation ignores the team size and uses the standard
// work-range implementation.  In the future maybe we should try and use
// hyperthreads in a useful way.  That would require:
//   -- interpreting the team-size differently, rather as the sacado size
//   -- determining the vector size of the architecture
//   -- laying out the threads differently to use hyperthreads across the
//      the sacado dimension
template< class FunctorType >
class ParallelFor< FunctorType , MPVectorWorkConfig< Serial > > :
  public ParallelFor< FunctorType , Kokkos::RangePolicy< Serial > > {
  typedef Kokkos::RangePolicy< Serial > Policy ;
public:
  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< Serial > & work_config ) :
    ParallelFor< FunctorType , Policy >( functor ,
                                         Policy( 0, work_config.range ) ) {}
};
#endif // defined(KOKKOS_ENABLE_SERIAL)

#if defined( KOKKOS_ENABLE_CUDA ) && defined( __CUDACC__ )

// Specialization of ParallelFor<> for MPVectorWorkConfig on Cuda
// Here we use threadIdx.x for each entry in the specified team-size
template< class FunctorType >
class ParallelFor< FunctorType , MPVectorWorkConfig< Cuda > > {
public:

  const FunctorType m_functor ;
  const MPVectorWorkConfig< Cuda > m_config;
  const Cuda::size_type m_work ;

  inline
  __device__
  void operator()(void) const
  {
    const Cuda::size_type work_stride = blockDim.y * gridDim.x ;

    for ( Cuda::size_type iwork = threadIdx.y + blockDim.y * blockIdx.x ;
          iwork < m_work ;
          iwork += work_stride ) {
      m_functor( iwork , threadIdx.x );
    }
  }

  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< Cuda > & work_config )
    : m_functor( functor ) ,
      m_config( work_config ) ,
      m_work( work_config.range )
  {
  }

  inline
  void execute() const
  {
    // To do:  query number of registers used by functor and adjust
    // nwarp accordingly to get maximum occupancy

    Cuda::size_type nwarp = 0;
    if (m_config.team > CudaTraits::WarpSize) {
      const Cuda::size_type warps_per_team =
        ( m_config.team + CudaTraits::WarpSize-1 ) / CudaTraits::WarpSize;
      nwarp = cuda_internal_maximum_warp_count() / warps_per_team;
    }
    else {
      const Cuda::size_type teams_per_warp =
        CudaTraits::WarpSize / m_config.team ;
      nwarp = cuda_internal_maximum_warp_count() * teams_per_warp;
    }
    const dim3 block( m_config.team , nwarp , 1 );

    Cuda::size_type nblock =
      std::min( (m_work + block.y - 1 ) / block.y ,
                cuda_internal_maximum_grid_count() );
    const dim3 grid( nblock , 1 , 1 );

    const Cuda::size_type shared = m_config.shared;
    CudaParallelLaunch< ParallelFor >( *this , grid , block , shared );
  }
};

#endif

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_ATOMIC_MP_VECTOR_HPP */
