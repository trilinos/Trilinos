// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  template< class ExecSpace, class Tag = void >
struct MPVectorWorkConfig {

  typedef MPVectorWorkConfig execution_policy ;
  typedef ExecSpace          execution_space ;
  typedef Tag                work_tag ;

  execution_space space_;
  size_t range;
  size_t team;
  size_t shared;


  /*! \brief in the provided execution space instance */
  MPVectorWorkConfig( const execution_space &space,
                      const size_t range_,
                      const size_t team_,
                      const size_t shared_ = 0 ) :
    space_(space), range(range_), team(team_), shared(shared_) {}

  /*! \brief in the default execution space instance */
  MPVectorWorkConfig( const size_t range_,
                      const size_t team_,
                      const size_t shared_ = 0 ) :
    MPVectorWorkConfig(execution_space(), range_, team_, shared_) {}

  ExecSpace space() const { return space_; }
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
template< class FunctorType, class Tag >
class ParallelFor< FunctorType , MPVectorWorkConfig< Threads, Tag > > :
  public ParallelFor< FunctorType , Kokkos::RangePolicy< Tag, Threads > > {
  typedef Kokkos::RangePolicy< Tag, Threads > Policy ;
public:
  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< Threads, Tag > & work_config ) :
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
template< class FunctorType, class Tag >
class ParallelFor< FunctorType , MPVectorWorkConfig< OpenMP, Tag > > :
    public ParallelFor< FunctorType , Kokkos::RangePolicy< Tag, OpenMP > > {
  typedef Kokkos::RangePolicy< Tag, OpenMP > Policy ;
public:
  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< OpenMP, Tag > & work_config ) :
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
template< class FunctorType, class Tag >
class ParallelFor< FunctorType , MPVectorWorkConfig< Serial, Tag > > :
  public ParallelFor< FunctorType , Kokkos::RangePolicy< Tag, Serial > > {
  typedef Kokkos::RangePolicy< Tag, Serial > Policy ;
public:
  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< Serial, Tag > & work_config ) :
    ParallelFor< FunctorType , Policy >( functor ,
                                         Policy( 0, work_config.range ) ) {}
};
#endif // defined(KOKKOS_ENABLE_SERIAL)

#if defined( KOKKOS_ENABLE_CUDA ) && defined( __CUDACC__ )

// Specialization of ParallelFor<> for MPVectorWorkConfig on Cuda
// Here we use threadIdx.x for each entry in the specified team-size
template< class FunctorType, class Tag >
class ParallelFor< FunctorType , MPVectorWorkConfig< Cuda, Tag > > {
public:

  typedef Kokkos::RangePolicy< Tag, Cuda > Policy;

  const FunctorType m_functor ;
  const MPVectorWorkConfig< Cuda, Tag > m_config;
  const Cuda::size_type m_work ;
  const Policy m_policy;

  template <class TagType>
  inline __device__
  typename std::enable_if<std::is_same<TagType, void>::value>::type
  exec_range(const Cuda::size_type i, Cuda::size_type j) const {
    m_functor(i, j);
  }

  template <class TagType>
  inline __device__
  typename std::enable_if<!std::is_same<TagType, void>::value>::type
  exec_range(const Cuda::size_type i, Cuda::size_type j) const {
    m_functor(TagType(), i, j);
  }

  Policy const& get_policy() const { return m_policy; }

  inline
  __device__
  void operator()(void) const
  {
    const Cuda::size_type work_stride = blockDim.y * gridDim.x ;

    for ( Cuda::size_type iwork = threadIdx.y + blockDim.y * blockIdx.x ;
          iwork < m_work ;
          iwork += work_stride ) {
      this->template exec_range<Tag>(iwork, threadIdx.x);
    }
  }

  ParallelFor( const FunctorType        & functor ,
               const MPVectorWorkConfig< Cuda, Tag > & work_config )
    : m_functor( functor ) ,
      m_config( work_config ) ,
      m_work( work_config.range ),
      m_policy()
  {
  }

  inline
  void execute() const
  {
    // To do:  query number of registers used by functor and adjust
    // nwarp accordingly to get maximum occupancy

    auto const maxWarpCount = std::min<unsigned>(
        m_policy.space().cuda_device_prop().maxThreadsPerBlock / CudaTraits::WarpSize,
        CudaTraits::WarpSize);

    Cuda::size_type nwarp = 0;
    if (m_config.team > CudaTraits::WarpSize) {
      const Cuda::size_type warps_per_team =
        ( m_config.team + CudaTraits::WarpSize-1 ) / CudaTraits::WarpSize;
      nwarp = maxWarpCount / warps_per_team;
    }
    else {
      const Cuda::size_type teams_per_warp =
        CudaTraits::WarpSize / m_config.team ;
      nwarp = maxWarpCount * teams_per_warp;
    }
    const dim3 block( m_config.team , nwarp , 1 );

    const Cuda::size_type maxGridSizeX = m_policy.space().cuda_device_prop().maxGridSize[0];
    Cuda::size_type nblock =
      std::min( (m_work + block.y - 1 ) / block.y , maxGridSizeX );
    const dim3 grid( nblock , 1 , 1 );

    const Cuda::size_type shared = m_config.shared;
    CudaParallelLaunch< ParallelFor >( *this , grid , block , shared , m_policy.space().impl_internal_space_instance() );
  }
};

#endif

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_ATOMIC_MP_VECTOR_HPP */
