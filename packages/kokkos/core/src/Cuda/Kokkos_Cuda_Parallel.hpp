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

#include <iostream>
#include <stdio.h>

#if defined( __CUDACC__ )

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_CudaExec.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_Internal.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename Type >
struct CudaJoinFunctor {
  typedef Type value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    volatile const value_type & input )
    { update += input ; }
};

class CudaTeamMember {
private:

  typedef Kokkos::Cuda                           execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  void                * m_team_reduce ;
  scratch_memory_space  m_team_shared ;
  int                   m_league_rank ;
  int                   m_league_size ;

public:

#if defined( __CUDA_ARCH__ )

  __device__ inline
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared ; }

  __device__ inline int league_rank() const { return m_league_rank ; }
  __device__ inline int league_size() const { return m_league_size ; }
  __device__ inline int team_rank() const { return threadIdx.y ; }
  __device__ inline int team_size() const { return blockDim.y ; }

  __device__ inline void team_barrier() const { __syncthreads(); }

  template< class JoinOp >
  __device__ inline
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & op ) const
    {
      typename JoinOp::value_type * const base_data = (typename JoinOp::value_type *) m_team_reduce ;

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if ( 0 == threadIdx.y ) { base_data[0] = 0 ; }

      base_data[ threadIdx.y ] = value ;

      Impl::cuda_intra_block_reduce_scan<false>( op , base_data );

      return base_data[ blockDim.y - 1 ];
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  __device__ inline Type team_scan( const Type & value , Type * const global_accum ) const
    {
      Type * const base_data = (Type *) m_team_reduce ;

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if ( 0 == threadIdx.y ) { base_data[0] = 0 ; }

      base_data[ threadIdx.y + 1 ] = value ;

      Impl::cuda_intra_block_reduce_scan<true>( Impl::CudaJoinFunctor<Type>() , base_data + 1 );

      if ( global_accum ) {
        if ( blockDim.y == threadIdx.y + 1 ) {
          base_data[ blockDim.y ] = atomic_fetch_add( global_accum , base_data[ blockDim.y ] );
        }
        __syncthreads(); // Wait for atomic
        base_data[ threadIdx.y ] += base_data[ blockDim.y ] ;
      }

      return base_data[ threadIdx.y ];
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  __device__ inline Type team_scan( const Type & value ) const
    { return this->template team_scan<Type>( value , 0 ); }

  //----------------------------------------
  // Private for the driver

  __device__ inline
  CudaTeamMember( void * shared
                , const int shared_begin
                , const int shared_size
                , const int arg_league_rank
                , const int arg_league_size )
    : m_team_reduce( shared )
    , m_team_shared( ((char *)shared) + shared_begin , shared_size )
    , m_league_rank( arg_league_rank ) 
    , m_league_size( arg_league_size ) 
    {}

#else

  const execution_space::scratch_memory_space & team_shmem() const ;

  int league_rank() const ;
  int league_size() const ;
  int team_rank() const ;
  int team_size() const ;

  void team_barrier() const ;

  template< class JoinOp >
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & op ) const ;

  template< typename Type >
  Type team_scan( const Type & value , Type * const global_accum ) const ;

  template< typename Type >
  Type team_scan( const Type & value ) const ;

  //----------------------------------------
  // Private for the driver

  CudaTeamMember( void * shared
                , const int shared_begin
                , const int shared_end
                , const int arg_league_rank
                , const int arg_league_size );

#endif /* #if ! defined( __CUDA_ARCH__ ) */

};

template<unsigned VECTOR_LENGTH>
class CudaTeamVectorMember {
private:

  typedef Kokkos::Cuda                           execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  void                * m_team_reduce ;
  scratch_memory_space  m_team_shared ;
  int                   m_league_rank ;
  int                   m_league_size ;

public:

  enum { VectorLength = VECTOR_LENGTH };

#if defined( __CUDA_ARCH__ )

  __device__ inline
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared ; }

  __device__ inline int league_rank() const { return m_league_rank ; }
  __device__ inline int league_size() const { return m_league_size ; }
  __device__ inline int team_rank() const { return threadIdx.y ; }
  __device__ inline int team_size() const { return blockDim.y ; }

  __device__ inline void team_barrier() const { __syncthreads(); }

  template< class JoinOp >
  __device__ inline
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & op ) const
    {
      typename JoinOp::value_type * const base_data = (typename JoinOp::value_type *) m_team_reduce ;

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if(threadIdx.y == 0)
        base_data[ threadIdx.y ] = value ;

      Impl::cuda_intra_block_reduce_scan<false>( op , base_data );

      return base_data[ blockDim.y - 1 ];
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  __device__ inline Type team_scan( const Type & value , Type * const global_accum ) const
    {
      Type * const base_data = (Type *) m_team_reduce ;

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if ( 0 == threadIdx.y ) { base_data[0] = 0 ; }

      base_data[ threadIdx.y + 1 ] = value ;

      Impl::cuda_intra_block_reduce_scan<true>( Impl::CudaJoinFunctor<Type>() , base_data + 1 );

      if ( global_accum ) {
        if ( blockDim.y == threadIdx.y + 1 ) {
          base_data[ blockDim.y ] = atomic_fetch_add( global_accum , base_data[ blockDim.y ] );
        }
        __syncthreads(); // Wait for atomic
        base_data[ threadIdx.y ] += base_data[ blockDim.y ] ;
      }

      return base_data[ threadIdx.y ];
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  __device__ inline Type team_scan( const Type & value ) const
    { return this->template team_scan<Type>( value , 0 ); }

#ifdef KOKKOS_HAVE_CXX11
  template< typename iType, class Operation, typename ValueType >
  __device__ inline void team_par_for(const iType n, const Operation & op ) const {
    for(int _i = threadIdx.y; _i < n; _i += team_size())
      op(_i);
  }

  template< typename iType, class Operation, typename ValueType >
  __device__ inline void team_par_reduce(const iType n, const Operation & op, ValueType& result
  ) const {

    ValueType val = ValueType();

    for(int _i = threadIdx.y; _i < n; _i += team_size())
      op(_i,val);

    __shared__ ValueType sh_val[blockDim.y];

    if( threadIdx.x == 0 ) sh_val[threadIdx.y] = val;
    for(int _i = 1; _i < blockDim.y; _i*=2) {
      if( ( threadIdx.x == 0) && ( threadIdx.y + _i<blockDim.y ) ) {
        sh_val[threadIdx.y] += sh_val[threadIdx.y+_i];
      }
      __threadfence_block();
    }

    result = sh_val[0];
  }

  template< typename iType, class Operation, typename ValueType, class JoinType >
  __device__ inline void team_par_reduce(const iType n, const Operation & op, ValueType& result, const JoinType & join
  ) const {

    ValueType val = ValueType();

    for(int _i = threadIdx.y; _i < n; _i += team_size())
      op(_i,val);

    __shared__ ValueType sh_val[blockDim.y];

    if( threadIdx.x == 0 ) sh_val[threadIdx.y] = val;
    for(int _i = 1; _i < blockDim.y; _i*=2) {
      if( ( threadIdx.x == 0) && ( threadIdx.y + _i<blockDim.y ) ) {
        join(sh_val[threadIdx.y], sh_val[threadIdx.y+_i]);
      }
      __threadfence_block();
    }

    result = sh_val[0];
  }

  template< typename iType, class Operation, typename ValueType >
  __device__ inline void team_par_scan(const iType n, const Operation & op, ValueType& scan_val ) const {

    scan_val = ValueType();

    iType loop_bound = ((n+blockDim.y-1)/blockDim.y) * blockDim.y;

    for(int _i = threadIdx.y; _i < loop_bound; _i += blockDim.y) {
      ValueType val = ValueType();
      if(_i<n)
        op(_i , val , false);

      __shared__ ValueType sh_val[blockDim.y];

      if( threadIdx.x==0 ) sh_val[threadIdx.y] = val;
      for(int _i = 1; _i < blockDim.y; _i++) {
        if( (threadIdx.x==0) && (threadIdx.y - _i >= 0) ) {
          sh_val[threadIdx.y] += sh_val[threadIdx.y-_i];
        }
        __threadfence_block();
      }
      val = scan_val + sh_val[threadIdx.y] - val;
      scan_val += sh_val[blockDim.y-1];

      if(_i<n)
        op(_i , val , true);
    }
  }

  template< class Operation >
  __device__ inline void vector_single(const Operation & op) const {
    if(threadIdx.x == 0)
      op();
  }

  template< typename iType, class Operation, typename ValueType >
  __device__ inline void vector_par_for(const iType n, const Operation & op ) const {
    for(int _i = threadIdx.x; _i < n; _i += VectorLength)
      op(_i);
  }

#if (__CUDA_ARCH__ >= 300)
  template< typename iType, class Operation, typename ValueType >
  __device__ inline void vector_par_reduce(const iType n, const Operation & op, ValueType& result
      , typename Impl::enable_if< (sizeof(ValueType) & 3u) == 0 >::type * = 0
    ) const {

    ValueType val = ValueType();

    for(int _i = threadIdx.x; _i < n; _i += VectorLength)
      op(_i,val);

    result = val;

    if (VectorLength > 1)
      result += shfl_down(result, 1,VectorLength);
    if (VectorLength > 2)
      result += shfl_down(result, 2,VectorLength);
    if (VectorLength > 4)
      result += shfl_down(result, 4,VectorLength);
    if (VectorLength > 8)
      result += shfl_down(result, 8,VectorLength);
    if (VectorLength > 16)
      result += shfl_down(result, 16,VectorLength);

    result = shfl(result,0,VectorLength);
  }
#endif

#if (__CUDA_ARCH__ >= 300)
  template< typename iType, class Operation, typename ValueType , class JoinType>
  __device__ inline void vector_par_reduce(const iType n, const Operation & op, ValueType& result , const JoinType & join
      , typename Impl::enable_if< (sizeof(ValueType) & 3u) == 0 >::type * = 0
  ) const {

    ValueType val = result;

    for(int _i = threadIdx.x; _i < n; _i += VectorLength) {
      ValueType tmp;
      op(_i,tmp);
      join(val,tmp);
    }

    result = val;

    if (VectorLength > 1)
      join(result, shfl_down(result, 1,VectorLength));
    if (VectorLength > 2)
      join(result, shfl_down(result, 2,VectorLength));
    if (VectorLength > 4)
      join(result, shfl_down(result, 4,VectorLength));
    if (VectorLength > 8)
      join(result, shfl_down(result, 8,VectorLength));
    if (VectorLength > 16)
      join(result, shfl_down(result,16,VectorLength));

    result = shfl(result,0,VectorLength);
  }
#endif



#if (__CUDA_ARCH__ >= 300)
  template< typename iType, class Operation, typename ValueType >
  __device__ inline void vector_par_scan(const iType n, const Operation & op, ValueType& scan_val = ValueType()
      , typename Impl::enable_if< (sizeof(ValueType) & 3u) == 0 >::type * = 0
    ) const {

    scan_val = ValueType();
    iType loop_bound = ((n+VectorLength-1)/VectorLength) * VectorLength;
    for(int _i = threadIdx.x; _i < loop_bound; _i += VectorLength) {
      ValueType val = ValueType();
      if(_i<n)
        op(_i , val , false);

      ValueType tmp = val;
      ValueType result_i = ValueType();

      if(threadIdx.x%VectorLength == 0)
        result_i = tmp;
      if (VectorLength > 1)
        tmp += shfl_up(tmp, 1,VectorLength);
      if (threadIdx.x%VectorLength == 1)
        result_i = tmp;
      if (VectorLength > 2)
        tmp += shfl_up(tmp, 2,VectorLength);
      if ((threadIdx.x%VectorLength >= 2) &&
          (threadIdx.x%VectorLength < 4))
        result_i = tmp;
      if (VectorLength > 4)
        tmp += shfl_up(tmp, 4,VectorLength);
      if ((threadIdx.x%VectorLength >= 4) &&
          (threadIdx.x%VectorLength < 8))
        result_i = tmp;
      if (VectorLength > 8)
        tmp += shfl_up(tmp, 8,VectorLength);
      if ((threadIdx.x%VectorLength >= 8) &&
          (threadIdx.x%VectorLength < 16))
        result_i = tmp;
      if (VectorLength > 16)
        tmp += shfl_up(tmp, 16,VectorLength);
      if (threadIdx.x%VectorLength >= 16)
        result_i = tmp;

      val = scan_val + result_i - val;
      scan_val += shfl(tmp,VectorLength-1,VectorLength);
      if(_i<n)
        op(_i , val , true);
    }
  }
#endif

#endif
  //----------------------------------------
  // Private for the driver

  __device__ inline
  CudaTeamVectorMember( void * shared
                , const int shared_begin
                , const int shared_size
                , const int arg_league_rank
                , const int arg_league_size )
    : m_team_reduce( shared )
    , m_team_shared( ((char *)shared) + shared_begin , shared_size )
    , m_league_rank( arg_league_rank )
    , m_league_size( arg_league_size )
    {}

#else

  const execution_space::scratch_memory_space & team_shmem() const ;

  int league_rank() const ;
  int league_size() const ;
  int team_rank() const ;
  int team_size() const ;

  void team_barrier() const ;

  template< typename Type >
  Type team_scan( const Type & value , Type * const global_accum ) const ;

  template< typename Type >
  Type team_scan( const Type & value ) const ;

  template< typename iType, class Operation >
  void team_par_for(const iType n, const Operation & op) const {}

  template< typename iType, class Operation, typename ValueType >
  void team_par_reduce(const iType n, const Operation & op, ValueType& result) const {}

  template< typename iType, class Operation, typename ValueType , class JoinType>
  void team_par_reduce(const iType n, const Operation & op, ValueType& result, const JoinType & join) const {}

  template< typename iType, class Operation, typename ValueType >
  void team_par_scan(const iType n, const Operation & op, ValueType& result) const {}

  template< class Operation >
  void vector_single(const Operation & op) const {}

  template< typename iType, class Operation >
  void vector_par_for(const iType n, const Operation & op) const {}

  template< typename iType, class Operation, typename ValueType >
  void vector_par_reduce(const iType n, const Operation & op, ValueType& result) const {}

  template< typename iType, class Operation, typename ValueType , class JoinType>
  void vector_par_reduce(const iType n, const Operation & op, ValueType& result, const JoinType & join) const {}

  template< typename iType, class Operation, typename ValueType >
  void vector_par_scan(const iType n, const Operation & op, ValueType& result = ValueType()) const {}

  //----------------------------------------
  // Private for the driver

  CudaTeamVectorMember( void * shared
                , const int shared_begin
                , const int shared_end
                , const int arg_league_rank
                , const int arg_league_size );

#endif /* #if ! defined( __CUDA_ARCH__ ) */

};

} // namespace Impl
} // namespace Kokkos


namespace Kokkos {

template< class WorkArgTag >
class TeamPolicy< Kokkos::Cuda , WorkArgTag > {
private:

  enum { MAX_WARP = 8 };

  const int m_league_size ;
  const int m_team_size ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Cuda               execution_space ; ///< Execution space

  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size , int team_size_request )
    : m_league_size( league_size )
    , m_team_size( std::min( team_size_request , int( MAX_WARP * Impl::CudaTraits::WarpSize ) ) )
    { }

  TeamPolicy( int league_size , int team_size_request )
    : m_league_size( league_size )
    , m_team_size( std::min( team_size_request , int( MAX_WARP * Impl::CudaTraits::WarpSize ) ) )
    { }

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & functor )
    {
      int n = MAX_WARP * Impl::CudaTraits::WarpSize ;

      for ( ; n ; n >>= 1 ) {
        const int shmem_size =
          /* for global reduce */ Impl::cuda_single_inter_block_reduce_scan_shmem<false>( functor , n )
          /* for team   reduce */ + ( n + 2 ) * sizeof(double)
          /* for team   shared */ + Impl::FunctorTeamShmemSize< FunctorType >::value( functor , n );

        if ( shmem_size < Impl::CudaTraits::SharedMemoryCapacity ) break ;
      }

      return n ;
    }

  typedef Kokkos::Impl::CudaTeamMember member_type ;
};

template< unsigned VectorLength, class WorkArgTag >
class TeamVectorPolicy< VectorLength, Kokkos::Cuda , WorkArgTag > {
private:

  enum { MAX_WARP = 8 };

  const int m_league_size ;
  const int m_team_size ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Cuda               execution_space ; ///< Execution space

  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamVectorPolicy( execution_space & , int league_size , int team_size_request )
    : m_league_size( league_size )
    , m_team_size( std::min( team_size_request , int( MAX_WARP * Impl::CudaTraits::WarpSize ) ) )
    { }

  TeamVectorPolicy( int league_size , int team_size_request )
    : m_league_size( league_size )
    , m_team_size( std::min( team_size_request , int( MAX_WARP * Impl::CudaTraits::WarpSize ) ) )
    { }

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & functor )
    {
      int n = MAX_WARP * Impl::CudaTraits::WarpSize / VectorLength;

      for ( ; n ; n >>= 1 ) {
        const int shmem_size =
          /* for global reduce */ Impl::cuda_single_inter_block_reduce_scan_shmem<false>( functor , n )
          /* for team   reduce */ + ( n + 2 ) * sizeof(double)
          /* for team   shared */ + Impl::FunctorTeamShmemSize< FunctorType >::value( functor , n );

        if ( shmem_size < Impl::CudaTraits::SharedMemoryCapacity ) break ;
      }

      return n ;
    }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType & functor, int reduce_flag = 1 )
    {
      if(reduce_flag)
        return Impl::cuda_get_opt_block_size< Impl::ParallelReduce<FunctorType,TeamVectorPolicy,Cuda> >(functor);
      else
        return Impl::cuda_get_opt_block_size< Impl::ParallelFor<FunctorType,TeamVectorPolicy,Cuda> >(functor);
    }

  typedef Kokkos::Impl::CudaTeamVectorMember<VectorLength> member_type ;
};


} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 >
                 , Kokkos::Cuda >
{
public:

  typedef FunctorType                                              Functor ;
  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 > Policy ;

private:


  const FunctorType  m_functor ;
  const Policy       m_policy ;  

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             )
    { functor( iwork ); }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< ! Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             )
    { functor( Tag() , iwork ); }

public:

  inline
  __device__
  void operator()(void) const
    {
      const typename Policy::member_type work_stride = blockDim.y * gridDim.x ;
      const typename Policy::member_type work_end    = m_policy.end();

      for ( typename Policy::member_type
              iwork =  m_policy.begin() + threadIdx.y + blockDim.y * blockIdx.x ;
              iwork <  work_end ;
              iwork += work_stride ) {
        ParallelFor::template driver< typename Policy::work_tag >( m_functor, iwork );
      }
    }

  ParallelFor( const FunctorType  & functor ,
               const Policy       & policy )
    : m_functor( functor )
    , m_policy(  policy )
    {
      const dim3 block(  1 , CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1);
      const dim3 grid( std::min( ( int( policy.end() - policy.begin() ) + block.y - 1 ) / block.y
                               , cuda_internal_maximum_grid_count() )
                     , 1 , 1);

      CudaParallelLaunch< ParallelFor >( *this , grid , block , 0 );
    }
};

template< class FunctorType >
class ParallelFor< FunctorType , Kokkos::TeamPolicy< Kokkos::Cuda , void > , Kokkos::Cuda >
{
public:
  typedef Kokkos::TeamPolicy< Kokkos::Cuda , void >   Policy ;
  typedef FunctorType                                 Functor ;

  typedef typename Policy::member_type                team_member ;
  typedef Cuda::size_type                             size_type ;

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType m_functor ;
  size_type         m_shmem_begin ;
  size_type         m_shmem_size ;
  size_type         m_league_size ;

  __device__ inline
  void operator()(void) const
  {
    // Iterate this block through the league
    for ( int league_rank = blockIdx.x ; league_rank < m_league_size ; league_rank += gridDim.x ) {

      const team_member member( kokkos_impl_cuda_shared_memory<void>()
                              , m_shmem_begin
                              , m_shmem_size
                              , league_rank
                              , m_league_size );

      m_functor( member );
    }
  }


  ParallelFor( const FunctorType  & functor 
             , const Policy       & policy 
             )
  : m_functor( functor )
  , m_shmem_begin( sizeof(double) * ( policy.team_size() + 2 ) )
  , m_shmem_size( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
  , m_league_size( policy.league_size() )
  {
    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_shmem_begin + m_shmem_size ;

    if ( CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
    }

    const dim3 grid( std::min( int(policy.league_size()) , int(cuda_internal_maximum_grid_count()) ) , 1 , 1 );
    const dim3 block( 1 , policy.team_size() , 1 );

    CudaParallelLaunch< ParallelFor >( *this, grid, block, shmem_size_total ); // copy to device and execute
  }
};

template< unsigned VectorLength, class FunctorType>
class ParallelFor< FunctorType , Kokkos::TeamVectorPolicy<VectorLength,Kokkos::Cuda,void> ,Kokkos::Cuda >
{

public:
  typedef TeamVectorPolicy<VectorLength,Kokkos::Cuda,void>  Policy ;
  typedef FunctorType  Functor ;

  typedef typename Policy::member_type                team_member ;
  typedef Cuda::size_type                             size_type ;

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType m_functor ;
  size_type         m_shmem_begin ;
  size_type         m_shmem_size ;
  size_type         m_league_size ;

  __device__ inline
  void operator()(void) const
  {
    // Iterate this block through the league
    for ( int league_rank = blockIdx.x ; league_rank < m_league_size ; league_rank += gridDim.x ) {

      const team_member member( kokkos_impl_cuda_shared_memory<void>()
                              , m_shmem_begin
                              , m_shmem_size
                              , league_rank
                              , m_league_size );

      m_functor( member );
    }
  }


  ParallelFor( const FunctorType  & functor
             , const Policy       & policy
             )
  : m_functor( functor )
  , m_shmem_begin( sizeof(double) * ( policy.team_size() + 2 ) )
  , m_shmem_size( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
  , m_league_size( policy.league_size() )
  {
    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_shmem_begin + m_shmem_size ;

    if ( CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
    }

    const dim3 grid( std::min( int(policy.league_size()) , int(cuda_internal_maximum_grid_count()) ) , 1 , 1 );
    const dim3 block( VectorLength, policy.team_size() , 1 );

    CudaParallelLaunch< ParallelFor >( *this, grid, block, shmem_size_total ); // copy to device and execute
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType 
                    , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 >
                    , Kokkos::Cuda
                    >
{
public:
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef FunctorType                         Functor ;
  typedef Kokkos::RangePolicy<Arg0,Arg1,Arg2> Policy ;
  typedef Cuda::size_type                     size_type ;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y == blockDim.z == 1

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type *       m_unified_space ;

  // Determine block size constrained by shared memory:
  static inline
  unsigned local_block_size( const FunctorType & f )
    {
      unsigned n = CudaTraits::WarpSize * 8 ;
      while ( n && CudaTraits::SharedMemoryCapacity < cuda_single_inter_block_reduce_scan_shmem<false>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value )
    { functor( iwork , value ); }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< ! Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value )
    { functor( Tag() , iwork , value ); }

  __device__ inline
  void operator()(void) const
  {
    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    {
      reference_type value =
        Reduce::init( m_functor , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

      // Number of blocks is bounded so that the reduction can be limited to two passes.
      // Each thread block is given an approximately equal amount of work to perform.
      // Accumulate the values for this block.
      // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

      const Policy range( m_policy , blockIdx.x , gridDim.x );

      for ( typename Policy::member_type iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
            iwork < iwork_end ; iwork += blockDim.y ) {
        ParallelReduce::template driver< typename Policy::work_tag >( m_functor , iwork , value );
      }
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false>(
           m_functor , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.y == 0 ) { Reduce::final( m_functor , shared ); }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor 
                , const Policy       & policy 
                , const HostViewType & result
                )
  : m_functor( functor )
  , m_policy(  policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  {
    const int block_size  = local_block_size( functor );
    const int block_count = std::min( int(block_size)
                                    , ( int(policy.end() - policy.begin()) + block_size - 1 ) / block_size
                                    );

    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * block_count );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
    m_unified_space = cuda_internal_scratch_unified( Reduce::value_size( functor ) );

    const dim3 grid( block_count , 1 , 1 );
    const dim3 block( 1 , block_size , 1 );
    const int shmem = cuda_single_inter_block_reduce_scan_shmem<false>( m_functor , block.y );

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem ); // copy to device and execute

    Cuda::fence();

    if ( result.ptr_on_device() ) {
      if ( m_unified_space ) {
        const int count = Reduce::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { result.ptr_on_device()[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = Reduce::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( result.ptr_on_device() , m_scratch_space , size );
      }
    }
  }
};

template< class FunctorType >
class ParallelReduce< FunctorType , Kokkos::TeamPolicy< Kokkos::Cuda , void > , Kokkos::Cuda >
{
public:
  typedef Kokkos::TeamPolicy< Kokkos::Cuda , void >   Policy ;
  typedef FunctorType                                 Functor ;
  typedef typename Policy::member_type                team_member ;
  typedef ReduceAdapter< FunctorType >                Reduce ;
  typedef typename Reduce::pointer_type               pointer_type ;
  typedef typename Reduce::reference_type             reference_type ;
  typedef Cuda::size_type                             size_type ;

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ global reduce space ]
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type *       m_unified_space ;
  size_type         m_team_begin ;
  size_type         m_shmem_begin ;
  size_type         m_shmem_size ;
  size_type         m_league_size ;


  __device__ inline
  void operator()(void) const
  {
    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    reference_type value =
      Reduce::init( m_functor , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

    // Iterate this block through the league
    for ( int league_rank = blockIdx.x ; league_rank < m_league_size ; league_rank += gridDim.x ) {

      const team_member member( kokkos_impl_cuda_shared_memory<char>() + m_team_begin
                              , m_shmem_begin
                              , m_shmem_size
                              , league_rank
                              , m_league_size );

      m_functor( member , value );
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false>(
           m_functor , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.y == 0 ) { Reduce::final( m_functor , shared ); }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }


  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor 
                , const Policy       & policy 
                , const HostViewType & result
                )
  : m_functor( functor )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_team_begin( cuda_single_inter_block_reduce_scan_shmem<false>( functor , policy.team_size() ) )
  , m_shmem_begin( sizeof(double) * ( policy.team_size() + 2 ) )
  , m_shmem_size( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
  , m_league_size( policy.league_size() )
  {
    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;
    const int not_power_of_two = 0 != ( policy.team_size() & ( policy.team_size() - 1 ) );

    if ( not_power_of_two ||  CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    const int block_count = std::min( policy.league_size() , policy.team_size() );

    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * block_count );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
    m_unified_space = cuda_internal_scratch_unified( Reduce::value_size( functor ) );

    const dim3 grid( block_count , 1 , 1 );
    const dim3 block( 1 , policy.team_size() , 1 );

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem_size_total ); // copy to device and execute

    Cuda::fence();

    if ( result.ptr_on_device() ) {
      if ( m_unified_space ) {
        const int count = Reduce::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { result.ptr_on_device()[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = Reduce::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( result.ptr_on_device() , m_scratch_space , size );
      }
    }
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 >
                  , Kokkos::Cuda
                  >
{
public:
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef FunctorType                         Functor ;
  typedef Kokkos::RangePolicy<Arg0,Arg1,Arg2> Policy ;
  typedef Cuda::size_type                     size_type ;

  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  // Determine block size constrained by shared memory:
  static inline
  unsigned local_block_size( const FunctorType & f )
    {
      // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
      // gridDim.x <= blockDim.y * blockDim.y
      //
      // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

      unsigned n = CudaTraits::WarpSize * 4 ;
      while ( n && CudaTraits::SharedMemoryCapacity < cuda_single_inter_block_reduce_scan_shmem<false>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
        size_type   m_final ;
  
  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value 
             , const bool     final )
    { functor( iwork , value , final ); }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< ! Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value
             , const bool     final )
    { functor( Tag() , iwork , value , final ); }

  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    size_type * const shared_value = kokkos_impl_cuda_shared_memory<size_type>() + word_count.value * threadIdx.y ;

    Reduce::init( m_functor , shared_value );

    // Number of blocks is bounded so that the reduction can be limited to two passes.
    // Each thread block is given an approximately equal amount of work to perform.
    // Accumulate the values for this block.
    // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

    const Policy range( m_policy , blockIdx.x , gridDim.x );

    for ( typename Policy::member_type iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
          iwork < iwork_end ; iwork += blockDim.y ) {
      ParallelScan::template driver< typename Policy::work_tag >
        ( m_functor , iwork , Reduce::reference( shared_value ) , false );
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups' totals.
    // Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i < gridDim.x
    cuda_single_inter_block_reduce_scan<true>( m_functor , blockIdx.x , gridDim.x , kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags );
  }

  //----------------------------------------

  __device__ inline
  void final(void) const
  {
    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_data   = kokkos_impl_cuda_shared_memory<size_type>();
    size_type * const shared_prefix = shared_data + word_count.value * threadIdx.y ;
    size_type * const shared_accum  = shared_data + word_count.value * ( blockDim.y + 1 );

    // Starting value for this thread block is the previous block's total.
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.y ) {
      Reduce::init( m_functor , shared_accum );
    }

    const Policy range( m_policy , blockIdx.x , gridDim.x );

    for ( typename Policy::member_type iwork_base = range.begin(); iwork_base < range.end() ; iwork_base += blockDim.y ) {

      const typename Policy::member_type iwork = iwork_base + threadIdx.y ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      Reduce::init( m_functor , shared_prefix + word_count.value );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i] ;
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < range.end() ) {
        ParallelScan::template driver< typename Policy::work_tag >
          ( m_functor , iwork , Reduce::reference( shared_prefix + word_count.value ) , false );
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true>( m_functor , Reduce::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.y ;
        for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

      // Call functor with exclusive scan value
      if ( iwork < range.end() ) {
        ParallelScan::template driver< typename Policy::work_tag >
          ( m_functor , iwork , Reduce::reference( shared_prefix ) , true );
      }
    }
  }

  //----------------------------------------

  __device__ inline
  void operator()(void) const
  {
    if ( ! m_final ) {
      initial();
    }
    else {
      final();
    }
  }

  ParallelScan( const FunctorType  & functor ,
                const Policy       & policy )
  : m_functor( functor )
  , m_policy( policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_final( false )
  {
    enum { GridMaxComputeCapability_2x = 0x0ffff };

    const int block_size = local_block_size( functor );

    const int grid_max = ( block_size * block_size ) < GridMaxComputeCapability_2x ?
                         ( block_size * block_size ) : GridMaxComputeCapability_2x ;

    // At most 'max_grid' blocks:
    const int nwork    = policy.end() - policy.begin();
    const int max_grid = std::min( int(grid_max) , int(( nwork + block_size - 1 ) / block_size ));

    // How much work per block:
    const int work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

    // How many block are really needed for this much work:
    const dim3 grid( ( nwork + work_per_block - 1 ) / work_per_block , 1 , 1 );
    const dim3 block( 1 , block_size , 1 );
    const int shmem = Reduce::value_size( functor ) * ( block_size + 2 );

    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * grid.x );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) * 1 );

    m_final = false ;
    CudaParallelLaunch< ParallelScan >( *this, grid, block, shmem ); // copy to device and execute

    m_final = true ;
    CudaParallelLaunch< ParallelScan >( *this, grid, block, shmem ); // copy to device and execute
  }

  void wait() const { Cuda::fence(); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */

