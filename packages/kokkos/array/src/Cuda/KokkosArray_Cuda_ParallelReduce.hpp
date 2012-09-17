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

#ifndef KOKKOSARRAY_CUDA_PARALLELREDUCE_HPP
#define KOKKOSARRAY_CUDA_PARALLELREDUCE_HPP

#if defined( __CUDACC__ )

#include <KokkosArray_ParallelReduce.hpp>

#include <stdlib.h>
#include <iostream>

#include <vector>
#include <stdexcept>

#include <Cuda/KokkosArray_Cuda_Parallel.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
// Maximize shared memory and minimize L1 cache:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferShared );
// For 2.0 capability: 48 KB shared and 16 KB L1
//----------------------------------------------------------------------------
// Must have consistent '__shared__' statement across all device kernels.
// Since there may be more than one kernel in a file then have to make this
// a simple array of words.
//----------------------------------------------------------------------------

template< class ReduceTraits >
struct CudaReduceShared {
  typedef          Cuda::size_type          size_type ;
  typedef typename ReduceTraits::value_type value_type ;

  enum { WarpSize          = Impl::CudaTraits::WarpSize };
  enum { WarpStride        = WarpSize + 1 };
  enum { HalfWarpSize      = Impl::CudaTraits::WarpSize >> 1 };
  enum { WarpIndexMask     = Impl::CudaTraits::WarpIndexMask };
  enum { WarpIndexShift    = Impl::CudaTraits::WarpIndexShift };
  enum { SharedMemoryBanks = Impl::CudaTraits::SharedMemoryBanks_13 };

  enum { ValueWordCount = ( sizeof(value_type) + sizeof(size_type) - 1 )
                          / sizeof(size_type) };

  /** \brief  If the reduction value occupies an
   *          exact multiple of shared memory banks
   *          then it must be padded to avoid bank conflicts.
   */
  enum { ValueWordStride = ValueWordCount +
          ( ValueWordCount % SharedMemoryBanks ? 0 : 2 ) };

  enum { WordsPerWarp       = ValueWordStride * WarpSize };
  enum { WordsPerWarpStride = ValueWordStride * WarpStride };

  //--------------------------------------------------------------------------

  static inline
  size_type warp_count_max()
  {
    const size_type maximum_shared_words =
      cuda_internal_maximum_shared_words();

    // Start with maximum number of warps per block:
    size_type warps_per_block = cuda_internal_maximum_warp_count();

    // Reduce number of warps to fit per-thread reduction data in shared memory
    while ( maximum_shared_words < warps_per_block * WordsPerWarpStride ) {
      warps_per_block >>= 1 ;
    }

    return warps_per_block ;
  }

  static inline
  size_type shmem_size( const size_t warp_count )
  {
    return sizeof(size_type) *
           ( ValueWordStride * ( WarpStride * warp_count - 1 ) + 1 );
  }


private:

  size_type * m_scratch_flags ;
  size_type * m_scratch_space ;
  size_type   m_thread_count ;
  size_type   m_warp_count ;
  size_type   m_shared_flag_offset ;
  size_type   m_block_begin ;
  size_type   m_block_count ;
  size_type   m_group_init_count ;
  size_type   m_group_init_offset ;
  size_type   m_group_init_skip ;
  size_type   m_group_init_use ;
  size_type   m_group_init_use_end ;
  size_type   m_group_init_full_end ;


public:

  CudaReduceShared( const size_type warp_count ,
                    const size_type block_begin ,
                    const size_type block_count )
  {
    m_warp_count         = warp_count ;
    m_thread_count       = m_warp_count * CudaTraits::WarpSize ;
    m_shared_flag_offset = ValueWordStride * ( WarpStride * m_warp_count - 1 );
    m_block_begin        = block_begin ;
    m_block_count        = block_count ;

    m_group_init_count  = 1 ;
    m_group_init_offset = 1 ;

    for ( size_type g =  m_thread_count ; g < block_count ;
                    g *= m_thread_count ) {
      m_group_init_count  =  g ;
      m_group_init_offset += g ;
    }

    m_group_init_skip = ( m_group_init_count * m_thread_count - block_count ) /
                        ( m_thread_count - 1 );

    m_group_init_use      = m_group_init_count - m_group_init_skip ;
    m_group_init_use_end  = block_count - m_group_init_skip ;
    m_group_init_full_end = ( m_group_init_use - 1 ) * m_thread_count ;

    // Global memory scratch space for inter-block reduction:

    const size_type n =
      ( m_group_init_offset + block_count ) * ValueWordStride ;

    m_scratch_flags = cuda_internal_scratch_flags( m_group_init_offset );
    m_scratch_space = cuda_internal_scratch_space( n );

#if 0
std::cout
  << "CudaReduceShared( warp_count " << warp_count
  << " , block_begin " << block_begin
  << " , block_count " << block_count
  << " )"
  << std::endl
  << "  groups[ 0 .. " << m_group_init_use - 1
  << " ) <= blocks[ 0 .. " << m_group_init_full_end << " ) full reduction group"
  << std::endl
  << "  groups[ " << m_group_init_use - 1
  << " .. " << m_group_init_use
  << " ) <= blocks[ " << m_group_init_full_end
  << " .. " << m_group_init_use_end << " ) partial reduction group"
  << std::endl
  << "  groups[ " << m_group_init_use
  << " .. " << m_group_init_count
  << " ) <= blocks[ " << m_group_init_use_end
  << " .. " << m_block_count << " ) skip one reduction cycle"
  << std::endl ;
#endif

  }

private:

  //--------------------------------------------------------------------------

  __device__
  static inline
  size_type shared_data_offset( const size_type tx , const size_type ty  )
    { return ValueWordStride * ( tx + WarpStride * ty ); }

  //--------------------------------------------------------------------------
  // Reduce intra-block contributions with
  //   ty = thread's warp index
  //   tx = thread's index within warp ty

  __device__ 
  inline
  void reduce_intrablock( const size_type tx , const size_type ty ) const
  {
    typedef const volatile value_type * cvvp ;
    typedef       volatile value_type * vvp ;

    extern __shared__ size_type shared_data[];

    // Phase A: Reduce within my warp:
    //          Warp's reads occur before joins and writes
    //          so there is no race condition.
    //          Declare shared data to be volatile to
    //          prevent compiler from introducing a race condition.
    //
    if ( tx < HalfWarpSize ) {
      enum { n1  = ValueWordStride * 1 };
      enum { n2  = ValueWordStride * 2 };
      enum { n4  = ValueWordStride * 4 };
      enum { n8  = ValueWordStride * 8 };
      enum { n16 = ValueWordStride * 16 };

      size_type * const data = shared_data + shared_data_offset(tx,ty) ;

      ReduceTraits::join( *((vvp) data), *((cvvp)( data + n16 )) );
      ReduceTraits::join( *((vvp) data), *((cvvp)( data +  n8 )) );
      ReduceTraits::join( *((vvp) data), *((cvvp)( data +  n4 )) );
      ReduceTraits::join( *((vvp) data), *((cvvp)( data +  n2 )) );
      ReduceTraits::join( *((vvp) data), *((cvvp)( data +  n1 )) );
    }

    // Phase B: Use a single warp to reduce results from each warp.
    //          This requires: m_warp_count <= WarpSize
    //

    __syncthreads();

    if ( 0 == ty && tx + 1 < m_warp_count ) {
      enum { n1  = WarpStride * ValueWordStride * 1 };
      enum { n2  = WarpStride * ValueWordStride * 2 };
      enum { n4  = WarpStride * ValueWordStride * 4 };
      enum { n8  = WarpStride * ValueWordStride * 8 };
      enum { n16 = WarpStride * ValueWordStride * 16 };

      size_type * const data = shared_data + shared_data_offset(0,tx);

      if ( tx + 2 < m_warp_count ) {
        if ( tx + 4 < m_warp_count ) {
          if ( tx + 8 < m_warp_count ) {
            if ( tx + 16 < m_warp_count ) {
              ReduceTraits::join( *((vvp) data) , *((cvvp)( data + n16 )) );
            }
            ReduceTraits::join( *((vvp) data) , *((cvvp)( data + n8 )) );
          }
          ReduceTraits::join( *((vvp) data) , *((cvvp)( data + n4 )) );
        }
        ReduceTraits::join( *((vvp) data) , *((cvvp)( data + n2 )) );
      }
      ReduceTraits::join( *((vvp) data) , *((cvvp)( data + n1 )) );
    }
  }

public:

  //--------------------------------------------------------------------------

  __device__
  static inline
  value_type & value( const size_type thread_of_block )
  {
    extern __shared__ size_type shared_data[];

    const size_type tidx = thread_of_block &  CudaTraits::WarpIndexMask ;
    const size_type tidy = thread_of_block >> CudaTraits::WarpIndexShift ;

    return * ((value_type*)( shared_data + shared_data_offset(tidx,tidy) ) );
  }

  //--------------------------------------------------------------------------
  // REQUIRE number of threads is a multiple of warp size

  __device__
  inline size_type operator()( const size_type thread_of_block ,
                               const size_type block_of_grid ) const
  {
    extern __shared__ size_type shared_data[];

    const size_type tidx = thread_of_block &   CudaTraits::WarpIndexMask ;
    const size_type tidy = thread_of_block >>  CudaTraits::WarpIndexShift ;
    const size_type widy = thread_of_block &  ~CudaTraits::WarpIndexMask ;

    size_type block_count   = m_block_count ;
    size_type block_id      = m_block_begin + block_of_grid ;
    size_type group_count   = m_group_init_count ;
    size_type group_offset  = m_group_init_offset ;
    size_type group_size    = m_thread_count ;

    //------------------------------------
    // Initialize block assignment to reduction groups
    //   Reduction groups                    Blocks
    // A: [ 0      .. g_full ) <= blocks [ 0          .. g_full_end )
    // B: [ g_full .. g_use  ) <= blocks [ g_full_end .. g_use_end )
    // C: [ g_use  .. g_max  ) <= blocks [ g_use_end  .. block_count )
    //
    // A: These are reduction groups of 'size_max' blocks.
    // B: This is one reduction group of '< size_max' blocks.
    // C: These blocks skip the first reduction step.

    if ( m_group_init_use_end <= block_id ) {
      // C: This block is a member of a one-block group,
      //    skip the first reduction step.
      group_offset  -= group_count ;
      block_id      -= m_group_init_use_end - m_group_init_use ;
      block_count    = group_count ;
      group_count   /= m_thread_count ;
    }
    else if ( m_group_init_full_end <= block_id ) {
      // B: This block is  member of a partially filled group
      group_size  = m_group_init_use_end - m_group_init_full_end ;
    }

    //------------------------------------

    for (;;) {
      // Reduce this block's thread's contributions to a single value:
      reduce_intrablock( tidx , tidy );

      // If one or less blocks then done.
      // Return if this thread has the final reduction value.
      if ( block_count <= 1 ) return 0 == thread_of_block ;

      //----------------------------------
      // This block is a member of group of fan-in blocks.

      const size_type group_id    = block_id / m_thread_count ;
      const size_type flag_offset = group_offset - group_count ;

      //----------------------------------
      // Coalesced global memory write of this block's contribution
      {
        size_type * const scratch =
           m_scratch_space + ( group_offset + block_id ) * ValueWordStride ;

        for ( size_type i = thread_of_block ;
                        i < ValueWordStride ;
                        i += m_thread_count ) {
          scratch[i] = shared_data[i] ;
        }

        __threadfence(); // Wait for write to complete
      }

      //----------------------------------
      // Warp #0 Thread #0 :
      // Check if this is the last block to finish in the group:

      if ( 0 == thread_of_block ) {
        // Need a unique location for this flag which takes into
        // account the block-fanin iteration loop

        size_type * const flag = m_scratch_flags + flag_offset + group_id ;

        // Inform the entire block of the last_block status:
        // atomicInc returns value prior to increment.
        shared_data[ m_shared_flag_offset ] =
          ( group_size == 1 + atomicInc(flag,group_size+1) );

        // Reset the flag for the next reduction operation
        if ( shared_data[ m_shared_flag_offset ] ) *flag = 0 ;
      }

      // All threads of block wait for last_block flag to be set.
      __syncthreads();

      // If block is finished then it does not have the final value
      if ( ! shared_data[ m_shared_flag_offset ] ) return 0 ;

      // Last block to complete performs this group's reduction
      //----------------------------------
      // Each warp reads own data.
      // A warp's shared memory is contiguous but
      // there is a gap between warps to avoid bank conflicts.

      {
        size_type * const shared = shared_data + shared_data_offset(0,tidy);

        size_type * const scratch =
          m_scratch_space + ValueWordStride * (
            group_offset       + // offset for this reduction cycle
            group_id * m_thread_count + // offset for this block
            widy );              // offset for this warp

        for ( size_type i = tidx ;
                        i < WordsPerWarp ;
                        i += CudaTraits::WarpSize ) {
          shared[i] = scratch[i] ;
        }
      }

      // If the group was short data then junk was read into it's value.
      // Initialize the thread's data.

      if ( group_size <= thread_of_block ) {
        ReduceTraits::init(
          *((value_type*) ( shared_data + shared_data_offset(tidx,tidy)) ) );
      }

      //------------------------
      // Next iteration...

      block_id      =  group_id ;
      block_count   =  group_count ;
      group_size    =  m_thread_count ;
      group_count  /=  m_thread_count ;
      group_offset  =  flag_offset ;
    }
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ReduceTraits , class FinalizeType >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Cuda >
{
public:

  typedef          Cuda                      device_type ;
  typedef          Cuda::size_type           size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  typedef CudaReduceShared< ReduceTraits > ReduceType ;

  //----------------------------------------------------------------------

  const FunctorType    m_work_functor ;
  const FinalizeType   m_work_finalize ;
  const ReduceType     m_reduce_shared ;
  const size_type      m_work_count ;

  //----------------------------------------------------------------------

public:

  static
  size_type warp_count_max()
    { return ReduceType::warp_count_max(); }

  ParallelReduce( const FunctorType  & functor ,
                  const FinalizeType & finalize ,
                  const size_type      work_count ,
                  const size_type      warp_count ,
                  const size_type      global_block_begin ,
                  const size_type      global_block_count )
    : m_work_functor(  functor )
    , m_work_finalize( finalize )
    , m_reduce_shared( warp_count , global_block_begin , global_block_count )
    , m_work_count(    work_count )
  {}

  //--------------------------------------------------------------------------

  inline
  __device__
  void execute_on_device() const
  {
    value_type & value = ReduceType::value( threadIdx.x );

    ReduceTraits::init( value );

    const size_type work_stride = blockDim.x * gridDim.x ;

    // Reduce to per-thread contributions
    for ( Cuda::size_type iwork = threadIdx.x + blockDim.x * blockIdx.x ;
          iwork < m_work_count ; iwork += work_stride ) {
      m_work_functor( iwork , value );
    }

    if ( m_reduce_shared( threadIdx.x , blockIdx.x ) ) {
      // Last block thread #0 has final value
      m_work_finalize( value );
    }
  }

  static
  void execute( const size_t         work_count ,
                const FunctorType  & functor ,
                const FinalizeType & finalize )
  {
    const size_type warp_count = warp_count_max();
    const dim3 block( Impl::CudaTraits::WarpSize * warp_count , 1 , 1 );

    const dim3 grid(
      std::min( cuda_internal_maximum_grid_count() ,
                size_type( ( work_count + block.x - 1 ) / block.x ) ) , 1 , 1 );

    const size_type shmem_size = ReduceType::shmem_size( warp_count );

    ParallelReduce driver( functor, finalize, work_count, warp_count, 0, grid.x );

    CudaParallelLaunch< ParallelReduce >
      ::execute( driver , grid , block , shmem_size );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType >
class FunctorAssignment< ValueType , Cuda >
{
public :
  typedef ValueType value_type ;

  value_type * m_host ;
  value_type * m_dev ;

  FunctorAssignment( value_type & value )
    : m_host( & value )
    , m_dev( (value_type *)( cuda_internal_scratch_space(0) +
                             CudaTraits::WarpSize ) )
    {}

  template< typename ValueTypeDev >
  __device__
  void operator()( const ValueTypeDev & value ) const
  { *m_dev = value ; }

  ~FunctorAssignment()
  {
    CudaMemorySpace
      ::copy_to_host_from_device( m_host , m_dev , sizeof(value_type) );
  }
};

template< typename ValueType >
class FunctorAssignment< View< ValueType , Cuda , Cuda > , Cuda >
{
public :

  typedef ValueType value_type ;
  typedef View< value_type , Cuda , Cuda > view_type ;

  view_type m_view ;

  FunctorAssignment( const view_type & view )
    : m_view( view )
    {}

  template< typename ValueTypeDev >
  __device__
  void operator()( const ValueTypeDev & value ) const
  { *m_view = value ; }
};

#if 1
template< class FunctorType , class ReduceTraits , typename ValueType >
class ParallelReduce< FunctorType , ReduceTraits ,
                      View< ValueType , Cuda , Cuda > , Cuda >
{
public:

  typedef View< ValueType , Cuda , Cuda >  view_type ;

  struct FunctorAssignment
  {
    view_type m_view ;

    FunctorAssignment( const view_type & view ) : m_view( view ) {}

    template< typename ValueTypeDev >
    __device__
    void operator()( const ValueTypeDev & value ) const
      { *m_view = value ; }
  };

  static
  void execute( const size_t         work_count ,
                const FunctorType  & functor ,
                const view_type    & view )
  {
    ParallelReduce< FunctorType , ReduceTraits, FunctorAssignment , Cuda >
      ::execute( work_count , functor , FunctorAssignment( view ) );
  }
};
#endif

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template < class FunctorType , class ReduceTraits , class FinalizeType >
class CudaMultiFunctorParallelReduceMember ;

template < class ReduceTraits , class FinalizeType >
class CudaMultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType> {
public:

  typedef Cuda::size_type size_type ;

  const size_type  m_work_count ;
  const size_type  m_warp_count ;
  const size_type  m_block_count ;

private:

  static inline
  size_type block_count( const size_type work_count ,
                         const size_type warp_count )
  {
    const size_type nt = Impl::CudaTraits::WarpSize * warp_count ;

    return std::min( cuda_internal_maximum_grid_count() ,
                     ( work_count + nt - 1 ) / nt );
  }

protected:

  explicit
  CudaMultiFunctorParallelReduceMember( const size_type work_count ,
                                        const size_type warp_count )
    : m_work_count( work_count )
    , m_warp_count( warp_count )
    , m_block_count( block_count( work_count , warp_count ) )
    {} 

public:

  virtual ~CudaMultiFunctorParallelReduceMember() {}

  virtual
  void execute( const FinalizeType & finalize ,
                const size_type      global_block_offset ,
                const size_type      global_block_count ) const = 0 ;
};

template < class FunctorType , class ReduceTraits , class FinalizeType >
class CudaMultiFunctorParallelReduceMember :
  public CudaMultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType>
{
public:
  typedef CudaMultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType> base_type ;
  typedef ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Cuda >  driver_type ;
  
  typedef Impl::CudaReduceShared< ReduceTraits > ReduceType ;

  typedef Cuda            device_type ;
  typedef Cuda::size_type size_type ;

  FunctorType  m_functor ;

  CudaMultiFunctorParallelReduceMember(
    const FunctorType  & functor ,
    const size_type      work_count ,
    const size_type      warp_count )
  : base_type( work_count , warp_count )
  , m_functor( functor )
  {}

  virtual
  void execute( const FinalizeType & finalize ,
                const size_type      global_block_begin ,
                const size_type      global_block_count ) const
  {
    const dim3 block( Impl::CudaTraits::WarpSize * base_type::m_warp_count , 1 , 1 );

    const dim3 grid( base_type::m_block_count , 1 , 1 );

    const size_type shmem_size =
      ReduceType::shmem_size( base_type::m_warp_count );

    driver_type  driver( m_functor, finalize,
                         base_type::m_work_count ,
                         base_type::m_warp_count ,
                         global_block_begin ,
                         global_block_count );

    cuda_parallel_launch_local_memory< driver_type ><<< grid , block , shmem_size >>>( driver );
  }
};

} // namespace Impl

template < class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduce< ReduceTraits , FinalizeType , Cuda > {
public:
  typedef  Cuda               device_type ;
  typedef  Cuda::size_type    size_type ;

private:

  typedef Impl::CudaReduceShared< ReduceTraits > ReduceType ;

  typedef Impl::CudaMultiFunctorParallelReduceMember< void , ReduceTraits , FinalizeType > MemberType ;
  typedef std::vector< MemberType * > MemberVector ;

  MemberVector m_member_functors ;
  FinalizeType m_finalize ;

public:

  MultiFunctorParallelReduce( const FinalizeType & finalize )
   : m_member_functors()
   , m_finalize( finalize )
   {}

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & f )
  {
    MemberType * member =
      new Impl::CudaMultiFunctorParallelReduceMember< FunctorType , ReduceTraits , FinalizeType >
        ( f , work_count , ReduceType::warp_count_max() );

    m_member_functors.push_back( member );
  }

  void execute()
  {
    // Dispatch to streams, second and subsequent dispatches set the 
    // recycling flag so previous reduction results will be read.

    // Each dispatch has a fixed block count, dispatch the functor
    // on up to all streams until the requested block count is met.
    // When the stream_count * blocks_per_stream < blocks_requested
    // then the functor will iterate within the dispatch.

    typename MemberVector::iterator m ;

    size_type global_block_count = 0 ;

    for ( m = m_member_functors.begin() ; m != m_member_functors.end() ; ++m ) {
      global_block_count += (*m)->m_block_count ;
    }

    size_type global_block_offset = 0 ;

    for ( m = m_member_functors.begin() ; m != m_member_functors.end() ; ++m ) {
      MemberType & member = **m ;
      member.execute( m_finalize , global_block_offset , global_block_count );
      global_block_offset += member.m_block_count ;
    }
  }
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#endif /* defined( __CUDACC__ ) */

#endif /* KOKKOSARRAY_CUDA_PARALLELREDUCE_HPP */

