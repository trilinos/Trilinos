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

#ifndef KOKKOSARRAY_CUDA_PARALLELREDUCE_HPP
#define KOKKOSARRAY_CUDA_PARALLELREDUCE_HPP

#if defined( __CUDACC__ )

#include <stdlib.h>
#include <iostream>

#include <vector>
#include <stdexcept>

#include <KokkosArray_ParallelReduce.hpp>
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

template< class ValueType >
struct CudaReduceShared_AnalyzeSize
{
  enum { WordSize          = sizeof(Cuda::size_type) };
  enum { ValueSize         = sizeof(ValueType) };
  enum { WarpSize          = Impl::CudaTraits::WarpSize };
  enum { WarpStride        = WarpSize + 1 };
  enum { SharedMemoryBanks = Impl::CudaTraits::SharedMemoryBanks_13 };
  enum { SharedMemorySize  = 0x04000 /* assume only 16k on Capability 1.3 */ }; 
  enum { SharedMemoryWords = SharedMemorySize / WordSize };

  enum { ValueWordCount = ( ValueSize + WordSize - 1 ) / WordSize };

  enum { ValueWordStride =
           ValueWordCount + ( ValueWordCount % SharedMemoryBanks ? 0 : 2 ) };

  enum { WordsPerWarpStride = ValueWordStride * WarpStride };

  /* Maximum of 8 ( 1 << 3 ) warps per block.
   * For simple dot products using 16 warps is slower than 8 warps.
   */
  enum { WarpCountShift =
           ( WordsPerWarpStride << 3 ) < SharedMemoryWords ? 3 : (
           ( WordsPerWarpStride << 2 ) < SharedMemoryWords ? 2 : (
           ( WordsPerWarpStride << 1 ) < SharedMemoryWords ? 1 : (
           ( WordsPerWarpStride      ) < SharedMemoryWords ? 0 : 256 ))) };

  enum { WarpCount = WarpCountShift <= 4 ? ( 1 << WarpCountShift ) : 0 };
};

template< class MemberType >
struct CudaReduceShared_AnalyzeSize< MemberType[] >
{ enum { WarpCount = 0 }; };

template< class MemberType >
struct CudaReduceShared_AnalyzeSize< MemberType* >
{ enum { WarpCount = 0 }; };

//----------------------------------------------------------------------------

template< class ValueOper ,
          class FinalizeFunctor ,
          unsigned WarpCountAnalysis =
            CudaReduceShared_AnalyzeSize<
              typename ValueOper::value_type >::WarpCount >
struct CudaReduceShared ;


template< class ValueOper , class FinalizeType , unsigned WarpCountAnalysis >
struct CudaReduceShared
{
  typedef
    CudaReduceShared_AnalyzeSize< typename ValueOper::value_type >
      StaticSizeAnalysis ;

  enum { WarpSize          = Impl::CudaTraits::WarpSize };
  enum { HalfWarpSize      = Impl::CudaTraits::WarpSize >> 1 };
  enum { WarpIndexMask     = Impl::CudaTraits::WarpIndexMask };
  enum { WarpIndexShift    = Impl::CudaTraits::WarpIndexShift };

  enum { WarpStride        = StaticSizeAnalysis::WarpStride };
  enum { ValueWordCount    = StaticSizeAnalysis::ValueWordCount };
  enum { ValueWordStride   = StaticSizeAnalysis::ValueWordStride };
  enum { WarpCountShift    = StaticSizeAnalysis::WarpCountShift };
  enum { WarpCount         = StaticSizeAnalysis::WarpCount };

  enum { ThreadCountShift  = WarpCountShift + WarpIndexShift };
  enum { ThreadCount       = 1 << ThreadCountShift };
  enum { SharedFlagOffset  = ValueWordStride * ( WarpStride * WarpCount - 1 ) };

  typedef Cuda::size_type size_type ;

  typedef ReduceOperator< ValueOper , FinalizeType >  reduce_operator ;

  //--------------------------------------------------------------------------

  inline static
  size_type value_size( const FinalizeType & )
    { return StaticSizeAnalysis::ValueSize ; }

  /*  How many words required to hold the size of value */
  static inline
  size_type word_count( const size_type ) { return ValueWordCount ; }

  /*  If the reduction value occupies an
   *  exact multiple of shared memory banks
   *  then it must be padded to avoid bank conflicts.
   */
  static inline
  size_type word_stride( const size_type value_size )
  { return ValueWordStride ; }

  static inline
  size_type warp_count_max( const FinalizeType & finalize )
  { return WarpCount ; }

  inline size_type shmem_size() const
  {
    return sizeof(size_type) *
           ( ValueWordCount * ( WarpStride * WarpCount - 1 ) + 1 );
  }

private:

  const reduce_operator  m_reduce ;

  const size_type   m_block_begin ;
  const size_type   m_block_count ;

  size_type   m_group_init_count ;
  size_type   m_group_init_offset ;
  size_type   m_group_init_skip ;
  size_type   m_group_init_use ;
  size_type   m_group_init_use_end ;
  size_type   m_group_init_full_end ;
  size_type * m_scratch_flags ;
  size_type * m_scratch_space ;

public:

  CudaReduceShared( const FinalizeType & finalize ,
                    const size_type block_begin ,
                    const size_type block_count )
  : m_reduce( finalize )
  , m_block_begin( block_begin )
  , m_block_count( block_count )
  {
    m_group_init_count  = 1 ;
    m_group_init_offset = 1 ;

    for ( size_type g =  ThreadCount ; g < block_count ;
                    g *= ThreadCount ) {
      m_group_init_count  =  g ;
      m_group_init_offset += g ;
    }

    m_group_init_skip = ( m_group_init_count * ThreadCount - block_count ) /
                        ( ThreadCount - 1 );

    m_group_init_use      = m_group_init_count - m_group_init_skip ;
    m_group_init_use_end  = block_count - m_group_init_skip ;
    m_group_init_full_end = ( m_group_init_use - 1 ) * ThreadCount ;

    // Global memory scratch space for inter-block reduction:

    const size_type n =
      ( m_group_init_offset + block_count ) * ValueWordStride ;

    m_scratch_flags = cuda_internal_scratch_flags( m_group_init_offset * sizeof(size_type) );
    m_scratch_space = cuda_internal_scratch_space( n * sizeof(size_type) );

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
  inline
  size_type shared_data_offset( const size_type tx , const size_type ty  ) const
    { return ValueWordStride * ( tx + WarpStride * ty ); }

  //--------------------------------------------------------------------------
  // Reduce intra-block contributions with
  //   ty = thread's warp index
  //   tx = thread's index within warp ty

  __device__ 
  inline
  void reduce_intrablock( const size_type tx , const size_type ty ) const
  {
    extern __shared__ size_type shared_data[];

    // Phase A: Reduce within my warp:
    //          Warp's reads occur before joins and writes
    //          so there is no race condition.
    //          Declare shared data to be volatile to
    //          prevent compiler from introducing a race condition.
    //
    if ( tx < HalfWarpSize ) {
      enum { n16 = ValueWordStride << 4 };
      enum { n08 = ValueWordStride << 3 };
      enum { n04 = ValueWordStride << 2 };
      enum { n02 = ValueWordStride << 1 };
      enum { n01 = ValueWordStride      };

      size_type * const data = shared_data + shared_data_offset(tx,ty) ;

      m_reduce.join( data, data + n16 );
      m_reduce.join( data, data + n08 );
      m_reduce.join( data, data + n04 );
      m_reduce.join( data, data + n02 );
      m_reduce.join( data, data + n01 );
    }

    // Phase B: Use a single warp to reduce results from each warp.
    //          This requires: m_warp_count <= WarpSize
    //

    __syncthreads();


    if ( 0 == ty && tx + 1 < WarpCount ) {
      enum { ValueWarpStride = WarpStride * ValueWordStride };
      enum { n16 = ValueWarpStride << 4 };
      enum { n08 = ValueWarpStride << 3 };
      enum { n04 = ValueWarpStride << 2 };
      enum { n02 = ValueWarpStride << 1 };
      enum { n01 = ValueWarpStride      };

      size_type * const data = shared_data + shared_data_offset(0,tx);

      if ( tx + 2 < WarpCount ) {
        if ( tx + 4 < WarpCount ) {
          if ( tx + 8 < WarpCount ) {

// Hard-wired WarpCount <= 16 so don't need this step
//
//          if ( tx + 16 < WarpCount ) {
//            m_reduce.join( data , data + n16 ); 
//          } 

            m_reduce.join( data , data + n08 );
          }
          m_reduce.join( data , data + n04 );
        }
        m_reduce.join( data , data + n02 );
      }
      m_reduce.join( data , data + n01 );
    }
  }

public:

  //--------------------------------------------------------------------------

  typedef typename reduce_operator::reference_type reference_type ;

  __device__
  inline
  reference_type value( const size_type thread_of_block ) const
  {
    extern __shared__ size_type shared_data[];

    const size_type tidx = thread_of_block &  CudaTraits::WarpIndexMask ;
    const size_type tidy = thread_of_block >> CudaTraits::WarpIndexShift ;

    return m_reduce.init( shared_data + shared_data_offset(tidx,tidy) );
  }

  //--------------------------------------------------------------------------
  // REQUIRE number of threads is a multiple of warp size

  __device__
  void operator()( const size_type thread_of_block ,
                   const size_type block_of_grid ) const
  {
    extern __shared__ size_type shared_data[];

    // tidx == which thread within the warp
    // tidy == which warp   within the block
    // widy == first thread of the warp within the block

    const size_type tidx = thread_of_block &   CudaTraits::WarpIndexMask ;
    const size_type tidy = thread_of_block >>  CudaTraits::WarpIndexShift ;
    const size_type widy = thread_of_block &  ~CudaTraits::WarpIndexMask ;

    size_type block_count   = m_block_count ;
    size_type block_id      = m_block_begin + block_of_grid ;
    size_type group_count   = m_group_init_count ;
    size_type group_offset  = m_group_init_offset ;
    size_type group_size    = ThreadCount ;

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
      group_count  >>= ThreadCountShift ;
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
      // Thread 0 has the final reduction value
      if ( block_count <= 1 ) {
        if ( 0 == thread_of_block ) {
          m_reduce.finalize( shared_data );
        }
        break ;
      }

      //----------------------------------
      // This block is a member of group of fan-in blocks.

      const size_type group_id    = block_id >> ThreadCountShift ;
      const size_type flag_offset = group_offset - group_count ;

      //----------------------------------
      // Coalesced global memory write of this block's contribution
      {
        size_type * const scratch =
           m_scratch_space + ( group_offset + block_id ) * ValueWordStride ;

        for ( size_type i = thread_of_block ;
                        i < ValueWordStride ;
                        i += ThreadCount ) {
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
        shared_data[ SharedFlagOffset ] =
          ( group_size == 1 + atomicInc(flag,group_size+1) );

        // Reset the flag for the next reduction operation
        if ( shared_data[ SharedFlagOffset ] ) *flag = 0 ;
      }

      // All threads of block wait for last_block flag to be set.
      __syncthreads();

      // If this is not the last block in the group then
      // this block is finished.
      if ( ! shared_data[ SharedFlagOffset ] ) break ;

      // Last block to complete performs this group's reduction
      //----------------------------------
      // Each warp reads its own data.
      // A warp's shared memory is contiguous but
      // there is a gap between warps to avoid bank conflicts.

      {
        enum { WordsPerWarp = ValueWordStride << CudaTraits::WarpIndexShift };

        size_type * const shared = shared_data + shared_data_offset(0,tidy);

        size_type * const scratch =
          m_scratch_space + ValueWordStride * (
            group_offset           + // offset for this reduction cycle
            ( group_id << ThreadCountShift ) + // offset for this block
            widy );                  // offset for this warp

        for ( size_type i = tidx ;
                        i < WordsPerWarp ;
                        i += CudaTraits::WarpSize ) {
          shared[i] = scratch[i] ;
        }
      }

      // If the group was short data then junk was read into it's value.
      // Initialize the thread's data.

      if ( group_size <= thread_of_block ) {
        m_reduce.init( shared_data + shared_data_offset(tidx,tidy) );
      }

      //------------------------
      // Next iteration...

      block_id       =  group_id ;
      block_count    =  group_count ;
      group_size     =  ThreadCount ;
      group_count  >>=  ThreadCountShift ;
      group_offset   =  flag_offset ;
    }

    return ;
  }
};

//----------------------------------------------------------------------------

template< class ValueOper , class FinalizeType >
struct CudaReduceShared<ValueOper,FinalizeType, 0 /* Cannot determine size */ >
{
  enum { WarpSize          = Impl::CudaTraits::WarpSize };
  enum { WarpStride        = WarpSize + 1 };
  enum { HalfWarpSize      = Impl::CudaTraits::WarpSize >> 1 };
  enum { WarpIndexMask     = Impl::CudaTraits::WarpIndexMask };
  enum { WarpIndexShift    = Impl::CudaTraits::WarpIndexShift };
  enum { SharedMemoryBanks = Impl::CudaTraits::SharedMemoryBanks_13 };

  typedef Cuda::size_type size_type ;

  typedef ReduceOperator< ValueOper , FinalizeType >  reduce_operator ;

  //--------------------------------------------------------------------------
  inline static
  size_type value_size( const FinalizeType & f )
  { return reduce_operator::value_size( f ); }

  /*  How many words required to hold the size of value */
  static inline
  size_type word_count( const size_type value_size )
  { return ( value_size + sizeof(size_type) - 1 ) / sizeof(size_type); }

  /*  If the reduction value occupies an
   *  exact multiple of shared memory banks
   *  then it must be padded to avoid bank conflicts.
   */
  static inline
  size_type word_stride( const size_type value_size )
  {
    const size_type n = word_count( value_size );
    return n + ( n % SharedMemoryBanks ? 0 : 2 );
  }

  static inline
  size_type warp_count_max( const FinalizeType & finalize )
  {
    const size_type value_size = reduce_operator::value_size( finalize );
    const size_type value_word_stride = word_stride( value_size );

    const size_type maximum_shared_words =
      cuda_internal_maximum_shared_words();

    const size_type words_per_warp_stride = value_word_stride * WarpStride ;

    // Start with maximum number of warps per block:
    size_type warps_per_block = cuda_internal_maximum_warp_count();

    // Reduce number of warps to fit per-thread reduction data in shared memory
    while ( maximum_shared_words < warps_per_block * words_per_warp_stride ) {
      warps_per_block >>= 1 ;
    }

    return warps_per_block ;
  }

  inline size_type shmem_size() const
  {
    return sizeof(size_type) *
           ( m_value_word_stride * ( WarpStride * m_warp_count - 1 ) + 1 );
  }

private:

  const reduce_operator  m_reduce ;

  const size_type   m_warp_count ;
  const size_type   m_thread_count ;
  const size_type   m_value_word_count ;
  const size_type   m_value_word_stride ;
  const size_type   m_shared_flag_offset ;
  const size_type   m_block_begin ;
  const size_type   m_block_count ;

  size_type   m_group_init_count ;
  size_type   m_group_init_offset ;
  size_type   m_group_init_skip ;
  size_type   m_group_init_use ;
  size_type   m_group_init_use_end ;
  size_type   m_group_init_full_end ;
  size_type * m_scratch_flags ;
  size_type * m_scratch_space ;

public:

  CudaReduceShared( const FinalizeType & finalize ,
                    const size_type block_begin ,
                    const size_type block_count )
  : m_reduce( finalize )
  , m_warp_count( warp_count_max( finalize ) )
  , m_thread_count( m_warp_count * CudaTraits::WarpSize )
  , m_value_word_count( word_count( m_reduce.value_size() ) )
  , m_value_word_stride( word_stride( m_reduce.value_size() ) )
  , m_shared_flag_offset( m_value_word_stride * ( WarpStride * m_warp_count - 1 ) )
  , m_block_begin( block_begin )
  , m_block_count( block_count )
  {
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
      ( m_group_init_offset + block_count ) * m_value_word_stride ;

    m_scratch_flags = cuda_internal_scratch_flags( m_group_init_offset * sizeof(size_type) );
    m_scratch_space = cuda_internal_scratch_space( n * sizeof(size_type) );

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
  inline
  size_type shared_data_offset( const size_type tx , const size_type ty  ) const
    { return m_value_word_stride * ( tx + WarpStride * ty ); }

  //--------------------------------------------------------------------------
  // Reduce intra-block contributions with
  //   ty = thread's warp index
  //   tx = thread's index within warp ty

  __device__ 
  inline
  void reduce_intrablock( const size_type tx , const size_type ty ) const
  {
    extern __shared__ size_type shared_data[];

    // Phase A: Reduce within my warp:
    //          Warp's reads occur before joins and writes
    //          so there is no race condition.
    //          Declare shared data to be volatile to
    //          prevent compiler from introducing a race condition.
    //
    if ( tx < HalfWarpSize ) {

      size_type * const data = shared_data + shared_data_offset(tx,ty) ;

      m_reduce.join( data, data + ( m_value_word_stride << 4 ) );
      m_reduce.join( data, data + ( m_value_word_stride << 3 ) );
      m_reduce.join( data, data + ( m_value_word_stride << 2 ) );
      m_reduce.join( data, data + ( m_value_word_stride << 1 ) );
      m_reduce.join( data, data + ( m_value_word_stride      ) );
    }

    // Phase B: Use a single warp to reduce results from each warp.
    //          This requires: m_warp_count <= WarpSize
    //

    __syncthreads();

    if ( 0 == ty && tx + 1 < m_warp_count ) {

      size_type * const data = shared_data + shared_data_offset(0,tx);

      const size_type s = WarpStride * m_value_word_stride ;

      if ( tx + 2 < m_warp_count ) {
        if ( tx + 4 < m_warp_count ) {
          if ( tx + 8 < m_warp_count ) {
            if ( tx + 16 < m_warp_count ) {
              m_reduce.join( data , data + ( s << 4 ) );
            }
            m_reduce.join( data , data + ( s << 3 ) );
          }
          m_reduce.join( data , data + ( s << 2 ) );
        }
        m_reduce.join( data , data + ( s << 1 ) );
      }
      m_reduce.join( data , data + s );
    }
  }

public:

  //--------------------------------------------------------------------------

  typedef typename reduce_operator::reference_type reference_type ;

  __device__
  inline
  reference_type value( const size_type thread_of_block ) const
  {
    extern __shared__ size_type shared_data[];

    const size_type tidx = thread_of_block &  CudaTraits::WarpIndexMask ;
    const size_type tidy = thread_of_block >> CudaTraits::WarpIndexShift ;

    return m_reduce.init( shared_data + shared_data_offset(tidx,tidy) );
  }

  //--------------------------------------------------------------------------
  // REQUIRE number of threads is a multiple of warp size

  __device__
  void operator()( const size_type thread_of_block ,
                   const size_type block_of_grid ) const
  {
    extern __shared__ size_type shared_data[];

    // tidx == which thread within the warp
    // tidy == which warp   within the block
    // widy == first thread of the warp within the block

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
      // Thread 0 has the final reduction value
      if ( block_count <= 1 ) {
        if ( 0 == thread_of_block ) {
          m_reduce.finalize( shared_data );
        }
        break ;
      }

      //----------------------------------
      // This block is a member of group of fan-in blocks.

      const size_type group_id    = block_id / m_thread_count ;
      const size_type flag_offset = group_offset - group_count ;

      //----------------------------------
      // Coalesced global memory write of this block's contribution
      {
        size_type * const scratch =
           m_scratch_space + ( group_offset + block_id ) * m_value_word_stride ;

        for ( size_type i = thread_of_block ;
                        i < m_value_word_stride ;
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

      // If this is not the last block in the group then
      // this block is finished.
      if ( ! shared_data[ m_shared_flag_offset ] ) break ;

      // Last block to complete performs this group's reduction
      //----------------------------------
      // Each warp reads own data.
      // A warp's shared memory is contiguous but
      // there is a gap between warps to avoid bank conflicts.

      {
        size_type * const shared = shared_data + shared_data_offset(0,tidy);

        size_type * const scratch =
          m_scratch_space + m_value_word_stride * (
            group_offset              + // offset for this reduction cycle
            group_id * m_thread_count + // offset for this block
            widy );                     // offset for this warp

        const size_type words_per_warp =
          m_value_word_stride << CudaTraits::WarpIndexShift ;

        for ( size_type i = tidx ;
                        i < words_per_warp ;
                        i += CudaTraits::WarpSize ) {
          shared[i] = scratch[i] ;
        }
      }

      // If the group was short data then junk was read into it's value.
      // Initialize the thread's data.

      if ( group_size <= thread_of_block ) {
        m_reduce.init( shared_data + shared_data_offset(tidx,tidy) );
      }

      //------------------------
      // Next iteration...

      block_id      =  group_id ;
      block_count   =  group_count ;
      group_size    =  m_thread_count ;
      group_count  /=  m_thread_count ;
      group_offset  =  flag_offset ;
    }

    return ;
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ValueOper , class FinalizeType >
class ParallelReduce< FunctorType , ValueOper , FinalizeType , Cuda >
{
public:

  typedef          Cuda                   device_type ;
  typedef          Cuda::size_type        size_type ;
  typedef typename ValueOper::value_type  value_type ;

  typedef CudaReduceShared< ValueOper , FinalizeType > ReduceType ;

  //----------------------------------------------------------------------

  const FunctorType    m_work_functor ;
  const ReduceType     m_reduce_shared ;
  const size_type      m_work_count ;

  //----------------------------------------------------------------------

  static
  size_type block_count( const FinalizeType & f ,
                         const size_type work_count )
  {
    const size_type nt =
      Impl::CudaTraits::WarpSize * ReduceType::warp_count_max( f );

    return std::min( nt , size_type( ( work_count + nt - 1 ) / nt ) );
  }

public:

  // For the multi-functor reduce:
  ParallelReduce( const FunctorType  & functor ,
                  const FinalizeType & finalize ,
                  const size_type      work_count ,
                  const size_type      global_block_begin ,
                  const size_type      global_block_count )
    : m_work_functor(  functor )
    , m_reduce_shared( finalize ,
                       global_block_begin ,
                       global_block_count )
    , m_work_count(    work_count )
  {}

  //--------------------------------------------------------------------------

  ParallelReduce( const size_type      work_count ,
                  const FunctorType  & functor ,
                  const FinalizeType & finalize )
    : m_work_functor(  functor )
    , m_reduce_shared( finalize , 0 , block_count(finalize,work_count) )
    , m_work_count(    work_count )
  {
    const size_type nw = ReduceType::warp_count_max( finalize );

    const dim3 block( Impl::CudaTraits::WarpSize * nw , 1, 1 );
    const dim3 grid( block_count(finalize,work_count) , 1 , 1 );
    const size_type shmem_size = m_reduce_shared.shmem_size();

    CudaParallelLaunch< ParallelReduce >( *this , grid , block , shmem_size );
  }

  //--------------------------------------------------------------------------

  inline
  __device__
  void operator()(void) const
  {
    typename ReduceType::reference_type value =
      m_reduce_shared.value( threadIdx.x );

    const size_type work_stride = blockDim.x * gridDim.x ;

    // Reduce to per-thread contributions
    for ( Cuda::size_type iwork = threadIdx.x + blockDim.x * blockIdx.x ;
          iwork < m_work_count ; iwork += work_stride ) {
      m_work_functor( iwork , value );
    }

    m_reduce_shared( threadIdx.x , blockIdx.x );
  }
};

/* Reduction into a view on the host */

template< class FunctorType , class ValueOper ,
          class ValueType , class LayoutType >
class ParallelReduce< FunctorType ,
                      ValueOper , 
                      View< ValueType , LayoutType , Host > ,
                      Cuda >
{
public:

  typedef typename FunctorType::value_type value_type ;

  typedef View< ValueType , LayoutType , Host > host_view_type ;

  ParallelReduce( const Cuda::size_type  work_count ,
                  const FunctorType    & functor ,
                  const host_view_type & host_view )
  {
    typedef typename
      StaticAssertSame< typename host_view_type::scalar_type[] , value_type >
        ::type ok_match_type ;

    typedef typename
      StaticAssert< host_view_type::Rank == 1 >::type ok_rank ;

    if ( functor.value_count != host_view.dimension_0() ) {
      std::ostringstream msg ;
      msg << "KokkosArray::parallel_reduce( <array_type> ) ERROR "
          << "given incompatible array lengths: functor.value_count("
          << functor.value_count
          << ") != view.dimension_0("
          << host_view.dimension_0() << ")" ;
      throw std::runtime_error(msg.str());
    }

    typedef Impl::ParallelReduceFunctorValue< value_type , Cuda >
      FinalizeType ;

    const FinalizeType finalize( functor.value_count );

    ParallelReduce< FunctorType , ValueOper , FinalizeType , Cuda >
      ( work_count , functor , finalize );

    finalize.result( host_view.ptr_on_device() );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename T >
struct CudaReduceResult {

  inline static
  T * device_pointer( const Cuda::size_type count )
  {
    T * ptr = (T*) cuda_internal_scratch_unified( sizeof(T) * count );
    if ( 0 == ptr ) {
      ptr = (T*) cuda_internal_scratch_space( sizeof(T) * count );
    }
    return ptr ;
  }

  inline static
  void copy_to_host( T host_pointer[] , const Cuda::size_type count )
  {
    Cuda::fence();
    T * ptr = (T*) cuda_internal_scratch_unified( sizeof(T) * count );
    if ( 0 != ptr ) {
      for ( Cuda::size_type i = 0 ; i < count ; ++i )
        host_pointer[i] = ptr[i] ;
    }
    else {
      ptr = (T*) cuda_internal_scratch_space( sizeof(T) * count );
      CudaMemorySpace
        ::copy_to_host_from_device( host_pointer , ptr , sizeof(T) * count );
    }
  }

  inline static
  T return_to_host()
  {
    Cuda::fence();
    T * ptr = (T*) cuda_internal_scratch_unified( sizeof(T) );
    if ( 0 != ptr ) {
      return *ptr ;
    }
    else {
      ptr = (T*) cuda_internal_scratch_space( sizeof(T) );
      T value ;
      CudaMemorySpace::copy_to_host_from_device( & value , ptr , sizeof(T) );
      return value ;
    }
  }
};

template< typename ValueType >
class ParallelReduceFunctorValue< ValueType , Cuda >
{
private:
  ValueType * const m_dev ; 
public:

  typedef ValueType value_type ;

  __device__
  inline
  void operator()( const value_type & value ) const
    { *m_dev = value ; }

  inline
  ParallelReduceFunctorValue()
    : m_dev( CudaReduceResult<value_type>::device_pointer(1) ) {}

  inline
  value_type result() const
  { return CudaReduceResult<value_type>::return_to_host(); }
};

template< typename MemberType >
class ParallelReduceFunctorValue< MemberType[] , Cuda >
{
private:
  MemberType * const m_dev ; 
public:

  typedef MemberType     value_type[] ;
  const Cuda::size_type  value_count ;

  __device__
  inline
  void operator()( const MemberType input[] ) const
    {
      for ( Cuda::size_type i = 0 ; i < value_count ; ++i )
        m_dev[i] = input[i] ;
    }

  explicit
  ParallelReduceFunctorValue( Cuda::size_type n )
    : m_dev( CudaReduceResult<MemberType>::device_pointer(n) )
    , value_count( n )
    {}

  void result( MemberType result[] ) const
  { CudaReduceResult<MemberType>::copy_to_host( result , value_count ); }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template < class FunctorType , class ValueOper , class FinalizeType >
class CudaMultiFunctorParallelReduceMember ;

template < class ValueOper , class FinalizeType >
class CudaMultiFunctorParallelReduceMember<void,ValueOper,FinalizeType> {
public:

  typedef Cuda::size_type size_type ;

  const size_type  m_work_count ;

protected:

  explicit
  CudaMultiFunctorParallelReduceMember( const size_type work_count )
    : m_work_count( work_count )
    {} 

public:

  virtual ~CudaMultiFunctorParallelReduceMember() {}

  virtual
  void apply( const FinalizeType & finalize ,
              const size_type      thread_count ,
              const size_type      block_count ,
              const size_type      global_block_offset ,
              const size_type      global_block_count ) const = 0 ;
};

template < class FunctorType , class ValueOper , class FinalizeType >
class CudaMultiFunctorParallelReduceMember :
  public CudaMultiFunctorParallelReduceMember<void,ValueOper,FinalizeType>
{
public:
  typedef CudaMultiFunctorParallelReduceMember<void,ValueOper,FinalizeType> base_type ;
  typedef ParallelReduce< FunctorType , ValueOper , FinalizeType , Cuda >  driver_type ;
  
  typedef Impl::CudaReduceShared< ValueOper , FinalizeType > ReduceType ;

  typedef Cuda            device_type ;
  typedef Cuda::size_type size_type ;

  FunctorType  m_functor ;

  CudaMultiFunctorParallelReduceMember(
    const FunctorType  & functor ,
    const size_type      work_count )
  : base_type( work_count )
  , m_functor( functor )
  {}

  virtual
  void apply( const FinalizeType & finalize ,
              const size_type      thread_count ,
              const size_type      block_count ,
              const size_type      global_block_begin ,
              const size_type      global_block_count ) const
  {
    const dim3 block( thread_count , 1 , 1 );
    const dim3 grid(  block_count ,  1 , 1 );

    driver_type  driver( m_functor, finalize,
                         base_type::m_work_count ,
                         global_block_begin ,
                         global_block_count );

    const size_type shmem_size = driver.m_reduce_shared.shmem_size();

    cuda_parallel_launch_local_memory< driver_type ><<< grid , block , shmem_size >>>( driver );
  }
};

} // namespace Impl

template < class ValueOper , class FinalizeType >
class MultiFunctorParallelReduce< ValueOper , FinalizeType , Cuda > {
public:
  typedef  Cuda               device_type ;
  typedef  Cuda::size_type    size_type ;

private:

  typedef Impl::CudaReduceShared< ValueOper , FinalizeType > ReduceType ;

  typedef Impl::CudaMultiFunctorParallelReduceMember< void , ValueOper , FinalizeType > MemberType ;
  typedef std::vector< MemberType * > MemberVector ;

  MemberVector m_member_functors ;
  FinalizeType m_finalize ;

  inline static 
  size_type block_count( const size_type warp_count ,
                         const size_type work_count )
  {
    const size_type nt = Impl::CudaTraits::WarpSize * warp_count ;

    return std::min( Impl::cuda_internal_maximum_grid_count() ,
                     ( work_count + nt - 1 ) / nt );
  }

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
      new Impl::CudaMultiFunctorParallelReduceMember< FunctorType , ValueOper , FinalizeType >
        ( f , work_count );

    m_member_functors.push_back( member );
  }

  void execute()
  {
    const size_type warp_count = ReduceType::warp_count_max( m_finalize );

    const size_type thread_count = warp_count * Impl::CudaTraits::WarpSize ;

    // Dispatch to streams, second and subsequent dispatches set the 
    // recycling flag so previous reduction results will be read.

    // Each dispatch has a fixed block count, dispatch the functor
    // on up to all streams until the requested block count is met.
    // When the stream_count * blocks_per_stream < blocks_requested
    // then the functor will iterate within the dispatch.

    typename MemberVector::iterator m ;

    size_type global_block_count = 0 ;

    for ( m = m_member_functors.begin() ; m != m_member_functors.end() ; ++m ) {
      global_block_count += block_count( warp_count , (*m)->m_work_count );
    }

    size_type global_block_offset = 0 ;

    for ( m = m_member_functors.begin() ; m != m_member_functors.end() ; ++m ) {
      MemberType & member = **m ;

      const size_type n = block_count( warp_count , (*m)->m_work_count );

      member.apply( m_finalize , thread_count , n ,
                    global_block_offset , global_block_count );

      global_block_offset += n ;
    }
  }
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#endif /* defined( __CUDACC__ ) */

#endif /* KOKKOSARRAY_CUDA_PARALLELREDUCE_HPP */

