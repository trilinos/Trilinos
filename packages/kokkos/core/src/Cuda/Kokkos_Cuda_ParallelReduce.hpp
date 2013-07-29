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

#ifndef KOKKOS_CUDA_PARALLELREDUCE_HPP
#define KOKKOS_CUDA_PARALLELREDUCE_HPP

#if defined( __CUDACC__ )

#include <stdlib.h>
#include <iostream>

#include <vector>

#include <Kokkos_ParallelReduce.hpp>
#include <Cuda/Kokkos_Cuda_Parallel.hpp>
#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
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
//----------------------------------------------------------------------------

template< class ValueType >
struct CudaReduceSharedSizes {

  typedef Cuda::size_type size_type ;

private:

  enum { WarpStride         = CudaTraits::WarpSize + 1 };
  enum { WordSize           = sizeof(Cuda::size_type) };
  enum { ValueSize          = sizeof(ValueType) };
  enum { ValueWordCount     = ( ValueSize + WordSize - 1 ) / WordSize };
  enum { ValueWordStride =
           ValueWordCount + ( ValueWordCount % CudaTraits::SharedMemoryBanks ? 0 : 2 ) };

  enum { WordsPerWarpStride = ValueWordStride * WarpStride };
  enum { SharedMemoryWords  = CudaTraits::SharedMemoryUsage / WordSize };

  /* Maximum of 8 ( 1 << 3 ) warps per block.
   * For simple dot products using 16 warps is slower than 8 warps.
   */
  enum { WarpCount =
           ( WordsPerWarpStride << 3 ) < SharedMemoryWords ? ( 1 << 3 ) : (
           ( WordsPerWarpStride << 2 ) < SharedMemoryWords ? ( 1 << 2 ) : (
           ( WordsPerWarpStride << 1 ) < SharedMemoryWords ? ( 1 << 1 ) : (
           ( WordsPerWarpStride      ) < SharedMemoryWords ? 1 : 0 ))) };

public:

  enum { value_word_count   = ValueWordCount };
  enum { value_word_stride  = ValueWordStride };
  enum { value_warp_stride  = WordsPerWarpStride };
  enum { warp_count         = WarpCount };
  enum { thread_count       = warp_count * CudaTraits::WarpSize };
  enum { shmem_size         = warp_count * WordSize * value_warp_stride };

  explicit inline
  CudaReduceSharedSizes( const size_type )
  { }

  __device__ inline
  size_type shared_data_offset( const size_type tx ,
                                const size_type ty  ) const
    { return value_word_stride * ( tx + WarpStride * ty ); }
};


template< class MemberType >
struct CudaReduceSharedSizes< MemberType[] > {

  typedef Cuda::size_type size_type ;

private:

  enum { WarpStride         = CudaTraits::WarpSize + 1 };
  enum { WordSize           = sizeof(size_type) };
  enum { SharedMemoryWords  = CudaTraits::SharedMemoryUsage / WordSize };

public:

  const size_type  value_word_count ;
  const size_type  value_word_stride ;
  const size_type  value_warp_stride ;
  const size_type  warp_count ;
  const size_type  thread_count ;
  const size_type  shmem_size ;

  explicit inline
  CudaReduceSharedSizes( const size_type value_size )
    : value_word_count(  ( value_size + WordSize - 1 ) / WordSize )
      /*  If the reduction value occupies an
       *  exact multiple of shared memory banks
       *  then it must be padded to avoid bank conflicts.
       */
    , value_word_stride( value_word_count ? value_word_count + ( value_word_count % CudaTraits::SharedMemoryBanks ? 0 : 2 ) : 0 )
    , value_warp_stride( value_word_stride * WarpStride )
    , warp_count( ( value_warp_stride == 0 ) ? 0 : (
                  ( value_warp_stride << 3 ) < SharedMemoryWords ? (1 << 3) : (
                  ( value_warp_stride << 2 ) < SharedMemoryWords ? (1 << 2) : (
                  ( value_warp_stride << 1 ) < SharedMemoryWords ? (1 << 1) : (
                  ( value_warp_stride      ) < SharedMemoryWords ? 1 : 0 )))))
    , thread_count( warp_count * CudaTraits::WarpSize )
    , shmem_size( warp_count * WordSize * value_warp_stride )
    {
      if ( value_size && ! warp_count ) {
        // Not enough memory for reduction
      }
    }

  __device__ inline
  size_type shared_data_offset( const size_type tx ,
                                const size_type ty  ) const
    { return value_word_stride * ( tx + WarpStride * ty ); }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class ValueOper , class FinalizeType >
struct CudaReduceShared
{
public:

  typedef CudaReduceSharedSizes< typename ValueOper::value_type > reduce_sizes ;
  typedef ReduceAdapter< ValueOper >  reduce_operator ;

  //--------------------------------------------------------------------------

  typedef Cuda::size_type size_type ;

  static inline
  size_type warp_count( const size_type arg_value_size )
  {
    const reduce_sizes data( arg_value_size );
    return data.warp_count ;
  }

  inline size_type warp_count() const
  { return m_data.warp_count ; }

  inline size_type shmem_size() const
  { return m_data.shmem_size ; }

  const reduce_operator  m_reduce ;
  const reduce_sizes     m_data ;

private:

  size_type   m_block_begin ;
  size_type   m_block_count ;
  size_type   m_group_init_count ;
  size_type   m_group_init_offset ;
  size_type   m_group_init_use ;
  size_type   m_group_init_use_end ;
  size_type   m_group_init_full_end ;
  size_type * m_scratch_flags ;
  size_type * m_scratch_space ;
  size_type * m_output_space ;

public:

  void assign_block_range( const size_type block_begin ,
                           const size_type block_count )
  {
    m_block_begin = block_begin ;
    m_block_count = block_count ;
    m_group_init_count  = 1 ;
    m_group_init_offset = 1 ;

    for ( size_type g =  m_data.thread_count ; g < block_count ;
                    g *= m_data.thread_count ) {
      m_group_init_count  =  g ;
      m_group_init_offset += g ;
    }

    const size_type group_init_skip =
      ( m_group_init_count * m_data.thread_count - block_count ) /
      ( m_data.thread_count - 1 );

    m_group_init_use      = m_group_init_count - group_init_skip ;
    m_group_init_use_end  = block_count - group_init_skip ;
    m_group_init_full_end = m_data.thread_count *
                            ( m_group_init_use ? m_group_init_use - 1 : 0 );

    // Global memory scratch space for inter-block reduction:

    const size_type n =
      ( m_group_init_offset + block_count ) * m_data.value_word_stride ;

    m_scratch_flags = cuda_internal_scratch_flags( m_group_init_offset * sizeof(size_type) );
    m_scratch_space = cuda_internal_scratch_space( n * sizeof(size_type) );
    m_output_space  = cuda_internal_scratch_unified( m_reduce.value_size() );

    if ( 0 == m_output_space ) { m_output_space = m_scratch_space ; }

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

  CudaReduceShared( const ValueOper & op )
  : m_reduce( op )
  , m_data( m_reduce.value_size() )
  {}

private:

  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  Reduce constributions within a warp.
  //  Warp's reads occur before joins and writes so there is no race condition.
  //  Declare shared data to be volatile to the prevent compiler
  //  from optimizing away read of updated data.
  //    thread's index within warp = tx < WarpSize
  __device__ inline
  void reduce_intra_warp( volatile size_type * const thread_data ,
                          const size_type tx ) const
  {
    enum { HalfWarpSize = 1 << 4 };

    if ( tx < HalfWarpSize ) {
      m_reduce.join( thread_data ,
                     thread_data + ( m_data.value_word_stride << 4 ) );

      m_reduce.join( thread_data ,
                     thread_data + ( m_data.value_word_stride << 3 ) );

      m_reduce.join( thread_data ,
                     thread_data + ( m_data.value_word_stride << 2 ) );

      m_reduce.join( thread_data ,
                     thread_data + ( m_data.value_word_stride << 1 ) );

      m_reduce.join( thread_data ,
                     thread_data + ( m_data.value_word_stride      ) );
    }
  }

  //--------------------------------------------------------------------------
  // Use a single warp to reduce results from each warp.
  // This requires: m_data.warp_count <= WarpSize
  // Only warp #0 should have 'warp_data'.

  __device__ 
  inline
  void reduce_inter_warp( volatile size_type * const warp_data ,
                          const size_type tx ) const
  {
    __syncthreads(); // Wait for all warps

    if ( warp_data && tx + 1 < m_data.warp_count ) {
      if ( tx + 2 < m_data.warp_count ) {
        if ( tx + 4 < m_data.warp_count ) {
          if ( tx + 8 < m_data.warp_count ) {

// Hard-wired m_data.warp_count <= 16 so don't need this step
//
//          if ( tx + 16 < m_data.warp_count ) {
//            m_reduce.join( warp_data ,
//                           warp_data + ( m_data.value_warp_stride << 4 ) ); 
//            __threadfence_block();
//          }

            m_reduce.join( warp_data ,
                           warp_data + ( m_data.value_warp_stride << 3 ) );
            __threadfence_block();
          }
          m_reduce.join( warp_data ,
                         warp_data + ( m_data.value_warp_stride << 2 ) );
          __threadfence_block();
        }
        m_reduce.join( warp_data ,
                       warp_data + ( m_data.value_warp_stride << 1 ) );
        __threadfence_block();
      }
      m_reduce.join( warp_data ,
                     warp_data + ( m_data.value_warp_stride ) );
    }
  }

  //--------------------------------------------------------------------------
  // Coalesced global memory write of this block's contribution
  // which has been computed by Warp #0.

  __device__ inline
  void write_block_value( const size_type block ,
                          const size_type thread_of_block ,
                          const size_type * const shared_data ) const
  {
    // If ( CudaTraits::WarpSize < ValueWordStride )
    // then multiple warps will so must synchronize threads.

    if ( CudaTraits::WarpSize < m_data.value_word_stride ) { __syncthreads(); }

    volatile size_type * const scratch =
      m_scratch_space + m_data.value_word_stride * block ;

    for ( size_type i = thread_of_block ;
                    i < m_data.value_word_stride ;
                    i += m_data.thread_count ) {
      scratch[i] = shared_data[i] ;
    }

    __threadfence(); // Wait for write to complete.
  }

  //--------------------------------------------------------------------------

public:

  //--------------------------------------------------------------------------

  typedef typename reduce_operator::reference_type reference_type ;

  __device__
  inline
  reference_type value( const size_type thread_of_block ) const
  {
    extern __shared__ size_type shared_data[];

    void * const data = shared_data +
                        m_data.shared_data_offset(
                          /* tidx */ thread_of_block &  CudaTraits::WarpIndexMask ,
                          /* tidy */ thread_of_block >> CudaTraits::WarpIndexShift );

    m_reduce.init( data );

    return m_reduce.reference( data );
  }

  //--------------------------------------------------------------------------

  __device__
  void operator()( const size_type thread_of_block ,
                   const size_type block_of_grid ) const
  {
    extern __shared__ size_type shared_data[];

    // tidx == which thread within the warp
    // tidy == which warp   within the block
    // wbeg == first thread of the warp within the block

    const size_type tidx = thread_of_block &   CudaTraits::WarpIndexMask ;
    const size_type tidy = thread_of_block >>  CudaTraits::WarpIndexShift ;
    const size_type wbeg = thread_of_block &  ~CudaTraits::WarpIndexMask ;

    size_type * const thread_data =
      shared_data + m_data.shared_data_offset(tidx,tidy);

    size_type * const inter_warp_data =
      tidy ? 0 : shared_data + m_data.shared_data_offset(0,tidx);

    size_type block_id      = m_block_begin + block_of_grid ;
    size_type group_count   = m_group_init_count ;
    size_type group_offset  = m_group_init_offset ;
    size_type group_size    = m_data.thread_count ;

    //------------------------------------
    // Initialize block assignment to reduction groups
    //   Reduction groups                    Blocks
    // A: [ 0      .. g_full ) <= blocks [ 0          .. g_full_end )
    // B: [ g_full .. g_use  ) <= blocks [ g_full_end .. g_use_end )
    // C: [ g_use  .. g_max  ) <= blocks [ g_use_end  .. m_block_count )
    //
    // A: These are reduction groups of 'size_max' blocks.
    // B: This is one reduction group of '< size_max' blocks.
    // C: These are reduction groups of 'size_max' blocks that
    //    skip the first reduction step.

    if ( m_group_init_use_end <= block_id ) {
      // C: This block is a member of a one-block group,
      //    skip the first reduction step.
      group_offset  -= group_count ;
      block_id      -= m_group_init_use_end - m_group_init_use ;
      group_count   /= m_data.thread_count ;
    }
    else if ( m_group_init_full_end <= block_id ) {
      // B: This block is  member of a partially filled group
      group_size = m_group_init_use_end - m_group_init_full_end ;
    }

    //------------------------------------

    for (;;) {

      reduce_intra_warp( thread_data , tidx );
      reduce_inter_warp( inter_warp_data , tidx );

      // If one or less blocks then done.
      // Thread 0 has the final reduction value
      if ( ! group_count ) {
#if 0
        if ( 0 == thread_of_block ) {
          m_reduce.final( shared_data );
        }
#endif
        for ( unsigned i = thread_of_block ; i < m_data.value_word_count ; ++i ) {
          m_output_space[i] = shared_data[i] ;
        }
        break ;
      }

      //----------------------------------
      // This block is a member of group of fan-in blocks.

      const size_type group_id    = block_id / m_data.thread_count ;
      const size_type flag_offset = group_offset - group_count ;

      //----------------------------------

      write_block_value( group_offset + block_id ,
                         thread_of_block , shared_data );

      //----------------------------------
      {
        // Need a unique location for this flag which takes into
        // account the block-fanin iteration loop
        //
        // Thread #0 increment the group write counter.
        // If the old value is 'group_size-1' then atomicInc
        // sets the value to zero and returns the old value.
        //
        // Inform the entire block of the last_block status.
        // If this is not the last block in the group then
        // this block is finished.

        const int not_last_block = thread_of_block ? 0 :
          ( 1 + atomicInc( m_scratch_flags + flag_offset + group_id , 
                           group_size - 1 ) ) < group_size ;

        if ( __syncthreads_or( not_last_block ) ) break ;
      }

      // Last block to complete performs this group's reduction
      //----------------------------------
      // Read data for group_size blocks.
      // This warp reads [ wbeg .. wbeg + WarpSize ).
      // A warp's shared memory is contiguous but
      // there is a gap between warps to avoid bank conflicts.

      m_reduce.init( thread_data );

      if ( wbeg < group_size ) {

        // This warp has data to read.

        size_type * const shared =
          shared_data + m_data.shared_data_offset(0,tidy);

        volatile size_type * const scratch =
          m_scratch_space + m_data.value_word_stride * (
            group_offset                       + // offset for reduction cycle
            ( group_id * m_data.thread_count ) + // offset for block
            wbeg );                              // offset for warp

        const size_type count = group_size - wbeg ;

        const size_type read_count = m_data.value_word_stride * (
          count < CudaTraits::WarpSize ? count : CudaTraits::WarpSize );

        for ( size_type i = tidx ;
                        i < read_count ;
                        i += CudaTraits::WarpSize ) {
          shared[i] = scratch[i] ;
        }
      }

      //------------------------
      // Next iteration...

      block_id       =  group_id ;
      group_size     =  m_data.thread_count ;
      group_count   /=  m_data.thread_count ;
      group_offset   =  flag_offset ;
    }

    return ;
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelReduce< FunctorType , WorkSpec , Cuda >
{
public:

  typedef          Cuda                   device_type ;
  typedef          Cuda::size_type        size_type ;
  typedef typename FunctorType::value_type  value_type ;
  typedef typename ReduceAdapter< FunctorType >::pointer_type pointer_type ;

  typedef CudaReduceShared< FunctorType , void > ReduceType ;

  //----------------------------------------------------------------------

        ReduceType     m_reduce_shared ;
  const size_type      m_work_count ;
  pointer_type         m_host_ptr ;

  //----------------------------------------------------------------------

public:

  // For the multi-functor reduce:
  ParallelReduce( const FunctorType  & functor ,
                  const size_type      work_count ,
                  const size_type      global_block_begin ,
                  const size_type      global_block_count )
    : m_reduce_shared( functor )
    , m_work_count(    work_count )
  {
    m_reduce_shared.assign_block_range( global_block_begin , global_block_count );
  }

  //--------------------------------------------------------------------------

  ParallelReduce( const size_type      work_count ,
                  const FunctorType  & functor ,
                  pointer_type         result )
    : m_reduce_shared( functor )
    , m_work_count(    work_count )
    , m_host_ptr(    result )
  {
    const size_type value_size = m_reduce_shared.m_reduce.value_size();
    const size_type nt = Impl::CudaTraits::WarpSize * ReduceType::warp_count( value_size );

    if ( 0 == nt ) return ;

    // One level:
    size_type nb = std::min( nt , size_type( ( work_count + nt - 1 ) / nt ) );

    // Two level:
    // size_type nb = std::min( nt * nt , size_type( ( work_count + nt - 1 ) / nt ) );

    m_reduce_shared.assign_block_range( 0 , nb );

    if ( m_work_count ) {
      const size_type nw = m_reduce_shared.warp_count();

      const dim3 block( Impl::CudaTraits::WarpSize * nw , 1, 1 );
      const dim3 grid( nb , 1 , 1 );
      const size_type shmem_size = m_reduce_shared.shmem_size();

      CudaParallelLaunch< ParallelReduce >( *this , grid , block , shmem_size );
    }
  }

  void wait() const
  {
    typedef typename ReduceType::reduce_operator::pointer_type ptr_type ;

    Cuda::fence();
    const unsigned size  = m_reduce_shared.m_reduce.value_size();
    const unsigned count = m_reduce_shared.m_reduce.value_count();
    ptr_type ptr = (ptr_type) cuda_internal_scratch_unified( size );
    if ( 0 != ptr ) {
      for ( Cuda::size_type i = 0 ; i < count ; ++i )
        m_host_ptr[i] = ptr[i] ;
    }
    else {
      ptr = (ptr_type) cuda_internal_scratch_space( size );
      DeepCopy<HostSpace,CudaSpace>( m_host_ptr , ptr , size );
    }
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
      m_reduce_shared.m_reduce.m_functor( iwork , value );
    }

    m_reduce_shared( threadIdx.x , blockIdx.x );
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
      DeepCopy<HostSpace,CudaSpace>( host_pointer , ptr , sizeof(T) * count );
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template < class FunctorType , class ValueOper >
class CudaMultiFunctorParallelReduceMember ;

template < class ValueOper >
class CudaMultiFunctorParallelReduceMember<void,ValueOper> {
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
  void apply( const size_type  thread_count ,
              const size_type  block_count ,
              const size_type  global_block_offset ,
              const size_type  global_block_count ) const = 0 ;
};

template < class FunctorType , class ValueOper >
class CudaMultiFunctorParallelReduceMember :
  public CudaMultiFunctorParallelReduceMember<void,ValueOper>
{
public:
  typedef CudaMultiFunctorParallelReduceMember<void,ValueOper> base_type ;
  typedef ParallelReduce< FunctorType , size_t , Cuda >  driver_type ;
  
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
  void apply( const size_type  thread_count ,
              const size_type  block_count ,
              const size_type  global_block_begin ,
              const size_type  global_block_count ) const
  {
    const dim3 block( thread_count , 1 , 1 );
    const dim3 grid(  block_count ,  1 , 1 );

    driver_type  driver( m_functor ,
                         base_type::m_work_count ,
                         global_block_begin ,
                         global_block_count );

    const size_type shmem_size = driver.m_reduce_shared.shmem_size();

    cuda_parallel_launch_local_memory< driver_type ><<< grid , block , shmem_size >>>( driver );
  }
};

} // namespace Impl

template < class ValueOper >
class MultiFunctorParallelReduce< ValueOper , Cuda > {
public:
  typedef  Cuda               device_type ;
  typedef  Cuda::size_type    size_type ;

private:

  typedef Impl::CudaReduceShared< ValueOper , void > ReduceType ;

  typedef Impl::CudaMultiFunctorParallelReduceMember< void , ValueOper > MemberType ;
  typedef std::vector< MemberType * > MemberVector ;

  MemberVector    m_member_functors ;

  Impl::ReduceAdapter< ValueOper > m_reduce ;

  inline static 
  size_type block_count( const size_type warp_count ,
                         const size_type work_count )
  {
    const size_type nt = Impl::CudaTraits::WarpSize * warp_count ;

    return std::min( Impl::cuda_internal_maximum_grid_count() ,
                     ( work_count + nt - 1 ) / nt );
  }

public:

  explicit MultiFunctorParallelReduce( const ValueOper & f )
    : m_member_functors()
    , m_reduce( f )
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
      new Impl::CudaMultiFunctorParallelReduceMember< FunctorType , ValueOper >( f , work_count );

    m_member_functors.push_back( member );
  }

  void execute()
  {
    if ( m_member_functors.size() ) {
      const size_type warp_count   = ReduceType::warp_count( m_reduce.value_size() );
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

        if ( (*m)->m_work_count ) {
          const size_type n = block_count( warp_count , (*m)->m_work_count );

          member.apply( thread_count , n ,
                        global_block_offset , global_block_count );

          global_block_offset += n ;
        }
      }
    }
  }

  void output( typename ReduceType::reduce_operator::pointer_type host_pointer ) const
  {
    typedef typename ReduceType::reduce_operator::pointer_type ptr_type ;

    Cuda::fence();
    ptr_type ptr = (ptr_type) Impl::cuda_internal_scratch_unified( m_reduce.value_size() );
    if ( 0 != ptr ) {
      for ( Cuda::size_type i = 0 ; i < m_reduce.value_count() ; ++i )
        host_pointer[i] = ptr[i] ;
    }
    else {
      ptr = (ptr_type) Impl::cuda_internal_scratch_space( m_reduce.value_size() );
      DeepCopy<HostSpace,CudaSpace>( host_pointer , ptr , m_reduce.value_size() );
    }
  }

  //--------------------------------------------------------------------------

};

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#endif /* defined( __CUDACC__ ) */

#endif /* KOKKOS_CUDA_PARALLELREDUCE_HPP */

