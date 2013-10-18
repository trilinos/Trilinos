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

#ifndef KOKKOS_CUDA_PARALLELSCAN_HPP
#define KOKKOS_CUDA_PARALLELSCAN_HPP

#if defined( __CUDACC__ )

#include <stdlib.h>
#include <utility>

#include <Kokkos_Parallel.hpp>
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
/*
 *  Algorithmic constraints:
 *   (a) blockDimSize is a power of two
 *   (b) blockDim.x == BlockSize == 1 << BlockSizeShift
 *   (c) blockDim.y == blockDim.z == 1
 *   (d) gridDim.x  <= BlockSize * BlockSize
 *   (e) gridDim.y  == gridDim.z == 1
 */
template< class FunctorType , unsigned WarpCount >
struct CudaReduceScan
{
  typedef Cuda::size_type               size_type ;
  typedef ReduceAdapter< FunctorType >  Reduce ;

  enum { BlockSize      = CudaTraits::WarpSize << power_of_two< WarpCount >::value };
  enum { BlockSizeShift = power_of_two< BlockSize >::value };
  enum { BlockSizeMask  = BlockSize - 1 };

  const integral_nonzero_constant< size_type , Reduce::StaticValueSize >> power_of_two<sizeof(size_type)>::value >
    word_count ;

  CudaReduceScan( const FunctorType & f )
    : word_count( Reduce::value_size( f ) >> power_of_two<sizeof(size_type)>::value ) {}

  __device__ inline
  size_type * intra_block_reduce( const FunctorType & functor , size_type * const base_data ) const
  {
#define BLOCK_REDUCE_STEP( S ) \
  if ( ! ( rtid & ((1<<(S+1))-1) ) ) { Reduce::join(functor,tdata,tdata - (word_count.value << S) ); }

    { // Intra-warp reduction:
      const unsigned rtid = threadIdx.x ^ BlockSizeMask ;
      size_type * const tdata = base_data + word_count.value * threadIdx.x ;
      BLOCK_REDUCE_STEP(0)
      BLOCK_REDUCE_STEP(1)
      BLOCK_REDUCE_STEP(2)
      BLOCK_REDUCE_STEP(3)
      BLOCK_REDUCE_STEP(4)
    }

    __syncthreads(); // Wait for all warps to reduce

    { // Inter-warp reduction by a single warp to avoid extra synchronizations
      const unsigned rtid = ( threadIdx.x ^ BlockSizeMask ) << CudaTraits::WarpIndexShift ;
      if ( rtid < BlockSize ) {
        enum { DO_5 = (1<<5) < BlockSizeMask };
        enum { DO_6 = (1<<6) < BlockSizeMask };
        enum { DO_7 = (1<<7) < BlockSizeMask };
        enum { DO_8 = (1<<8) < BlockSizeMask };
        size_type * const tdata = base_data + word_count.value * ( rtid ^ BlockSizeMask );
        if ( DO_5 ) {                        BLOCK_REDUCE_STEP(5) }
        if ( DO_6 ) { __threadfence_block(); BLOCK_REDUCE_STEP(6) }
        if ( DO_7 ) { __threadfence_block(); BLOCK_REDUCE_STEP(7) }
        if ( DO_8 ) { __threadfence_block(); BLOCK_REDUCE_STEP(8) }
      }
    }

    __syncthreads(); // Wait for final value to be available

    return base_data + word_count.value * BlockSizeMask ; // Location of grand total

#undef BLOCK_REDUCE_STEP
  }

  __device__
  size_type * intra_block_reduce_scan( const FunctorType & functor , size_type * const data ) const
  {
#define BLOCK_SCAN_STEP( S ) \
  if ( n == (1<<S) ) { Reduce::join( functor , tdata , tdata - (word_count.value << S) ); }

    size_type * const reduce_total = intra_block_reduce( functor , data );

    { // Intra warp scan by a single warp to avoid extra synchronization:

      const unsigned rtid = ( threadIdx.x ^ BlockSizeMask ) << CudaTraits::WarpIndexShift ;

      if ( rtid < BlockSize ) {
        size_type * const tdata = data + word_count.value * ( rtid ^ BlockSizeMask );

        int n = ( rtid &  32 ) ?  32 : (
                ( rtid &  64 ) ?  64 : (
                ( rtid & 128 ) ? 128 : (
                ( rtid & 256 ) ? 256 : 0 )));

        if ( ! ( rtid + n < BlockSize ) ) n = 0 ;

        BLOCK_SCAN_STEP(8)
        BLOCK_SCAN_STEP(7)
        BLOCK_SCAN_STEP(6)
        BLOCK_SCAN_STEP(5)
      }
    }

    __syncthreads() ;

    {
      const unsigned rtid = threadIdx.x ^ BlockSizeMask ;

      size_type * const tdata = data + word_count.value * threadIdx.x ;

      int n = ( rtid &  1 ) ?  1 : (
              ( rtid &  2 ) ?  2 : (
              ( rtid &  4 ) ?  4 : (
              ( rtid &  8 ) ?  8 : (
              ( rtid & 16 ) ? 16 : 0 ))));

      if ( ! ( rtid + n < BlockSize ) ) n = 0 ;

      BLOCK_SCAN_STEP(4) __threadfence_block();
      BLOCK_SCAN_STEP(3) __threadfence_block();
      BLOCK_SCAN_STEP(2) __threadfence_block();
      BLOCK_SCAN_STEP(1) __threadfence_block();
      BLOCK_SCAN_STEP(0)
    }

    return reduce_total ;

#undef BLOCK_SCAN_STEP
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelScan< FunctorType , WorkSpec , Cuda >
{
public:
  // Algorithmic constraints:
  //  (a) blockDimSize is a power of two
  //  (b) blockDim.x == BlockSize == 1 << BlockSizeShift
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.x * blockDim.x
  //  (d) gridDim.y  == gridDim.z == 1

  // blockDim.x must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
  // gridDim.x < blockDim.x * blockDim.x
  //
  // Choose 8 warps, which enables latency hiding and limits gridDim.x to 65k

  enum { WarpCount = 8 };

  typedef CudaReduceScan< FunctorType , WarpCount > ReduceScan ;
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef Cuda::size_type                     size_type ;

  enum { GridMaxComputeCapability_2x = 0x0ffff };
  enum { GridMax = ( ReduceScan::BlockSize * ReduceScan::BlockSize ) < GridMaxComputeCapability_2x
                 ? ( ReduceScan::BlockSize * ReduceScan::BlockSize ) : GridMaxComputeCapability_2x };

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  const size_type   m_work ;
        size_type   m_work_per_block ;
        size_type   m_final ;
  const ReduceScan  m_scan ;
  
  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    extern __shared__ size_type shared_data[];

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_value = shared_data + m_scan.word_count.value * ( threadIdx.x + 1 );

    if ( 0 == threadIdx.x ) { m_functor.init( Reduce::reference( shared_data ) ); }

    m_functor.init( Reduce::reference( shared_value ) );

    {
      // Number of blocks is bounded so that the reduction can be limited to two passes.
      // Each thread block is given an approximately equal amount of work to perform.
      // Accumulate the values for this block.
      // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

      const size_type iwork_beg = blockIdx.x * m_work_per_block ;
      const size_type iwork_end = iwork_beg + m_work_per_block < m_work 
                                ? iwork_beg + m_work_per_block : m_work ;

      for ( size_type iwork = threadIdx.x + iwork_beg ; iwork < iwork_end ; iwork += ReduceScan::BlockSize ) {
        m_functor( iwork , Reduce::reference( shared_value ) , false );
      }
    }

    {
      // Reduction for exclusive scan startes at location[1]
      size_type * const block_total = m_scan.intra_block_reduce( m_functor , shared_data + m_scan.word_count.value );

      // Reduce the accumulation for the entire block and write to global scratch space:
      size_type * const global = m_scratch_space + m_scan.word_count.value * blockIdx.x ;

      // Write out reduction total for this block
      for ( size_type i = threadIdx.x ; i < m_scan.word_count.value ; i += ReduceScan::BlockSize ) { global[i] = block_total[i] ; }
    }

    // if gridDim.x <= BlockSize then one pass else two pass
    // first pass global values are contiguous,
    // second pass global values from last entry of each group of blocks

    for ( size_type valueTotal = gridDim.x , valueId = blockIdx.x , first_pass_contiguous = 1 ; 1 < valueTotal ; ) {

      // This block contributed to the range [ value_begin .. value_begin + value_count )
      // The last block to contribute to this range is responsible for scanning this range.

      const size_type value_begin = unsigned(valueId) & ~unsigned(ReduceScan::BlockSizeMask);
      const size_type value_count = ( ( value_begin + ReduceScan::BlockSize ) < valueTotal )
                                  ? ReduceScan::BlockSize : valueTotal - value_begin ;

      // This block is a member of a group which reduces to the following value
      valueId >>= ReduceScan::BlockSizeShift ;

      // How many values to reduce after this reduction pass?
      valueTotal = ( valueTotal + ReduceScan::BlockSizeMask ) >> ReduceScan::BlockSizeShift ;

      {
        // Contributing blocks note that their contribution has been completed via an atomic-increment flag
        size_type * const flag = m_scratch_flags + ( 1 < valueTotal ? 1 + valueId : 0 );

        // If this block is not the last block to contribute to this group then the block is done.
        const int not_last_block =
          __syncthreads_or( threadIdx.x ? 0 : ( 1 + atomicInc( flag , value_count - 1 ) < value_count ) );

        if ( not_last_block ) break ;
      }

      {
        // thread reads or initializes its data
        // First pass scans a contiguous span of global data
        // Second pass (if any) scans the last value of each block ( value_begin == 0 )

        size_type * const global = m_scratch_space + m_scan.word_count.value *
          ( first_pass_contiguous ? ( value_begin + threadIdx.x )
                                  : ( ( ( threadIdx.x + 1 ) << ReduceScan::BlockSizeShift ) - 1 ) );
        Reduce::copy( m_functor , shared_value , ( threadIdx.x < value_count ? global : shared_data ) );

        m_scan.intra_block_reduce_scan( m_functor , shared_data + m_scan.word_count.value );

        if ( threadIdx.x < value_count ) {
          Reduce::copy( m_functor , global , shared_value );
        }
      }

      first_pass_contiguous = 0 ;
    }
  }

  //----------------------------------------

  __device__ inline
  void final(void) const
  {
    extern __shared__ size_type shared_data[];

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_prefix = shared_data + m_scan.word_count.value * threadIdx.x ;
    size_type * const shared_accum  = shared_data + m_scan.word_count.value * ( ReduceScan::BlockSize + 1 );

    // Starting value for this thread block:
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + m_scan.word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.x ; i < m_scan.word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.x ) {
      m_functor.init( Reduce::reference( shared_accum ) );
    }

    // The first pass performed an inclusive scan with each group of block
    // and an inlusive scan across the total of each group of blocks.
    // If the exclusive scan value for this block is not a group-total
    // then must sum the prior group's inclusive scan total.

    if ( ( 0 == threadIdx.x ) &&
         ( ReduceScan::BlockSize < blockIdx.x ) &&
         ( blockIdx.x & ReduceScan::BlockSizeMask ) ) {
        /* Not the first group of blocks AND Not the global reduction block */

      m_functor.join( Reduce::reference( shared_accum ) ,
                      Reduce::reference( m_scratch_space + m_scan.word_count.value *
                                         ( ( blockIdx.x & ~ReduceScan::BlockSizeMask ) - 1 ) ) );
    }

          unsigned iwork_beg = blockIdx.x * m_work_per_block ;
    const unsigned iwork_end = iwork_beg + m_work_per_block ;

    for ( ; iwork_beg < iwork_end ; iwork_beg += ReduceScan::BlockSize ) {

      const unsigned iwork = threadIdx.x + iwork_beg ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      m_functor.init( Reduce::reference( shared_prefix + m_scan.word_count.value ) );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.x ; i < m_scan.word_count.value ; ++i ) {
        shared_data[i + m_scan.word_count.value] = shared_data[i] = shared_accum[i] ;
      }

      if ( CudaTraits::WarpSize < m_scan.word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < m_work ) { m_functor( iwork , Reduce::reference( shared_prefix + m_scan.word_count.value ) , false ); }

      // Scan block values into locations shared_data[1..ReduceScan::BlockSize]
      size_type * const block_total = m_scan.intra_block_reduce_scan( m_functor , shared_data + m_scan.word_count.value );

      for ( unsigned i = threadIdx.x ; i < m_scan.word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }

      // Call functor with exclusive scan value
      if ( iwork < m_work ) { m_functor( iwork , Reduce::reference( shared_prefix ) , true ); }
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
                const size_t         nwork )
  : m_functor( functor )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_work( nwork )
  , m_work_per_block( 0 )
  , m_final( false )
  , m_scan( functor )
  {
    // At most 'max_grid' blocks:
    const int max_grid = std::min( int(GridMax) , int(( nwork + ReduceScan::BlockSizeMask ) / ReduceScan::BlockSize ));

    // How much work per block:
    m_work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

    // How many block are really needed for this much work:
    const dim3 grid( ( nwork + m_work_per_block - 1 ) / m_work_per_block , 1 , 1 );
    const dim3 block( ReduceScan::BlockSize , 1 , 1 );
    const int shmem = Reduce::value_size( functor ) * ( ReduceScan::BlockSize + 2 );

    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * grid.x );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) * ( 1 + grid.x ) );

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

#endif /* KOKKOS_CUDA_PARALLELSCAN_HPP */

