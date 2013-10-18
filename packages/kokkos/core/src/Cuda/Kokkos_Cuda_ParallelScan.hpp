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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

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
 *   (a) blockDim.x is a power of two
 *   (b) blockDim.x <= 512
 *   (c) blockDim.y == blockDim.z == 1
 *   (d) gridDim.x  <= BlockDim.x * BlockDim.x
 *   (e) gridDim.y  == gridDim.z == 1
 */
template< bool DoScan , class FunctorType >
__device__
void cuda_intra_block_reduce_scan( const FunctorType & functor ,
                                   const typename ReduceAdapter< FunctorType >::pointer_type base_data )
{
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  const unsigned value_count   = Reduce::value_count( functor );
  const unsigned BlockSizeMask = blockDim.x - 1 ;

  // Must have power of two thread count

  if ( BlockSizeMask & blockDim.x ) { cuda_abort("Cuda::cuda_intra_block_scan requires power-of-two blockDim"); }

#define BLOCK_REDUCE_STEP( R , TD , S )  \
  if ( ! ( R & ((1<<(S+1))-1) ) ) \
    { functor.join( Reduce::reference(TD) , Reduce::reference(TD - (value_count<<S))); }

#define BLOCK_SCAN_STEP( TD , N , S )  \
  if ( N == (1<<S) ) \
    { functor.join( Reduce::reference(TD) , Reduce::reference(TD - (value_count<<S))); }

  const unsigned     rtid_intra = threadIdx.x ^ BlockSizeMask ;
  const pointer_type tdata_intra = base_data + value_count * threadIdx.x ;

  { // Intra-warp reduction:
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,0)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,1)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,2)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,3)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,4)
  }

  __syncthreads(); // Wait for all warps to reduce

  { // Inter-warp reduce-scan by a single warp to avoid extra synchronizations
    const unsigned rtid_inter = ( threadIdx.x ^ BlockSizeMask ) << CudaTraits::WarpIndexShift ;

    if ( rtid_inter < blockDim.x ) {

      const pointer_type tdata_inter = base_data + value_count * ( rtid_inter ^ BlockSizeMask );

      if ( (1<<5) < BlockSizeMask ) {                        BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,5) }
      if ( (1<<6) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,6) }
      if ( (1<<7) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,7) }
      if ( (1<<8) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,8) }

      if ( DoScan ) {

        int n = ( rtid_inter &  32 ) ?  32 : (
                ( rtid_inter &  64 ) ?  64 : (
                ( rtid_inter & 128 ) ? 128 : (
                ( rtid_inter & 256 ) ? 256 : 0 )));

        if ( ! ( rtid_inter + n < blockDim.x ) ) n = 0 ;

        BLOCK_SCAN_STEP(tdata_inter,n,8)
        BLOCK_SCAN_STEP(tdata_inter,n,7)
        BLOCK_SCAN_STEP(tdata_inter,n,6)
        BLOCK_SCAN_STEP(tdata_inter,n,5)
      }
    }
  }

  __syncthreads(); // Wait for inter-warp reduce-scan to complete

  if ( DoScan ) {
    int n = ( rtid_intra &  1 ) ?  1 : (
            ( rtid_intra &  2 ) ?  2 : (
            ( rtid_intra &  4 ) ?  4 : (
            ( rtid_intra &  8 ) ?  8 : (
            ( rtid_intra & 16 ) ? 16 : 0 ))));

    if ( ! ( rtid_intra + n < blockDim.x ) ) n = 0 ;

    BLOCK_SCAN_STEP(tdata_intra,n,4) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,3) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,2) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,1) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,0)
  }

#undef BLOCK_SCAN_STEP
#undef BLOCK_REDUCE_STEP
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelScan< FunctorType , WorkSpec , Cuda >
{
public:
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef Cuda::size_type                     size_type ;

  // Algorithmic constraints:
  //  (a) blockSize is a power of two
  //  (b) blockDim.x == BlockSize == 1 << BlockSizeShift
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.x * blockDim.x
  //  (d) gridDim.y  == gridDim.z == 1

  // blockDim.x must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
  // gridDim.x <= blockDim.x * blockDim.x
  //
  // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

  enum { WarpCount      = 4 };
  enum { BlockSize      = CudaTraits::WarpSize << power_of_two< WarpCount >::value };
  enum { BlockSizeShift = power_of_two< BlockSize >::value };
  enum { BlockSizeMask  = BlockSize - 1 };

  enum { GridMaxComputeCapability_2x = 0x0ffff };
  enum { GridMax = ( BlockSize * BlockSize ) < GridMaxComputeCapability_2x
                 ? ( BlockSize * BlockSize ) : GridMaxComputeCapability_2x };

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  const size_type   m_work ;
        size_type   m_work_per_block ;
        size_type   m_final ;
  
  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    extern __shared__ size_type shared_data[];

    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_count( m_functor ) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_value = shared_data + word_count.value * ( threadIdx.x + 1 );

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

      for ( size_type iwork = threadIdx.x + iwork_beg ; iwork < iwork_end ; iwork += BlockSize ) {
        m_functor( iwork , Reduce::reference( shared_value ) , false );
      }
    }

    {
      cuda_intra_block_reduce_scan<false>( m_functor , pointer_type(shared_data + word_count.value) );

      size_type * const block_total = shared_data + word_count.value * BlockSize ;

      // Reduce the accumulation for the entire block and write to global scratch space:
      size_type * const global = m_scratch_space + word_count.value * blockIdx.x ;

      // Write out reduction total for this block
      for ( size_type i = threadIdx.x ; i < word_count.value ; i += BlockSize ) { global[i] = block_total[i] ; }
    }

    // if gridDim.x <= BlockSize then one pass else two pass
    // first pass global values are contiguous,
    // second pass global values from last entry of each group of blocks

    for ( size_type valueTotal = gridDim.x , valueId = blockIdx.x , first_pass_contiguous = 1 ; 1 < valueTotal ; ) {

      // This block contributed to the range [ value_begin .. value_begin + value_count )
      // The last block to contribute to this range is responsible for scanning this range.

      const size_type value_begin = unsigned(valueId) & ~unsigned(BlockSizeMask);
      const size_type value_count = ( ( value_begin + BlockSize ) < valueTotal )
                                  ? BlockSize : valueTotal - value_begin ;

      // This block is a member of a group which reduces to the following value
      valueId >>= BlockSizeShift ;

      // How many values to reduce after this reduction pass?
      valueTotal = ( valueTotal + BlockSizeMask ) >> BlockSizeShift ;

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

        size_type * const global = m_scratch_space + word_count.value *
          ( first_pass_contiguous ? ( value_begin + threadIdx.x )
                                  : ( ( ( threadIdx.x + 1 ) << BlockSizeShift ) - 1 ) );
        Reduce::copy( m_functor , shared_value , ( threadIdx.x < value_count ? global : shared_data ) );

        cuda_intra_block_reduce_scan<true>( m_functor , pointer_type(shared_data + word_count.value) );

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

    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_count( m_functor ) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_prefix = shared_data + word_count.value * threadIdx.x ;
    size_type * const shared_accum  = shared_data + word_count.value * ( BlockSize + 1 );

    // Starting value for this thread block:
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.x ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.x ) {
      m_functor.init( Reduce::reference( shared_accum ) );
    }

    // The first pass performed an inclusive scan with each group of block
    // and an inclusive scan across the total of each group of blocks.
    // If the exclusive scan value for this block is not a group-total
    // then must sum the prior group's inclusive scan total.

    if ( ( 0 == threadIdx.x ) &&
         ( BlockSize < blockIdx.x ) &&
         ( blockIdx.x & BlockSizeMask ) ) {
        /* Not the first group of blocks AND Not the global reduction block */

      m_functor.join( Reduce::reference( shared_accum ) ,
                      Reduce::reference( m_scratch_space + word_count.value *
                                         ( ( blockIdx.x & ~BlockSizeMask ) - 1 ) ) );
    }

          unsigned iwork_beg = blockIdx.x * m_work_per_block ;
    const unsigned iwork_end = iwork_beg + m_work_per_block ;

    for ( ; iwork_beg < iwork_end ; iwork_beg += BlockSize ) {

      const unsigned iwork = threadIdx.x + iwork_beg ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      m_functor.init( Reduce::reference( shared_prefix + word_count.value ) );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.x ; i < word_count.value ; ++i ) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i] ;
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < m_work ) { m_functor( iwork , Reduce::reference( shared_prefix + word_count.value ) , false ); }

      // Scan block values into locations shared_data[1..BlockSize]
      cuda_intra_block_reduce_scan<true>( m_functor , Reduce::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.x ;
        for ( unsigned i = threadIdx.x ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

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
  {
    // At most 'max_grid' blocks:
    const int max_grid = std::min( int(GridMax) , int(( nwork + BlockSizeMask ) / BlockSize ));

    // How much work per block:
    m_work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

    // How many block are really needed for this much work:
    const dim3 grid( ( nwork + m_work_per_block - 1 ) / m_work_per_block , 1 , 1 );
    const dim3 block( BlockSize , 1 , 1 );
    const int shmem = Reduce::value_size( functor ) * ( BlockSize + 2 );

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

#if defined( __CUDA_ARCH__ )

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

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template< typename TypeLocal , typename TypeGlobal >
__device__ inline TypeGlobal Cuda::team_scan( const TypeLocal & value , TypeGlobal * const global_accum )
{
  enum { BlockSizeMax = 512 };

  __shared__ TypeGlobal base_data[ BlockSizeMax + 1 ];

  __syncthreads(); // Don't write in to shared data until all threads have entered this function

  if ( 0 == threadIdx.x ) { base_data[0] = 0 ; }

  base_data[ threadIdx.x + 1 ] = value ;

  Impl::cuda_intra_block_reduce_scan<true>( Impl::CudaJoinFunctor<TypeGlobal>() , base_data + 1 );

  if ( global_accum ) {
    if ( blockDim.x == threadIdx.x + 1 ) {
      base_data[ blockDim.x ] = atomic_fetch_add( global_accum , base_data[ blockDim.x ] );
    }
    __syncthreads(); // Wait for atomic
    base_data[ threadIdx.x ] += base_data[ blockDim.x ] ;
  }

  return base_data[ threadIdx.x ];
}

template< typename Type >
__device__ inline Type Cuda::team_scan( const Type & value )
{ return team_scan( value , (Type*) 0 ); }

} // namespace Kokkos

#else

namespace Kokkos {

template< typename Type > inline Type Cuda::team_scan( const Type & ) { return 0 ; }

template< typename TypeLocal , typename TypeGlobal >
inline TypeGlobal Cuda::team_scan( const TypeLocal & , TypeGlobal * const ) { return 0 ; }

} // namespace Kokkos

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined( __CUDACC__ ) */

#endif /* KOKKOS_CUDA_PARALLELSCAN_HPP */

