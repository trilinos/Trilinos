/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_DEVICECUDA_PARALLELREDUCE_HPP
#define KOKKOS_DEVICECUDA_PARALLELREDUCE_HPP

#include <iostream>

#include <DeviceCuda/Kokkos_DeviceCuda_ParallelDriver.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Must have consistent '__shared__' statement across all device kernels.
// Since there may be more than one kernel in a file then have to make this
// a simple array of words.

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class DriverType >
__device__
void cuda_reduce_shared( const DeviceCuda::size_type used_warp_count )
{
  typedef          DeviceCuda::size_type  size_type ;
  typedef typename DriverType::value_type value_type ;

  typedef volatile value_type * vvp ;
  typedef volatile const value_type * cvvp ;

  enum { ValueWordStride = DriverType::ValueWordStride };
  enum { WarpStride      = DriverType::WarpStride };
  enum { HalfWarpSize    = Impl::DeviceCudaTraits::WarpSize >> 1 };

  extern __shared__ size_type shared_data[];

  // threadIdx.x == index within warp [ 0 .. WarpSize - 1 ]
  // threadIdx.y == which warp        [ 0 .. used_warp_count - 1 ]

  // Phase A: Reduce within my warp:
  //          Warp's reads occur before joins and writes
  //          so there is no race condition.
  //          Declare shared data to be volatile to
  //          prevent compiler from introducing a race condition.
  //
  if ( threadIdx.y < used_warp_count && threadIdx.x < HalfWarpSize ) {
    enum { n1  = ValueWordStride * 1 };
    enum { n2  = ValueWordStride * 2 };
    enum { n4  = ValueWordStride * 4 };
    enum { n8  = ValueWordStride * 8 };
    enum { n16 = ValueWordStride * 16 };

    size_type * const data = shared_data + DriverType::shared_data_offset();

    DriverType::join( *((vvp) data), *((cvvp)( data + n16 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n8 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n4 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n2 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n1 )) );
  }

  // Phase B: Use a single warp to reduce results from each warp.
  //          This requires: used_warp_count <= WarpSize
  //

  __syncthreads();

  if ( 0 == threadIdx.y && threadIdx.x + 1 < used_warp_count ) {
    enum { n1  = WarpStride * ValueWordStride * 1 };
    enum { n2  = WarpStride * ValueWordStride * 2 };
    enum { n4  = WarpStride * ValueWordStride * 4 };
    enum { n8  = WarpStride * ValueWordStride * 8 };
    enum { n16 = WarpStride * ValueWordStride * 16 };

    size_type * const data = shared_data + DriverType::shared_data_offset( 0 , threadIdx.x );

    if ( threadIdx.x + 2 < used_warp_count ) {
      if ( threadIdx.x + 4 < used_warp_count ) {
        if ( threadIdx.x + 8 < used_warp_count ) {
          if ( threadIdx.x + 16 < used_warp_count ) {
            DriverType::join( *((vvp) data) , *((cvvp)( data + n16 )) );
          }
          DriverType::join( *((vvp) data) , *((cvvp)( data + n8 )) );
        }
        DriverType::join( *((vvp) data) , *((cvvp)( data + n4 )) );
      }
      DriverType::join( *((vvp) data) , *((cvvp)( data + n2 )) );
    }
    DriverType::join( *((vvp) data) , *((cvvp)( data + n1 )) );
  }
}

//----------------------------------------------------------------------------

template< class DriverType >
__device__
void cuda_reduce_global( const DriverType & driver )
{
  typedef          DeviceCuda::size_type  size_type ;
  typedef typename DriverType::value_type value_type ;

  enum { WarpSize       = Impl::DeviceCudaTraits::WarpSize };
  enum { WarpIndexMask  = Impl::DeviceCudaTraits::WarpIndexMask };
  enum { WarpIndexShift = Impl::DeviceCudaTraits::WarpIndexShift };
  enum { ValueWordCount = DriverType::ValueWordCount };
  enum { WordsPerWarp   = DriverType::WordsPerWarp };

  extern __shared__ size_type shared_data[];

  const size_type thread_of_block = threadIdx.x + blockDim.x * threadIdx.y ;
  size_type i , j ;

  // Phase A: Output block's results to global memory
  //          and then input results into last block.
  //
  // ** REQUIRED: driver.m_grid_size <= blockDim.x * blockDim.y **

  if ( thread_of_block < ValueWordCount ) {
    // Coalesced global memory write, wait for write to complete
    // Write by multiblock_id and
    // Read  by ( threadIdx.x + blockDim.x * threadIdx.y ) == multiblock_id
    // Determine correct location into scratch memory.

    const size_type multiblock_id = driver.m_grid_offset + blockIdx.x ;
    const size_type thread_stride = blockDim.x * blockDim.y ;

    size_type * const scratch = driver.m_scratch_space +
      DriverType::shared_data_offset(
        multiblock_id &  WarpIndexMask  /* for threadIdx.x */ ,
        multiblock_id >> WarpIndexShift /* for threadIdx.y */ );

    for ( i = thread_of_block ; i < ValueWordCount ; i += thread_stride ) {
      scratch[i] = shared_data[i] ;
    }

    __threadfence(); // Wait for write to complete
  }

  // Check if this is the last block to finish:

  if ( 0 == thread_of_block ) {
    // atomicInc returns value prior to increment.

    shared_data[ DriverType::shared_flag_offset() ] =
      driver.m_grid_size == 1 + atomicInc( driver.m_scratch_flag , driver.m_grid_size + 1 );
  }
  __syncthreads();

  // The last block to finish reduces the results from all blocks.

  if ( shared_data[ DriverType::shared_flag_offset() ] ) {

    // Let the last thread of the last warp reinitialize the flag
    // for the next reduce operation.
    // This warp is least likely to be doing any more work.

    if ( blockDim.y == threadIdx.y + 1 &&
         blockDim.x == threadIdx.x + 1 ) {
      *(driver.m_scratch_flag) = 0 ;
    }

    // Each warp does a coalesced read of its own data.

    if ( threadIdx.y < driver.m_scratch_warp ) {

      // Coalesced global memory read for this warp's data.

      i = DriverType::shared_data_offset( 0 , threadIdx.y );
      j = i + WordsPerWarp ;

      if ( driver.m_scratch_upper < j ) {
        j = driver.m_scratch_upper ;

        // Only partial data will be read by this warp
        // so initialize the values before reading.
        DriverType::init( *((value_type *)( shared_data + DriverType::shared_data_offset() )) );
      }

      for ( i += threadIdx.x ; i < j ; i += WarpSize ) {
        shared_data[i] = driver.m_scratch_space[i] ;
      }
    }

    // Phase B: Reduce these contributions
    cuda_reduce_shared< DriverType >( driver.m_scratch_warp );
  }
}

//----------------------------------------------------------------------------

template< class DriverType >
__global__

#if KOKKOS_DEVICE_CUDA_USE_CONSTANT_MEMORY 

static void cuda_parallel_reduce()
{
  // The driver functor has been copied to constant memory

  const DriverType & driver =
    *((const DriverType *) kokkos_device_cuda_constant_memory_buffer );

#else

static void cuda_parallel_reduce( const DriverType driver )
{

#endif

  typedef          DeviceCuda::size_type  size_type ;
  typedef typename DriverType::value_type value_type ;

  extern __shared__ size_type shared_data[];

  value_type & value =
    *( (value_type *)( shared_data + DriverType::shared_data_offset() ) );

  DriverType::init( value );

  // Phase 1: Reduce to per-thread contributions
  {
    const size_type work_count  = driver.m_work_count ;
    const size_type work_stride = driver.m_work_stride ;

    size_type iwork =
      threadIdx.x + blockDim.x * ( threadIdx.y + blockDim.y * blockIdx.x );

    for ( ; iwork < work_count ; iwork += work_stride ) {
      driver.m_work_functor( iwork , value );
    }
  }

  // Phase 2: Reduce this block's thread's contributions
  //          to a single reduction value.
  cuda_reduce_shared< DriverType >( blockDim.y );

  // Phase 3: Reduce contributions from multiple blocks

  int last_block = 1 == driver.m_grid_size ;

  if ( ! last_block ) {

    cuda_reduce_global< DriverType >( driver );

    last_block = shared_data[ DriverType::shared_flag_offset() ];
  }

  // Phase 4: Thread #0 of last block performs serial finalization
  if ( last_block && 0 == threadIdx.x && 0 == threadIdx.y ) {
    driver.m_work_finalize( value );
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class FinalizeType >
class ParallelReduce< FunctorType , FinalizeType , DeviceCuda > {
public:
  typedef          DeviceCuda               device_type ;
  typedef          DeviceCuda::size_type    size_type ;
  typedef typename FunctorType::value_type  value_type ;

  enum { WarpSize   = Impl::DeviceCudaTraits::WarpSize };
  enum { WarpStride = WarpSize + 1 };
  enum { WarpIndexMask  = Impl::DeviceCudaTraits::WarpIndexMask };
  enum { WarpIndexShift = Impl::DeviceCudaTraits::WarpIndexShift };

  enum { SharedMemoryBanks = Impl::DeviceCudaTraits::SharedMemoryBanks_13 };

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

  //----------------------------------------------------------------------

  const FunctorType  m_work_functor ;
  const FinalizeType m_work_finalize ;
  const size_type    m_work_count ;
  const size_type    m_work_stride ;
  const size_type    m_grid_size ;
  const size_type    m_grid_offset ; /// to be used for multi-functor / multi-block reduction 

  // Scratch space for multi-block reduction
  // m_scratch_warp  == number of warps required
  // m_scratch_upper == upper bound of reduction scratch space used.

  size_type * const m_scratch_space ;
  size_type * const m_scratch_flag ;
  const size_type   m_scratch_warp ;
  const size_type   m_scratch_upper ;
  
  //----------------------------------------------------------------------

  static inline
  __device__
  size_type shared_data_offset()
  { return ValueWordStride * ( threadIdx.x + WarpStride * threadIdx.y ); }

  static inline
  __device__
  size_type shared_data_offset( size_type x , size_type y )
  { return ValueWordStride * ( x + WarpStride * y ); }

  static inline
  __device__
  size_type shared_flag_offset()
  { return ValueWordStride * ( WarpStride * blockDim.y - 1 ); }

  static inline
  __device__
  void join( volatile       value_type & update ,
             volatile const value_type & input )
    { FunctorType::join( update , input ); }

  static inline
  __device__
  void init( value_type & update )
    { FunctorType::init( update ); }

  //----------------------------------------------------------------------

private:

  ParallelReduce( const FunctorType  & functor ,
                  const FinalizeType & finalize ,
                  const size_type      work_count ,
                  const size_type      work_stride ,
                  const size_type      grid_size )
    : m_work_functor(  functor )
    , m_work_finalize( finalize )
    , m_work_count(    work_count )
    , m_work_stride(   work_stride )
    , m_grid_size(     grid_size )
    , m_grid_offset(   0 )
    , m_scratch_space( device_type::reduce_multiblock_scratch_space() )
    , m_scratch_flag(  device_type::reduce_multiblock_scratch_flag() )
    , m_scratch_warp( ( grid_size >> WarpIndexShift ) +
                      ( grid_size &  WarpIndexMask ? 1 : 0 ) )
    , m_scratch_upper( ValueWordStride * ( grid_size + m_scratch_warp - 1 ) )
  {}

public:
  //----------------------------------------------------------------------

  static
  void execute( const size_t         work_count ,
                const FunctorType  & functor ,
                const FinalizeType & finalize )
  {
    typedef ParallelReduce< FunctorType , FinalizeType , DeviceCuda > self_type ;

    const size_type maximum_shared_words = device_type::maximum_shared_words();

    dim3 block( Impl::DeviceCudaTraits::WarpSize , 
                device_type::maximum_warp_count() , 1 );

    while ( maximum_shared_words < block.y * WordsPerWarpStride ) {
      block.y >>= 1 ;
    }

    dim3 grid( 1 , 1 , 1 );

    if ( 0 == work_count ) { // Avoid infinite loop.
      block.y = 1 ;
    }
    else if ( work_count <= WarpSize * block.y ) { // Need at most one block
      while ( work_count <= WarpSize * ( block.y >> 1 ) ) { block.y >>= 1 ; }
    }
    else {
      const size_t threads_per_block = WarpSize * block.y ;

      grid.x = ( work_count + threads_per_block - 1 ) / threads_per_block ;

      // At most one block per thread so that the final reduction
      // operation can process one reduction value per thread.
      if ( grid.x > threads_per_block ) { grid.x = threads_per_block ; }
    }

    const size_type shmem_size =
      sizeof(size_type) * ( ValueWordStride * (WarpStride * block.y - 1) + 1 );

    device_type::set_dispatch_functor();

    self_type driver( functor , finalize , work_count , block.x * block.y * grid.x , grid.x );

    device_type::clear_dispatch_functor();

#if KOKKOS_DEVICE_CUDA_USE_CONSTANT_MEMORY 

    // Copy functor to constant memory on the device
    cudaMemcpyToSymbol( kokkos_device_cuda_constant_memory_buffer , & driver , sizeof(self_type) );

    // Invoke the driver function on the device
    cuda_parallel_reduce< self_type ><<< grid , block , shmem_size >>>();

#else

    cuda_parallel_reduce< self_type ><<< grid , block , shmem_size >>>( driver );

#endif
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelReduce< FunctorType , void , DeviceCuda > 
{
public:
  typedef typename FunctorType::value_type     value_type ;
  typedef ValueView< value_type , DeviceCuda > view_type ;

  void execute( const size_t        work_count ,
                const FunctorType & work_functor ,
                      value_type  & result )
  {
    view_type tmp =
      create_labeled_value< value_type , DeviceCuda >(
        std::string("parallel_reduce_temporary_result") );

    ParallelReduce< FunctorType , view_type , DeviceCuda >
      ::execute( work_count , work_functor , tmp );

    deep_copy( result , tmp );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

#include <Kokkos_DeviceClear_macros.hpp>

#endif /* KOKKOS_DEVICECUDA_PARALLELREDUCE_HPP */

