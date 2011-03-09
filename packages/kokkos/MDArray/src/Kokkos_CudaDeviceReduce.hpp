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

#ifndef KOKKOS_CUDADEVICEREDUCE_HPP
#define KOKKOS_CUDADEVICEREDUCE_HPP

#include <Kokkos_ParallelReduce.hpp>
#include <Kokkos_ValueView.hpp>
#include <Kokkos_CudaDevice.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

typedef unsigned int CudaWordType ;

template< typename ValueType >
union CudaSharedMemoryReduceType {
  // Pad the data size to an odd number of words to provide the proper
  // coalesced shared memory access.
  enum { WordCount =
        ( sizeof(ValueType) + sizeof(CudaWordType) - 1 ) / sizeof(CudaWordType) };

  ValueType    value ;
  unsigned int storage[ WordCount + ( WordCount % 32 ? 0 : 1 ) ];
};

//----------------------------------------------------------------------------

template< class ReduceOperators >
__device__
void reduce_shared_on_cuda(
  CudaSharedMemoryReduceType< typename ReduceOperators::value_type >
    * const shared_local )
{
  for ( unsigned int j = blockDim.x ; j ; ) {
    j >>= 1 ;

#ifndef __DEVICE_EMULATION__
    // Wait for contributing thread from a different half-warp
    if ( warpSize < j ) { __syncthreads(); }
#else
    __syncthreads();
#endif

    if ( threadIdx.x < j ) {
      ReduceOperators::join( shared_local->value , shared_local[j].value );
    }
  }
}

//----------------------------------------------------------------------------

template< typename ValueType >
inline
__device__
CudaSharedMemoryReduceType< ValueType > *
cuda_local_thread_reduce_shared_memory()
{
  typedef CudaSharedMemoryReduceType< ValueType > reduce_type ;

  // Compiler insists that all 'extern __shared__'
  // statements be identical, cannot template the type.

  extern __shared__ CudaWordType shared_memory[] ;

  // Now cast to the actual type:
  return ( (reduce_type *) shared_memory ) + threadIdx.x ;
}

//----------------------------------------------------------------------------
// Single block, reduce contributions

template< class FunctorType , class FinalizeType >
__global__
void run_reduce_operator_on_cuda(
  CudaSharedMemoryReduceType< typename FunctorType::value_type >
    * const block_result ,
  const FinalizeType finalize )
{
  typedef typename FunctorType::value_type         value_type ;
  typedef CudaSharedMemoryReduceType< value_type > reduce_type ;

  // Now cast to the actual type:
  reduce_type * const shared_local =
    cuda_local_thread_reduce_shared_memory< value_type >();

  shared_local->value = block_result[ threadIdx.x ].value ;

  reduce_shared_on_cuda< FunctorType >( shared_local );

  finalize();
}

//----------------------------------------------------------------------------

template< class FunctorType , class FinalizeType >
__global__
void run_reduce_functor_on_cuda(
  const size_t       work_count ,
  const FunctorType  functor ,
  const FinalizeType finalize )
{
  typedef typename FunctorType::value_type         value_type ;
  typedef CudaSharedMemoryReduceType< value_type > reduce_type ;

  reduce_type * const shared_local =
    cuda_local_thread_reduce_shared_memory< value_type >();

  FunctorType::init( shared_local->value );

  const size_t work_stride = blockDim.x * gridDim.x ;

  for ( size_t iwork = threadIdx.x + blockDim.x * blockIdx.x ;
        iwork < work_count ; iwork += work_stride ) {
    functor( iwork , shared_local->value );
  }

  reduce_shared_on_cuda< FunctorType >( shared_local );

  finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** \brief  Finalize via application provided serial functor */
template< typename ValueType , class FinalizeFunctor >
struct CudaParallelReduceFinalizeFunctor {

  FinalizeFunctor serial_finalize ;

  ~CudaParallelReduceFinalizeFunctor() {}

  CudaParallelReduceFinalizeFunctor( const FinalizeFunctor & rhs )
    : serial_finalize( rhs ) {}

  CudaParallelReduceFinalizeFunctor( const CudaParallelReduceFinalizeFunctor & rhs )
    : serial_finalize( rhs.serial_finalize ) {}

  __device__
  void operator()() const
  {
    typedef CudaSharedMemoryReduceType< ValueType > reduce_type ;

    reduce_type * const shared_local =
      cuda_local_thread_reduce_shared_memory< ValueType >();

    if ( gridDim.x == 1 && threadIdx.x == 0 ) {
      serial_finalize( shared_local->value );
    }
  }

private:
  CudaParallelReduceFinalizeFunctor();
  CudaParallelReduceFinalizeFunctor & operator =
    ( const CudaParallelReduceFinalizeFunctor & );
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct CudaParallelReduceFinalizeBlock {

  typedef CudaSharedMemoryReduceType< ValueType > reduce_type ;

  reduce_type * const block_value ;

  CudaParallelReduceFinalizeBlock( reduce_type * ptr )
    : block_value( ptr ) {}

  ~CudaParallelReduceFinalizeBlock() {}

  __device__
  void operator()() const
  {
    reduce_type * const shared_local =
      cuda_local_thread_reduce_shared_memory< ValueType >();

    reduce_type * const output = block_value + blockIdx.x ;

    // Output a single 'value' for the whole block.
    // Have multiple threads in the block work together to output
    // this thread-block's reduction value.
    // reduce_type is a union { value , storage[N] }

    for ( int i = threadIdx.x ; i < reduce_type::WordCount ; i += blockDim.x ) {
      output->storage[i] = shared_local->storage[i];
    }
  }

private:
  CudaParallelReduceFinalizeBlock();
  CudaParallelReduceFinalizeBlock & operator =
    ( const CudaParallelReduceFinalizeBlock & );
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct CudaParallelReduceFinalizeValue {

  typedef CudaSharedMemoryReduceType< ValueType > reduce_type ;

  ValueType * const value ;

  CudaParallelReduceFinalizeValue( ValueType * ptr )
    : value( ptr ) {}

  ~CudaParallelReduceFinalizeValue() {}

  __device__
  void operator()() const
  {
    reduce_type * const shared_local =
      cuda_local_thread_reduce_shared_memory< ValueType >();

    // Output a single 'value' for the whole block.
    // Have multiple threads in the block work together to output
    // this thread-block's reduction value.
    // reduce_type is a union { value , storage[N] }

    CudaWordType * const output = (CudaWordType *) value ;

    for ( int i = threadIdx.x ; i < reduce_type::WordCount ; i += blockDim.x ) {
      output[i] = shared_local->storage[i];
    }
  }

private:
  CudaParallelReduceFinalizeValue();
  CudaParallelReduceFinalizeValue & operator =
    ( const CudaParallelReduceFinalizeValue & );
};

//----------------------------------------------------------------------------

template< class FunctorType , class FinalizeType >
void cuda_parallel_reduce( const size_t work_count , const FunctorType & functor , const FinalizeType & finalize )
{
  typedef size_t                                        size_type ;
  typedef typename FunctorType::value_type              value_type ;
  typedef CudaSharedMemoryReduceType< value_type >      reduce_type ;
  typedef CudaParallelReduceFinalizeBlock< value_type > finalize_block_type ;

  // Size of result for shared memory
  const size_type reduce_size = sizeof(reduce_type);

  const size_type max_thread_count = CudaDevice::reduction_thread_max( reduce_size );

  if ( work_count < max_thread_count ) {

    // Small amount of work, use a single thread block

    size_type thread_count = max_thread_count ;

    // Reduce thread count until nearly every thread will have work
    while ( work_count <= ( thread_count >> 1 ) ) { thread_count >>= 1 ; }

    run_reduce_functor_on_cuda< FunctorType , FinalizeType >
      <<< 1, thread_count, reduce_size * thread_count >>>
      ( work_count , functor , finalize );
  }
  else {

    // Large amount of work, use multiple thread blocks
    // requiring partial reductions from each thread block
    // with a final reduction of the partial values.

    size_type block_count = CudaDevice::block_count_max();

    // Block count must be less than or equal to max_thread_count
    // so that the final reduction has one value per block.

    if ( max_thread_count < block_count ) { block_count = max_thread_count ; }

    // Reduce block count until nearly every block will have work.
    while ( work_count <= max_thread_count * ( block_count >> 1 ) ) { block_count >>= 1 ; }

    // Block reductions are stored in global memory
    finalize_block_type finalize_block(
      (reduce_type *) CudaDevice::allocate_memory( reduce_size, block_count, std::string() ) );

    run_reduce_functor_on_cuda< FunctorType , finalize_block_type >
      <<< block_count , max_thread_count , reduce_size * max_thread_count >>>
      ( work_count , functor , finalize_block );

    // Reduce block reduction partial values to a single value
    run_reduce_operator_on_cuda< FunctorType , FinalizeType >
      <<< 1 , block_count , reduce_size * block_count >>>
      ( finalize_block.block_value , finalize );

    CudaDevice::deallocate_memory( finalize_block.block_value );
  }
};

//----------------------------------------------------------------------------
/** \brief  Return reduce value */
template< class FunctorType >
struct ParallelReduce< FunctorType , void , CudaDevice >
{
  typedef typename FunctorType::value_type         value_type ;

  static value_type run( const size_t work_count , const FunctorType & functor )
  {
    typedef CudaParallelReduceFinalizeValue< value_type > finalize_functor ;

    value_type tmp ;

    ValueView< value_type , CudaDevice >
      view( create_value< value_type , CudaDevice >() );

    finalize_functor finalize( view.address_on_device() );

    cuda_parallel_reduce( work_count , functor , finalize );

    deep_copy( tmp , view );

    return tmp ;
  }
};

//----------------------------------------------------------------------------
/** \brief  Place reduce value in value view  */
template< class FunctorType >
struct ParallelReduce< FunctorType ,
                       ValueView< typename FunctorType::value_type , CudaDevice > ,
                       CudaDevice >
{
  typedef typename FunctorType::value_type      value_type ;
  typedef ValueView< value_type , CudaDevice >  view_type ;

  static void run( const size_t work_count , const FunctorType & functor , const view_type & view )
  {
    typedef CudaParallelReduceFinalizeValue< value_type > finalize_functor ;

    finalize_functor finalize( view.address_on_device() );

    cuda_parallel_reduce( work_count , functor , finalize );
  }
};

//----------------------------------------------------------------------------
/** \brief  Process reduce value via a finalize functor */
template< class FunctorType , class FinalizeType >
struct ParallelReduce< FunctorType , FinalizeType , CudaDevice >
{
  typedef typename FunctorType::value_type value_type ;

  static void run( const size_t work_count , const FunctorType & functor , const FinalizeType & finalize )
  {
    CudaParallelReduceFinalizeFunctor< value_type , FinalizeType >
      tmp( finalize );

    cuda_parallel_reduce( work_count , functor , tmp );
  }
};

//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_CUDADEVICEREDUCE_HPP */

