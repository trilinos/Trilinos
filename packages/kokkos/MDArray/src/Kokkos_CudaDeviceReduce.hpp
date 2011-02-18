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

#include <Kokkos_CudaDevice.hpp>

namespace Kokkos {

template< class FunctorType , class DeviceType > struct ParallelReduce ;

template< class FunctorType >
inline
void parallel_reduce( const FunctorType & functor ,
                      const typename FunctorType::reduce_type & result )
{
  ParallelReduce< FunctorType , typename FunctorType::device_type >::run( functor , result );
}

//----------------------------------------------------------------------------

template< class ReduceOperators >
__device__
void reduce_shared_on_cuda(
  const typename ReduceOperators::reduce_type & result )
{
  typedef typename ReduceOperators::reduce_type       reduce_type ;
  typedef typename ReduceOperators::reduce_value_type reduce_value_type ;
  
  reduce_value_type * const ptr = result.address_on_device();

  reduce_type source ; source.assign_on_device( result );

  for ( unsigned int j = blockDim.x ; j ; ) {
    j >>= 1 ;

    // Wait for contributing thread from a different half-warp
    if ( warpSize < j ) { __syncthreads(); }

    if ( threadIdx.x < j ) {

      // View for source data

      source.assign_on_device( ptr + j );

      ReduceOperators::join( result , source );
    }
  }
}

template< typename Scalar >
__device__
Scalar * get_shared_values_on_cuda()
{
  extern __shared__ Scalar shared_values[];
  return & shared_values[0] ;
}

template< class FunctorType >
__global__
void run_reduce_functor_on_cuda(
  const          FunctorType              functor ,
  const typename FunctorType::reduce_type block_result ,
  const unsigned int                      reduce_count )
{
  typedef typename FunctorType::reduce_type       reduce_type ;
  typedef typename FunctorType::reduce_value_type reduce_value_type ;

  reduce_value_type * const shared_values = get_shared_values_on_cuda<reduce_value_type>();

  reduce_type local_result ;

  local_result.assign_on_device( block_result );
  local_result.assign_on_device( shared_values + threadIdx.x , blockDim.x );

  FunctorType::init( local_result );

  const unsigned int work_count  = functor.work_count();
  const unsigned int work_stride = blockDim.x * gridDim.x ;

  for ( unsigned int iwork = threadIdx.x + blockDim.x * blockIdx.x ;
        iwork < work_count ; iwork += work_stride ) {
    functor( iwork , local_result );
  }

  reduce_shared_on_cuda< FunctorType >( local_result );

  // Output partial sum per block: offset = blockIdx , stride = gridDim.x
  // Copy from local result: stride = blockDim.x

  {
    reduce_value_type * const output =
      block_result.address_on_device() + blockIdx.x ;

    for ( int i = threadIdx.x ; i < reduce_count ; i += blockDim.x ) {
      output[ i * gridDim.x ] = shared_values[ i * blockDim.x ];
    }
  }
}

// Single block, reduce contributions
template< class ReduceOperators >
__global__
void run_reduce_operator_on_cuda(
  const typename ReduceOperators::reduce_type result ,
  const typename ReduceOperators::reduce_type block_result ,
  const unsigned int reduce_count )
{
  typedef typename ReduceOperators::reduce_type reduce_type ;
  typedef typename reduce_type::value_type      reduce_value_type ;

  reduce_value_type * const shared_values = get_shared_values_on_cuda<reduce_value_type>();

  {
    const reduce_value_type * const input_values = block_result.address_on_device();

    const unsigned int count = reduce_count * blockDim.x ;

    for ( unsigned int i = threadIdx.x ; i < count ; i += blockDim.x ) {
      shared_values[i] = input_values[i] ;
    }
  }

  reduce_type local_result ;

  local_result.assign_on_device( block_result );
  local_result.assign_on_device( shared_values + threadIdx.x , blockDim.x );

  reduce_shared_on_cuda< ReduceOperators >( local_result );

  { // Copy reduced result to output.
    reduce_value_type * const output = result.address_on_device();

    for ( int i = threadIdx.x ; i < reduce_count ; i += blockDim.x ) {
      output[i] = shared_values[ i * blockDim.x ] ;
    }
  }
}

template< class FunctorType >
struct ParallelReduce< FunctorType , CudaDevice > {
  typedef typename FunctorType::device_type       device_type ;
  typedef typename FunctorType::reduce_type       reduce_type ;
  typedef typename FunctorType::reduce_value_type reduce_value_type ;
  typedef typename device_type::size_type         size_type ;


  static void run( const FunctorType & functor ,
                   const typename FunctorType::reduce_type & result )
  {
    size_type reduce_count = 1 ;
    for ( size_type i = 0 ; i < result.rank() ; ++i ) {
      reduce_count *= result.dimension(i);
    }
    // Size of result for shared memory
    const size_type work_count  = functor.work_count();
    const size_type reduce_size = sizeof(reduce_value_type) * reduce_count ;
    size_type thread_count = device_type::reduction_thread_max( reduce_size );

    if ( work_count < thread_count ) {

      // Small amount of work, use a single thread block
      // Reduce thread_count below work_count so that
      // the number of reduction terms does not exceed the
      // partial contributions.

      while ( work_count < thread_count ) { thread_count >>= 1 ; }

      const size_type shmem_size = thread_count * reduce_size ;

      run_reduce_functor_on_cuda<FunctorType> <<< 1, thread_count, shmem_size >>>
        ( functor , result , reduce_count );
    }
    else {

      // Large amount of work, use multiple thread blocks
      // requiring partial reductions from each thread block
      // with a final reduction of the partial values.

      size_type block_count_max = device_type::block_count_max();
      size_type block_count = 1 ;
      while ( ( block_count << 1 ) <= block_count_max ) { block_count <<= 1 ; }

      if ( thread_count < block_count ) { block_count = thread_count ; }

      const size_type shmem_size = thread_count * reduce_size ;

      reduce_type block_result ;

      const std::string block_result_name("reduce_scratch_array");

      switch( result.rank() ) {
      case 1 :
        block_result = device_type::template create_scratch_mdarray<reduce_value_type>(
                         result.dimension(0) ,
                         block_result_name , block_count );
        break ;
      case 2 :
        block_result = device_type::template create_scratch_mdarray<reduce_value_type>(
                         result.dimension(0) ,
                         result.dimension(1) ,
                         block_result_name , block_count );
        break ;
      case 3 :
        block_result = device_type::template create_scratch_mdarray<reduce_value_type>(
                         result.dimension(0) ,
                         result.dimension(1) ,
                         result.dimension(2) ,
                         block_result_name , block_count );
        break ;
      case 4 :
        block_result = device_type::template create_scratch_mdarray<reduce_value_type>(
                         result.dimension(0) ,
                         result.dimension(1) ,
                         result.dimension(2) ,
                         result.dimension(3) ,
                         block_result_name , block_count );
        break ;
      default:
        break ; // throw std::runtime_error(...);
      }

      run_reduce_functor_on_cuda<FunctorType> <<< block_count , thread_count , shmem_size >>>
        ( functor , block_result , reduce_count );

      run_reduce_operator_on_cuda<FunctorType>
        <<< 1 , block_count , shmem_size >>>( result , block_result , reduce_count );
    }

    cudaThreadSynchronize();
  }
};


}

#endif /* KOKKOS_CUDADEVICEREDUCE_HPP */

