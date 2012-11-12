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

#ifndef KOKKOSARRAY_CUDA_SYMMETRICDIAGONALSPEC_HPP
#define KOKKOSARRAY_CUDA_SYMMETRICDIAGONALSPEC_HPP

#include <stdexcept>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template<>
class Multiply< SymmetricDiagonalSpec< KokkosArray::Cuda > , void , void >
{
public:
  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;
  typedef SymmetricDiagonalSpec< KokkosArray::Cuda >  block_type ;

  __host__
  static dim3 thread_block( const block_type & block )
  {
    const int d = block.dimension();
    const int y = ( cuda_internal_maximum_warp_count() * CudaTraits::WarpSize ) / d ;

    if ( 0 == y ) {
      throw std::runtime_error( std::string("KokkosArray::Impl::Multiply< SymmetricDiagonalSpec<Cuda> > ERROR: block too large") );
    }

    // dimension X #diagonals to concurrently process
    return dim3( d , std::min( y , ( 1 + d ) / 2 ) , 1 );
  }

  template< typename VectorValue >
  __host__
  static size_type shmem_size( const block_type & block )
  {
    const dim3 d = thread_block( block );

    return sizeof(VectorValue) * d.x * d.y ;
  }

  __host__
  static size_type matrix_size( const block_type & block )
    { return block.matrix_size(); }

  // Required: blockDim.x == block.dimension()
  // Required: blockDim.y > 1
  // Computing: Y[ threadIdx.x ]
  template< typename MatrixValue , typename VectorValue >
  __device__
  static VectorValue apply( const block_type  & block ,
                            const MatrixValue * const a ,
                            const VectorValue * const x )
  {
    const int  dimension = block.dimension();
    const int  dim_half  = ( dimension + 1 ) >> 1 ;

    VectorValue * const shX = kokkos_impl_cuda_shared_memory<VectorValue>();

    int ia = -1 ;
    int ix = -1 ;

    VectorValue y = 0 ;

    if ( 0 == threadIdx.y ) {
      // Load vector and multiply the main diagonal (first diagonal)

      shX[ threadIdx.x ] = x[ threadIdx.x ]; // Load 'x'

      y = shX[ threadIdx.x ] * a[ threadIdx.x ];
    }

    __syncthreads(); // Wait for load of 'x'

    if ( 0 == threadIdx.y && ! ( dimension & 01 ) ) {

      // If even number of diagonals then the last diagonal is half-length
      // and only used once.

      ix = threadIdx.x ;
      ia = threadIdx.x + dim_half * dimension ;

      if ( threadIdx.x < dim_half ) {
        ix += dim_half ;
      }
      else {
        ix -= dim_half ;
        ia -= dim_half ;
      }
      y += shX[ ix ] * a[ ia ];
    }

    //------------------------------------

    const int A_stride = blockDim.y * dimension ;

    int d = 1 + threadIdx.y ;

    const MatrixValue * A = a + d * dimension ;

    for ( ; d < dim_half ; d += blockDim.y , A += A_stride ) {

      ix = threadIdx.x + d ; if ( dimension <= ix ) ix -= dimension ;
      ia = threadIdx.x - d ; if ( ia < 0 ) ia += dimension ;

      // A 'threadIdx.y' group accesses the matrix diagonal
      // A[ 0 .. dimension - 1 ] twice in the following statement.
      // Spatial and temporal locality are excellent for L1 cache reuse
      // for the second access.

      y += shX[ ix ] * A[ threadIdx.x ] +
           shX[ ia ] * A[ ia ];
    }

    if ( 0 < threadIdx.y ) {
      shX[ threadIdx.x + threadIdx.y * dimension ] = y ;
    }

    __syncthreads(); // Wait for load of contributions to 'y'

    for ( ix = 1 ; ix < blockDim.y ; ++ix ) {
      y += shX[ threadIdx.x + ix * dimension ];
    }

    return y ;
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_CUDA_SYMMETRICDIAGONALSPEC_HPP */
