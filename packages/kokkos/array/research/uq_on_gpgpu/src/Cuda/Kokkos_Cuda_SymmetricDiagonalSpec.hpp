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

#ifndef KOKKOS_CUDA_SYMMETRICDIAGONALSPEC_HPP
#define KOKKOS_CUDA_SYMMETRICDIAGONALSPEC_HPP

#include <Cuda/Kokkos_Cuda_Parallel.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class Multiply< SymmetricDiagonalSpec< Cuda > , void , void >
{
public:
  typedef Cuda                           device_type ;
  typedef device_type::size_type         size_type ;
  typedef SymmetricDiagonalSpec< Cuda >  block_type ;

  __host__
  static dim3 thread_block( const block_type & block )
  {
    const int y = ( cuda_internal_maximum_warp_count() * CudaTraits::WarpSize )
                / block.dimension();

    if ( y < 2 ) {
      throw std::runtime_error( std::string("Kokkos::Impl::Multiply< SymmetricDiagonalSpec<Cuda> > ERROR: block too large") );
    }
    return dim3( block.dimension() , y , 1 );
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

    if ( 0 == threadIdx.y ) {
      // Multiply the main diagonal (first diagonal)
      ix = threadIdx.x ;
      ia = threadIdx.x ;
      shX[ threadIdx.x ] = x[ threadIdx.x ]; // Load 'x'
    }
    else if ( blockDim.y == 1 + threadIdx.y && ! ( dimension & 01 ) ) {
      // If even number of diagonals then the last diagonal is half-length
      ix = threadIdx.x ;
      ia = threadIdx.x + dim_half * dimension ;

      if ( threadIdx.x < dim_half ) {
        ix += dim_half ;
      }
      else {
        ix -= dim_half ;
        ia -= dim_half ;
      }
    }

    __syncthreads(); // Wait for load of 'x'

    VectorValue y = ( ia < 0 ) ? 0 : shX[ ix ] * a[ ia ];

    for ( int d = 1 + threadIdx.y ; d < dim_half ; d += blockDim.y ) {

      ix = threadIdx.x + d ; if ( dimension <= ix ) ix -= dimension ;
      ia = threadIdx.x - d ; if ( ia < 0 ) ia += dimension ;

      const MatrixValue * const A = a + d * dimension ;

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
} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_SYMMETRICDIAGONALSPEC_HPP */
