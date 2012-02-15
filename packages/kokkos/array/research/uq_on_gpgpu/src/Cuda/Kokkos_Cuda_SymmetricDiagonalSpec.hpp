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
  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;
  typedef SymmetricDiagonalSpec< Cuda > block_type ;

  template< typename VectorValue >
  __host__
  static size_type shmem_size( const block_type & block )
  {
    return 2 * sizeof(VectorValue) *
           Impl::CudaTraits::warp_align( block.dimension() );
  }

  __host__
  static size_type matrix_size( const block_type & block )
    { return block.matrix_size(); }

  // Required: block.dimension() == blockDim.x
  template< typename MatrixValue , typename VectorValue >
  __device__
  static void apply( const block_type  & block ,
                     const MatrixValue *       a ,
                     const VectorValue * const x ,
                           VectorValue & y )
  {
    VectorValue * const sh = kokkos_impl_cuda_shared_memory<VectorValue>();
    const size_type dimension     = block.dimension();
    const size_type shared_offset = Impl::CudaTraits::warp_align( dimension );

    // Number of full diagonals
    // If even number of diagonals then the last diagonal is half length
    const size_type dim_half = ( dimension + 1 ) >> 1 ;

    // Coalesced memory load of input vector
    sh[ threadIdx.x ] = x[ threadIdx.x ];

    // Multiply the main diagonal (first diagonal)
    y += a[ threadIdx.x ] * sh[ threadIdx.x ] ;

    // Multiply remaining full diagionals, each diagonal is accessed twice
    size_type kx  = threadIdx.x ;
    size_type kxr = threadIdx.x ;

    for ( size_type d = 1 ; d < dim_half ; ++d ) {

      a += dimension ; // next diagonal

      sh[ threadIdx.x + shared_offset ] = a[ threadIdx.x ];

      ++kx ;
      if ( dimension == kx ) kx = 0 ;
      if ( 0 == kxr ) kxr = dimension ;
      --kxr ;

      __syncthreads(); // Wait for matrix diagonal to load

      y += sh[ threadIdx.x + shared_offset ] * sh[ kx ];
      y += sh[ kxr         + shared_offset ] * sh[ kxr ];
    }

    // If even number of diagonals then the last diagonal is half-length
    if ( ! ( dimension & 01 ) ) {

      a += dimension ; // next diagonal

      ++kx ;
      if ( dimension == kx ) kx = 0 ;

      if ( threadIdx.x < dim_half ) {
        sh[ threadIdx.x + shared_offset ] = a[ threadIdx.x ];
        kxr = threadIdx.x ;
      }
      else {
        kxr = threadIdx.x - dim_half ;
      }

      __syncthreads(); // Wait for matrix diagonal to load

      y += sh[ kxr + shared_offset ] * sh[ kx ];
    }
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_SYMMETRICDIAGONALSPEC_HPP */
