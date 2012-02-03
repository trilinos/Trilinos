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
//----------------------------------------------------------------------------

namespace Kokkos {

template<>
class SymmetricDiagonalSpec< Cuda > {
public:
  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  inline
  __device__ __host__
  size_type dimension() const { return m_dimension ; }

  inline
  __device__ __host__
  size_type size() const
    { return ( m_dimension * ( m_dimension + 1 ) ) >> 1 ; }

  __device__ __host__
  size_type offset( const size_type row , const size_type column ) const
    {
      const int diag_count = 1 + ( m_dimension >> 1 );
      const int diag = (int) column - (int) row ;

      size_type offset = 0 ;

      if ( ( 0 <= diag && diag < diag_count ) || ( diag <= - diag_count ) ) {
        offset = row + m_dimension * ( ( m_dimension + diag ) % m_dimension );
      }
      else {
        offset = column + m_dimension * ( ( m_dimension - diag ) % m_dimension );
      }

      return offset ;
    }

  //--------------------------------------------------------------------------

  template< typename VectorValue >
  __host__
  size_type shmem_size() const
    { return 2 * sizeof(VectorValue) * Impl::CudaTraits::warp_align( m_dimension ); }
  

  // Required: m_dimension == blockDim.x
  template< typename MatrixValue , typename VectorValue >
  __device__
  void multiply( const MatrixValue *       a ,
                 const VectorValue * const x ,
                       VectorValue & y ) const
  {
    VectorValue * const sh = kokkos_impl_cuda_shared_memory<VectorValue>();

    const size_type shared_offset = Impl::CudaTraits::warp_align( m_dimension );

    // Number of full diagonals
    // If even number of diagonals then the last diagonal is half length
    const size_type dim_half = ( m_dimension + 1 ) >> 1 ;

    // Coalesced memory load of input vector
    sh[ threadIdx.x ] = x[ threadIdx.x ];

    // Multiply the main diagonal (first diagonal)
    y += a[ threadIdx.x ] * sh[ threadIdx.x ] ;

    // Multiply remaining full diagionals, each diagonal is accessed twice
    size_type kx  = threadIdx.x ;
    size_type kxr = threadIdx.x ;

    for ( size_type d = 1 ; d < dim_half ; ++d ) {

      a += m_dimension ; // next diagonal

      sh[ threadIdx.x + shared_offset ] = a[ threadIdx.x ];

      ++kx ;
      if ( m_dimension == kx ) kx = 0 ;
      if ( 0 == kxr ) kxr = m_dimension ;
      --kxr ;

      __syncthreads(); // Wait for matrix diagonal to load

      y += sh[ threadIdx.x + shared_offset ] * sh[ kx ];
      y += sh[ kxr         + shared_offset ] * sh[ kxr ];
    }

    // If even number of diagonals then the last diagonal is half-length
    if ( ! ( m_dimension & 01 ) ) {

      a += m_dimension ; // next diagonal

      ++kx ;
      if ( m_dimension == kx ) kx = 0 ;

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

  //--------------------------------------------------------------------------

  SymmetricDiagonalSpec()
    : m_dimension( 0 ) {}

  SymmetricDiagonalSpec( const SymmetricDiagonalSpec & rhs )
    : m_dimension( rhs.m_dimension ) {}

  SymmetricDiagonalSpec & operator =
    ( const SymmetricDiagonalSpec & rhs )
      { m_dimension = rhs.m_dimension ; return *this ; }

  explicit
  SymmetricDiagonalSpec( const size_type dim )
    : m_dimension( dim ) {}

private:
  size_type m_dimension ;
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_SYMMETRICDIAGONALSPEC_HPP */
