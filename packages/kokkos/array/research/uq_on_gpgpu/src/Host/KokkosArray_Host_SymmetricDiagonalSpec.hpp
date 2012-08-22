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


#ifndef KOKKOSARRAY_HOST_SYMMETRICDIAGONALSPEC_HPP
#define KOKKOSARRAY_HOST_SYMMETRICDIAGONALSPEC_HPP

namespace KokkosArray {
namespace Impl {

template<>
class Multiply< SymmetricDiagonalSpec< Host > , void , void > {
public:
  typedef Host                    device_type ;
  typedef device_type::size_type  size_type ;
  typedef SymmetricDiagonalSpec< device_type > block_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const block_type  & block ,
                     const MatrixValue *       a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type dimension = block.dimension();
    const size_type dim_half  = ( dimension + 1 ) >> 1 ;

    // Multiply the main diagonal (first diagonal)
    for ( size_type j = 0 ; j < dimension ; ++j ) {
      y[j] += a[j] * x[j] ; // Contiguous access
    }

    // Multiply remaining full diagionals, each diagonal is accessed twice
    for ( size_type d = 1 ; d < dim_half ; ++d ) {
      size_type kx  = d ;
      size_type kxr = dimension - d ;

      a += dimension ; // next diagonal

      for ( size_type j = 0 ; j < dimension ; ++j ) {
        y[j] += a[j] * x[kx] + a[kxr] * x[kxr]; // Contiguous access
        if ( dimension == ++kx )  kx = 0 ;
        if ( dimension == ++kxr ) kxr = 0 ;
      }
    }

    // If even number of diagonals then the last diagonal is half-length
    if ( ! ( dimension & 01 ) ) {
      size_type kx = dim_half ;

      a += dimension ; // next diagonal

      for ( size_type j = 0 ; j < dim_half ; ++j , ++kx ) {
        y[j]  += a[j] * x[kx] ; // Contiguous access
        y[kx] += a[j] * x[j] ;  // Contiguous access
      }
    }
  }

  static size_type matrix_size( const block_type & block )
    { return block.matrix_size(); }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_HOST_SYMMETRICDIAGONALSPEC_HPP */

