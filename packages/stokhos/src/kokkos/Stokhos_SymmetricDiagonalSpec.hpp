// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SYMMETRIC_DIAGONAL_SPEC_HPP
#define STOKHOS_SYMMETRIC_DIAGONAL_SPEC_HPP

#include "Kokkos_CrsArray.hpp"

namespace Stokhos {

/** \brief  Symmetric diagonal storage for a dense matrix.
 *
 *  Block storage size = dimension * ( dimension + 1 ) / 2
 *
 *  Given block_dim then total_diagonal_count = 1 + dimension / 2
 *
 *  If dimension is even then the last diagonal is only half length.
 *
 *  { a11 , a22 , a33 , a44 , a55 , ... }
 *  { a12 , a23 , a34 , a45 , a56 , ... }
 *  { a13 , a24 , a35 , a46 , a57 , ... }
 *
 */
template< class ExecutionSpace >
class SymmetricDiagonalSpec {
public:

  typedef unsigned size_type;

  /** \brief  Dimension of vector block */
  KOKKOS_INLINE_FUNCTION
  unsigned dimension() const { return m_dimension ; }

  /** \brief  Storage location for the (row,column) entry */
  KOKKOS_INLINE_FUNCTION
  unsigned matrix_offset( const unsigned row , const unsigned column ) const
    {
      const int diag_count = 1 + ( m_dimension >> 1 );
      const int diag = (int) column - (int) row ;

      unsigned offset = 0 ;

      if ( ( 0 <= diag && diag < diag_count ) || ( diag <= - diag_count ) ) {
        offset = row + m_dimension * ( ( m_dimension + diag ) % m_dimension );
      }
      else {
        offset = column + m_dimension * ( ( m_dimension - diag ) % m_dimension );
      }

      return offset ;
    }

  /** \brief  Storage size for block coefficients */
  KOKKOS_INLINE_FUNCTION
  unsigned matrix_size() const
    { return ( m_dimension * ( m_dimension + 1 ) ) >> 1 ; }

  SymmetricDiagonalSpec()
    : m_dimension( 0 ) {}

  SymmetricDiagonalSpec( const SymmetricDiagonalSpec & rhs )
    : m_dimension( rhs.m_dimension ) {}

  SymmetricDiagonalSpec & operator =
    ( const SymmetricDiagonalSpec & rhs )
      { m_dimension = rhs.m_dimension ; return *this ; }

  explicit
  SymmetricDiagonalSpec( const unsigned dim )
    : m_dimension( dim ) {}

private:
  unsigned m_dimension ;
};

template < typename Device >
class BlockMultiply< SymmetricDiagonalSpec< Device > > {
public:
  typedef Device execution_space ;
  typedef typename execution_space::size_type size_type ;
  typedef SymmetricDiagonalSpec< execution_space > block_type ;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
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

  KOKKOS_INLINE_FUNCTION
  static size_type matrix_size( const block_type & block )
    { return block.matrix_size(); }
};

} // namespace Stokhos

#endif /* #ifndef STOKHOS_SYMMETRIC_DIAGONAL_SPEC_HPP */
