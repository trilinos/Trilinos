/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_SYMMETRICDIAGONALSPEC_HPP
#define KOKKOSARRAY_SYMMETRICDIAGONALSPEC_HPP

#include <KokkosArray_CrsArray.hpp>

namespace KokkosArray {

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
template< class DeviceType >
class SymmetricDiagonalSpec {
public:

  /** \brief  Dimension of vector block */
  KOKKOSARRAY_INLINE_FUNCTION
  unsigned dimension() const { return m_dimension ; }

  /** \brief  Storage location for the (row,column) entry */
  KOKKOSARRAY_INLINE_FUNCTION
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
  KOKKOSARRAY_INLINE_FUNCTION
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

} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_SYMMETRICDIAGONALSPEC_HPP */


