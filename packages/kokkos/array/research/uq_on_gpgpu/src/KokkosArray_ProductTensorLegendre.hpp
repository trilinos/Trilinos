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

#ifndef KOKKOSARRAY_PRODUCTTENSORLEGENDRE_HPP
#define KOKKOSARRAY_PRODUCTTENSORLEGENDRE_HPP

#include <cmath>
#include <utility>
#include <vector>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_ProductTensorIndex.hpp>
#include <KokkosArray_LegendrePolynomial.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Triple product tensor for the Legendre polynomial bases.
 *
 *  Hard-coded to an upper bound of P = 7.
 */
class TripleProductTensorLegendre {
private:
  enum { N = 8 /* MaximumPolyDegree + 1 */ };
  enum { NONZERO_COUNT = 36 };

  unsigned char m_map[ N * N * N ];
  float m_terms[ NONZERO_COUNT ];

  static inline
  unsigned offset( const unsigned i , const unsigned j , const unsigned k ) 
  { return ( i << 6 ) | ( j << 3 ) | k ; }

public:

  enum { MaximumPolyDegree = 7 };

  TripleProductTensorLegendre();

  /** \brief  Value of <i,j,k> tensor entry */
  KOKKOSARRAY_INLINE_FUNCTION
  float operator()( const unsigned i , const unsigned j , const unsigned k ) const
  {
    return m_terms[ m_map[ offset(i,j,k) ] ];
  }

  /** \brief  Value of Product_v < I[v] , J[v] , K[v] > */
  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  double operator()( const unsigned n ,
                     const iType * const I ,
                     const iType * const J ,
                     const iType * const K ) const
  {
    double val = 1 ;
#if 1
    for ( unsigned iv = 0 ; iv < n ; ++iv ) {
      val *= m_terms[ m_map[ offset(I[iv],J[iv],K[iv]) ] ];
    }
    return val ;
#else
    unsigned m ;
    for ( unsigned iv = 0 ;
          iv < n &&
          0 != ( m = m_map[ offset(I[iv],J[iv],K[iv]) ] ) ; ++iv ) {
      val *= m_terms[ m ];
    }
#endif
    return val ;
  }

  /** \brief  Value of Product_v < I[v] , J[v] , K[v] > is non-zero. */
  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  bool is_non_zero( const unsigned n ,
                    const iType * const I ,
                    const iType * const J ,
                    const iType * const K ) const
  {
    unsigned iv = 0 ;
    for ( ; iv < n && m_map[ offset(I[iv],J[iv],K[iv]) ] ; ++iv );
    return n == iv ;
  }
};

//----------------------------------------------------------------------------
/** \brief  Generate the bases for a set of variables of a given
 *          polynomial degree and maximum polynomial degree.
 *
 *  The combinatorial-Legendre bases are ordered according to
 *  the number of non-zeros in the triple product associated with
 *  the bases.  Ordering is first by diagonal count (largest to smallest)
 *  and then by off-diagonal count (largest to smallest).
 */

class TripleProductTensorLegendreCombinatorialEvaluation {
private:

  const TripleProductTensorLegendre  m_tensor ;
  std::vector< unsigned char >       m_bases_map ;
  unsigned                           m_variable_count ;
  unsigned                           m_bases_count ;
  unsigned                           m_bases_count_with_diagonal ;

public:

  TripleProductTensorLegendreCombinatorialEvaluation(
    const unsigned                  maximum_polynomial_degree ,
    const std::vector< unsigned > & variable_polynomial_degree );

  inline
  unsigned bases_count() const
  { return m_bases_count ; }

  inline
  unsigned bases_count_with_diagonal() const
  { return m_bases_count_with_diagonal ; }

  inline
  unsigned variable_count() const
  { return m_variable_count ; }

  /** \brief  The variable_count Legendre polynomial indices
   *          for this basis.
   */
  inline
  const unsigned char * operator[]( const unsigned i ) const
  { return & m_bases_map[ i * m_variable_count ]; }

  /** \brief  The triple product value for <i,j,k> */
  inline
  double operator()( const unsigned i , const unsigned j , const unsigned k ) const
  {
    return m_tensor( m_variable_count ,
                     & m_bases_map[ i * m_variable_count ] ,
                     & m_bases_map[ j * m_variable_count ] ,
                     & m_bases_map[ k * m_variable_count ] );
  }

  /** \brief  The triple product value for <i,j,k> is nonzero. */
  inline
  bool is_non_zero( const unsigned i , const unsigned j , const unsigned k ) const
  {
    return m_tensor.is_non_zero( m_variable_count ,
                                 & m_bases_map[ i * m_variable_count ] ,
                                 & m_bases_map[ j * m_variable_count ] ,
                                 & m_bases_map[ k * m_variable_count ] );
  }
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_PRODUCTTENSORLEGENDRE_HPP */


