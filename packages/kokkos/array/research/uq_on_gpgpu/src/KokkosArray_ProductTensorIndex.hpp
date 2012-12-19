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

#ifndef KOKKOSARRAY_PRODUCTTENSORINDEX_HPP
#define KOKKOSARRAY_PRODUCTTENSORINDEX_HPP

#include <KokkosArray_Macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

/** \brief  Use symmetry to compress the product tensor index space.
 *
 *  By symmetry the coordinates of a particular index can be
 *  arbitrarily reordered.  Compression due to symmetry results
 *  in { D! / ( (D-Rank)! * Rank! ) } unique entries in the index space.
 *
 *  Indices are lexigraphically ordered.
 *  coord(0) >= coord(1) >= ... >= coord(Rank-1)
 */
template< unsigned Rank >
class ProductTensorIndex {
public:

  ProductTensorIndex();
  ProductTensorIndex( const ProductTensorIndex & );
  ProductTensorIndex & operator = ( const ProductTensorIndex & );

  explicit ProductTensorIndex( unsigned offset );

  // ProductTensorIndex( unsigned coord0 , unsigned coord1 , ... );

  ProductTensorIndex & operator ++ ();

  /** \brief  Coordinate 'c' where 0 <= c < Rank */
  unsigned coord( unsigned c ) const ;

  /** \brief  Offset of this entry in the index space */
  unsigned offset() const ;

  bool operator == ( const ProductTensorIndex & ) const ;
  bool operator != ( const ProductTensorIndex & ) const ;
  bool operator <  ( const ProductTensorIndex & ) const ;
  bool operator <= ( const ProductTensorIndex & ) const ;
  bool operator >  ( const ProductTensorIndex & ) const ;
  bool operator >= ( const ProductTensorIndex & ) const ;
};


template<>
class ProductTensorIndex< 3 > {
private:

  unsigned m_coord[4] ;

public:

  KOKKOSARRAY_INLINE_FUNCTION
  ~ProductTensorIndex() {}

  KOKKOSARRAY_INLINE_FUNCTION
  ProductTensorIndex()
  {
    m_coord[0] = 0 ;
    m_coord[1] = 0 ;
    m_coord[2] = 0 ;
    m_coord[3] = 0 ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  ProductTensorIndex( const ProductTensorIndex & rhs )
  {
    m_coord[0] = rhs.m_coord[0] ;
    m_coord[1] = rhs.m_coord[1] ;
    m_coord[2] = rhs.m_coord[2] ;
    m_coord[3] = rhs.m_coord[3] ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  ProductTensorIndex & operator = ( const ProductTensorIndex & rhs )
  {
    m_coord[0] = rhs.m_coord[0] ;
    m_coord[1] = rhs.m_coord[1] ;
    m_coord[2] = rhs.m_coord[2] ;
    m_coord[3] = rhs.m_coord[3] ;
    return *this ;
  }

private:

  template< unsigned I , unsigned J >
  KOKKOSARRAY_INLINE_FUNCTION
  void order()
  {
    if ( m_coord[I] < m_coord[J] ) {
      m_coord[3] = m_coord[I] ;
      m_coord[I] = m_coord[J] ;
      m_coord[J] = m_coord[3] ;
    }
  }

public:

  /** \brief  Construct with the dense symmetric compression map */
  template< typename iType , typename jType , typename kType >
  KOKKOSARRAY_INLINE_FUNCTION
  ProductTensorIndex( const iType & argI ,
                      const jType & argJ ,
                      const kType & argK )
  {
    m_coord[0] = argI ;
    m_coord[1] = argJ ;
    m_coord[2] = argK ;

    // Sort the indices I >= J >= K
    order<0,2>();
    order<1,2>();
    order<0,1>();

    m_coord[3] = ( m_coord[0] * ( m_coord[0] + 1 ) * ( m_coord[0] + 2 ) ) / 6 +
                 ( m_coord[1] * ( m_coord[1] + 1 ) ) / 2 + m_coord[2] ;
  }

  /*------------------------------------------------------------------------*/
  /** \brief  Initialize with dense symmetric compression map offset */
  KOKKOSARRAY_INLINE_FUNCTION
  explicit ProductTensorIndex( unsigned dense_offset )
  {
    unsigned long i , j ;

    m_coord[3] = dense_offset ;

    {
      const unsigned long n6 = 6 * dense_offset ;
      // Upper and lower bound for 'i'.  Arbitrary start at '8'.
      for ( j = 8 ; j * ( j + 1 ) * ( j + 2 ) < n6 ; j <<= 1 );
      for ( i = j >> 1 ; n6 < i * ( i + 1 ) * ( i + 2 ) ; i >>= 1 );
      // Binary search [i,j)
      j -= i ;
      while ( 0 < j ) {
        const unsigned long half = j >> 1 ;
        const unsigned long mid  = i + half + 1 ;
        if ( n6 < mid * ( mid + 1 ) * ( mid + 2 ) ) {
          j = half ;
        }
        else {
          i = mid ;
          j -= half + 1 ;
        }
      }
    }

    m_coord[0] = unsigned( i );

    dense_offset -= unsigned( ( i * ( i + 1 ) * ( i + 2 ) ) / 6 );

    {
      const unsigned long n2 = 2 * dense_offset ;

      // Binary search [0,i)
      j = 0 ;
      while ( 0 < i ) {
        const unsigned long half = i >> 1 ;
        const unsigned long mid  = j + half + 1 ;
        if ( n2 < mid * ( mid + 1 ) ) {
          i = half ;
        }
        else {
          j = mid ;
          i -= half + 1 ;
        }
      }
    }

    m_coord[1] = unsigned( j );
    m_coord[2] = dense_offset - ( j * ( j + 1 ) ) / 2 ;
  }

  /*------------------------------------------------------------------------*/

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned coord( unsigned c ) const { return m_coord[c]; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned offset() const { return m_coord[3]; }

  /*------------------------------------------------------------------------*/
  /** \brief  Increment in the dense map ordering */
  KOKKOSARRAY_INLINE_FUNCTION
  ProductTensorIndex & operator++()
  {
    if      ( m_coord[1] > m_coord[2] ) { ++m_coord[2] ; }
    else if ( m_coord[0] > m_coord[1] ) { ++m_coord[1] ; m_coord[2] = 0 ; }
    else { ++m_coord[0] ; m_coord[1] = 0 ; m_coord[2] = 0 ; }
    ++m_coord[3] ;
    return *this ;
  }
  /*------------------------------------------------------------------------*/

  KOKKOSARRAY_INLINE_FUNCTION
  bool operator == ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] == rhs.m_coord[3] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool operator < ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] < rhs.m_coord[3] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool operator > ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] > rhs.m_coord[3] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool operator != ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] != rhs.m_coord[3] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool operator <= ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] <= rhs.m_coord[3] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool operator >= ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] >= rhs.m_coord[3] ; }

  /*------------------------------------------------------------------------*/
};

} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_PRODUCTTENSORINDEX_HPP */


