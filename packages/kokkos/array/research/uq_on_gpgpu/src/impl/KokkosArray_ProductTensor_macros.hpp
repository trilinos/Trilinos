/*
//@HEADER
// ************************************************************************
// 
//           Kokkos: Node API and Parallel Node Kernels
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <impl/KokkosArray_ProductTensor_macros.hpp> without macros defined"

#else

namespace KokkosArray {

//----------------------------------------------------------------------------

template<>
class ProductTensorIndex< 3 , KOKKOSARRAY_MACRO_DEVICE > {
public:

  typedef KOKKOSARRAY_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  ~ProductTensorIndex() {}

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  ProductTensorIndex()
  {
    m_coord[0] = 0 ;
    m_coord[1] = 0 ;
    m_coord[2] = 0 ;
    m_coord[3] = 0 ;
  }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  ProductTensorIndex( const ProductTensorIndex & rhs )
  {
    m_coord[0] = rhs.m_coord[0] ;
    m_coord[1] = rhs.m_coord[1] ;
    m_coord[2] = rhs.m_coord[2] ;
    m_coord[3] = rhs.m_coord[3] ;
  }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
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
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
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
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
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

  template< class D >
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  ProductTensorIndex( const ProductTensorIndex<3,D> & rhs )
  {
    m_coord[0] = rhs.coord(0) ;
    m_coord[1] = rhs.coord(1) ;
    m_coord[2] = rhs.coord(2) ;
    m_coord[3] = rhs.offset();
  }

  template< class D >
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  ProductTensorIndex & operator = ( const ProductTensorIndex<3,D> & rhs )
  {
    m_coord[0] = rhs.coord(0) ;
    m_coord[1] = rhs.coord(1) ;
    m_coord[2] = rhs.coord(2) ;
    m_coord[3] = rhs.offset();
    return *this ;
  }

  /*------------------------------------------------------------------------*/
  /** \brief  Initialize with dense symmetric compression map offset */
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  explicit ProductTensorIndex( size_type dense_offset )
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

    m_coord[0] = size_type( i );

    dense_offset -= size_type( ( i * ( i + 1 ) * ( i + 2 ) ) / 6 );

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

    m_coord[1] = size_type( j );
    m_coord[2] = dense_offset - ( j * ( j + 1 ) ) / 2 ;
  }

  /*------------------------------------------------------------------------*/

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type coord( size_type c ) const { return m_coord[c]; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset() const { return m_coord[3]; }

  /*------------------------------------------------------------------------*/
  /** \brief  Increment in the dense map ordering */
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  ProductTensorIndex & operator++()
  {
    if      ( m_coord[1] > m_coord[2] ) { ++m_coord[2] ; }
    else if ( m_coord[0] > m_coord[1] ) { ++m_coord[1] ; m_coord[2] = 0 ; }
    else { ++m_coord[0] ; m_coord[1] = 0 ; m_coord[2] = 0 ; }
    ++m_coord[3] ;
    return *this ;
  }
  /*------------------------------------------------------------------------*/

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] == rhs.m_coord[3] ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator < ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] < rhs.m_coord[3] ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator > ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] > rhs.m_coord[3] ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] != rhs.m_coord[3] ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator <= ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] <= rhs.m_coord[3] ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator >= ( const ProductTensorIndex & rhs ) const
  { return m_coord[3] >= rhs.m_coord[3] ; }

  /*------------------------------------------------------------------------*/

private:

  size_type m_coord[4] ;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType >
class SparseProductTensor< 3 , ValueType , KOKKOSARRAY_MACRO_DEVICE > {
public:

  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;
  typedef ValueType               value_type ;

private:

  typedef KokkosArray::View< size_type[][3] , device_type >  map_type ;
  typedef KokkosArray::View< value_type[] ,   device_type >  vec_type ;

public:

  inline
  ~SparseProductTensor() {}

  inline
  SparseProductTensor() : m_coord() , m_value() , m_dimen(0) {}

  inline
  SparseProductTensor( const SparseProductTensor & rhs )
  : m_coord( rhs.m_coord ) , m_value( rhs.m_value ) , m_dimen( rhs.m_dimen ) {}

  inline
  SparseProductTensor & operator = ( const SparseProductTensor & rhs )
  {
    m_coord = rhs.m_coord ;
    m_value = rhs.m_value ;
    m_dimen = rhs.m_dimen ;
    return *this ;
  }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension() const { return m_dimen ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_count() const
  { return m_value.dimension_0(); }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  size_type coord( const size_type entry , const size_type c ) const
  { return m_coord( entry , c ); }

private:

  map_type   m_coord ;
  vec_type   m_value ;
  size_type  m_dimen ;

  template< class T , class I >
  friend class Impl::CreateSparseProductTensor ;
};

//----------------------------------------------------------------------------

template< typename ValueType >
class CrsProductTensor< 3 , ValueType , KOKKOSARRAY_MACRO_DEVICE > {
public:

  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;
  typedef ValueType               value_type ;

private:

  typedef KokkosArray::View< value_type[] , device_type >  vec_type ;

public:

  inline
  ~CrsProductTensor() {}

  inline
  CrsProductTensor() : m_coord() , m_value() , m_entry_max(0) {}

  inline
  CrsProductTensor( const CrsProductTensor & rhs )
  : m_coord( rhs.m_coord ) , m_value( rhs.m_value ) , m_entry_max( rhs.m_entry_max ) {}

  inline
  CrsProductTensor & operator = ( const CrsProductTensor & rhs )
  {
    m_coord = rhs.m_coord ;
    m_value = rhs.m_value ;
    m_entry_max = rhs.m_entry_max ;
    return *this ;
  }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension() const { return m_coord.row_map.dimension(0) - 1 ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_count() const
  { return m_coord.entries.dimension(0); }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type entry_maximum() const
  { return m_entry_max ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  size_type entry_begin( size_type i ) const
  { return m_coord.row_map[i]; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  size_type entry_end( size_type i ) const
  { return m_coord.row_map[i+1]; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  size_type coord( const size_type entry , const size_type c ) const
  { return m_coord.entries( entry , c ); }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

private:

  KokkosArray::CrsArray< size_type[2] , device_type >   m_coord ;
  KokkosArray::View< value_type[] , device_type >  m_value ;
  size_type                                        m_entry_max ;

  template< class T , class I >
  friend class Impl::CreateSparseProductTensor ;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace KokkosArray

#endif

