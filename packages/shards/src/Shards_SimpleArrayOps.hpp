// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Shards_SimpleArrayOps_hpp
#define Shards_SimpleArrayOps_hpp

namespace shards {

/** \ingroup shards_package
 *  \defgroup  shards_package_simple_array_ops Simple Compile-time Sized Array Operations
 *
 *  \brief  Simple operations such as Copy<N>, Sum<N>, Max<N>, etc.
 *
 *  \author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  \{
 */

/** \brief  Copy into an array.
 *  \tparam  n  Number of members to copy.
 */
template< unsigned n , unsigned i = 0 >
struct Copy {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] = src[0..N-1] */
  template<typename T>
  inline Copy( T * const dst , const T * const src )
    { dst[I] = src[I] ; Copy<N-1,I+1>(dst,src); }

  /** \brief  dst[0..N-1] = src */
  template<typename T>
  inline Copy( T * const dst , const T src )
    { dst[I] = src ; Copy<N-1,I+1>(dst,src); }

  Copy() {}
};

/** \brief  Sum into an array.
 *  \tparam  n  Number of members to sum.
 */
template< unsigned n , unsigned i = 0 >
struct Sum {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] += src[0..N-1] */
  template<typename T>
  inline Sum( T * const dst , const T * const src )
    { dst[I] += src[I] ; Sum<N-1,I+1>(dst,src); }

  /** \brief  dst[0..N-1] += a * src[0..N-1] */
  template<typename T>
  inline Sum( T * const dst , const T a , const T * const src )
    { dst[I] += a * src[I] ; Sum<N-1,I+1>(dst,a,src); }

  Sum() {}
};

/** \brief   Scale into an array.
 *  \tparam  n  Number of members to multiply.
 */
template< unsigned n , unsigned i = 0 >
struct Prod {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] *= src[0..N-1] */
  template<typename T>
  inline Prod( T * const dst , const T * const src )
    { dst[I] *= src[I] ; Prod<N-1,I+1>(dst,src); }
  Prod() {}
};

/** \brief   Bitwise-or into an array.
 *  \tparam  n  Number of members to or.
 */
template< unsigned n , unsigned i = 0 >
struct BitOr {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] |= src[0..N-1] */
  template<typename T>
  inline BitOr( T * const dst , const T * const src )
    { dst[I] |= src[I] ; BitOr<N-1,I+1>(dst,src); }

  BitOr() {}
};

/** \brief   Bitwise-and into an array.
 *  \tparam  n  Number of members to and.
 */
template< unsigned n , unsigned i = 0 >
struct BitAnd {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] &= src[0..N-1] */
  template<typename T>
  inline BitAnd( T * const dst , const T * const src )
    { dst[I] &= src[I] ; BitAnd<N-1,I+1>(dst,src); }

  BitAnd() {}
};

/** \brief   Take maximum value of each member of two arrays.
 *  \tparam  n  Number of members to iterate.
 */
template< unsigned n , unsigned i = 0 >
struct Max {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] = max( dst[0..N-1] , src[0..N-1] ) */
  template<typename T>
  inline Max( T * const dst , const T * const src )
    { if ( dst[I] < src[I] ) { dst[I] = src[I] ; } Max<N-1,I+1>(dst,src); }

  Max() {}
};

/** \brief   Take minimum value of each member of two arrays.
 *  \tparam  n  Number of members to iterate.
 */
template< unsigned n , unsigned i = 0 >
struct Min {
  enum { N = n , I = i };

  /** \brief  dst[0..N-1] = min( dst[0..N-1] , src[0..N-1] ) */
  template<typename T>
  inline Min( T * const dst , const T * const src )
    { if ( src[I] < dst[I] ) { dst[I] = src[I] ; } Min<N-1,I+1>(dst,src); }

  Min() {}
};

/** \brief  Inner product of two arrays
 *  \tparam  n  Number of members to iterate.
 */
template< unsigned n , unsigned i = 0 >
struct InnerProduct {
  enum { N = n , I = i };

  /** \brief  value += sum[ k = 0..N-1 ]( x[k] * y[k] ) */
  template<typename T>
  inline InnerProduct( T & value , const T * const x , const T * const y )
    { value += x[I] * y[I] ; InnerProduct<N-1,I+1>( value , x , y ); }

  InnerProduct() {}
};

/** \brief  Lexicographical comparison of two arrays.
 *  \tparam  n  Number of members to iterate.
 */
template< unsigned n , unsigned i = 0 >
struct Compare {
  enum { N = n , I = i };

  /** \brief  All members are equal */
  template<typename T>
  inline static bool equal( const T * const x , const T * const y )
    { return x[I] == y[I] && Compare<N-1,I+1>::equal(x,y); }

  /** \brief  All members are not equal */
  template<typename T>
  inline static bool not_equal( const T * const x , const T * const y )
    { return x[I] != y[I] || Compare<N-1,I+1>::not_equal(x,y); }

  /** \brief  First non-equal members satisfy x[k] < y[k] */
  template<typename T>
  inline static bool less( const T * const x , const T * const y )
    {
      return x[I] != y[I] ? x[I] < y[I] : Compare<N-1,I+1>::less(x,y);
    }

  /** \brief  First non-equal members satisfies x[k] <= y[k] */
  template<typename T>
  inline static bool less_equal( const T * const x , const T * const y )
    {
      return x[I] != y[I] ? x[I] < y[I] : Compare<N-1,I+1>::less_equal(x,y);
    }

  /** \brief  First non-equal members satisfies x[k] > y[k] */
  template<typename T>
  inline static bool greater( const T * const x , const T * const y )
    {
      return x[I] != y[I] ? x[I] > y[I] : Compare<N-1,I+1>::greater(x,y);
    }

  /** \brief  First non-equal members satisfies x[k] >= y[k] */
  template<typename T>
  inline static bool greater_equal( const T * const x , const T * const y )
    {
      return x[I] != y[I] ? x[I] > y[I] : Compare<N-1,I+1>::greater_equal(x,y);
    }

  Compare() {}
};

/** \} */

//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

template<unsigned i>
struct Copy<0,i> {
  enum { N = 0 };
  Copy() {}
  template<typename T> inline Copy( T * const , const T * const ) {}
  template<typename T> inline Copy( T * const , const T ) {}
};

template<unsigned i>
struct Sum<0,i> {
  enum { N = 0 };
  Sum() {}
  template<typename T> inline Sum( T * const , const T * const ) {}
  template<typename T> inline Sum( T * const , const T , const T * const ) {}
};

template<unsigned i>
struct Prod<0,i> {
  enum { N = 0 };
  Prod() {}
  template<typename T> inline Prod( T * const , const T * const ) {}
};

template<unsigned i>
struct Max<0,i> {
  enum { N = 0 };
  Max() {}
  template<typename T> inline Max( T * const , const T * const ) {}
};

template<unsigned i>
struct Min<0,i> {
  enum { N = 0 };
  Min() {}
  template<typename T> inline Min( T * const , const T * const ) {}
};

template<unsigned i>
struct BitOr<0,i> {
  enum { N = 0 };
  BitOr() {}
  template<typename T> inline BitOr( T * const , const T * const ) {}
};

template<unsigned i>
struct BitAnd<0,i> {
  enum { N = 0 };
  BitAnd() {}
  template<typename T> inline BitAnd( T * const , const T * const ) {}
};

template<unsigned i>
struct InnerProduct<0,i> {
  enum { N = 0 };
  InnerProduct() {}
  template<typename T>
  inline InnerProduct( T & , const T * const , const T * const ) {}
};

template<unsigned i>
struct Compare<0,i> {
  enum { N = 0 };
  Compare() {}

  template<typename T>
  inline static bool equal( const T * const , const T * const )
    { return true ; }

  template<typename T>
  inline static bool not_equal( const T * const , const T * const )
    { return false ; }

  template<typename T>
  inline static bool less( const T * const , const T * const )
    { return false ; }

  template<typename T>
  inline static bool less_equal( const T * const , const T * const )
    { return true ; }

  template<typename T>
  inline static bool greater( const T * const , const T * const )
    { return false ; }

  template<typename T>
  inline static bool greater_equal( const T * const , const T * const )
    { return true ; }
};


#endif /* DOXYGEN_COMPILE */

} // namespace shards

#endif

