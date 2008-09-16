/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/

#ifndef util_SimpleArrayOps_hpp
#define util_SimpleArrayOps_hpp

namespace PHX {

/** \defgroup util_simple_array_ops  Simple C-Array Operations
 *  \author H. Carter Edwards  <hcedwar@sandia.gov>
 *  \date   November 2006
 *
 *  \brief  Basic operations for compile-time fixed length arrays
 *
 *  Example:  Sum<5>(x,y)   results in  x[i] += y[i]  for i=0..4
 *
 *  \{
 */

template< unsigned N , typename T > void Copy( T * , const T & );
template< unsigned N , typename T > void Copy( T * , const T * );
template< unsigned N , typename T > void Sum(  T * , const T * );
template< unsigned N , typename T > void Sum(  T * , const T & , const T * );
template< unsigned N , typename T > void Prod( T * , const T * );
template< unsigned N , typename T > void Min(  T * , const T * );
template< unsigned N , typename T > void Max(  T * , const T * );

template< unsigned N , typename T > void BitOr(  T * , const T * );
template< unsigned N , typename T > void BitAnd( T * , const T * );

template< unsigned N , typename T > T InnerProduct( const T * , const T * );

template< unsigned N , typename T > bool Equal(        const T * , const T * );
template< unsigned N , typename T > bool NotEqual(     const T * , const T * );
template< unsigned N , typename T > bool Less(         const T * , const T * );
template< unsigned N , typename T > bool LessEqual(    const T * , const T * );
template< unsigned N , typename T > bool Greater(      const T * , const T * );
template< unsigned N , typename T > bool GreaterEqual( const T * , const T * );

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace PHX {
namespace Op {

template< typename T, unsigned n , unsigned i = 0 > struct Copy ;
template< typename T, unsigned n , unsigned i = 0 > struct Sum ;
template< typename T, unsigned n , unsigned i = 0 > struct Prod ;
template< typename T, unsigned n , unsigned i = 0 > struct Max ;
template< typename T, unsigned n , unsigned i = 0 > struct Min ;
template< typename T, unsigned n , unsigned i = 0 > struct BitOr ;
template< typename T, unsigned n , unsigned i = 0 > struct BitAnd ;

template< typename T, unsigned n, unsigned i = 0 > struct InnerProduct ;
template< typename T, unsigned n, unsigned i = 0 > struct Equal ;
template< typename T, unsigned n, unsigned i = 0 > struct NotEqual ;
template< typename T, unsigned n, unsigned i = 0 > struct Less ;
template< typename T, unsigned n, unsigned i = 0 > struct LessEqual ;
template< typename T, unsigned n, unsigned i = 0 > struct Greater ;
template< typename T, unsigned n, unsigned i = 0 > struct GreaterEqual ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< typename T , unsigned n >
struct Copy<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}
  static void op( T * , const T & ) {}

  Copy() : ptr(0) {}
  Copy( T * p ) : ptr(p) {}
  Copy( const Copy & rhs ) : ptr( rhs.ptr ) {}
  Copy & operator = ( const Copy & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct Sum<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}
  static void op( T * , const T & , const T * ) {}

  Sum() : ptr(0) {}
  Sum( T * p ) : ptr(p) {}
  Sum( const Sum & rhs ) : ptr( rhs.ptr ) {}
  Sum & operator = ( const Sum & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct Prod<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}

  Prod() : ptr(0) {}
  Prod( T * p ) : ptr(p) {}
  Prod( const Prod & rhs ) : ptr( rhs.ptr ) {}
  Prod & operator = ( const Prod & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct Max<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}

  Max() : ptr(0) {}
  Max( T * p ) : ptr(p) {}
  Max( const Max & rhs ) : ptr( rhs.ptr ) {}
  Max & operator = ( const Max & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct Min<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}

  Min() : ptr(0) {}
  Min( T * p ) : ptr(p) {}
  Min( const Min & rhs ) : ptr( rhs.ptr ) {}
  Min & operator = ( const Min & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct BitOr<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}

  BitOr() : ptr(0) {}
  BitOr( T * p ) : ptr(p) {}
  BitOr( const BitOr & rhs ) : ptr( rhs.ptr ) {}
  BitOr & operator = ( const BitOr & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct BitAnd<T,n,n> {
  typedef T type ;
  enum { N = 0 };
  static void op( T * , const T * ) {}

  BitAnd() : ptr(0) {}
  BitAnd( T * p ) : ptr(p) {}
  BitAnd( const BitAnd & rhs ) : ptr( rhs.ptr ) {}
  BitAnd & operator = ( const BitAnd & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n >
struct Equal<T,n,n> {
  static bool op( const T * , const T * ) { return true ; }
};

template< typename T , unsigned n >
struct NotEqual<T,n,n> {
  static bool op( const T * , const T * ) { return false ; }
};

template< typename T , unsigned n >
struct InnerProduct<T,n,n> {
  static T op( const T * , const T * ) { return 0 ; }
};

template< typename T , unsigned n >
struct Less<T,n,n> {
  static bool op( const T * , const T * ) { return false ; }
};

template< typename T , unsigned n >
struct LessEqual<T,n,n> {
  static bool op( const T * , const T * ) { return true ; }
};

template< typename T , unsigned n >
struct Greater<T,n,n> {
  static bool op( const T * , const T * ) { return false ; }
};

template< typename T , unsigned n >
struct GreaterEqual<T,n,n> {
  static bool op( const T * , const T * ) { return true ; }
};

//----------------------------------------------------------------------

template< typename T , unsigned n , unsigned i >
struct Copy {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { dst[i] = src[i] ; Op::Copy<T,N,i+1>::op(dst,src); }

  static void op( T * dst , const T & src )
    { dst[i] = src ; Op::Copy<T,N,i+1>::op(dst,src); }

  Copy() : ptr(0) {}
  Copy( T * p ) : ptr(p) {}
  Copy( const Copy & rhs ) : ptr( rhs.ptr ) {}
  Copy & operator = ( const Copy & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n , unsigned i >
struct Sum {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { dst[i] += src[i] ; Op::Sum<T,N,i+1>::op(dst,src); }

  static void op( T * dst , const T & a , const T * src )
    { dst[i] += a * src[i] ; Op::Sum<T,N,i+1>::op(dst,a,src); }

  Sum() : ptr(0) {}
  Sum( T * p ) : ptr(p) {}
  Sum( const Sum & rhs ) : ptr( rhs.ptr ) {}
  Sum & operator = ( const Sum & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n , unsigned i >
struct Prod {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { dst[i] *= src[i] ; Op::Prod<T,N,i+1>::op(dst,src); }

  Prod() : ptr(0) {}
  Prod( T * p ) : ptr(p) {}
  Prod( const Prod & rhs ) : ptr( rhs.ptr ) {}
  Prod & operator = ( const Prod & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n , unsigned i >
struct BitOr {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { dst[i] |= src[i] ; Op::BitOr<T,N,i+1>::op(dst,src); }

  BitOr() : ptr(0) {}
  BitOr( T * p ) : ptr(p) {}
  BitOr( const BitOr & rhs ) : ptr( rhs.ptr ) {}
  BitOr & operator = ( const BitOr & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n , unsigned i >
struct BitAnd {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { dst[i] |= src[i] ; Op::BitAnd<T,N,i+1>::op(dst,src); }

  BitAnd() : ptr(0) {}
  BitAnd( T * p ) : ptr(p) {}
  BitAnd( const BitAnd & rhs ) : ptr( rhs.ptr ) {}
  BitAnd & operator = ( const BitAnd & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n , unsigned i >
struct Max {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { if ( dst[i] < src[i] ) { dst[i] = src[i] ; }
      Op::Max<T,N,i+1>::op(dst,src); }

  Max() : ptr(0) {}
  Max( T * p ) : ptr(p) {}
  Max( const Max & rhs ) : ptr( rhs.ptr ) {}
  Max & operator = ( const Max & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};

template< typename T , unsigned n , unsigned i >
struct Min {
  typedef T type ;
  enum { N = n };

  static void op( T * dst , const T * src )
    { if ( src[i] < dst[i] ) { dst[i] = src[i] ; }
      Op::Min<T,N,i+1>::op(dst,src); }

  Min() : ptr(0) {}
  Min( T * p ) : ptr(p) {}
  Min( const Min & rhs ) : ptr( rhs.ptr ) {}
  Min & operator = ( const Min & rhs ) { ptr = rhs.ptr ; return *this ; }

  T * ptr ;
};


template< typename T , unsigned n , unsigned i >
struct InnerProduct {
  static T op( const T * x , const T * y )
    { return x[i] * y[i] + Op::InnerProduct<T,n,i+1>::op( x , y ); }
};

template< typename T , unsigned n , unsigned i >
struct Equal {
  static bool op( const T * x , const T * y )
    { return x[i] == y[i] && Op::Equal<T,n,i+1>::op(x,y); }
};

template< typename T , unsigned n , unsigned i >
struct NotEqual {
  static bool op( const T * x , const T * y )
    { return x[i] != y[i] || Op::NotEqual<T,n,i+1>::op(x,y); }
};

template< typename T , unsigned n , unsigned i >
struct Less {
  static bool op( const T * const lhs , const T * const rhs )
    { return lhs[i] != rhs[i] ? lhs[i] < rhs[i]
                              : Op::Less<T,n,i+1>::op(lhs,rhs); }
};

template< typename T , unsigned n , unsigned i >
struct LessEqual {
  static bool op( const T * const lhs , const T * const rhs )
    { return lhs[i] != rhs[i] ? lhs[i] < rhs[i]
                              : Op::LessEqual<T,n,i+1>::op(lhs,rhs); }
};

template< typename T , unsigned n , unsigned i >
struct Greater {
  static bool op( const T * const lhs , const T * const rhs )
    { return lhs[i] != rhs[i] ? lhs[i] > rhs[i]
                              : Op::Greater<T,n,i+1>::op(lhs,rhs); }
};

template< typename T , unsigned n , unsigned i >
struct GreaterEqual {
  static bool op( const T * const lhs , const T * const rhs )
    { return lhs[i] != rhs[i] ? lhs[i] > rhs[i]
                              : Op::GreaterEqual<T,n,i+1>::op(lhs,rhs); }
};

} // namespace Op
} // namespace PHX

//-----------------------------------

namespace PHX {

template< unsigned N , typename T >
void Copy( T * dst , const T * src )
{ Op::Copy<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void Copy( T * dst , const T & src )
{ Op::Copy<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void Sum( T * dst , const T * src )
{ Op::Sum<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void Sum( T * dst , const T & a , const T * src )
{ Op::Sum<T,N,0>::op( dst , a , src ); }

template< unsigned N , typename T >
void Prod( T * dst , const T * src )
{ Op::Prod<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void Max( T * dst , const T * src )
{ Op::Max<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void Min( T * dst , const T * src )
{ Op::Min<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void BitOr( T * dst , const T * src )
{ Op::BitOr<T,N,0>::op( dst , src ); }

template< unsigned N , typename T >
void BitAnd( T * dst , const T * src )
{ Op::BitAnd<T,N,0>::op( dst , src ); }

//-----------------------------------

template< unsigned N , typename T >
Op::Copy<T,N> Copy( T * dst )
{ return Op::Copy<T,N,0>( dst ); }

template< unsigned N , typename T >
Op::Sum<T,N> Sum( T * dst )
{ return Op::Sum<T,N,0>( dst ); }

template< unsigned N , typename T >
Op::Prod<T,N> Prod( T * dst )
{ return Op::Prod<T,N,0>( dst ); }

template< unsigned N , typename T >
Op::Max<T,N> Max( T * dst )
{ return Op::Max<T,N,0>( dst ); }

template< unsigned N , typename T >
Op::Min<T,N> Min( T * dst )
{ return Op::Min<T,N,0>( dst ); }

template< unsigned N , typename T >
Op::BitOr<T,N> BitOr( T * dst )
{ return Op::BitOr<T,N,0>( dst ); }

template< unsigned N , typename T >
Op::BitAnd<T,N> BitAnd( T * dst )
{ return Op::BitAnd<T,N,0>( dst ); }

//-----------------------------------

template< unsigned N , typename T >
T InnerProduct( const T * x , const T * y )
{ return Op::InnerProduct<T,N,0>::op( x , y ); }

template< unsigned N , typename T >
bool Equal( const T * x , const T * y )
{ return Op::Equal<T,N,0>::op( x , y ); }

template< unsigned N , typename T >
bool NotEqual( const T * x , const T * y )
{ return Op::NotEqual<T,N,0>::op( x , y ); }

template< unsigned N , typename T >
bool Less( const T * x , const T * y )
{ return Op::Less<T,N,0>::op( x , y ); }

template< unsigned N , typename T >
bool LessEqual( const T * x , const T * y )
{ return Op::LessEqual<T,N,0>::op( x , y ); }

template< unsigned N , typename T >
bool Greater( const T * x , const T * y )
{ return Op::Greater<T,N,0>::op( x , y ); }

template< unsigned N , typename T >
bool GreaterEqual( const T * x , const T * y )
{ return Op::GreaterEqual<T,N,0>::op( x , y ); }

} // namespace PHX

//-----------------------------------

/** \} */

#endif

