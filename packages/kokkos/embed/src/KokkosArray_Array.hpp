/*
//@HEADER
// ************************************************************************
//
//                             KokkosArray
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_ARRAY_HPP
#define KOKKOSARRAY_ARRAY_HPP

#include <ostream>
#include <KokkosArray_Macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename T , unsigned N , class Proxy = void >
class Array ;

//----------------------------------------------------------------------------
/** \brief  Print value of an array expression. */
template< class T , unsigned N , class Proxy >
std::ostream & operator <<
  ( std::ostream & s , const KokkosArray::Array<T,N,Proxy> & a )
{
  s << "{" ;
  for ( unsigned i = 0 ; i < a.size() ; ++i ) { s << " " << a[i] ; }
  s << " }" ;
  return s ;
}

//----------------------------------------------------------------------------
/** \brief  Print expression tree of an array expression. */
template< typename T , unsigned N , class Proxy >
void print_expression( std::ostream & s ,
                       const Array<T,N,Proxy> & exp )
{
  s << "Array<T," << N << ", " ;
  exp.print_expression( s );
  s << " >" << std::endl ;
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

#if defined( KOKKOSARRAY_EXPRESSION_CHECK )

void array_check_bounds_throw( const unsigned , const unsigned );

inline
void array_check_bounds( const unsigned i , const unsigned count )
  { array_check_bounds_throw( i , N ); }

#else

KOKKOSARRAY_INLINE_FUNCTION
void array_check_bounds( const unsigned , const unsigned ) {}

#endif

//----------------------------------------------------------------------------
/** \brief  A static size array */

template< typename T , unsigned N >
class Array< T , N , void >
{
private:

  T elems[ N ];

  template< class , unsigned , class > friend class Array ;

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return N ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const volatile { return N ; }

  //------------------------------------

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  value_type & operator[]( const iType & i )
    { array_check_bounds(i,N); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile value_type & operator[]( const iType & i ) volatile
    { array_check_bounds(i,N); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,N); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const volatile value_type & operator[]( const iType & i ) const volatile
    { array_check_bounds(i,N); return elems[i] ; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array() {}

  //------------------------------------
  // Copy constructors:

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; }

  /** \brief  Take the value once before performing assignment */
  KOKKOSARRAY_INLINE_FUNCTION
  Array( const value_type rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs ; }

  template< typename rT , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array< rT , N , P > & rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; }

  //------------------------------------
  // Assignment operators for non-volatile type:

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; return *this ; }

  /** \brief  Take the value once before performing assignment */
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const value_type rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs ; return *this ; }

  template< typename rT , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array< rT , N , P > & rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; return *this ; }

  template< typename rT , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const volatile Array< rT , N , P > & rhs )
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; return *this ; }

  //------------------------------------
  // Assignment operators for volatile type:

  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator = ( const Array & rhs ) volatile
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; return *this ; }

  /** \brief  Take the value once before performing assignment */
  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator = ( const value_type rhs ) volatile
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs ; return *this ; }

  template< typename rT , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator =
    ( const Array< rT , N , P > & rhs ) volatile
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; return *this ; }

  template< typename rT , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator =
    ( const volatile Array< rT , N , P > & rhs ) volatile
  { for (unsigned i = 0 ; i < N ; ++i) elems[i] = rhs[i] ; return *this ; }

  //------------------------------------

  void print_expression( std::ostream & s ) const
  { s << "[" << (void*) elems << "]" ; }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename T >
class Array< T , 0 , void >
{
private:

  T * const elems ;
  const unsigned count ;

  template< class , unsigned , class > friend class Array ;

public:

  typedef T      value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count ; }

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const volatile { return count ; }

  //------------------------------------

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  value_type & operator[]( const iType & i )
    { array_check_bounds(i,count); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile value_type & operator[]( const iType & i ) volatile
    { array_check_bounds(i,count); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const volatile value_type & operator[]( const iType & i ) const volatile
    { array_check_bounds(i,count); return elems[i] ; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array() : elems(0), count(0) {}

  KOKKOSARRAY_INLINE_FUNCTION
  ~Array() { if ( elems ) delete[] elems ; }

  //------------------------------------
  // Copy constructors:

  inline
  Array( const Array & rhs )
    : elems( new value_type[ rhs.count ] )
    , count( rhs.count )
    { for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ; }

  template< typename rT , unsigned rN , class P >
  inline
  Array( const Array< rT , rN , P > & rhs )
    : elems( rhs.size() ? new value_type[ rhs.size() ] : 0 )
    , count( rhs.size() )
    { for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const value_type rhs , const unsigned N )
    : elems( N ? new value_type[N] : 0 )
    , count(N)
    { for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs ; }

private:
  /*  Only allow this constructor if meaningful semantics can be defined */
  KOKKOSARRAY_INLINE_FUNCTION
  Array( const value_type rhs );
public:

  //------------------------------------
  // Assignment operators for non-volatile type:

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
    {
      array_check_bounds( count , rhs.count + 1 );
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ;
      return *this ;
    }

  /** \brief  Take the value once before performing assignment */
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const value_type rhs )
    {
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs ;
      return *this ;
    }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array< rT , rN , P > & rhs )
    {
      array_check_bounds( count , rhs.size() + 1 );
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ;
      return *this ;
    }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const volatile Array< rT , rN , P > & rhs )
    {
      array_check_bounds( count , rhs.size() + 1 );
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ;
      return *this ;
    }

  //------------------------------------
  // Assignment operators for volatile type:

  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator = ( const Array & rhs ) volatile
    {
      array_check_bounds( count , rhs.count + 1 );
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ;
      return *this ;
    }

  /** \brief  Take the value once before performing assignment */
  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator = ( const value_type rhs ) volatile
    {
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs ;
      return *this ;
    }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator = ( const Array< rT , rN , P > & rhs ) volatile
    {
      array_check_bounds( count , rhs.size() + 1 );
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ;
      return *this ;
    }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile Array & operator =
    ( const volatile Array< rT , rN , P > & rhs ) volatile
    {
      array_check_bounds( count , rhs.size() + 1 );
      for (unsigned i = 0 ; i < count ; ++i) elems[i] = rhs[i] ;
      return *this ;
    }
  
  //------------------------------------

  void print_expression( std::ostream & s ) const
  { s << "[" << (void*) elems << " : " << count << "]" ; }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_ArrayProxy.hpp>
#include <impl/KokkosArray_ArrayExp.hpp>

#endif /* #ifndef KOKKOSARRAY_ARRAY_HPP */

