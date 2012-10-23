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

#ifndef KOKKOSARRAY_ARRAYPROXY_HPP
#define KOKKOSARRAY_ARRAYPROXY_HPP

#include <ostream>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_Array.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

struct ArrayProxyContiguous {};
struct ArrayProxyStrided {};
struct ArrayProxyValue {};

//----------------------------------------------------------------------------
// Proxy for a proxy-array is itself:

template< class ProxyType >
struct ArrayProxy { typedef ProxyType type ; };

// Proxy for an array value is a pointer to that data with
// the cv-qualifier carried on the proxy type.

template<>
struct ArrayProxy< void >
{ typedef const ArrayProxyContiguous type ; };

template<>
struct ArrayProxy< volatile void >
{ typedef const volatile ArrayProxyContiguous type ; };

//----------------------------------------------------------------------------

template< unsigned Count >
struct ArrayCount
{
  static const unsigned N = Count ;

  KOKKOSARRAY_INLINE_FUNCTION
  explicit ArrayCount( unsigned ) {}
};

template<>
struct ArrayCount< 0u >
{
  const unsigned N ;

  KOKKOSARRAY_INLINE_FUNCTION
  explicit ArrayCount( unsigned c ) : N( c ) {}
};

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename Type , unsigned N >
class Array< Type , N , ArrayProxyValue >
{
private:

  const Type value ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const          Type & arg ) : value( arg ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const volatile Type & arg ) : value( arg ) {}

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const Type & operator[]( const iType & ) const { return value ; }

  void print_expression( std::ostream & s ) const
    { s << value ; }
};

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

// four cv-qualified possibilities:
//
//   Array<T,C,                 ArrayProxyContiguous >
//   Array<T,C,       volatile  ArrayProxyContiguous >
//   Array<T,C, const           ArrayProxyContiguous >
//   Array<T,C, const volatile  ArrayProxyContiguous >
//

/** \brief  Static array reference to a contiguous
 *          to a chunk of memory.
 */
template< typename T , unsigned ArgN >
class Array< T , ArgN , ArrayProxyContiguous >
{
private:

  T * const elems ;

  ArrayCount<ArgN> count ;

  Array();

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  value_type & operator[]( const iType & i )
    { array_check_bounds(i,count.N); return elems[i]; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i]; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( value_type * rhs )
    : elems( rhs ), count( ArgN ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( value_type * rhs , unsigned c )
    : elems( rhs ), count(c) {}

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs[i] ;
    return *this ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const value_type rhs )
  { for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs ; return *this ; }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array< rT , rN , P > & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs[i] ;
    return *this ;
  }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const volatile Array< rT , rN , P > & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs[i] ;
    return *this ;
  }

  //------------------------------------

  void print_expression( std::ostream & s ) const
  { s << "[" << (void*) elems << "]" ; }
};

template< typename T , unsigned ArgN >
class Array< T , ArgN , volatile ArrayProxyContiguous >
{
private:

  volatile T * const elems ;

  ArrayCount<ArgN> count ;

  Array();

public:

  typedef T  value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile value_type & operator[]( const iType & i )
    { array_check_bounds(i,count.N); return elems[i] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const volatile value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i] ; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.size() ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( volatile value_type * rhs )
    : elems( rhs ), count( ArgN ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( volatile value_type * rhs , unsigned c )
    : elems( rhs ), count( c ) {}

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs[i] ;
    return *this ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const value_type rhs )
  { for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs ; return *this ; }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array< rT , rN , P > & rhs )
  { for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs[i] ; return *this ; }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const volatile Array< rT , rN , P > & rhs )
  { for (unsigned i = 0 ; i < count.N ; ++i) elems[i] = rhs[i] ; return *this ; }

  //------------------------------------

  void print_expression( std::ostream & s ) const
  { s << "[V " << (void*) elems << "]" ; }
};

template< typename T , unsigned ArgN >
class Array< T , ArgN , const ArrayProxyContiguous >
{
private:

  const T * const elems ;

  ArrayCount<ArgN> count ;

  Array();
  Array & operator = ( const Array & );

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( const value_type * rhs ) : elems( rhs ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const value_type * rhs , unsigned c )
    : elems( rhs ), count(c) {}

  template< unsigned rN >
  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( const Array<T,rN,void> & rhs )
    : elems( rhs.elems ), count( rhs.size() ) {}

  //------------------------------------

  void print_expression( std::ostream & s ) const
    { s << "[C " << (void*) elems << "]" ; }
};

template< typename T , unsigned ArgN >
class Array< T , ArgN , const volatile ArrayProxyContiguous >
{
private:

  const volatile T * const elems ;

  ArrayCount<ArgN> count ;

  Array();
  Array & operator = ( const Array & );

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const volatile value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( const volatile value_type * rhs )
   : elems( rhs ), count( ArgN ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( const volatile value_type * rhs , unsigned c )
   : elems( rhs ), count( c ) {}

  template< unsigned rN >
  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  Array( const volatile Array<T,rN,void> & rhs )
    : elems( rhs.elems ), count( rhs.size() ) {}

  //------------------------------------

  void print_expression( std::ostream & s ) const
  { s << "[CV " << (void*) elems << "]" ; }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

// four cv-qualified possibilities:
//
//   Array<T,C,                 ArrayProxyStrided >
//   Array<T,C,       volatile  ArrayProxyStrided >
//   Array<T,C, const           ArrayProxyStrided >
//   Array<T,C, const volatile  ArrayProxyStrided >
//

template< typename T , unsigned ArgN >
class Array< T , ArgN , ArrayProxyStrided >
{
private:

  T * const elems ;
  ArrayCount<ArgN>  count ;
  const unsigned    stride ;

  Array();

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  //------------------------------------

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  value_type & operator[]( const iType & i )
    { array_check_bounds(i,count.N); return elems[i*stride] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i*stride] ; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N ), stride( rhs.stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( value_type * rhs_pointer , unsigned rhs_stride )
    : elems( rhs_pointer ), count( ArgN ), stride( rhs_stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( value_type * rhs_pointer , unsigned rhs_stride , unsigned rhs_count )
    : elems( rhs_pointer ), count( rhs_count ), stride( rhs_stride ) {}

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs[i] ;
    return *this ;
  }


  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const value_type rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs ;
    return *this ;
  }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array< rT , rN , P > & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs[i] ;
    return *this ;
  }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const volatile Array< rT , rN , P > & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs[i] ;
    return *this ;
  }

  void print_expression( std::ostream & s ) const
  { s << "[" << (void*) elems << " : " << stride << "]" ; }
};

template< typename T , unsigned ArgN >
class Array< T , ArgN , volatile ArrayProxyStrided >
{
private:

  volatile T * const elems ;
  ArrayCount<ArgN>   count ;
  const unsigned     stride ;

  Array();

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  //------------------------------------

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  volatile value_type & operator[]( const iType & i )
    { array_check_bounds(i,count.N); return elems[i*stride] ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const volatile value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i*stride] ; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N ), stride( rhs.stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( volatile value_type * rhs_pointer , unsigned rhs_stride )
    : elems( rhs_pointer ), count( ArgN ), stride( rhs_stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( volatile value_type * rhs_pointer , unsigned rhs_stride , unsigned rhs_count )
    : elems( rhs_pointer ), count( rhs_count ), stride( rhs_stride ) {}

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs[i] ;
    return *this ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const value_type rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs ;
    return *this ;
  }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const Array< rT , rN , P > & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs[i] ;
    return *this ;
  }

  template< typename rT , unsigned rN , class P >
  KOKKOSARRAY_INLINE_FUNCTION
  Array & operator = ( const volatile Array< rT , rN , P > & rhs )
  {
    for (unsigned i = 0 ; i < count.N ; ++i) elems[i*stride] = rhs[i] ;
    return *this ;
  }

  void print_expression( std::ostream & s ) const
  { s << "[V " << (void*) elems << " : " << stride << "]" ; }
};

template< typename T , unsigned ArgN >
class Array< T , ArgN , const ArrayProxyStrided >
{
private:

  const T * const  elems ;
  ArrayCount<ArgN> count ;
  const unsigned   stride ;

  Array();

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i*stride] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N), stride( rhs.stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const value_type * rhs , const unsigned stride )
    : elems( rhs ), count( ArgN ), stride( stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const value_type * rhs_pointer , unsigned rhs_stride , unsigned rhs_count )
    : elems( rhs_pointer ), count( rhs_count ), stride( rhs_stride ) {}

  void print_expression( std::ostream & s ) const
  { s << "[C " << (void*) elems << " : " << stride << "]" ; }
};

template< typename T , unsigned ArgN >
class Array< T , ArgN , const volatile ArrayProxyStrided >
{
private:

  const volatile T * const elems ;
  ArrayCount<ArgN>         count ;
  const unsigned           stride ;

  Array();

public:

  typedef T value_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return count.N ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  const volatile value_type & operator[]( const iType & i ) const
    { array_check_bounds(i,count.N); return elems[i*stride] ; }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array & rhs )
    : elems( rhs.elems ), count( rhs.count.N ), stride( rhs.stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const volatile value_type * rhs , const unsigned stride )
    : elems( rhs ), count( ArgN ), stride( stride ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const volatile value_type * rhs_pointer , unsigned rhs_stride , unsigned rhs_count )
    : elems( rhs_pointer ), count( rhs_count ), stride( rhs_stride ) {}

  void print_expression( std::ostream & s ) const
  { s << "[CV " << (void*) elems << " : " << stride << "]" ; }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_ARRAYPROXY_HPP */

