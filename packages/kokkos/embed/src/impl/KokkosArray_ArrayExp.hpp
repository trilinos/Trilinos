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

#ifndef KOKKOSARRAY_ARRAYEXP_HPP
#define KOKKOSARRAY_ARRAYEXP_HPP

#include <math.h>
#include <stdlib.h>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class >
struct ERROR__cannot_assign_value_to_constant_expression ;

//----------------------------------------------------------------------------

template< typename TypeRHS , typename TypeLHS >
struct BinaryExpressionType ;

template< typename T >
struct BinaryExpressionType<T,T> { typedef T type ; };

#define KOKKOSARRAY_BINARYEXPRESSIONTYPE( TD , TS ) \
template<> struct BinaryExpressionType<TD,TS> { typedef TD type ; }; \
template<> struct BinaryExpressionType<TS,TD> { typedef TD type ; }; \

KOKKOSARRAY_BINARYEXPRESSIONTYPE(double,float)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(double,long)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(double,int)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(double,unsigned long)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(double,unsigned int)

KOKKOSARRAY_BINARYEXPRESSIONTYPE(float,long)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(float,int)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(float,unsigned long)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(float,unsigned int)

KOKKOSARRAY_BINARYEXPRESSIONTYPE(long,int)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(long,unsigned long)
KOKKOSARRAY_BINARYEXPRESSIONTYPE(long,unsigned int)

KOKKOSARRAY_BINARYEXPRESSIONTYPE(int,unsigned int)

#undef KOKKOSARRAY_BINARYEXPRESSIONTYPE

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Assigment operators:

#define KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR  +=
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR

#define KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR  -=
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR

#define KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR  *=
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR

#define KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR  /=
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS  ArrayProxyMinus
#define KOKKOSARRAY_ARRAY_UNARY_OPERATOR        -
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS  ArrayProxyPlus
#define KOKKOSARRAY_ARRAY_UNARY_OPERATOR        +
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS  ArrayProxyNot
#define KOKKOSARRAY_ARRAY_UNARY_OPERATOR        !
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS  ArrayProxyAddition
#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR        +
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS

#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS  ArrayProxySubtraction
#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR        -
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS

#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS  ArrayProxyMultiplication
#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR        *
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS

#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS  ArrayProxyDivision
#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR        /
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS

#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS  ArrayProxyModulo
#define KOKKOSARRAY_ARRAY_BINARY_OPERATOR        %
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR
#undef  KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename Type >
struct ArrayWeakOrdering
{
  enum Result { EQUAL , LESS , GREATER , NOT_ORDERED };

  template< class ArrayLHS , class ArrayRHS >
  KOKKOSARRAY_INLINE_FUNCTION static
  Result compare_array_array( const ArrayLHS & lhs , const ArrayRHS & rhs )
  {
    Result result = EQUAL ; // Assume equal for 0 == Count

    for ( unsigned i = 0 ; i < lhs.size() && result != NOT_ORDERED ; ++i ) {
      const Type L = lhs[i] ;
      const Type R = rhs[i] ;

      if      ( L < R ) { result = result != GREATER ? LESS : NOT_ORDERED ; }
      else if ( R < L ) { result = result != LESS ? GREATER : NOT_ORDERED ; }
    }

    return result ;
  }

  template< class ArrayLHS >
  KOKKOSARRAY_INLINE_FUNCTION static
  Result compare_array_value( const ArrayLHS & lhs , const Type R )
  {
    Result result = EQUAL ; // Assume equal for 0 == Count

    for ( unsigned i = 0 ; i < lhs.size() && result != NOT_ORDERED ; ++i ) {
      const Type L = lhs[i] ;

      if      ( L < R ) { result = result != GREATER ? LESS : NOT_ORDERED ; }
      else if ( R < L ) { result = result != LESS ? GREATER : NOT_ORDERED ; }
    }

    return result ;
  }

  template< class ArrayRHS >
  KOKKOSARRAY_INLINE_FUNCTION static
  Result compare_value_array( const Type L , const ArrayRHS & rhs )
  {
    Result result = EQUAL ; // Assume equal for 0 == Count

    for ( unsigned i = 0 ; i < rhs.size() && result != NOT_ORDERED ; ++i ) {
      const Type R = rhs[i] ;

      if      ( L < R ) { result = result != GREATER ? LESS : NOT_ORDERED ; }
      else if ( R < L ) { result = result != LESS ? GREATER : NOT_ORDERED ; }
    }

    return result ;
  }

  KOKKOSARRAY_INLINE_FUNCTION static
  bool equal( const Result result ) { return EQUAL == result ; }

  KOKKOSARRAY_INLINE_FUNCTION static
  bool not_equal( const Result result ) { return EQUAL != result ; }

  KOKKOSARRAY_INLINE_FUNCTION static
  bool less( const Result result ) { return LESS == result ; }

  KOKKOSARRAY_INLINE_FUNCTION static
  bool less_or_equal( const Result result )
    { return LESS == result || EQUAL == result ; }

  KOKKOSARRAY_INLINE_FUNCTION static
  bool greater( const Result result )
    { return GREATER == result ; }

  KOKKOSARRAY_INLINE_FUNCTION static
  bool greater_or_equal( const Result result )
    { return GREATER == result || EQUAL == result ; }
};

} // namespace Impl
} // namespace KokkosArray

#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR       ==
#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL  equal
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR

#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR       !=
#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL  not_equal
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR

#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR       <
#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL  less
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR

#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR       >
#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL  greater
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR

#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR       <=
#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL  less_or_equal
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR

#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR       >=
#define KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL  greater_or_equal
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyAbs
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         abs
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::abs
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyFabs
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         fabs
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::fabs
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxySqrt
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         sqrt
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::sqrt
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyExp
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         exp
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::exp
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyLog
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         log
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::log
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyLog10
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         log10
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::log10
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxySin
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         sin
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::sin
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyCos
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         cos
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::cos
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyTan
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         tan
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::tan
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyAcos
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         acos
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::acos
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyAsin
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         asin
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::asin
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS   ArrayProxyAtan
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION         atan
#define KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER  ::atan
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename T >
KOKKOSARRAY_INLINE_FUNCTION
const T & min( const T & lhs , const T & rhs )
{ return lhs < rhs ? lhs : rhs ; }

template< typename T >
KOKKOSARRAY_INLINE_FUNCTION
const T & max( const T & lhs , const T & rhs )
{ return lhs < rhs ? rhs : lhs ; }

}

#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS   ArrayProxyMin
#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION         min
#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER  KokkosArray::min
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS   ArrayProxyMax
#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION         max
#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER  KokkosArray::max
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS

#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS   ArrayProxyAtan2
#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION         atan2
#define KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER  ::atan2
#include <impl/KokkosArray_ArrayExp_macros.hpp>
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION
#undef  KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_ARRAY_HPP */

