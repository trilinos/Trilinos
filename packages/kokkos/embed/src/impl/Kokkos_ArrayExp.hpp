/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
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

#ifndef KOKKOS_COREEXP_HPP
#define KOKKOS_COREEXP_HPP

#include <math.h>
#include <stdlib.h>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class >
struct ERROR__cannot_assign_value_to_constant_expression ;

//----------------------------------------------------------------------------

template< typename TypeRHS , typename TypeLHS >
struct BinaryExpressionType ;

template< typename T >
struct BinaryExpressionType<T,T> { typedef T type ; };

#define KOKKOS_BINARYEXPRESSIONTYPE( TD , TS ) \
template<> struct BinaryExpressionType<TD,TS> { typedef TD type ; }; \
template<> struct BinaryExpressionType<TS,TD> { typedef TD type ; }; \

KOKKOS_BINARYEXPRESSIONTYPE(double,float)
KOKKOS_BINARYEXPRESSIONTYPE(double,long)
KOKKOS_BINARYEXPRESSIONTYPE(double,int)
KOKKOS_BINARYEXPRESSIONTYPE(double,unsigned long)
KOKKOS_BINARYEXPRESSIONTYPE(double,unsigned int)

KOKKOS_BINARYEXPRESSIONTYPE(float,long)
KOKKOS_BINARYEXPRESSIONTYPE(float,int)
KOKKOS_BINARYEXPRESSIONTYPE(float,unsigned long)
KOKKOS_BINARYEXPRESSIONTYPE(float,unsigned int)

KOKKOS_BINARYEXPRESSIONTYPE(long,int)
KOKKOS_BINARYEXPRESSIONTYPE(long,unsigned long)
KOKKOS_BINARYEXPRESSIONTYPE(long,unsigned int)

KOKKOS_BINARYEXPRESSIONTYPE(int,unsigned int)

#undef KOKKOS_BINARYEXPRESSIONTYPE

//----------------------------------------------------------------------------

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Assigment operators:

#define KOKKOS_CORE_ASSIGN_OPERATOR_CLASS ArrayAssignAddition
#define KOKKOS_CORE_ASSIGN_OPERATOR  +=
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef KOKKOS_CORE_ASSIGN_OPERATOR
#undef KOKKOS_CORE_ASSIGN_OPERATOR_CLASS

#define KOKKOS_CORE_ASSIGN_OPERATOR_CLASS ArrayAssignSubtraction
#define KOKKOS_CORE_ASSIGN_OPERATOR  -=
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef KOKKOS_CORE_ASSIGN_OPERATOR
#undef KOKKOS_CORE_ASSIGN_OPERATOR_CLASS

#define KOKKOS_CORE_ASSIGN_OPERATOR_CLASS ArrayAssignMultiplication
#define KOKKOS_CORE_ASSIGN_OPERATOR  *=
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef KOKKOS_CORE_ASSIGN_OPERATOR
#undef KOKKOS_CORE_ASSIGN_OPERATOR_CLASS

#define KOKKOS_CORE_ASSIGN_OPERATOR_CLASS ArrayAssignDivision
#define KOKKOS_CORE_ASSIGN_OPERATOR  /=
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef KOKKOS_CORE_ASSIGN_OPERATOR
#undef KOKKOS_CORE_ASSIGN_OPERATOR_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOS_CORE_UNARY_OPERATOR_CLASS  ArrayProxyMinus
#define KOKKOS_CORE_UNARY_OPERATOR        -
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_OPERATOR
#undef  KOKKOS_CORE_UNARY_OPERATOR_CLASS

#define KOKKOS_CORE_UNARY_OPERATOR_CLASS  ArrayProxyPlus
#define KOKKOS_CORE_UNARY_OPERATOR        +
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_OPERATOR
#undef  KOKKOS_CORE_UNARY_OPERATOR_CLASS

#define KOKKOS_CORE_UNARY_OPERATOR_CLASS  ArrayProxyNot
#define KOKKOS_CORE_UNARY_OPERATOR        !
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_OPERATOR
#undef  KOKKOS_CORE_UNARY_OPERATOR_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOS_CORE_BINARY_OPERATOR_CLASS  ArrayProxyAddition
#define KOKKOS_CORE_BINARY_OPERATOR        +
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_OPERATOR
#undef  KOKKOS_CORE_BINARY_OPERATOR_CLASS

#define KOKKOS_CORE_BINARY_OPERATOR_CLASS  ArrayProxySubtraction
#define KOKKOS_CORE_BINARY_OPERATOR        -
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_OPERATOR
#undef  KOKKOS_CORE_BINARY_OPERATOR_CLASS

#define KOKKOS_CORE_BINARY_OPERATOR_CLASS  ArrayProxyMultiplication
#define KOKKOS_CORE_BINARY_OPERATOR        *
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_OPERATOR
#undef  KOKKOS_CORE_BINARY_OPERATOR_CLASS

#define KOKKOS_CORE_BINARY_OPERATOR_CLASS  ArrayProxyDivision
#define KOKKOS_CORE_BINARY_OPERATOR        /
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_OPERATOR
#undef  KOKKOS_CORE_BINARY_OPERATOR_CLASS

#define KOKKOS_CORE_BINARY_OPERATOR_CLASS  ArrayProxyModulo
#define KOKKOS_CORE_BINARY_OPERATOR        %
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_OPERATOR
#undef  KOKKOS_CORE_BINARY_OPERATOR_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename Type >
struct ArrayWeakOrdering
{
  enum Result { EQUAL , LESS , GREATER , NOT_ORDERED };

  template< class ArrayLHS , class ArrayRHS >
  KOKKOS_INLINE_FUNCTION static
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
  KOKKOS_INLINE_FUNCTION static
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
  KOKKOS_INLINE_FUNCTION static
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

  KOKKOS_INLINE_FUNCTION static
  bool equal( const Result result ) { return EQUAL == result ; }

  KOKKOS_INLINE_FUNCTION static
  bool not_equal( const Result result ) { return EQUAL != result ; }

  KOKKOS_INLINE_FUNCTION static
  bool less( const Result result ) { return LESS == result ; }

  KOKKOS_INLINE_FUNCTION static
  bool less_or_equal( const Result result )
    { return LESS == result || EQUAL == result ; }

  KOKKOS_INLINE_FUNCTION static
  bool greater( const Result result )
    { return GREATER == result ; }

  KOKKOS_INLINE_FUNCTION static
  bool greater_or_equal( const Result result )
    { return GREATER == result || EQUAL == result ; }
};

} // namespace Impl
} // namespace Kokkos

#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR       ==
#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL  equal
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR

#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR       !=
#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL  not_equal
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR

#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR       <
#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL  less
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR

#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR       >
#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL  greater
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR

#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR       <=
#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL  less_or_equal
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR

#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR       >=
#define KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL  greater_or_equal
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL
#undef  KOKKOS_CORE_WEAK_ORDERING_OPERATOR

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyAbs
#define KOKKOS_CORE_UNARY_FUNCTION         abs
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::abs
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyFabs
#define KOKKOS_CORE_UNARY_FUNCTION         fabs
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::fabs
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxySqrt
#define KOKKOS_CORE_UNARY_FUNCTION         sqrt
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::sqrt
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyExp
#define KOKKOS_CORE_UNARY_FUNCTION         exp
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::exp
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyLog
#define KOKKOS_CORE_UNARY_FUNCTION         log
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::log
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyLog10
#define KOKKOS_CORE_UNARY_FUNCTION         log10
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::log10
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxySin
#define KOKKOS_CORE_UNARY_FUNCTION         sin
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::sin
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyCos
#define KOKKOS_CORE_UNARY_FUNCTION         cos
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::cos
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyTan
#define KOKKOS_CORE_UNARY_FUNCTION         tan
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::tan
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyAcos
#define KOKKOS_CORE_UNARY_FUNCTION         acos
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::acos
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyAsin
#define KOKKOS_CORE_UNARY_FUNCTION         asin
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::asin
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

#define KOKKOS_CORE_UNARY_FUNCTION_CLASS   ArrayProxyAtan
#define KOKKOS_CORE_UNARY_FUNCTION         atan
#define KOKKOS_CORE_UNARY_FUNCTION_MEMBER  ::atan
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_UNARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_UNARY_FUNCTION
#undef  KOKKOS_CORE_UNARY_FUNCTION_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< typename T >
KOKKOS_INLINE_FUNCTION
const T & min( const T & lhs , const T & rhs )
{ return lhs < rhs ? lhs : rhs ; }

template< typename T >
KOKKOS_INLINE_FUNCTION
const T & max( const T & lhs , const T & rhs )
{ return lhs < rhs ? rhs : lhs ; }

}

#define KOKKOS_CORE_BINARY_FUNCTION_CLASS   ArrayProxyMin
#define KOKKOS_CORE_BINARY_FUNCTION         min
#define KOKKOS_CORE_BINARY_FUNCTION_MEMBER  Kokkos::min
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_BINARY_FUNCTION
#undef  KOKKOS_CORE_BINARY_FUNCTION_CLASS

#define KOKKOS_CORE_BINARY_FUNCTION_CLASS   ArrayProxyMax
#define KOKKOS_CORE_BINARY_FUNCTION         max
#define KOKKOS_CORE_BINARY_FUNCTION_MEMBER  Kokkos::max
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_BINARY_FUNCTION
#undef  KOKKOS_CORE_BINARY_FUNCTION_CLASS

#define KOKKOS_CORE_BINARY_FUNCTION_CLASS   ArrayProxyAtan2
#define KOKKOS_CORE_BINARY_FUNCTION         atan2
#define KOKKOS_CORE_BINARY_FUNCTION_MEMBER  ::atan2
#include <impl/Kokkos_ArrayExp_macros.hpp>
#undef  KOKKOS_CORE_BINARY_FUNCTION_MEMBER
#undef  KOKKOS_CORE_BINARY_FUNCTION
#undef  KOKKOS_CORE_BINARY_FUNCTION_CLASS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CORE_HPP */

