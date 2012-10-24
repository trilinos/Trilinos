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


#define PRINT_ARRAY_PROXY_NAME_BASE( S ) # S
#define PRINT_ARRAY_PROXY_NAME( S ) PRINT_ARRAY_PROXY_NAME_BASE( S )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR )

// ArrayExp OP= cv-ArrayExp
// ArrayExp OP= scalar

namespace KokkosArray {

template< typename T, unsigned N, class ProxyLHS, class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< T , N , ProxyLHS > &
operator KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR
  ( Array< T , N , ProxyLHS > & lhs ,
    Array< T , N , ProxyRHS > & rhs )
{
  for ( unsigned i = 0 ; i < lhs.size() ; ++i ) {
    lhs[i] KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR rhs[i] ;
  }
  return lhs ;
}

template< typename T , unsigned N , class Proxy , class RHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< T , N , Proxy > &
operator KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR
  ( Array< T , N , Proxy > & lhs , const RHS rhs )
{
  for ( unsigned i = 0 ; i < lhs.size() ; ++i ) {
    lhs[i] KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR rhs ;
  }
  return lhs ;
}

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_ARRAY_ASSIGN_OPERATOR ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS ) && \
    defined( KOKKOSARRAY_ARRAY_UNARY_OPERATOR )

namespace KokkosArray {

template< class > struct KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS {};

template< typename TypeRHS , unsigned N , class ProxyRHS >
class Array< TypeRHS , N ,
             KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS< ProxyRHS > >
{
private:

  const Array< TypeRHS , N , ProxyRHS > rhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return rhs.size(); }

  template< class ArgType >
  KOKKOSARRAY_INLINE_FUNCTION explicit
  Array( const ArgType & arg ) : rhs( arg ) {}

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  TypeRHS operator[]( const iType & i ) const
    {
      array_check_bounds(i,N);
      return KOKKOSARRAY_ARRAY_UNARY_OPERATOR rhs[i] ;
    }

  void print_expression( std::ostream & s ) const
  {
    s << "( " << PRINT_ARRAY_PROXY_NAME( KOKKOSARRAY_ARRAY_UNARY_OPERATOR )
      << " " ;
    rhs.print_expression(s);
    s << " )" ;
  }
};

template< typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< TypeRHS , N ,
       KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS
       < typename ArrayProxy< ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_UNARY_OPERATOR
  ( const Array< TypeRHS , N , ProxyRHS > & a )
{
  return 
    Array< TypeRHS , N ,
           KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS
           < typename ArrayProxy< ProxyRHS >::type > >( a );
}

template< typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< TypeRHS , N ,
       KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS
       < typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_UNARY_OPERATOR
  ( const volatile Array< TypeRHS , N , ProxyRHS > & a )
{
  return 
    Array< TypeRHS , N ,
           KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS
           < typename ArrayProxy< volatile ProxyRHS >::type > >( a );
}

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_ARRAY_UNARY_OPERATOR_CLASS ) && \
              defined( KOKKOSARRAY_ARRAY_UNARY_OPERATOR ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS ) && \
    defined( KOKKOSARRAY_ARRAY_BINARY_OPERATOR )

//
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]
// ArrayExp <= scalar      OP cv-ArrayExp   [2]
// ArrayExp <= cv-ArrayExp OP scalar        [2]
//

namespace KokkosArray {

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
struct KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS {};

//----------------------------------------------------------------------------
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]

template< typename T , unsigned N ,
          typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
class Array< T , N ,
             KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
             <TypeLHS,CountLHS,ProxyLHS, TypeRHS,CountRHS,ProxyRHS> >
{
private:

  const Array< TypeLHS , CountLHS , ProxyLHS > lhs ;
  const Array< TypeRHS , CountRHS , ProxyRHS > rhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const
    { return lhs.size() < rhs.size() ? lhs.size() : rhs.size(); }

  template< class ArgLHS , class ArgRHS >
  KOKKOSARRAY_INLINE_FUNCTION
  Array( const ArgLHS & arg_lhs , const ArgRHS & arg_rhs )
    : lhs( arg_lhs ), rhs( arg_rhs ) {}

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { return lhs[i] KOKKOSARRAY_ARRAY_BINARY_OPERATOR rhs[i] ; }

  void print_expression( std::ostream & s ) const
  {
    s << "(" ;
    lhs.print_expression(s);
    s << " " << PRINT_ARRAY_PROXY_NAME( KOKKOSARRAY_ARRAY_BINARY_OPERATOR )
      << " " ;
    rhs.print_expression(s);
    s << ")" ;
  }

  template< class X >
  void operator = ( const X & )
  { ERROR__cannot_assign_value_to_constant_expression< Array >(); }
};

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy<          ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const          Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy<          ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy<          ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const          Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy<          ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// scalar OP cv-array

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , N , ArrayProxyValue ,
         TypeRHS , N , typename ArrayProxy< ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const TypeLHS & lhs , const Array< TypeRHS , N , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , N , ArrayProxyValue ,
             TypeRHS , N , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , N , ArrayProxyValue ,
         TypeRHS , N , typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const TypeLHS & lhs ,
    const volatile Array< TypeRHS , N , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , N , ArrayProxyValue ,
             TypeRHS , N , typename ArrayProxy< volatile ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// cv-array OP scalar

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
       KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
       < TypeLHS , N , typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , N , ArrayProxyValue > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const Array< TypeLHS , N , ProxyLHS > & lhs , const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , N , typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , N , ArrayProxyValue > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOSARRAY_INLINE_FUNCTION
class Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
             KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
             < TypeLHS , N , typename ArrayProxy< volatile ProxyLHS >::type ,
               TypeRHS , N , ArrayProxyValue > >
operator KOKKOSARRAY_ARRAY_BINARY_OPERATOR
  ( const volatile Array< TypeLHS , N , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS
           < TypeLHS , N , typename ArrayProxy< volatile ProxyLHS >::type ,
             TypeRHS , N , ArrayProxyValue > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_ARRAY_BINARY_OPERATOR_CLASS ) && \
              defined( KOKKOSARRAY_ARRAY_BINARY_OPERATOR ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS ) && \
    defined( KOKKOSARRAY_ARRAY_BINARY_FUNCTION ) && \
    defined( KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER )

//
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]
// ArrayExp <= scalar      OP cv-ArrayExp   [2]
// ArrayExp <= cv-ArrayExp OP scalar        [2]
//

namespace KokkosArray {

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
struct KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS {};

//----------------------------------------------------------------------------
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]

template< typename T , unsigned N ,
          typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
class Array< T , N ,
             KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
               < TypeLHS , CountLHS , ProxyLHS ,
                 TypeRHS , CountRHS , ProxyRHS > >
{
private:

  const Array< TypeLHS , CountLHS , ProxyLHS > lhs ;
  const Array< TypeRHS , CountRHS , ProxyRHS > rhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const
    { return lhs.size() < rhs.size() ? lhs.size() : rhs.size(); }

  template< class ArgLHS , class ArgRHS >
  KOKKOSARRAY_INLINE_FUNCTION
  Array( const ArgLHS & arg_lhs , const ArgRHS & arg_rhs )
    : lhs( arg_lhs ), rhs( arg_rhs ) {}

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { return KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER( lhs[i] , rhs[i] ); }

  void print_expression( std::ostream & s ) const
  {
    s << PRINT_ARRAY_PROXY_NAME( KOKKOSARRAY_ARRAY_BINARY_FUNCTION ) << "( " ; 
    lhs.print_expression(s);
    s  << " , " ;
    rhs.print_expression(s);
    s << " )" ;
  }

  template< class X >
  void operator = ( const X & )
  { ERROR__cannot_assign_value_to_constant_expression< Array >(); }
};


template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS ,          typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const          Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS ,          typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS ,          typename ArrayProxy< ProxyRHS >::type > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const          Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS ,          typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// scalar OP cv-array

template< typename TypeLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
       < TypeLHS , CountRHS , ArrayProxyValue ,
         TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const TypeLHS & lhs ,
    const Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
           < TypeLHS , CountRHS , ArrayProxyValue ,
             TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
       < TypeLHS , CountRHS , ArrayProxyValue ,
         TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const TypeLHS & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
           < TypeLHS , CountRHS , ArrayProxyValue ,
             TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// cv-array OP scalar

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , CountLHS , ArrayProxyValue > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , CountLHS , ArrayProxyValue > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS >
KOKKOSARRAY_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
       KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
       < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , CountLHS , ArrayProxyValue > >
KOKKOSARRAY_ARRAY_BINARY_FUNCTION
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
           KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS
           < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , CountLHS , ArrayProxyValue > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_ARRAY_BINARY_FUNCTION_CLASS ) && \
              defined( KOKKOSARRAY_ARRAY_BINARY_FUNCTION ) && \
              defined( KOKKOSARRAY_ARRAY_BINARY_FUNCTION_MEMBER ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR ) && \
    defined( KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL )

//
// bool <= cv-ArrayExp OP cv-ArrayExp   [4]
// bool <= scalar      OP cv-ArrayExp   [2]
// bool <= cv-ArrayExp OP scalar        [2]
//

namespace KokkosArray {

//----------------------------------------------------------------------------
// bool <= cv-ArrayExp OP cv-ArrayExp   [4]

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const volatile Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const          Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const volatile Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const          Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

//----------------------------------------------------------------------------
// scalar OP cv-array

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const TypeLHS & lhs , const Array< TypeRHS , N , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_value_array( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const TypeLHS & lhs ,
    const volatile Array< TypeRHS , N , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_value_array( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

//----------------------------------------------------------------------------
// cv-array OP scalar

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const Array< TypeLHS , N , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_value( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR
  ( const volatile Array< TypeLHS , N , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_value( lhs , rhs );

  return WeakOrdering::KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL( result );
}

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR ) && \
              defined( KOKKOSARRAY_ARRAY_WEAK_ORDERING_OPERATOR_EVAL ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS ) && \
    defined( KOKKOSARRAY_ARRAY_UNARY_FUNCTION ) && \
    defined( KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER )

namespace KokkosArray {

template< class > struct KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS {};

template< typename T , unsigned N , class LHS >
class Array< T , N , KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS< LHS > >
{
private:

  typedef typename ArrayProxy< LHS >::type ProxyLHS ;

  const Array<T,N,ProxyLHS > lhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned size() const { return lhs.size(); }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { 
      array_check_bounds(i,N);
      return KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER( lhs[i] );
    }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const Array<T,N,LHS> & arg ) : lhs( arg ) {}

  void print_expression( std::ostream & s ) const
  {
    s << PRINT_ARRAY_PROXY_NAME( KOKKOSARRAY_ARRAY_UNARY_FUNCTION );
    s << "(" ;
    lhs.print_expression(s);
    s << ")" ;
  }
};

template< typename T , unsigned N , class LHS >
class Array< T , N ,
      KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS< volatile LHS > >
{
private:

  typedef typename ArrayProxy< volatile LHS >::type ProxyLHS ;

  const Array<T,N,ProxyLHS > lhs ;

  Array();
  Array & operator = ( const Array & );

public:

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { 
      array_check_bounds(i,N);
      return KOKKOSARRAY_ARRAY_UNARY_FUNCTION_MEMBER( lhs[i] );
    }

  KOKKOSARRAY_INLINE_FUNCTION
  Array( const volatile Array<T,N,LHS> & arg ) : lhs( arg ) {}

  void print_expression( std::ostream & s ) const
  {
    s << PRINT_ARRAY_PROXY_NAME( KOKKOSARRAY_ARRAY_UNARY_FUNCTION );
    s << "(" ;
    lhs.print_expression(s);
    s << ")" ;
  }
};

template< typename T , unsigned N , class Proxy >
KOKKOSARRAY_INLINE_FUNCTION
Array< T , N ,
       KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS < Proxy > >
KOKKOSARRAY_ARRAY_UNARY_FUNCTION
  ( const Array< T , N , Proxy > & a )
{
  return 
    Array< T , N ,
           KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS < Proxy > >( a );
}

template< typename T , unsigned N , class Proxy >
KOKKOSARRAY_INLINE_FUNCTION
Array< T , N ,
       KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS < volatile Proxy > >
KOKKOSARRAY_ARRAY_UNARY_FUNCTION
  ( const volatile Array< T , N , Proxy > & a )
{
  return 
    Array< T , N ,
           KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS < volatile Proxy > >( a );
}

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_ARRAY_UNARY_FUNCTION_CLASS ) && \
              defined( KOKKOSARRAY_ARRAY_UNARY_FUNCTION ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#undef PRINT_ARRAY_PROXY_NAME
#undef PRINT_ARRAY_PROXY_NAME_BASE



