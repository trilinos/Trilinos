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


#define PRINT_ARRAY_PROXY_NAME_BASE( S ) # S
#define PRINT_ARRAY_PROXY_NAME( S ) PRINT_ARRAY_PROXY_NAME_BASE( S )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_CORE_ASSIGN_OPERATOR_CLASS ) && \
    defined( KOKKOS_CORE_ASSIGN_OPERATOR )

namespace Kokkos {

template< class LHS , class RHS >
struct KOKKOS_CORE_ASSIGN_OPERATOR_CLASS {};

//----------------------------------------------------------------------------
// ArrayExp OP= cv-ArrayExp
// ArrayExp OP= scalar

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
class KOKKOS_CORE_ASSIGN_OPERATOR_CLASS
  < Array< TypeLHS,CountLHS,ProxyLHS > ,
    Array< TypeRHS,CountRHS,ProxyRHS > >
{
public:
  KOKKOS_INLINE_FUNCTION
  KOKKOS_CORE_ASSIGN_OPERATOR_CLASS
    (       Array<TypeLHS,CountLHS,ProxyLHS > & lhs ,
      const Array<TypeRHS,CountRHS,ProxyRHS > & rhs )
  {
    for ( unsigned i = 0 ; i < lhs.size() ; ++i ) {
      lhs[i] KOKKOS_CORE_ASSIGN_OPERATOR rhs[i] ;
    }
  }
};

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS , class RHS >
class KOKKOS_CORE_ASSIGN_OPERATOR_CLASS
  < Array< TypeLHS,CountLHS,ProxyLHS > , RHS >
{
public:
  KOKKOS_INLINE_FUNCTION
  KOKKOS_CORE_ASSIGN_OPERATOR_CLASS
    ( Array<TypeLHS,CountLHS,ProxyLHS > & lhs , const RHS rhs )
  {
    for ( unsigned i = 0 ; i < lhs.size() ; ++i ) {
      lhs[i] KOKKOS_CORE_ASSIGN_OPERATOR rhs ;
    }
  }
};

template< typename T , unsigned N , class Proxy , class RHS >
KOKKOS_INLINE_FUNCTION
Array< T , N , Proxy > &
operator KOKKOS_CORE_ASSIGN_OPERATOR
  ( Array< T , N , Proxy > & lhs , const RHS & rhs )
{
  KOKKOS_CORE_ASSIGN_OPERATOR_CLASS < Array< T , N , Proxy > , RHS >( lhs , rhs );
  return lhs ;
}

} /* namespace Kokkos */

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_CORE_UNARY_OPERATOR_CLASS ) && \
    defined( KOKKOS_CORE_UNARY_OPERATOR )

namespace Kokkos {

template< class > struct KOKKOS_CORE_UNARY_OPERATOR_CLASS {};

template< typename TypeRHS , unsigned N , class ProxyRHS >
class Array< TypeRHS , N ,
             KOKKOS_CORE_UNARY_OPERATOR_CLASS< ProxyRHS > >
{
private:

  const Array< TypeRHS , N , ProxyRHS > rhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOS_INLINE_FUNCTION
  unsigned size() const { return rhs.size(); }

  template< class ArgType >
  KOKKOS_INLINE_FUNCTION explicit
  Array( const ArgType & arg ) : rhs( arg ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  TypeRHS operator[]( const iType & i ) const
    {
      array_check_bounds(i,rhs.size());
      return KOKKOS_CORE_UNARY_OPERATOR rhs[i] ;
    }

  void print_expression( std::ostream & s ) const
  {
    s << "( " << PRINT_ARRAY_PROXY_NAME( KOKKOS_CORE_UNARY_OPERATOR )
      << " " ;
    rhs.print_expression(s);
    s << " )" ;
  }
};

template< typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< TypeRHS , N ,
       KOKKOS_CORE_UNARY_OPERATOR_CLASS
       < typename ArrayProxy< ProxyRHS >::type > >
operator KOKKOS_CORE_UNARY_OPERATOR
  ( const Array< TypeRHS , N , ProxyRHS > & a )
{
  return 
    Array< TypeRHS , N ,
           KOKKOS_CORE_UNARY_OPERATOR_CLASS
           < typename ArrayProxy< ProxyRHS >::type > >( a );
}

template< typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< TypeRHS , N ,
       KOKKOS_CORE_UNARY_OPERATOR_CLASS
       < typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOS_CORE_UNARY_OPERATOR
  ( const volatile Array< TypeRHS , N , ProxyRHS > & a )
{
  return 
    Array< TypeRHS , N ,
           KOKKOS_CORE_UNARY_OPERATOR_CLASS
           < typename ArrayProxy< volatile ProxyRHS >::type > >( a );
}

} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_CORE_UNARY_OPERATOR_CLASS ) && \
              defined( KOKKOS_CORE_UNARY_OPERATOR ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_CORE_BINARY_OPERATOR_CLASS ) && \
    defined( KOKKOS_CORE_BINARY_OPERATOR )

//
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]
// ArrayExp <= scalar      OP cv-ArrayExp   [2]
// ArrayExp <= cv-ArrayExp OP scalar        [2]
//

namespace Kokkos {

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
struct KOKKOS_CORE_BINARY_OPERATOR_CLASS {};

//----------------------------------------------------------------------------
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]

template< typename T , unsigned N ,
          typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
class Array< T , N ,
             KOKKOS_CORE_BINARY_OPERATOR_CLASS
             <TypeLHS,CountLHS,ProxyLHS, TypeRHS,CountRHS,ProxyRHS> >
{
private:

  const Array< TypeLHS , CountLHS , ProxyLHS > lhs ;
  const Array< TypeRHS , CountRHS , ProxyRHS > rhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOS_INLINE_FUNCTION
  unsigned size() const
    { return lhs.size() < rhs.size() ? lhs.size() : rhs.size(); }

  template< class ArgLHS , class ArgRHS >
  KOKKOS_INLINE_FUNCTION
  Array( const ArgLHS & arg_lhs , const ArgRHS & arg_rhs )
    : lhs( arg_lhs ), rhs( arg_rhs ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { return lhs[i] KOKKOS_CORE_BINARY_OPERATOR rhs[i] ; }

  void print_expression( std::ostream & s ) const
  {
    s << "(" ;
    lhs.print_expression(s);
    s << " " << PRINT_ARRAY_PROXY_NAME( KOKKOS_CORE_BINARY_OPERATOR )
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
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy<          ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const          Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy<          ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy<          ProxyRHS >::type > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const          Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy<          ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
         TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< volatile ProxyLHS >::type ,
             TypeRHS , CountRHS , typename ArrayProxy< volatile ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// scalar OP cv-array

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , N , ArrayProxyValue ,
         TypeRHS , N , typename ArrayProxy< ProxyRHS >::type > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const TypeLHS & lhs , const Array< TypeRHS , N , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , N , ArrayProxyValue ,
             TypeRHS , N , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , N , ArrayProxyValue ,
         TypeRHS , N , typename ArrayProxy< volatile ProxyRHS >::type > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const TypeLHS & lhs ,
    const volatile Array< TypeRHS , N , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , N , ArrayProxyValue ,
             TypeRHS , N , typename ArrayProxy< volatile ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// cv-array OP scalar

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
       KOKKOS_CORE_BINARY_OPERATOR_CLASS
       < TypeLHS , N , typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , N , ArrayProxyValue > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const Array< TypeLHS , N , ProxyLHS > & lhs , const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , N , typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , N , ArrayProxyValue > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOS_INLINE_FUNCTION
class Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
             KOKKOS_CORE_BINARY_OPERATOR_CLASS
             < TypeLHS , N , typename ArrayProxy< volatile ProxyLHS >::type ,
               TypeRHS , N , ArrayProxyValue > >
operator KOKKOS_CORE_BINARY_OPERATOR
  ( const volatile Array< TypeLHS , N , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , N ,
           KOKKOS_CORE_BINARY_OPERATOR_CLASS
           < TypeLHS , N , typename ArrayProxy< volatile ProxyLHS >::type ,
             TypeRHS , N , ArrayProxyValue > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------

} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_CORE_BINARY_OPERATOR_CLASS ) && \
              defined( KOKKOS_CORE_BINARY_OPERATOR ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_CORE_BINARY_FUNCTION_CLASS ) && \
    defined( KOKKOS_CORE_BINARY_FUNCTION ) && \
    defined( KOKKOS_CORE_BINARY_FUNCTION_MEMBER )

//
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]
// ArrayExp <= scalar      OP cv-ArrayExp   [2]
// ArrayExp <= cv-ArrayExp OP scalar        [2]
//

namespace Kokkos {

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
struct KOKKOS_CORE_BINARY_FUNCTION_CLASS {};

//----------------------------------------------------------------------------
// ArrayExp <= cv-ArrayExp OP cv-ArrayExp   [4]

template< typename T , unsigned N ,
          typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
class Array< T , N ,
             KOKKOS_CORE_BINARY_FUNCTION_CLASS
               < TypeLHS , CountLHS , ProxyLHS ,
                 TypeRHS , CountRHS , ProxyRHS > >
{
private:

  const Array< TypeLHS , CountLHS , ProxyLHS > lhs ;
  const Array< TypeRHS , CountRHS , ProxyRHS > rhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOS_INLINE_FUNCTION
  unsigned size() const
    { return lhs.size() < rhs.size() ? lhs.size() : rhs.size(); }

  template< class ArgLHS , class ArgRHS >
  KOKKOS_INLINE_FUNCTION
  Array( const ArgLHS & arg_lhs , const ArgRHS & arg_rhs )
    : lhs( arg_lhs ), rhs( arg_rhs ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { return KOKKOS_CORE_BINARY_FUNCTION_MEMBER( lhs[i] , rhs[i] ); }

  void print_expression( std::ostream & s ) const
  {
    s << PRINT_ARRAY_PROXY_NAME( KOKKOS_CORE_BINARY_FUNCTION ) << "( " ; 
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
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS ,          typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const          Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS ,          typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
       ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
         < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
           TypeRHS , CountRHS ,          typename ArrayProxy< ProxyRHS >::type > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const          Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type ,
           ( CountLHS < CountRHS ? CountLHS : CountRHS ) ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
             < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
               TypeRHS , CountRHS ,          typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// scalar OP cv-array

template< typename TypeLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
       < TypeLHS , CountRHS , ArrayProxyValue ,
         TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const TypeLHS & lhs ,
    const Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
           < TypeLHS , CountRHS , ArrayProxyValue ,
             TypeRHS , CountRHS , typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

template< typename TypeLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
       < TypeLHS , CountRHS , ArrayProxyValue ,
         TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const TypeLHS & lhs ,
    const volatile Array< TypeRHS , CountRHS , ProxyRHS > & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountRHS ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
           < TypeLHS , CountRHS , ArrayProxyValue ,
             TypeRHS , CountRHS , volatile typename ArrayProxy< ProxyRHS >::type > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------
// cv-array OP scalar

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
       < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , CountLHS , ArrayProxyValue > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
           < TypeLHS , CountLHS , typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , CountLHS , ArrayProxyValue > >
      ( lhs , rhs );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS >
KOKKOS_INLINE_FUNCTION
Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
       KOKKOS_CORE_BINARY_FUNCTION_CLASS
       < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
         TypeRHS , CountLHS , ArrayProxyValue > >
KOKKOS_CORE_BINARY_FUNCTION
  ( const volatile Array< TypeLHS , CountLHS , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  return 
    Array< typename BinaryExpressionType<TypeLHS,TypeRHS>::type , CountLHS ,
           KOKKOS_CORE_BINARY_FUNCTION_CLASS
           < TypeLHS , CountLHS , volatile typename ArrayProxy< ProxyLHS >::type ,
             TypeRHS , CountLHS , ArrayProxyValue > >
      ( lhs , rhs );
}

//----------------------------------------------------------------------------

} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_CORE_BINARY_FUNCTION_CLASS ) && \
              defined( KOKKOS_CORE_BINARY_FUNCTION ) && \
              defined( KOKKOS_CORE_BINARY_FUNCTION_MEMBER ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_CORE_WEAK_ORDERING_OPERATOR ) && \
    defined( KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL )

//
// bool <= cv-ArrayExp OP cv-ArrayExp   [4]
// bool <= scalar      OP cv-ArrayExp   [2]
// bool <= cv-ArrayExp OP scalar        [2]
//

namespace Kokkos {

//----------------------------------------------------------------------------
// bool <= cv-ArrayExp OP cv-ArrayExp   [4]

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const volatile Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const          Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const volatile Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned CountLHS , class ProxyLHS ,
          typename TypeRHS , unsigned CountRHS , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const volatile Array< TypeRHS , CountLHS , ProxyLHS > & lhs ,
    const          Array< TypeLHS , CountRHS , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_array( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

//----------------------------------------------------------------------------
// scalar OP cv-array

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const TypeLHS & lhs , const Array< TypeRHS , N , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_value_array( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS ,
          typename TypeRHS , unsigned N , class ProxyRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const TypeLHS & lhs ,
    const volatile Array< TypeRHS , N , ProxyRHS > & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_value_array( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

//----------------------------------------------------------------------------
// cv-array OP scalar

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const Array< TypeLHS , N , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_value( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

template< typename TypeLHS , unsigned N , class ProxyLHS ,
          typename TypeRHS >
KOKKOS_INLINE_FUNCTION
bool operator KOKKOS_CORE_WEAK_ORDERING_OPERATOR
  ( const volatile Array< TypeLHS , N , ProxyLHS > & lhs ,
    const TypeRHS & rhs )
{
  typedef typename BinaryExpressionType<TypeLHS,TypeRHS>::type T ;

  typedef Impl::ArrayWeakOrdering<T> WeakOrdering ;

  const typename WeakOrdering::Result result =
    WeakOrdering::compare_array_value( lhs , rhs );

  return WeakOrdering::KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL( result );
}

//----------------------------------------------------------------------------

} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_CORE_WEAK_ORDERING_OPERATOR ) && \
              defined( KOKKOS_CORE_WEAK_ORDERING_OPERATOR_EVAL ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_CORE_UNARY_FUNCTION_CLASS ) && \
    defined( KOKKOS_CORE_UNARY_FUNCTION ) && \
    defined( KOKKOS_CORE_UNARY_FUNCTION_MEMBER )

namespace Kokkos {

template< class > struct KOKKOS_CORE_UNARY_FUNCTION_CLASS {};

template< typename T , unsigned N , class LHS >
class Array< T , N , KOKKOS_CORE_UNARY_FUNCTION_CLASS< LHS > >
{
private:

  typedef typename ArrayProxy< LHS >::type ProxyLHS ;

  const Array<T,N,ProxyLHS > lhs ;

  Array();
  Array & operator = ( const Array & );

public:

  KOKKOS_INLINE_FUNCTION
  unsigned size() const { return lhs.size(); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { 
      array_check_bounds(i,lhs.size());
      return KOKKOS_CORE_UNARY_FUNCTION_MEMBER( lhs[i] );
    }

  KOKKOS_INLINE_FUNCTION
  Array( const Array<T,N,LHS> & arg ) : lhs( arg ) {}

  void print_expression( std::ostream & s ) const
  {
    s << PRINT_ARRAY_PROXY_NAME( KOKKOS_CORE_UNARY_FUNCTION );
    s << "(" ;
    lhs.print_expression(s);
    s << ")" ;
  }
};

template< typename T , unsigned N , class LHS >
class Array< T , N ,
      KOKKOS_CORE_UNARY_FUNCTION_CLASS< volatile LHS > >
{
private:

  typedef typename ArrayProxy< volatile LHS >::type ProxyLHS ;

  const Array<T,N,ProxyLHS > lhs ;

  Array();
  Array & operator = ( const Array & );

public:

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  T operator[]( const iType & i ) const
    { 
      array_check_bounds(i,lhs.size());
      return KOKKOS_CORE_UNARY_FUNCTION_MEMBER( lhs[i] );
    }

  KOKKOS_INLINE_FUNCTION
  Array( const volatile Array<T,N,LHS> & arg ) : lhs( arg ) {}

  void print_expression( std::ostream & s ) const
  {
    s << PRINT_ARRAY_PROXY_NAME( KOKKOS_CORE_UNARY_FUNCTION );
    s << "(" ;
    lhs.print_expression(s);
    s << ")" ;
  }
};

template< typename T , unsigned N , class Proxy >
KOKKOS_INLINE_FUNCTION
Array< T , N ,
       KOKKOS_CORE_UNARY_FUNCTION_CLASS < Proxy > >
KOKKOS_CORE_UNARY_FUNCTION
  ( const Array< T , N , Proxy > & a )
{
  return 
    Array< T , N ,
           KOKKOS_CORE_UNARY_FUNCTION_CLASS < Proxy > >( a );
}

template< typename T , unsigned N , class Proxy >
KOKKOS_INLINE_FUNCTION
Array< T , N ,
       KOKKOS_CORE_UNARY_FUNCTION_CLASS < volatile Proxy > >
KOKKOS_CORE_UNARY_FUNCTION
  ( const volatile Array< T , N , Proxy > & a )
{
  return 
    Array< T , N ,
           KOKKOS_CORE_UNARY_FUNCTION_CLASS < volatile Proxy > >( a );
}

} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_CORE_UNARY_FUNCTION_CLASS ) && \
              defined( KOKKOS_CORE_UNARY_FUNCTION ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#undef PRINT_ARRAY_PROXY_NAME
#undef PRINT_ARRAY_PROXY_NAME_BASE



