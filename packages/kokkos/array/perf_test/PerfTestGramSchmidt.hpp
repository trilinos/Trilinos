/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef UNIT_TEST_GRAM_SCHMIDT_OPERATORS
#define UNIT_TEST_GRAM_SCHMIDT_OPERATORS

#include <cmath>

namespace KokkosArray {

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct MultiVectorScale ;

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct MultiVectorYPAX ;

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct DotSingle ;

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct Dot ;

}

#endif /* #ifndef UNIT_TEST_GRAM_SCHMIDT_OPERATORS */

//----------------------------------------------------------------------------
// Types and static methods for summation

namespace KokkosArray {

template< typename Scalar >
struct DotSingle<Scalar, KOKKOS_MACRO_DEVICE , 1 > {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;
  typedef KokkosArray::MultiVector<value_type,device_type> vector_type ;

  vector_type X ;

  DotSingle( const vector_type & argX ) : X( argX ) {}

  DotSingle( const vector_type & argX ,
             const size_type     argI ) : X( argX , argI ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & update ) const
  {
    const Scalar value = X(iwork);
    update += value * value ;
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update , const volatile value_type & source )
  { update += source ; }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
};

template< typename Scalar >
struct Dot< Scalar, KOKKOS_MACRO_DEVICE , 1 > {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;
  typedef KokkosArray::MultiVector<value_type,device_type> vector_type ;

  vector_type X , Y ;

  Dot( const vector_type & argX ,
       const vector_type & argY )
  : X( argX ) , Y( argY ) {}

  Dot( const vector_type argX , const size_type argI ,
       const vector_type argY , const size_type argJ )
  : X( argX , argI ), Y( argY , argJ ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & update ) const
  { update += X(iwork) * Y(iwork); }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update , const volatile value_type & source )
  { update += source ; }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
};

//----------------------------------------------------------------------------
// Y(:) *= S ;
template< typename Scalar >
struct MultiVectorScale<Scalar, KOKKOS_MACRO_DEVICE ,1> {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;

  KokkosArray::MultiVector<Scalar,device_type> Y ;
  KokkosArray::Value<Scalar,device_type>     S ;

  MultiVectorScale(
    const KokkosArray::MultiVector<Scalar,device_type> & argY ,
    const KokkosArray::Value<Scalar,device_type>       & argS )
    : Y( argY ), S( argS ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) *= *S ; }
};

//----------------------------------------------------------------------------
// Y(:) += S * X(:)
template< typename Scalar >
struct MultiVectorYPAX<Scalar, KOKKOS_MACRO_DEVICE ,1> {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;

  KokkosArray::MultiVector<Scalar,device_type> Y , X ;
  KokkosArray::Value<Scalar,device_type>     A ;

  MultiVectorYPAX( const KokkosArray::MultiVector<Scalar,device_type> & argY ,
                   const KokkosArray::Value<Scalar,device_type>     & argA ,
                   const KokkosArray::MultiVector<Scalar,device_type> & argX )
   : Y( argY ), X( argX ), A( argA ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) += X(iwork) * *A ; }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct ModifiedGramSchmidt ;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar >
struct ModifiedGramSchmidt< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE                     device_type ;
  typedef device_type::size_type                  size_type ;
  typedef KokkosArray::MultiVector<Scalar,device_type> multivector_type ;
  typedef KokkosArray::Value<Scalar,device_type>       Value ;

  typedef KokkosArray::DotSingle< Scalar , device_type , 1 > DotSingle ;
  typedef KokkosArray::Dot< Scalar , device_type , 1 > Dot ;
  typedef KokkosArray::MultiVectorScale< Scalar , device_type , 1 > Scale ;
  typedef KokkosArray::MultiVectorYPAX<  Scalar , device_type , 1 > YPAX ;

  // Reduction   : result = dot( Q(:,j) , Q(:,j) );
  // PostProcess : R(j,j) = result ; inv = 1 / result ;
  struct InvNorm2 {
    multivector_type R ;
    Value       inv ;
    size_type   j ;

    InvNorm2( const multivector_type & argR ,
              const Value       & argInv ,
              size_type argJ )
      : R( argR )
      , inv( argInv )
      , j( argJ )
      {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( Scalar & result ) const
    {
      const Scalar value = sqrt( result );
      R(j,j) = value ;
      if ( 0 < value ) {
        *inv = 1.0 / value ;
      }
      else {
        *inv = 0 ;
      }
    }
  };

  // PostProcess : tmp = - ( R(j,k) = result );
  struct DotM {

    multivector_type R ;
    Value       tmp ;
    size_type j , k ;

    DotM( const multivector_type & argR ,
          const Value       & argTmp ,
          size_type argJ ,
          size_type argK )
      : R( argR )
      , tmp( argTmp )
      , j( argJ )
      , k( argK )
      {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( Scalar & result ) const
    {
       R(j,k) = result ;
       *tmp   = - result ;
    }
  };

  multivector_type Q ;
  multivector_type R ;
  double seconds ;

  ModifiedGramSchmidt( const typename multivector_type::HostMirror & A )
  : Q( KokkosArray::create_multivector<multivector_type>( A.length(), A.count()) )
  , R( KokkosArray::create_multivector<multivector_type>( A.count() , A.count()) )
  {
    const size_type N = A.length();
    Value tmp = KokkosArray::create_value<Value>();

    KokkosArray::deep_copy( Q , A );

    KokkosArray::Impl::Timer timer ;

    for ( size_type j = 0 ; j < Q.count() ; ++j ) {
      // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      // PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
      KokkosArray::parallel_reduce( N , DotSingle( Q , j ),
                                   InvNorm2( R , tmp , j ) );

      // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
      KokkosArray::parallel_for( N , Scale( multivector_type( Q , j ) , tmp ) );

      for ( size_t k = j + 1 ; k < Q.count() ; ++k ) {
        // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
        // PostProcess : tmp = - R(j,k);
        KokkosArray::parallel_reduce( N , Dot( Q , j , Q , k ) ,
                                     DotM( R , tmp , j , k ) );

        // Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
        KokkosArray::parallel_for( N , YPAX( multivector_type( Q , k ) , tmp , multivector_type( Q , j ) ) );
      }
    }

    device_type::fence();

    seconds = timer.seconds();
  }

  //--------------------------------------------------------------------------

  static double test( const size_t length ,
                      const size_t count )
  {
    typedef typename multivector_type::HostMirror HostMultiVector ;

    // Create and fill A on the host

    HostMultiVector A( KokkosArray::create_multivector<HostMultiVector>( "A" , length , count ) );

    for ( size_type j = 0 ; j < A.count() ; ++j ) {
      for ( size_type i = 0 ; i < A.length() ; ++i ) {
        A(i,j) = ( i + 1 ) * ( j + 1 );
      }
    }

    multivector_type Q ;
    multivector_type R ;

    double seconds ;
    {
      ModifiedGramSchmidt factorization( A );

      Q = factorization.Q ; // Save a view
      R = factorization.R ; // Save a view

      seconds = factorization.seconds ;
    }
    // A = Q * R

    return seconds ;
  }

};


