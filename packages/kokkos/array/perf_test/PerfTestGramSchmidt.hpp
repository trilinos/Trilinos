/*--------------------------------------------------------------------*/
/*    Copyright 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef UNIT_TEST_GRAM_SCHMIDT_OPERATORS
#define UNIT_TEST_GRAM_SCHMIDT_OPERATORS

namespace Kokkos {

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

namespace Kokkos {

template< typename Scalar >
struct DotSingle<Scalar, KOKKOS_MACRO_DEVICE , 1 > {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;
  typedef Kokkos::MultiVectorView<value_type,device_type> vector_type ;

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
  static void join( value_type & update , const value_type & source )
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
  typedef Kokkos::MultiVectorView<value_type,device_type> vector_type ;

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
  static void join( value_type & update , const value_type & source )
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

  Kokkos::MultiVectorView<Scalar,device_type> Y ;
  Kokkos::ValueView<Scalar,device_type>     S ;

  MultiVectorScale(
    const Kokkos::MultiVectorView<Scalar,device_type> & argY ,
    const Kokkos::ValueView<Scalar,device_type>       & argS )
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

  Kokkos::MultiVectorView<Scalar,device_type> Y , X ;
  Kokkos::ValueView<Scalar,device_type>     A ;

  MultiVectorYPAX( const Kokkos::MultiVectorView<Scalar,device_type> & argY ,
                   const Kokkos::ValueView<Scalar,device_type>     & argA ,
                   const Kokkos::MultiVectorView<Scalar,device_type> & argX )
   : Y( argY ), X( argX ), A( argA ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) += X(iwork) * *A ; }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct ModifiedGramSchmidt ;

template< typename Scalar , class DeviceType >
struct MultiVectorTestFill ;

//----------------------------------------------------------------------------

template< typename Scalar >
struct MultiVectorTestFill< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  Kokkos::MultiVectorView< Scalar , device_type > A ;

  MultiVectorTestFill( const Kokkos::MultiVectorView< Scalar , device_type > & argA )
    : A( argA ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  {
    for ( size_type j = 0 ; j < A.count() ; ++j ) {
      A(iwork,j) = ( iwork + 1 ) * ( j + 1 );
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar >
struct ModifiedGramSchmidt< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE                         device_type ;
  typedef device_type::size_type                      size_type ;
  typedef Kokkos::MultiVectorView<Scalar,device_type> MultiVector ;
  typedef Kokkos::ValueView<Scalar,device_type>       Value ;

  typedef Kokkos::DotSingle< Scalar , device_type , 1 > DotSingle ;
  typedef Kokkos::Dot< Scalar , device_type , 1 > Dot ;
  typedef Kokkos::MultiVectorScale< Scalar , device_type , 1 > Scale ;
  typedef Kokkos::MultiVectorYPAX<  Scalar , device_type , 1 > YPAX ;

  // Reduction   : result = dot( Q(:,j) , Q(:,j) );
  // PostProcess : R(j,j) = result ; inv = 1 / result ;
  struct InvNorm2 {
    MultiVector R ;
    Value       inv ;
    size_type   j ;

    InvNorm2( const MultiVector & argR ,
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

    MultiVector R ;
    Value       tmp ;
    size_type j , k ;

    DotM( const MultiVector & argR ,
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

  MultiVector Q ;
  MultiVector R ;

  ModifiedGramSchmidt( const Kokkos::MultiVectorView<Scalar,Kokkos::DeviceHost> & A )
    : Q( Kokkos::create_multivector<Scalar,device_type>( A.length() , A.count() ) )
    , R( Kokkos::create_multivector<Scalar,device_type>( A.count() , A.count() ) )
  {
    const size_type N = A.length();
    Value tmp = Kokkos::create_value<Scalar,device_type>();

    Kokkos::deep_copy( Q , A );

    for ( size_type j = 0 ; j < Q.count() ; ++j ) {
      // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      // PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
      Kokkos::parallel_reduce( N , DotSingle( Q , j ),
                                   InvNorm2( R , tmp , j ) );

      // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
      Kokkos::parallel_for( N , Scale( MultiVector( Q , j ) , tmp ) );

      for ( size_t k = j + 1 ; k < Q.count() ; ++k ) {
        // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
        // PostProcess : tmp = - R(j,k);
        Kokkos::parallel_reduce( N , Dot( Q , j , Q , k ) ,
                                     DotM( R , tmp , j , k ) );

        // Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
        Kokkos::parallel_for( N , YPAX( MultiVector( Q , k ) , tmp , MultiVector( Q , j ) ) );
      }
    }
  }

  //--------------------------------------------------------------------------

  static void test( const size_t length ,
                    const size_t count )
  {
    typedef Kokkos::MultiVectorView<Scalar,Kokkos::DeviceHost> HostMultiVector ;

    // Create and fill A on the host

    HostMultiVector A( Kokkos::create_labeled_multivector<HostMultiVector>( "A" , length , count ) );

    Kokkos::parallel_for( length , MultiVectorTestFill<Scalar,device_type>( A ) );

    MultiVector Q ;
    MultiVector R ;

    {
      ModifiedGramSchmidt factorization( A );

      Q = factorization.Q ; // Save a view
      R = factorization.R ; // Save a view
    }
    // A = Q * R
  }

};


