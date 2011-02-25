/*--------------------------------------------------------------------*/
/*    Copyright 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

namespace Kokkos {

template< typename Scalar , class DeviceMap , unsigned N = 1 >
struct MVectorScale ;

template< typename Scalar , class DeviceMap , unsigned N = 1 >
struct MVectorYPAX ;

template< typename Scalar , class Device , unsigned N = 1 >
struct ReduceSum ;

//----------------------------------------------------------------------------
// Types and static methods for summation

template< typename Scalar , class Device >
struct ReduceSum<Scalar,Device,1> {
  typedef Device                             device_type ;
  typedef Kokkos::MDArrayView<Scalar,Device> reduce_type ;

  KOKKOS_DEVICE_FUNCTION
  static void join( const reduce_type & result ,
                    const reduce_type & source )
  { result(0) += source(0); }

  KOKKOS_DEVICE_FUNCTION
  static void init( const reduce_type & result )
  { result(0) = 0 ; }
};

//----------------------------------------------------------------------------
// Y(:) *= S ;
template< typename Scalar , class DeviceMap >
struct MVectorScale<Scalar,DeviceMap,1> {
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;

  Kokkos::MVectorView<Scalar,DeviceMap>   Y ;
  Kokkos::MDArrayView<Scalar,device_type> S ;

  MVectorScale( const Kokkos::MVectorView<Scalar,DeviceMap>   & argY ,
                const Kokkos::MDArrayView<Scalar,device_type> & argS )
    : Y( argY ), S( argS ) {}

  KOKKOS_DEVICE_FUNCTION
  size_type work_count() const { return Y.length(); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) *= S(0); }
};

//----------------------------------------------------------------------------
// Y(:) += S * X(:)
template< typename Scalar , class DeviceMap >
struct MVectorYPAX<Scalar,DeviceMap,1> {
  // Required API typedef:
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;

  Kokkos::MVectorView<Scalar,DeviceMap>   Y , X ;
  Kokkos::MDArrayView<Scalar,device_type> A ;

  MVectorYPAX( const Kokkos::MVectorView<Scalar,DeviceMap>   & argY ,
               const Kokkos::MDArrayView<Scalar,device_type> & argA ,
               const Kokkos::MVectorView<Scalar,DeviceMap>   & argX )
   : Y( argY ), X( argX ), A( argA ) {}

  KOKKOS_DEVICE_FUNCTION
  size_type work_count() const { return Y.length(); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) += X(iwork) * A(0); }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class DeviceMap >
struct ModifiedGramSchmidt {
  // My convenience typedef:
  typedef typename DeviceMap::device_type        Device ;
  typedef Kokkos::MVectorView<Scalar,DeviceMap>  MVector ;
  typedef Kokkos::MDArrayView<Scalar,Device>     MDArray ;
  typedef typename Device::size_type             size_type ;
  typedef MDArray                                reduce_type ;

  typedef Kokkos::MVectorScale< Scalar , DeviceMap , 1 > Scale ;
  typedef Kokkos::MVectorYPAX<  Scalar , DeviceMap , 1 > YPAX ;

  // Reduction   : result = dot( Q(:,j) , Q(:,j) );
  // PostProcess : R(j,j) = result ; inv = 1 / result ;
  struct InvNorm2 : public Kokkos::ReduceSum<Scalar,Device,1> {

    MVector QJ ;
    MDArray R ;
    MDArray inv ;
    size_type j ;

    InvNorm2( const MVector & argQ ,
              const MDArray & argR ,
              const MDArray & argInv ,
              size_type argJ )
      : QJ( argQ , argJ )
      , R( argR )
      , inv( argInv )
      , j( argJ )
      {}

    KOKKOS_DEVICE_FUNCTION
    size_type work_count() const { return QJ.length(); }

    KOKKOS_DEVICE_FUNCTION
    void operator()( size_type iwork , const reduce_type & result ) const
    {
      const Scalar value = QJ(iwork);
      result(0) += value * value ;
    }

    KOKKOS_DEVICE_FUNCTION
    void post_process( const reduce_type & result ) const
    {
      const Scalar value = Kokkos::sqrt( result(0) );
      R(j,j) = value ;
      if ( 0 < value ) {
        inv(0) = 1.0 / value ;
      }
      else {
        inv(0) = 0 ;
      }
    }
  };

  // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
  // PostProcess : tmp = - R(j,k);
  struct DotM : public Kokkos::ReduceSum<Scalar,Device,1> {

    MVector QJ , QK ;
    MDArray R ;
    MDArray tmp ;
    size_type j , k ;

    DotM( const MVector & argQ ,
          const MDArray & argR ,
          const MDArray & argTmp ,
          size_type argJ ,
          size_type argK )
      : QJ( argQ , argJ )
      , QK( argQ , argK )
      , R( argR )
      , tmp( argTmp )
      , j( argJ )
      , k( argK )
      {}

    KOKKOS_DEVICE_FUNCTION
    size_type work_count() const { return QJ.length(); }

    KOKKOS_DEVICE_FUNCTION
    void operator()( size_type iwork , const reduce_type & result ) const
    { result(0) += QJ(iwork) * QK(iwork); }

    KOKKOS_DEVICE_FUNCTION
    void post_process( const reduce_type & result ) const
    {
       R(j,k) = result(0);
       tmp(0) = - result(0);
    }
  };

  template< class HostMap >
  ModifiedGramSchmidt( const Kokkos::MVectorView<Scalar,HostMap> & A )
    : Q( Kokkos::create_mvector<Scalar,DeviceMap>( A.length() , A.count() ) )
    , R( Device::template create_mdarray<Scalar>( A.count() , A.count() ) )
    {
      MDArray tmp( Device::template create_mdarray<Scalar>( 1 ) );

      Kokkos::deep_copy( Q , A );

      for ( size_type j = 0 ; j < Q.count() ; ++j ) {
        // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
        // PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
        Kokkos::parallel_reduce( InvNorm2( Q , R , tmp , j ) , tmp );

        // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
        Kokkos::parallel_for( Scale( MVector( Q , j ) , tmp ) );

        for ( size_t k = j + 1 ; k < Q.count() ; ++k ) {
          // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
          // PostProcess : tmp = - R(j,k);
          Kokkos::parallel_reduce( DotM( Q , R , tmp , j , k ) , tmp );

          // Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
          Kokkos::parallel_for( YPAX( MVector( Q , k ) , tmp , MVector( Q , j ) ) );
        }
      }
    }

  MVector Q ;
  MDArray R ;
};

//----------------------------------------------------------------------------

template< typename Scalar , class DeviceMap >
struct MVectorTestFill {
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;

  Kokkos::MVectorView< Scalar , DeviceMap > A ;

  MVectorTestFill( const Kokkos::MVectorView< Scalar , DeviceMap > & argA )
    : A( argA ) {}

  KOKKOS_DEVICE_FUNCTION
  size_type work_count() const { return A.length(); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  {
    for ( size_type j = 0 ; j < A.count() ; ++j ) {
      A(iwork,j) = ( iwork + 1 ) * ( j + 1 );
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <Kokkos_HostMap.hpp>

template< typename Scalar , class DeviceMap >
void test_modified_gram_schmidt( const size_t length ,
                                 const size_t count )
{
  typedef Kokkos::HostMap                 HostMap ;
  typedef typename DeviceMap::device_type Device ;

  typedef Kokkos::MVectorView<Scalar,DeviceMap> MVector ;
  typedef Kokkos::MDArrayView<Scalar,Device>    MDArray ;

  // Create and fill A on the host

  MVector A( Kokkos::create_labeled_mvector<Scalar,HostMap>( length , count , "A" ) );

  Kokkos::parallel_for( MVectorTestFill<Scalar,DeviceMap>( A ) );

  MVector Q ;
  MDArray R ;

  {
    ModifiedGramSchmidt<Scalar,DeviceMap> factorization( A );

    Q = factorization.Q ; // Save a view
    R = factorization.R ; // Save a view
  }
  // A = Q * R
}



