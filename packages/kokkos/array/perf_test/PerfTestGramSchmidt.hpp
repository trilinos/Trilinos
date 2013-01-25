/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <cmath>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct ModifiedGramSchmidt
{
  typedef DeviceType  device_type ;
  typedef typename device_type::size_type  size_type ;

  typedef KokkosArray::View< Scalar** ,
                             KokkosArray::LayoutLeft ,
                             device_type > multivector_type ;

  typedef KokkosArray::View< Scalar* ,
                             KokkosArray::LayoutLeft ,
                             device_type > vector_type ;

  typedef KokkosArray::View< Scalar ,
                             KokkosArray::LayoutLeft ,
                             device_type > value_view ;

  // Reduction   : result = dot( Q(:,j) , Q(:,j) );
  // PostProcess : R(j,j) = result ; inv = 1 / result ;
  struct InvNorm2 {
    typedef Scalar value_type ;
    value_view  Rjj ;
    value_view  inv ;

    InvNorm2( const value_view & argR ,
              const value_view & argInv )
      : Rjj( argR )
      , inv( argInv )
      {}

    KOKKOSARRAY_INLINE_FUNCTION
    void operator()( const value_type & result ) const
    {
      const Scalar value = sqrt( result );
      *Rjj = value ;
      *inv = ( 0 < value ) ? 1.0 / value : 0 ;
    }
  };

  // PostProcess : tmp = - ( R(j,k) = result );
  struct DotM {
    typedef Scalar value_type ;

    value_view  Rjk ;
    value_view  tmp ;

    DotM( const value_view & argR ,
          const value_view & argTmp )
      : Rjk( argR )
      , tmp( argTmp )
      {}

    KOKKOSARRAY_INLINE_FUNCTION
    void operator()( const value_type & result ) const
    {
       *Rjk  = result ;
       *tmp  = - result ;
    }
  };

  multivector_type Q ;
  multivector_type R ;

  static double factorization( const multivector_type Q ,
                               const multivector_type R )
  {
    const size_type count  = Q.dimension_1();
    value_view tmp("tmp");
    value_view one("one");

    KokkosArray::deep_copy( one , (Scalar) 1 );

    KokkosArray::Impl::Timer timer ;

    for ( size_type j = 0 ; j < count ; ++j ) {
      // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      // PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
      const vector_type Qj( Q , j );
      const value_view  Rjj( R , j , j );

      KokkosArray::dot( Qj , InvNorm2( Rjj , tmp  ) );

      // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
      KokkosArray::scale( tmp , Qj );

      for ( size_t k = j + 1 ; k < count ; ++k ) {
        const vector_type Qk( Q , k );
        const value_view  Rjk( R , j , k );

        // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
        // PostProcess : tmp = - R(j,k);
        KokkosArray::dot( Qj , Qk , DotM( Rjk , tmp ) );

        // Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
        KokkosArray::axpby( tmp , Qj , one , Qk );
      }
    }

    device_type::fence();

    return timer.seconds();
  }

  //--------------------------------------------------------------------------

  static double test( const size_t length ,
                      const size_t count ,
                      const size_t iter = 1 )
  {
    multivector_type Q( "Q" , length , count );
    multivector_type R( "R" , count , count );

    typename multivector_type::HostMirror A =
      KokkosArray::create_mirror( Q );

    // Create and fill A on the host

    for ( size_type j = 0 ; j < count ; ++j ) {
      for ( size_type i = 0 ; i < length ; ++i ) {
        A(i,j) = ( i + 1 ) * ( j + 1 );
      }
    }

    double dt_min = 0 ;

    for ( size_t i = 0 ; i < iter ; ++i ) {

      KokkosArray::deep_copy( Q , A );

      // A = Q * R

      const double dt = factorization( Q , R );

      if ( 0 == i ) dt_min = dt ;
      else dt_min = dt < dt_min ? dt : dt_min ;
    }

    return dt_min ;
  }
};


