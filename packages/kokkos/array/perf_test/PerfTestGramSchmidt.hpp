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

#include <cmath>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct ModifiedGramSchmidt ;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar >
struct ModifiedGramSchmidt< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  typedef KokkosArray::View< Scalar[][0] ,
                             KokkosArray::LayoutLeft ,
                             device_type > multivector_type ;

  typedef KokkosArray::View< Scalar[] ,
                             KokkosArray::LayoutLeft ,
                             device_type > vector_type ;

  typedef KokkosArray::View< Scalar ,
                             KokkosArray::LayoutLeft ,
                             device_type > value_type ;

  // Reduction   : result = dot( Q(:,j) , Q(:,j) );
  // PostProcess : R(j,j) = result ; inv = 1 / result ;
  struct InvNorm2 {
    value_type  Rjj ;
    value_type  inv ;

    InvNorm2( const value_type & argR ,
              const value_type & argInv )
      : Rjj( argR )
      , inv( argInv )
      {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( Scalar & result ) const
    {
      const Scalar value = sqrt( result );
      *Rjj = value ;
      *inv = ( 0 < value ) ? 1.0 / value : 0 ;
    }
  };

  // PostProcess : tmp = - ( R(j,k) = result );
  struct DotM {

    value_type  Rjk ;
    value_type  tmp ;

    DotM( const value_type & argR ,
          const value_type & argTmp )
      : Rjk( argR )
      , tmp( argTmp )
      {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( Scalar & result ) const
    {
       *Rjk  = result ;
       *tmp  = - result ;
    }
  };

  multivector_type Q ;
  multivector_type R ;
  double seconds ;

  static double factorization( const multivector_type Q ,
                               const multivector_type R )
  {
    const size_type count  = Q.dimension_1();
    value_type tmp("tmp");
    value_type one("one");

    KokkosArray::deep_copy( one , (Scalar) 1 );

    KokkosArray::Impl::Timer timer ;

    for ( size_type j = 0 ; j < count ; ++j ) {
      // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      // PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
      vector_type Qj  = KokkosArray::view( Q , j );
      value_type Rjj = KokkosArray::view( R , j , j );

      KokkosArray::dot( Qj , InvNorm2( Rjj , tmp  ) );

      // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
      KokkosArray::scale( tmp , Qj );

      for ( size_t k = j + 1 ; k < count ; ++k ) {
        vector_type Qk  = KokkosArray::view( Q , k );
        value_type Rjk = KokkosArray::view( R , j , k );

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
                      const size_t count )
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

    KokkosArray::deep_copy( Q , A );

    // A = Q * R

    return factorization( Q , R );
  }

};


