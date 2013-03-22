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

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <utility>

#include <KokkosArray_Host.hpp>

#include <KokkosArray_LegendrePolynomial.hpp>
#include <KokkosArray_ProductTensorLegendre.hpp>
#include <KokkosArray_CrsProductTensorLegendre.hpp>

namespace unit_test {
namespace {

//----------------------------------------------------------------------------
// Verify integration cubature rules for manufactured polynomial

double poly( const double x , const double c , const int p , const int d )
{
  double y = 0 ;
  if ( 0 <= p ) {
    if ( d <= p ) y = 1 ;
    for ( int i = 0 ; i < p - d ; ++i ) { y *= ( x - c ); }
    for ( int i = 0 ; i < d ; ++i ) { y *= p - i ; }
  }
  return y ;
}

void gauss_legendre_poly( const unsigned ng , const unsigned p , const double c )
{
  double pt[ 1024 ] , wgt[ 1024 ] ;

  const double y_correct = poly(  1.0 , c , p+1 , 0 ) -
                           poly( -1.0 , c , p+1 , 0 );

  KokkosArray::Impl::gauss_legendre( ng , pt , wgt );

  double y = 0 ;

  for ( int i = 0 ; i < (int) ng ; ++i ) {
    y += wgt[i] * poly( pt[i] , c , p+1 , 1 );
  }

  // Allow 'tol' roundoff error
  const double tol = ng * 3.0e-16 ;
  const double e = std::abs( y - y_correct );
  const double r = std::abs( y_correct );
  const bool ok = tol < r ? e / r < tol : e < tol ;

  if ( ! ok ) {
    std::cout << "Test Integral[-1,1]{ " << p + 1
              << " (x-" << c << ")^" << p << " } = "
              << y_correct
              << " : Failed ng(" << ng 
              << ") error = " << y - y_correct ;

    if ( tol < r ) {
      std::cout << " , ratio = " << e / r ;
    }
    std::cout << std::endl ;
  }
}

//----------------------------------------------------------------------------

template< unsigned P , bool Print >
void normalized_legendre_polynomial_bases()
{
  //const double sqrt2 = std::sqrt( 2.0 );
  double pt[32] , wgt[32] ;
  double value[P+1];
  double inner_product[P+1][P+1];
  double integral[P+1];

  KokkosArray::NormalizedLegendrePolynomialBases<P> bases ;

  const unsigned p_integral = 1 + P * P / 2 ;

  if ( 32 < p_integral ) {
    std::cout << "Cannot perform inner product on polynomial degree " << P << std::endl ;
  }

  const unsigned ng_end = p_integral <= 20 ? 20 : 32 ;
  unsigned ng = p_integral <= 20 ? p_integral : 32 ;

  for ( ; ng <= ng_end ; ++ng ) {

    const double tol = ng * 3.0e-16 ;

    KokkosArray::Impl::gauss_legendre( ng , pt , wgt );

    // Inner product of all function pairs:
    for ( unsigned j = 0 ; j <= P ; ++j ) {
      integral[j] = 0 ;
    }
    for ( unsigned j = 0 ; j <= P ; ++j ) {
    for ( unsigned k = 0 ; k <= P ; ++k ) {
      inner_product[j][k] = 0 ;
    }}

    for ( unsigned i = 0 ; i < ng ; ++i ) {

      bases.evaluate( P , pt[i] , value );

      for ( unsigned j = 0 ; j <= P ; ++j ) {
        integral[j] += 0.5 * wgt[i] * value[j] ;
      }
      for ( unsigned j = 0 ; j <= P ; ++j ) {
      for ( unsigned k = 0 ; k <= P ; ++k ) {
        inner_product[j][k] += 0.5 * wgt[i] * value[j] * value[k] ;
      }}
    }

    // Verify orthonormal
    for ( unsigned j = 0 ; j <= P ; ++j ) {
    for ( unsigned k = 0 ; k <= P ; ++k ) {
      const double e = std::abs( inner_product[j][k] - ( j == k ? 1 : 0 ) );

      if ( tol < e ) {
        std::cout << "inner_product failed "
                  << "<" << j << "," << k << "> ng(" << ng << ") = " << e
                  << std::endl ;
      }
    }}

    // Verify integrals
    for ( unsigned j = 0 ; j <= P ; ++j ) {
      const double e = std::abs( integral[j] - ( 0 == j ? 1 : 0 ) );

      if ( Print ) {
        std::cout << "  integral <" << j << "> = " << integral[j] << std::endl ;
      }
 
      if ( tol < e ) {
        std::cout << "integral failed "
                  << "<" << j << "> ng(" << ng << ") = " << e
                  << std::endl ;
      }
    }
  }
}

//----------------------------------------------------------------------------

template< unsigned P , bool Print >
void product_tensor_integration()
{
  const KokkosArray::NormalizedLegendrePolynomialBases<P> bases ;
  const KokkosArray::GaussLegendre<6*P> gauss ;
  // const double tol = gauss.N * 3.0e-16 ;
  const double zero_tol  = 1.0e-14 ;
  const double float_tol = 3.0e-8 ;

  double integral[ P + 1 ][ P + 1 ][ P + 1 ];
  double scale = 0 ;

  for ( unsigned i = 0 ; i <= P ; ++i ) {
  for ( unsigned j = 0 ; j <= P ; ++j ) {
  for ( unsigned k = 0 ; k <= P ; ++k ) {
    integral[i][j][k] = 0 ;
  }}}

  for ( unsigned ig = 0 ; ig < gauss.N ; ++ig ) {
    double value[ P + 1 ];

    bases.evaluate( P , gauss.points[ig] , value );

    // Integral of <0,i,i> should equal integral <0>

    scale += gauss.weights[ig] * value[0] ;
    for ( unsigned i = 0 ; i <= P ; ++i ) {
    for ( unsigned j = 0 ; j <= P ; ++j ) {
    for ( unsigned k = 0 ; k <= P ; ++k ) {
      integral[i][j][k] += 0.5 * gauss.weights[ig] * value[i] * value[j] * value[k] ;
    }}}
  }

  enum { DegreeTensor = 7 };
  KokkosArray::TripleProductTensorLegendre tensor ;

  for ( unsigned i = 0 ; i <= P ; ++i ) {
  for ( unsigned j = 0 ; j <= P ; ++j ) {
  for ( unsigned k = 0 ; k <= P ; ++k ) {
    if ( i <= DegreeTensor && j <= DegreeTensor && k <= DegreeTensor ) {
      if ( float_tol < std::fabs( integral[i][j][k] - tensor(i,j,k) ) ) {
        std::cout << "product_tensor_integration failed "
                  << "  <" << i << "," << j << "," << k << "> integral( "
                  << integral[i][j][k]
                  << " ) = tensor( " << tensor(i,j,k)
                  << " ) = " << integral[i][j][k] - tensor(i,j,k)
                  << std::endl ;
      }
    }
  }}}

  if ( Print ) {

    unsigned nonzero_count = 0 ;

    for ( unsigned i = 0 ; i <= P ; ++i ) {
    for ( unsigned j = i ; j <= P ; ++j ) {
    for ( unsigned k = j ; k <= P ; ++k ) {
      if ( zero_tol < std::fabs( integral[i][j][k] ) ) {
        ++nonzero_count ;
      }
    }}}

    std::cout << "unit_test::product_tensor_integration<" << P
              << "> nonzero_count = " << nonzero_count
              << std::endl ;

    for ( unsigned i = 0 ; i <= P ; ++i ) {
    for ( unsigned j = i ; j <= P ; ++j ) {
    for ( unsigned k = j ; k <= P ; ++k ) {
      if ( zero_tol < std::fabs( integral[i][j][k] ) ) {
        std::cout << "  <" << i << "," << j << "," << k << "> = "
                  << integral[i][j][k] << std::endl ;
      }
    }}}
  }
}

//----------------------------------------------------------------------------

template< bool Print >
void product_tensor_bases( const unsigned nvar ,
                           const unsigned poly_deg ,
                           const unsigned poly_deg_max )
{
  std::vector<unsigned> vdeg( nvar , poly_deg );

  KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation 
    combinatorial( vdeg , poly_deg_max );

  //------------------------------------
  // Triple products of the bases:

  if ( Print ) {
    std::cout << "triple_product_tensor_evaluation row nonzeros:" << std::endl ;
  }

  unsigned count_diag = 0 ;
  unsigned count_off_diag = 0 ;

  for ( unsigned i = 0 ; i < combinatorial.bases_count() ; ++i ) {
    const unsigned prev_count_diag     = count_diag ;
    const unsigned prev_count_off_diag = count_off_diag ;
    count_diag = 0 ;
    count_off_diag = 0 ;

    for ( unsigned j = 0 ; j < combinatorial.bases_count() ; ++j ) {

      if ( combinatorial.is_non_zero(i,j,j) ) ++count_diag ;

      for ( unsigned k = j + 1 ; k < combinatorial.bases_count() ; ++k ) {
        if ( combinatorial.is_non_zero(i,j,k) ) count_off_diag += 1 ;
      }
    }

    if ( Print ) {
      std::cout << "  {" ;
      for ( unsigned j = 0 ; j < combinatorial.variable_count() ; ++j ) {
        if ( j ) std::cout << "," ;
        std::cout << (unsigned) combinatorial[i][j] ;
      }
      std::cout << "} " ;
      std::cout << " diag_count = " << count_diag
                << " , off_diag_count = " << count_off_diag
                << std::endl ;
    }

    if ( i ) {
      const unsigned prev_count_total = prev_count_diag + prev_count_off_diag ;
      const unsigned count_total      = count_diag + count_off_diag ;

      const bool ok = prev_count_total  != count_total ? prev_count_total > count_total : (
                      prev_count_diag   != count_diag  ? prev_count_diag  > count_diag : true );

      if ( ! ok ) {
        std::cout << "product_tensor_bases FAILED : bad ordering at " << i << std::endl ;
      }
    }
  }
}


} // namespace
} // namespace unit_test

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace unit_test {

//----------------------------------------------------------------------------

void product_tensor_legendre()
{
  //--------------------------------
  // Verify correct integration from [-1,1]
  // of manufactured polynomial function
  //     f(x) = ( x - c )^p
  // For ranges of 'p' and 'c'.

  for ( int j = 0 ; j <= 32 ; ++j ) {
    const double c = double(j-16) / 8.0 ;
    for ( unsigned p = 0 ; p <= 32 ; ++p ) {
      for ( unsigned ng = 1 + p / 2 ; ng <= 20 ; ++ng ) {
        gauss_legendre_poly( ng , p , c );
      }
      gauss_legendre_poly(   32 , p , c );
      gauss_legendre_poly(   64 , p , c );
      gauss_legendre_poly(  128 , p , c );
      gauss_legendre_poly(  256 , p , c );
      gauss_legendre_poly(  512 , p , c );
      gauss_legendre_poly( 1024 , p , c );
    }
  }
  //--------------------------------
  // Verify orthonormality of bases.
  
  normalized_legendre_polynomial_bases<2,false>();
  normalized_legendre_polynomial_bases<5,false>();
  normalized_legendre_polynomial_bases<7,false>();

  //--------------------------------

  product_tensor_integration<3,false>();
  product_tensor_integration<4,false>();
  product_tensor_integration<5,false>();
  product_tensor_integration<6,false>();
  product_tensor_integration<7,false>();
  product_tensor_integration<12,false>();

  //--------------------------------

  {
    const unsigned nvar = 6 ;
    const unsigned pdeg = 5 ;
    const unsigned pmax = 5 ;
    product_tensor_bases<false>( nvar , pdeg , pmax );
  }
}

} // name unit_test

