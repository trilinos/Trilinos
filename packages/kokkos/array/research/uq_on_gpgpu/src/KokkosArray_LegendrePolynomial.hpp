/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_LEGENDREPOLYNOMIALS_HPP
#define KOKKOSARRAY_LEGENDREPOLYNOMIALS_HPP

#include <cmath>
#include <KokkosArray_Macros.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Evaluate normalized Legendre polynomial bases.
 *
 *
 *  P_{0}(x) = 1 ;
 *  P_{1}(x) = x ;
 *  P_{k}(x) = ((2*k-1)/k) * x * P_{k-1}(x) - ((k-1)/k) * P_{k-2}(x) ;
 *
 *  Normalized over the norm
 *    < P_{k}(x) , P_{n}(x) > = L2 norm on the interval [-1,1]
 *                            = kd_{k,n} * 2 / ( 2 * k + 1 ) 
 *
 */
template< unsigned MaximumDegree >
class NormalizedLegendrePolynomialBases {
public:

  template< typename Scalar >
  KOKKOSARRAY_INLINE_FUNCTION
  void evaluate( unsigned N , Scalar x , Scalar value[] ) const
  {
    double vkm1 = 0 ;
    double vk   = 1 ;

    value[0] = vk * m_norm[0];

    for ( unsigned k = 1 ; k <= N ; ++k ) {
      const double vkm2 = vkm1 ;
      vkm1 = vk ;
      vk = ( vkm1 * x * (2*k-1) - vkm2 * (k-1) ) / double(k);
      value[k] = vk * m_norm[k];
    }
  }

  NormalizedLegendrePolynomialBases()
  {
    for ( unsigned k = 0 ; k <= MaximumDegree ; ++k ) {
      m_norm[k] = std::sqrt( 1.0 + 2.0*k );
    }
  }

private:

  double m_norm[ MaximumDegree + 1 ];

};

//----------------------------------------------------------------------------

namespace Impl {
unsigned guass_legendre_available( unsigned n );
void     gauss_legendre( unsigned n , double point[] , double weight[] );
}

//----------------------------------------------------------------------------
/** \brief  Gauss-Legendre numerical integration points and weights
 *          for exact integration (within round-off) of
 *          polynomial degree 'PolyDegree' on the interval [-1,1].
 */
template< unsigned PolyDegree >
class GaussLegendre {
private:
  GaussLegendre( const GaussLegendre & );
  GaussLegendre operator = ( const GaussLegendre & );

  enum { Desired = unsigned( 1 + PolyDegree / 2 ) };

public:

  /** \brief  N = number of points and weights for available 
   *          integration rules that satisfy the desired degree.
   */

  enum { N = 1024 < Desired ?    0 : (
              512 < Desired ? 1024 : (
              256 < Desired ?  512 : (
              128 < Desired ?  256 : (
               64 < Desired ?  128 : (
               32 < Desired ?   64 : (
               20 < Desired ?   32 : Desired )))))) };

  ~GaussLegendre() {}
  GaussLegendre() { Impl::gauss_legendre( N , points , weights ); }

  double points[N] ;
  double weights[N] ;
};

} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_LEGENDREPOLYNOMIALS_HPP */


