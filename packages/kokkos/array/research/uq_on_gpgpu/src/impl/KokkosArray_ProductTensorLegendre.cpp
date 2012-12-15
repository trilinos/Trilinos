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

#include <sstream>
#include <stdexcept>

#include <KokkosArray_ProductTensorLegendre.hpp>
#include <KokkosArray_LegendrePolynomial.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief
 *
 *
 */
TripleProductTensorLegendre::TripleProductTensorLegendre()
{
  enum { IntegDegree = 3 * MaximumPolyDegree };

  const NormalizedLegendrePolynomialBases< MaximumPolyDegree > bases ;
  const GaussLegendre< IntegDegree > gauss ;
  const double tol = 1.0e-14 ;

  double tmp[ N * N * N ];

  for ( unsigned i = 0 ; i < N*N*N ; ++i ) tmp[i] = 0 ;
  for ( unsigned i = 0 ; i < N*N*N ; ++i ) m_map[i] = 0 ;

  for ( unsigned ig = 0 ; ig < gauss.N ; ++ig ) {
    double value[N] ;
    bases.evaluate( MaximumPolyDegree , gauss.points[ig] , value );

    unsigned entry_count = 0 ;
    for ( unsigned i = 0 ; i <= MaximumPolyDegree ; ++i ) {
    for ( unsigned j = i ; j <= MaximumPolyDegree ; ++j ) {
    for ( unsigned k = j ; k <= MaximumPolyDegree ; ++k ) {
      tmp[ offset(i,j,k) ] += gauss.weights[ig] * value[i] * value[j] * value[k] ;
    }}}
  }

  m_terms[0] = 0 ; // For all zero terms
  m_terms[1] = tmp[ offset(0,0,0) ]; // For all <0,i,i> terms

  for ( unsigned i = 0 ; i <= MaximumPolyDegree ; ++i ) {
    m_map[ offset(0,i,i) ] = 1 ;
    m_map[ offset(i,0,i) ] = 1 ;
    m_map[ offset(i,i,0) ] = 1 ;
  }

  unsigned count = 2 ;

  for ( unsigned i = 1 ; i <= MaximumPolyDegree ; ++i ) {
  for ( unsigned j = i ; j <= MaximumPolyDegree ; ++j ) {
  for ( unsigned k = j ; k <= MaximumPolyDegree ; ++k ) {
    if ( tol < std::abs( tmp[ offset(i,j,k) ] ) ) {

      if ( NONZERO_COUNT <= count ) {
        std::ostringstream msg ;
        msg << "KokkosArray::TripleProductTensorLegendre FATAL ERROR: "
            << "non-zero terms exceeded previously determined value "
            << NONZERO_COUNT ;
        throw std::runtime_error( msg.str() );
      }

      m_terms[ count ] = tmp[ offset(i,j,k) ];

      m_map[ offset(i,j,k) ] = count ;
      m_map[ offset(i,k,j) ] = count ;
      m_map[ offset(j,k,i) ] = count ;
      m_map[ offset(j,i,k) ] = count ;
      m_map[ offset(k,i,j) ] = count ;
      m_map[ offset(k,j,i) ] = count ;
      ++count ;
    }
  }}}
}

//----------------------------------------------------------------------------

}

