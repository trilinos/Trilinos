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

#include <iostream>

#include <sstream>
#include <stdexcept>
#include <algorithm>

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

namespace {

void product_tensor_bases( const std::vector<unsigned> & poly_deg_var ,
                           const unsigned poly_deg_total ,
                           const unsigned poly_deg_target ,
                           const unsigned ivar ,
                           std::vector<unsigned> & iv ,
                           std::vector<unsigned char> & bases_map )
{
  const unsigned nvar = poly_deg_var.size();
  const unsigned max  = std::min( poly_deg_var[ivar] , std::min( poly_deg_total , poly_deg_target ) );

  for ( int p = max ; 0 <= p ; --p ) {
    iv[ivar] = p ;

    if ( ivar + 1 < nvar ) {
      product_tensor_bases( poly_deg_var , poly_deg_total - p , poly_deg_target , ivar + 1 , iv ,
                            bases_map );
    }
    else if ( poly_deg_total == p ) {

      // No variables exceeded poly_deg_target
      // At least one variable must be equal to poly_deg_target

      bool ok = false ;

      for ( unsigned j = 0 ; j < nvar && ! ( ok = poly_deg_target == iv[j] ) ; ++j );

      if ( ok ) {
        bases_map.insert( bases_map.end() , iv.begin() , iv.end() );
      }
    }
  }
}

struct ProductTensorRowCount {
  unsigned count_diag ;
  unsigned count_off_diag ;
  unsigned basis ;

  ProductTensorRowCount()
  : count_diag(0)
  , count_off_diag(0)
  , basis(0)
  {}
};

struct CompareProductTensorRowCount {

  bool operator()( const ProductTensorRowCount & lhs ,
                   const ProductTensorRowCount & rhs ) const
  {
    return
           lhs.count_diag     != rhs.count_diag     ? lhs.count_diag     > rhs.count_diag : (
           lhs.count_off_diag != rhs.count_off_diag ? lhs.count_off_diag > rhs.count_off_diag : (
           lhs.basis < rhs.basis ));

  }
};

}

TripleProductTensorLegendreCombinatorialEvaluation::
TripleProductTensorLegendreCombinatorialEvaluation(
  const unsigned                  maximum_polynomial_degree ,
  const std::vector< unsigned > & variable_polynomial_degree )
  : m_tensor()
  , m_bases_map()
  , m_variable_count( variable_polynomial_degree.size() )
  , m_bases_count(0)
  , m_bases_count_with_diagonal(0)
{
  // Generate a combinatorial bases map to start the ordering process.
  // This ordering will be used when two rows' have equal non-zero density.

  std::vector< unsigned char > tmp_bases_map ;

  {
    std::vector<unsigned> iv( m_variable_count , 0u );
    for ( unsigned d = 0 ; d <= maximum_polynomial_degree ; ++d ) {
      for ( unsigned p = 0 ; p <= d ; ++p ) {
        product_tensor_bases( variable_polynomial_degree , d , p , 0 , iv , tmp_bases_map );
      }
    }
  }

  m_bases_count = tmp_bases_map.size() / m_variable_count ;

  std::vector< ProductTensorRowCount > row_count( m_bases_count );

  // Sorting largest to smallest...

  for ( unsigned i = 0 ; i < m_bases_count ; ++i ) {
    row_count[i].basis = i ;

    // <i,i,i>
    if ( m_tensor.is_non_zero( m_variable_count ,
                               & tmp_bases_map[i*m_variable_count] ,
                               & tmp_bases_map[i*m_variable_count] ,
                               & tmp_bases_map[i*m_variable_count] ) ) {
      row_count[i].count_diag += 1 ;
    }

    for ( unsigned j = i + 1 ; j < m_bases_count ; ++j ) {

      // <i,j,j> , i < j < m_bases_count
      if ( m_tensor.is_non_zero( m_variable_count ,
                                 & tmp_bases_map[i*m_variable_count] ,
                                 & tmp_bases_map[j*m_variable_count] ,
                                 & tmp_bases_map[j*m_variable_count] ) ) {
        row_count[i].count_diag += 1 ;
        row_count[j].count_off_diag += 2 ;
      }

      // <i,i,j> , i < j < m_bases_count
      if ( m_tensor.is_non_zero( m_variable_count ,
                                 & tmp_bases_map[i*m_variable_count] ,
                                 & tmp_bases_map[i*m_variable_count] ,
                                 & tmp_bases_map[j*m_variable_count] ) ) {
        row_count[j].count_diag += 1 ;
        row_count[i].count_off_diag += 2 ;
      }

      // <i,j,k> , i < j < k < m_bases_count
      for ( unsigned k = j + 1 ; k < m_bases_count ; ++k ) {
        if ( m_tensor.is_non_zero( m_variable_count ,
                                   & tmp_bases_map[i*m_variable_count] ,
                                   & tmp_bases_map[j*m_variable_count] ,
                                   & tmp_bases_map[k*m_variable_count] ) ) {
          row_count[i].count_off_diag += 2 ;
          row_count[j].count_off_diag += 2 ;
          row_count[k].count_off_diag += 2 ;
        }
      }
    }
  }

  std::sort( row_count.begin() , row_count.end() , CompareProductTensorRowCount() );

  // Count diagonals

  while ( m_bases_count_with_diagonal < m_bases_count &&
          row_count[ m_bases_count_with_diagonal ].count_diag ) {
    ++m_bases_count_with_diagonal ;
  }

  // Generate map:

  m_bases_map.resize( tmp_bases_map.size() );

  for ( unsigned i = 0 ; i < m_bases_count ; ++i ) {

    const unsigned ib = row_count[i].basis ;

    for ( unsigned j = 0 ; j < m_variable_count ; ++j ) {
      m_bases_map[ i * m_variable_count + j ] =
        tmp_bases_map[ ib * m_variable_count + j ];
    }

#if 0

    std::cout << "  {" ;
    for ( unsigned j = 0 ; j < variable_count() ; ++j ) {
      if ( j ) std::cout << "," ;
      std::cout << (unsigned) m_bases_map[i*m_variable_count+j] ;
    }
    std::cout << "} " ;
    std::cout << " diag_count = " << row_count[i].count_diag
              << " , off_diag_count = " << row_count[i].count_off_diag
              << " , old_index = " << ib
              << std::endl ;

#endif

  }
}

}

