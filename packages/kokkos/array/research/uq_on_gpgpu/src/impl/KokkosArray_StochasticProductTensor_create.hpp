/*
//@HEADER
// ************************************************************************
// 
//           Kokkos: Node API and Parallel Node Kernels
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_STOCHASTICPRODUCTTENSORCREATE_HPP
#define KOKKOSARRAY_STOCHASTICPRODUCTTENSORCREATE_HPP

#include <iostream>
#include <cmath>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType , class PolynomialType ,
          class Device , typename IntType ,
          template< unsigned , typename , class > class TensorType >
class CreateSparseProductTensor<
  StochasticProductTensor< ValueType , PolynomialType , Device , TensorType > ,
  std::vector< IntType > >
{
public:

  typedef Device device_type ;
  typedef typename device_type::size_type size_type ;
  typedef StochasticProductTensor< ValueType ,
                                   PolynomialType ,
                                   device_type ,
                                   TensorType > type ;

  typedef std::vector< IntType > input_type ;

private:

  typedef ProductTensorIndex< 3 , device_type > product_tensor_index_type ;
  typedef typename type::tensor_type tensor_type ;
  typedef std::map< product_tensor_index_type , ValueType >  tensor_map_type ;

  PolynomialType      m_poly ;
  tensor_map_type     m_tensor_map ;
  std::vector< int >  m_variable_degree ;
  std::vector< int >  m_bases_index ;
  std::vector< int >  m_bases_map ; // [ m_variable_count X m_bases_count ]

  int  m_variable_count ;
  int  m_maximum_degree ;
  int  m_bases_count ;

  CreateSparseProductTensor( const input_type & input )
    :  m_poly()
    ,  m_tensor_map()
    ,  m_variable_degree( input.size() )
    ,  m_bases_index( input.size() , int(0) )
    ,  m_bases_map()
    ,  m_variable_count( input.size() )
    ,  m_maximum_degree(0)
    ,  m_bases_count(0)
    {
      for ( size_type i = 0 ; i < input.size() ; ++i ) {
        m_variable_degree[i] = input[i];
        m_maximum_degree = std::max( m_maximum_degree , m_variable_degree[i] );
      }

      // Maximum degree of combinatorial bases:
      // const int maximum_basis_degree = m_variable_count * m_maximum_degree ;

      // Truncate to minimum of 2 and maximum individual degree
      const int maximum_basis_degree = std::max( 2 , m_maximum_degree );

      for ( int p = 0 ; p <= maximum_basis_degree ; ++p ) {
        for ( int max_p = 0 ; max_p <= p ; ++max_p ) {
          populate_bases( 0 , p , max_p );
        }
      }

      populate_tensor();
    }

  //--------------------------------------------------------------------------

  void populate_bases( const int ivar ,
                       const int poly_degree ,
                       const int target_poly_degree )
  {
    for ( int p = std::min( m_variable_degree[ivar] , poly_degree ) ;
          0 <= p ; --p ) {
      m_bases_index[ ivar ] = p ;

      if ( ( ivar + 1 ) < m_variable_count ) {
        // More variables to map:
        populate_bases( (ivar + 1) , (poly_degree - p) , target_poly_degree );
      }
      else if ( poly_degree == p ) {
        // All variables are specified and the
        // complete product satisfies the target polynomial degree.

        // No variables exceeding target_poly_degree and at least
        // one variable equal to target_poly_degree

        bool ok_limit = true ;
        bool ok_equal = false ;

        for ( int kvar = 0 ; kvar < m_variable_count && ok_limit ; ++kvar ) {
          ok_limit = m_bases_index[kvar] <= target_poly_degree ;
          if ( target_poly_degree == m_bases_index[kvar] ) ok_equal = true ;
        }

        if ( ok_limit && ok_equal ) {
          ++m_bases_count ;

          // Insert this variable bases combination into map:

          m_bases_map.insert( m_bases_map.end() ,
                              m_bases_index.begin() , m_bases_index.end() );  
        }
      }
    }
  }

  void populate_tensor()
  {
    const double  tolerance( 1e-12 );
    const double  scaling( std::sqrt( (double)( 1 << m_variable_count )));

    // Integration rule required for integration of triple product of 
    // maximum degree polynomial

    const int rule_count = 1 + ( 3 * m_maximum_degree ) / 2 ;

    std::vector< double >  integral_sum( m_variable_count );
    std::vector< double >  poly_values(  m_maximum_degree + 1 );
    std::vector< double >  rule_points(  rule_count );
    std::vector< double >  rule_weights( rule_count );

    gauss_legendre( rule_count , & rule_points[0] , & rule_weights[0] );

    // Triple product this new basis with all previous bases.
    // Only compute unique values from symmetry: I >= J >= K.

    for ( int i = 0 ; i < m_bases_count ; ++i ) {
      const int iMap = i * m_variable_count ;

      for ( int j = 0 ; j <= i ; ++j ) {
        const int jMap = j * m_variable_count ;

        for ( int k = 0 ; k <= j ; ++k ) {
          const int kMap = k * m_variable_count ;

          for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
            integral_sum[iv] = 0 ;
          }

          for ( int r = 0 ; r < rule_count; ++r ) {

            m_poly.evaluate( m_maximum_degree ,
                             rule_points[r] , & poly_values[0] );

            for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
              integral_sum[iv] +=
                 rule_weights[r] *
                 poly_values[ m_bases_map[iMap+iv] ] *
                 poly_values[ m_bases_map[jMap+iv] ] *
                 poly_values[ m_bases_map[kMap+iv] ] ;
            }
          }

          double value = scaling ;
          for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
            value *= integral_sum[iv] ;
          }

          if ( tolerance < value ) {
            // A non-zero entry for the product tensor 
            m_tensor_map[ product_tensor_index_type(i,j,k) ] = value ;
          }
        }
      }
    }
  }

  void print_bases() const
  {
    std::cout.precision(8);

    std::cout << "KokkosArray::Impl::CreateProductTensor(" ;
    for ( int i = 0 ; i < m_variable_count ; ++i ) {
      std::cout << " " << m_variable_degree[i] ;
    }
    std::cout << " )" << std::endl ;

    std::cout << "Bases = " << std::endl ;
    {
      size_t k = 0 ;
      for ( int i = 0 ; i < m_bases_count ; ++i ) {
        std::cout << "  bases["  << i << "] = (" ;
        for ( int j = 0 ; j < m_variable_count ; ++j , ++k ) {
          std::cout << " " << m_bases_map[k] ;
        }
        std::cout << " )" << std::endl ;
      }
    }
  }

  void print_tensor() const
  {
    std::cout.precision(8);

    std::cout << "KokkosArray::Impl::CreateProductTensor(" ;
    std::cout << "Tensor (" << m_tensor_map.size() << ") = " << std::endl ;
    for ( typename tensor_map_type::const_iterator
          n =  m_tensor_map.begin() ;
          n != m_tensor_map.end() ; ++n ) {
      const int i = (*n).first.coord(0);
      const int j = (*n).first.coord(1);
      const int k = (*n).first.coord(2);

      std::cout << "  <" << i << "," << j << "," << k
                << "> = <" ;

      for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
        if ( 0 == iv ) std::cout << "(" ;
        else std::cout << "," ;
        std::cout << m_bases_map[iv+i*m_variable_count] ;
      }
      std::cout << ")," ;
      for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
        if ( 0 == iv ) std::cout << "(" ;
        else std::cout << "," ;
        std::cout << m_bases_map[iv+j*m_variable_count] ;
      }
      std::cout << ")," ;
      for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
        if ( 0 == iv ) std::cout << "(" ;
        else std::cout << "," ;
        std::cout << m_bases_map[iv+k*m_variable_count] ;
      }
      std::cout << ")> = " << (*n).second << std::endl ;
    }
  }

public:

  static
  type create( const input_type & variable_polynomial_degree )
  {
    // Manufactor all required data
    const CreateSparseProductTensor work( variable_polynomial_degree );

    // work.print_bases();
    // work.print_tensor();

    type tmp ;

    // Allocate and transfer data to the device-resident object.

    // SparseProductTensor {
    //   sparse-product-tensor evaluations
    //   variable-count
    //   polynomial-bases
    //   polydegree( bases-index , variable-index )
    //   scaling( bases-index )
    // };

    typedef View< size_type** , device_type > int_array_type ;
    typedef typename int_array_type  ::HostMirror host_int_array_type ;

    tmp.m_degree_map =
      int_array_type( "stochastic_tensor_degree_map" ,
                      work.m_bases_count + 1 ,
                      work.m_variable_count );

    tmp.m_variable  = work.m_variable_count ;

    {
      // If possible the mirror uses a view.
      host_int_array_type degree_map = create_mirror_view( tmp.m_degree_map );

      for ( int j = 0 ; j < work.m_variable_count ; ++j ) {
        degree_map(0,j) = work.m_variable_degree[j];
      }
      for ( int j = 0 ; j < work.m_variable_count ; ++j ) {
        for ( int i = 0 ; i < work.m_bases_count ; ++i ) {
          degree_map(i+1,j) = work.m_bases_map[j + i * work.m_variable_count ];
        }
      }
      deep_copy( tmp.m_degree_map , degree_map );
    }

    tmp.m_tensor = CreateSparseProductTensor< tensor_type , tensor_map_type >
                     ::create( work.m_tensor_map );

    return tmp ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_STOCHASTICPRODUCTTENSORCREATE_HPP */


