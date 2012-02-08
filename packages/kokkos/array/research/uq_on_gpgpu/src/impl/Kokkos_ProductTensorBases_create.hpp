/*
//@HEADER
// ************************************************************************
// 
//                         Kokkos Array
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

#ifndef KOKKOS_LEGENDREPOLYNOMIALSCREATE_HPP
#define KOKKOS_LEGENDREPOLYNOMIALSCREATE_HPP

#include <iostream>
#include <cmath>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType , class PolynomialType , class Device , typename IntType >
class CreateProductTensorFromBases<
  ProductTensorFromBases< ValueType , PolynomialType , Device > ,
  std::vector< IntType > >
{
private:

  typedef Device device_type ;
  typedef typename device_type::size_type size_type ;
  typedef ProductTensorFromBases< ValueType ,
                                  PolynomialType ,
                                  device_type > type ;

  typedef std::vector< IntType > input_type ;

  typedef ProductTensorIndex< 3 , device_type > product_tensor_index_type ;
  typedef std::map< product_tensor_index_type , ValueType >  tensor_map_type ;

  enum { is_host_memory =
           SameType< typename device_type::memory_space , Host >::value };

  PolynomialType  m_poly ;

  tensor_map_type         m_tensor_map ;
  std::vector< int >      m_variable_degree ;
  std::vector< int >      m_bases_index ;
  std::vector< int >      m_bases_map ; // [ m_variable_count X m_bases_count ]
  std::vector< double >   m_bases_scaling ;
  std::vector< double >   m_bases_values ;
  std::vector< double >   m_poly_values ;
  std::vector< double >   m_rule_points ;
  std::vector< double >   m_rule_weights ;

  int  m_variable_count ;
  int  m_maximum_degree ;
  int  m_bases_count ;
  int  m_rule_count ;
  const double m_tolerance ;

  CreateProductTensorFromBases( const input_type & input )
    :  m_poly()
    ,  m_tensor_map()
    ,  m_variable_degree( input.size() )
    ,  m_bases_index( input.size() , int(0) )
    ,  m_bases_map()
    ,  m_bases_scaling()
    ,  m_bases_values()
    ,  m_poly_values()
    ,  m_rule_points()
    ,  m_rule_weights()
    ,  m_variable_count( input.size() )
    ,  m_maximum_degree(0)
    ,  m_bases_count(0)
    ,  m_rule_count(0)
    , m_tolerance( 1e-14 )
    {
      for ( size_type i = 0 ; i < input.size() ; ++i ) {
        m_variable_degree[i] = input[i];
        m_maximum_degree = std::max( m_maximum_degree , m_variable_degree[i] );
      }

      // Maximum degree of combinatorial bases
      // const int maximum_basis_degree = m_variable_count * m_maximum_degree ;
      const int maximum_basis_degree = m_maximum_degree ;

      // Integration rule required for maximum degree triple product
      m_rule_count = 1 + ( 3 * maximum_basis_degree ) / 2 ;

      m_rule_points.resize(  m_rule_count );
      m_rule_weights.resize( m_rule_count );

      gauss_legendre( m_rule_count , & m_rule_points[0] , & m_rule_weights[0] );

      for ( int p = 0 ; p <= maximum_basis_degree ; ++p ) {
        for ( int max_p = 0 ; max_p <= p ; ++max_p ) {
          populate_bases( 0 , p , max_p );
        }
      }

      populate_tensor();
    }

  //--------------------------------------------------------------------------

  void evaluate_polynomials( double x )
  {
     m_poly_values.resize( m_maximum_degree + 1 );

     m_poly.evaluate( m_maximum_degree , x , & m_poly_values[0] );
  }

  void evaluate_bases( double x )
  {
    m_bases_values.resize( m_bases_count );

    evaluate_polynomials( x );

    int iIndex = 0 ;
    for ( int iBasis = 0 ; iBasis < m_bases_count ; ++iBasis ) {
      double value = m_bases_scaling[ iBasis ];

      // Variable's bases degree stored in:
      // m_bases_map[ iIndex .. iIndex + m_variable_count - 1 ]

      for ( int iVar = 0 ; iVar < m_variable_count ; ++iVar , ++iIndex ) {
        value *= m_poly_values[ m_bases_map[ iIndex ] ] ;
      }

      m_bases_values[iBasis] = value ;
    }
  }

  // A candidate basis in 'm_bases_index'.
  // Check for non-zero norm and save the scale factor for normalization.
  void candidate_basis()
  {
    // Determine if a non-zero inner product

    double sum = 0 ;

    for ( int r = 0 ; r < m_rule_count ; ++r ) {
      evaluate_polynomials( m_rule_points[r] );
      double tmp = 1 ;
      for ( int iv = 0 ; iv < m_variable_count ; ++iv ) {
        tmp *= m_poly_values[ m_bases_index[iv] ];
      }
      sum += m_rule_weights[r] * tmp * tmp ;
    }

    if ( m_tolerance < sum ) {
      // Non-zero inner product

      ++m_bases_count ;

      // Save the scaling for normalization

      sum = 1.0 / std::sqrt( sum );
      m_bases_scaling.push_back( sum );

      // Insert this variable bases combination into map:

      m_bases_map.insert( m_bases_map.end() ,
                          m_bases_index.begin() , m_bases_index.end() );  
    }
  }

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
          candidate_basis();
        }
      }
    }
  }

  void populate_tensor()
  {

    // Triple product this new basis with all previous bases.
    // Only compute unique values from symmetry: I >= J >= K.

    for ( int i = 0 ; i < m_bases_count ; ++i ) {
      for ( int j = 0 ; j <= i ; ++j ) {
        for ( int k = 0 ; k <= j ; ++k ) {
          double sum = 0 ;

          for ( int r = 0 ; r < m_rule_count; ++r ) {

            evaluate_bases( m_rule_points[r] );

            sum += m_rule_weights[r] * m_bases_values[i] *
                                       m_bases_values[j] *
                                       m_bases_values[k] ;
          }

          if ( m_tolerance < sum ) {
            // A non-zero entry for the product tensor 
            m_tensor_map[ product_tensor_index_type(i,j,k) ] = sum ;
          }
        }
      }
    }
  }

  void print() const
  {
    std::cout.precision(8);

    std::cout << "Kokkos::Impl::CreateProductTensorFromBases(" ;
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
        std::cout << " ) = " << m_bases_scaling[i] << std::endl ;
      }
    }


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
    const CreateProductTensorFromBases work( variable_polynomial_degree );

    work.print();

    type tmp ;

    // Allocate and transfer data to the device-resident object.

    // SparseProductTensor {
    //   sparse-product-tensor evaluations
    //   variable-count
    //   polynomial-bases
    //   polydegree( bases-index , variable-index )
    //   scaling( bases-index )
    // };

    typedef MDArray< size_type , device_type > int_array_type ;
    typedef Impl::MemoryView< ValueType , device_type > value_array_type ;

    typedef typename int_array_type  ::HostMirror host_int_array_type ;
    typedef typename value_array_type::HostMirror host_value_array_type ;
    typedef typename value_array_type::HostMirror host_value_array_type ;

    tmp.m_degree_map =
      create_mdarray< int_array_type >( work.m_bases_count + 1 ,
                                        work.m_variable_count );

    tmp.m_coord =
      create_mdarray< int_array_type >( work.m_tensor_map.size() , 3 );

    tmp.m_value.allocate( work.m_tensor_map.size() , std::string() );

    tmp.m_variable  = work.m_variable_count ;
    tmp.m_dimension = work.m_bases_count ;
    tmp.m_tensor_count = work.m_tensor_map.size() ;

    {
      host_int_array_type degree_map =
        Impl::CreateMirror< int_array_type , is_host_memory >
            ::create( tmp.m_degree_map );

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

    {
      host_int_array_type coord =
        Impl::CreateMirror< int_array_type , is_host_memory >
            ::create( tmp.m_coord );

      host_value_array_type value =
        Impl::CreateMirror< value_array_type , is_host_memory >
            ::create( tmp.m_value , work.m_tensor_map.size() );

      int j = 0 ;
      for ( typename tensor_map_type::const_iterator
            n = work.m_tensor_map.begin();
            n != work.m_tensor_map.begin(); ++n , ++j ) {
        coord(j,0) = (*n).first.coord(0);
        coord(j,1) = (*n).first.coord(1);
        coord(j,2) = (*n).first.coord(2);
        value[j]   = (*n).second ;
      }
      deep_copy( tmp.m_coord , coord );
      DeepCopy< value_array_type , host_value_array_type >
        ::run( tmp.m_value , value , work.m_tensor_map.size() );
    }

    return tmp ;
  }
};


} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_LEGENDREPOLYNOMIALSCREATE_HPP */


