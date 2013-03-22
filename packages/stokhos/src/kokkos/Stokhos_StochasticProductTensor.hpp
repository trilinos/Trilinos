// @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_STOCHASTICPRODUCTTENSOR_HPP
#define STOKHOS_STOCHASTICPRODUCTTENSOR_HPP

#include <ostream>

#include "KokkosArray_View.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace Stokhos {

//----------------------------------------------------------------------------

/** \brief  Bases defined by combinatorial product of polynomial bases.
 *
 *  Bases: \prod_{j=0}^{N-1} P_k(x) \forall j and k \in M(j)
 *  Where:  P_k      is a polynomial of degree k
 *  Where: <P_a,P_b> is the the integral on [-1,1]
 *  Where: <P_a,P_b> is the Kronecker delta \delta_{a,b}
 *                   thus the polynomials are normalized
 *                   with respect to this inner product.
 *
 *  Where: N = the number of variables expanded via polynomial bases
 *  Where: M(j) = the degree of a particular variable
 *
 *  Where: \psi_I(x) = is one basis function and I is a multi-index
 *                     of rank N, denoting one function from each 
 *                     variable's polynomial bases.
 *
 *  Were: <\psi_I,\psi_J,\psi_K> is the integral on [-1,1]
 *
 *  The bases space is sparse due to orthogonality within the
 *  expansion.
 */
template< typename ValueType , typename TensorType, class Device >
class StochasticProductTensor {
public:

  typedef Device                                      device_type ;
  typedef typename device_type::size_type             size_type ;
  typedef ValueType                                   value_type ;
  typedef TensorType                                  tensor_type ;

private:

  tensor_type                                     m_tensor ;
  KokkosArray::View< size_type** , device_type >  m_degree_map ;
  size_type                                       m_variable ;

public:

  inline
  ~StochasticProductTensor() {}

  inline
  StochasticProductTensor()
    : m_tensor()
    , m_degree_map()
    , m_variable(0)
    {}

  inline
  StochasticProductTensor( const StochasticProductTensor & rhs )
    : m_tensor(       rhs.m_tensor )
    , m_degree_map(   rhs.m_degree_map )
    , m_variable(     rhs.m_variable )
    {}

  inline
  StochasticProductTensor & operator = ( const StochasticProductTensor & rhs )
  {
    m_tensor       = rhs.m_tensor ;
    m_degree_map   = rhs.m_degree_map ;
    m_variable     = rhs.m_variable ;
    return *this ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  const tensor_type & tensor() const { return m_tensor ; }

  /** \brief  Dimension: number of bases and
   *          length of the vector block (and tensor).
   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_tensor.dimension(); }

  /** \brief  How many variables are being expanded. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type variable_count() const { return m_variable ; }

  /** \brief  Polynomial degree of a given variable */
  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_type variable_degree( const iType & iVariable ) const
    { return m_degree_map( 0 , iVariable ); }

  /** \brief  Basis function 'iBasis' is the product of
   *          'variable_count()' polynomials.  Return the
   *          polynomial degree of component 'iVariable'.
   */
  template< typename iType , typename jType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_type bases_degree( const iType & iBasis , const jType & iVariable ) const
    { return m_degree_map( iBasis + 1 , iVariable ); }

  void print( std::ostream & s ) const
  {
    for ( unsigned i = 1 ; i < m_degree_map.dimension_0() ; ++i ) {
      s << "  bases[" << i - 1 << "] (" ;
      for ( unsigned j = 0 ; j < m_degree_map.dimension_1() ; ++j ) {
        s << " " << m_degree_map(i,j);
      }
      s << " )" << std::endl ;
    }
  }

  template <typename OrdinalType>
  static StochasticProductTensor 
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
	  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk)
  {
    StochasticProductTensor spt ;

    // Allocate and transfer data to the device-resident object.

    typedef KokkosArray::View< size_type** , device_type > int_array_type ;
    typedef typename int_array_type::HostMirror host_int_array_type ;

    OrdinalType basis_sz = basis.size();
    OrdinalType basis_dim = basis.dimension();
    Stokhos::MultiIndex<OrdinalType> max_orders = basis.getMaxOrders();

    spt.m_degree_map =
      int_array_type( "stochastic_tensor_degree_map" ,
                      basis_sz + 1 ,
		      basis_dim );

    spt.m_variable  = basis_dim ;

    // Build degree_map
    host_int_array_type degree_map = 
      KokkosArray::create_mirror_view( spt.m_degree_map );
    for ( OrdinalType j = 0 ; j < basis_dim ; ++j )
      degree_map(0,j) = max_orders[j];
    for ( OrdinalType i = 0 ; i < basis_sz ; ++i ) {
      const Stokhos::MultiIndex<OrdinalType>& term = basis.term(i);
      for ( OrdinalType j = 0 ; j < basis_dim ; ++j ) {
	degree_map(i+1,j) = term[j];
      }
    }
    KokkosArray::deep_copy( spt.m_degree_map , degree_map );

    // Build 3 tensor
    spt.m_tensor = tensor_type::create( basis, Cijk );

    return spt ;
  }
};

template<  typename TensorType, typename OrdinalType , typename ValueType >
StochasticProductTensor<ValueType, TensorType, typename TensorType::device_type>
create_stochastic_product_tensor( 
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk)
{
  typedef typename TensorType::device_type Device;
  return StochasticProductTensor<ValueType, TensorType, Device>::create( basis,
									 Cijk);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



} // namespace Stokhos

#endif /* #ifndef STOKHOS_STOCHASTICPRODUCTTENSOR_HPP */


