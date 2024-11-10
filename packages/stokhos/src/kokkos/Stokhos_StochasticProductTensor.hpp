// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STOCHASTICPRODUCTTENSOR_HPP
#define STOKHOS_STOCHASTICPRODUCTTENSOR_HPP

#include <ostream>

#include "Kokkos_Core.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Teuchos_ParameterList.hpp"

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

  typedef Device                          execution_space ;
  typedef ValueType                       value_type ;
  typedef TensorType                      tensor_type ;
  typedef typename tensor_type::size_type size_type ;

private:

  tensor_type                                m_tensor ;
  Kokkos::View< size_type** , execution_space >  m_degree_map ;
  size_type                                  m_variable ;

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

  KOKKOS_INLINE_FUNCTION
  const tensor_type & tensor() const { return m_tensor ; }

  /** \brief  Dimension: number of bases and
   *          length of the vector block (and tensor).
   */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_tensor.dimension(); }

   /** \brief  Aligned dimension: length of the vector block properly aligned.
   */
  KOKKOS_INLINE_FUNCTION
  size_type aligned_dimension() const {
    const bool is_cuda =
#if defined( KOKKOS_ENABLE_CUDA )
      std::is_same<execution_space,Kokkos::Cuda>::value;
#else
      false ;
#endif
    const size_type AlignBytes = is_cuda ? 128 : 64;
    const size_type NumAlign = AlignBytes/sizeof(value_type);
    return (dimension() + NumAlign-1) & ~(NumAlign-1);
  }

  /** \brief  How many variables are being expanded. */
  KOKKOS_INLINE_FUNCTION
  size_type variable_count() const { return m_variable ; }

  /** \brief  Polynomial degree of a given variable */
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_type variable_degree( const iType & iVariable ) const
    { return m_degree_map( 0 , iVariable ); }

  /** \brief  Basis function 'iBasis' is the product of
   *          'variable_count()' polynomials.  Return the
   *          polynomial degree of component 'iVariable'.
   */
  template< typename iType , typename jType >
  KOKKOS_INLINE_FUNCTION
  size_type bases_degree( const iType & iBasis , const jType & iVariable ) const
    { return m_degree_map( iBasis + 1 , iVariable ); }

  void print( std::ostream & s ) const
  {
    for ( unsigned i = 1 ; i < m_degree_map.extent(0) ; ++i ) {
      s << "  bases[" << i - 1 << "] (" ;
      for ( unsigned j = 0 ; j < m_degree_map.extent(1) ; ++j ) {
        s << " " << m_degree_map(i,j);
      }
      s << " )" << std::endl ;
    }
  }

  template <typename OrdinalType, typename CijkType>
  static StochasticProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const CijkType& Cijk,
          const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    StochasticProductTensor spt ;

    // Allocate and transfer data to the device-resident object.

    typedef Kokkos::View< size_type** , execution_space > int_array_type ;
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
      Kokkos::create_mirror_view( spt.m_degree_map );
    for ( OrdinalType j = 0 ; j < basis_dim ; ++j )
      degree_map(0,j) = max_orders[j];
    for ( OrdinalType i = 0 ; i < basis_sz ; ++i ) {
      const Stokhos::MultiIndex<OrdinalType>& term = basis.term(i);
      for ( OrdinalType j = 0 ; j < basis_dim ; ++j ) {
        degree_map(i+1,j) = term[j];
      }
    }
    Kokkos::deep_copy( spt.m_degree_map , degree_map );

    // Build 3 tensor
    spt.m_tensor = tensor_type::create( basis, Cijk, params );

    return spt ;
  }
};

template<  typename TensorType, typename OrdinalType , typename ValueType, typename CijkType >
StochasticProductTensor<ValueType, TensorType, typename TensorType::execution_space>
create_stochastic_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const CijkType& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList())
{
  typedef typename TensorType::execution_space Device;
  return StochasticProductTensor<ValueType, TensorType, Device>::create(
    basis, Cijk, params);
}

template < typename ValueType , typename Device, class TensorType >
class BlockMultiply< StochasticProductTensor< ValueType, TensorType, Device > >
{
public:
  typedef Device execution_space ;
  typedef typename execution_space::size_type size_type ;
  typedef StochasticProductTensor< ValueType, TensorType, execution_space > block_type ;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const block_type  & block ,
                     const MatrixValue *       a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    typedef BlockMultiply< typename block_type::tensor_type > tensor_multiply ;

    tensor_multiply::apply( block.tensor() , a , x , y );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



} // namespace Stokhos

#endif /* #ifndef STOKHOS_STOCHASTICPRODUCTTENSOR_HPP */
