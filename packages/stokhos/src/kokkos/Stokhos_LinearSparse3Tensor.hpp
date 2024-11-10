// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_LINEAR_SPARSE_3_TENSOR_HPP
#define STOKHOS_LINEAR_SPARSE_3_TENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_TinyVec.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 */
template< typename ValueType , class ExecutionSpace , int BlockSize >
class LinearSparse3Tensor {
public:

  typedef ExecutionSpace                       execution_space ;
  typedef typename execution_space::size_type  size_type ;
  typedef ValueType                        value_type ;

  static const int block_size = BlockSize;

private:

  typedef Kokkos::View< value_type[], execution_space > value_array_type ;

  value_array_type   m_value ;
  size_type          m_dim ;
  size_type          m_aligned_dim ;
  size_type          m_nnz ;
  size_type          m_flops ;
  bool               m_symmetric ;

public:

  inline
  ~LinearSparse3Tensor() {}

  inline
  LinearSparse3Tensor() :
    m_value() ,
    m_dim() ,
    m_aligned_dim(),
    m_nnz(0) ,
    m_flops(0) ,
    m_symmetric(false) {}

  inline
  LinearSparse3Tensor( const LinearSparse3Tensor & rhs ) :
    m_value( rhs.m_value ) ,
    m_dim( rhs.m_dim ),
    m_aligned_dim( rhs.m_aligned_dim ),
    m_nnz( rhs.m_nnz ) ,
    m_flops( rhs.m_flops ) ,
    m_symmetric( rhs.m_symmetric ) {}

  inline
  LinearSparse3Tensor & operator = ( const LinearSparse3Tensor & rhs )
  {
    m_value = rhs.m_value ;
    m_dim = rhs.m_dim ;
    m_aligned_dim = rhs.m_aligned_dim;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    m_symmetric = rhs.m_symmetric;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_dim ; }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type aligned_dimension() const { return m_aligned_dim ; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const
  { return m_value.extent(0); }

  /** \brief Is tensor built from symmetric PDFs. */
   KOKKOS_INLINE_FUNCTION
   bool symmetric() const
   { return m_symmetric; }

  /** \brief  Value for entry 'entry' */
  KOKKOS_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const
  { return m_flops; }

  template <typename OrdinalType>
  static LinearSparse3Tensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params)
  {
    const bool symmetric = params.get<bool>("Symmetric");

    // Allocate tensor data -- currently assuming isotropic
    const size_type dim = basis.size();
    LinearSparse3Tensor tensor ;
    tensor.m_dim = dim;
    tensor.m_aligned_dim = dim;
    if (tensor.m_aligned_dim % block_size)
      tensor.m_aligned_dim += block_size - tensor.m_aligned_dim % block_size;
    tensor.m_symmetric = symmetric;
    tensor.m_nnz = symmetric ? 2 : 3 ;
    tensor.m_value = value_array_type( "value" , tensor.m_nnz );

    // Create mirror, is a view if is host memory
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );

    // Get Cijk values
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases = basis.getCoordinateBases();
    Teuchos::RCP< Stokhos::Dense3Tensor<OrdinalType,ValueType> > cijk =
      bases[0]->computeTripleProductTensor();
    // For non-isotropic, need to take products of these over basis components
    host_value(0) = (*cijk)(0,0,0);
    host_value(1) = (*cijk)(0,1,1);
    if (!symmetric)
      host_value(2) = (*cijk)(1,1,1);

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_value , host_value );

    tensor.m_flops = 8*dim;
    if (!symmetric)
      tensor.m_flops += 2*dim ;

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType , int BlockSize >
LinearSparse3Tensor<ValueType, Device,BlockSize>
create_linear_sparse_3_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params)
{
  return LinearSparse3Tensor<ValueType, Device, BlockSize>::create(
    basis, Cijk, params );
}

template < typename ValueType, typename Device, int BlockSize >
class BlockMultiply< LinearSparse3Tensor< ValueType , Device , BlockSize > >
{
public:

  typedef typename Device::size_type size_type ;
  typedef LinearSparse3Tensor< ValueType , Device , BlockSize > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type block_size = tensor_type::block_size;
    typedef TinyVec<ValueType,block_size,true> TV;
    const size_type dim = tensor.dimension();

    const ValueType c0 = tensor.value(0);
    const ValueType c1 = tensor.value(1);
    const ValueType a0 = a[0];
    const ValueType x0 = x[0];

    if (block_size > 1) {

      TV vc0(c0), vc1(c1), va0(a0), vx0(x0), vy0;
      TV ai, ai2, xi, yi;

      const MatrixValue *aa = a;
      const VectorValue *xx = x;
      VectorValue *yy = y;
      vy0.zero();

      const size_type nBlock = dim / block_size;
      const size_type iEnd = nBlock * block_size;

      if (tensor.symmetric()) {

        size_type i=0;
        for ( ; i < iEnd; i+=block_size,aa+=block_size,xx+=block_size,yy+=block_size) {
          ai.aligned_load(aa);
          ai2 = ai;
          xi.aligned_load(xx);
          yi.aligned_load(yy);

          // y[i] += c1*(a0*xi + ai*x0);
          ai.times_equal(vx0);
          ai2.times_equal(xi);
          xi.times_equal(va0);
          xi.plus_equal(ai);
          xi.times_equal(vc1);
          yi.plus_equal(xi);
          yi.aligned_scatter(yy);

          // y0  += c1*ai*xi;
          ai2.times_equal(vc1);
          vy0.plus_equal(ai2);
        }
        ValueType y0 = vy0.sum();

        // Do remaining entries with a scalar loop
        for ( ; i < dim; ++i) {
          const ValueType ai = *aa++;
          const ValueType xi = *xx++;
          *yy++ += c1*(a0*xi + ai*x0);
          y0  += c1*ai*xi;
        }
        y[0] += y0 + (c0-3.0*c1)*a0*x0;
      }
      else {

        const ValueType c2 = tensor.value(2);
        TV vc2(c2);
        size_type i=0;
        for ( ; i < iEnd; i+=block_size,aa+=block_size,xx+=block_size,yy+=block_size) {
          ai.aligned_load(aa);
          ai2 = ai;
          xi.aligned_load(xx);
          yi.aligned_load(yy);

          // y[i] += c1*(a0*xi + ai*x0) + c2*aixi;
          ai.times_equal(vx0);
          ai2.times_equal(xi);
          xi.times_equal(va0);
          xi.plus_equal(ai);
          xi.times_equal(vc1);
          yi.plus_equal(xi);
          ai = ai2;
          ai.times_equal(vc2);
          yi.plus_equal(ai);
          yi.aligned_scatter(yy);

          // y0  += c1*aixi;
          ai2.times_equal(vc1);
          vy0.plus_equal(ai2);
        }
        ValueType y0 = vy0.sum();

        // Do remaining entries with a scalar loop
        for ( ; i < dim; ++i) {
          const ValueType ai = *aa++;
          const ValueType xi = *xx++;
          const ValueType aixi = ai*xi;
          *yy++ += c1*(a0*xi + ai*x0) + c2*aixi;
          y0  += c1*aixi;
        }
        y[0] += y0 + (c0-3.0*c1-c2)*a0*x0;

      }

    }

    else {

      ValueType y0 = c0*a0*x0;

      if (tensor.symmetric()) {

        for ( size_type i = 1; i < dim; ++i) {
          const ValueType ai = a[i];
          const ValueType xi = x[i];
          y[i] += c1*(a0*xi + ai*x0);
          y0  += c1*ai*xi;
        }
        y[0] += y0;

      }
      else {

        const ValueType c2 = tensor.value(2);
        for ( size_type i = 1; i < dim; ++i) {
          const ValueType ai = a[i];
          const ValueType xi = x[i];
          const ValueType aixi = ai*xi;
          y[i] += c1*(a0*xi + ai*x0) + c2*aixi;
          y0  += c1*aixi;
        }
        y[0] += y0;

      }

    }

  }

  KOKKOS_INLINE_FUNCTION
  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  KOKKOS_INLINE_FUNCTION
  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_LINEAR_SPARSE_3_TENSOR_HPP */
