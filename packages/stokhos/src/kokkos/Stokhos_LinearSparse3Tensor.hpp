// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_LINEAR_SPARSE_3_TENSOR_HPP
#define STOKHOS_LINEAR_SPARSE_3_TENSOR_HPP

#include "KokkosArray_View.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 */
template< typename ValueType , class DeviceType , int BlockSize = 1 >
class LinearSparse3Tensor {
public:

  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

  static const int block_size = BlockSize;

private:

  typedef KokkosArray::View< value_type[], device_type > value_array_type ;

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
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_dim ; }

  /** \brief  Dimension of the tensor. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type aligned_dimension() const { return m_aligned_dim ; }

  /** \brief  Number of sparse entries. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_value.dimension_0(); }

  /** \brief Is tensor built from symmetric PDFs. */
   KOKKOSARRAY_INLINE_FUNCTION
   bool symmetric() const 
   { return m_symmetric; }

  /** \brief  Value for entry 'entry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  /** \brief Number of non-zero's */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_non_zeros() const
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_flops() const
  { return m_flops; }

  template <typename OrdinalType>
  static LinearSparse3Tensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis ,
          const bool symmetric )
  {
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
      host_value = KokkosArray::create_mirror_view( tensor.m_value );

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
    KokkosArray::deep_copy( tensor.m_value , host_value );

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
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
{
  return LinearSparse3Tensor<ValueType, Device, BlockSize>::create( basis, Cijk );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_LINEAR_SPARSE_3_TENSOR_HPP */
