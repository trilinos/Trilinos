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

#ifndef STOKHOS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP
#define STOKHOS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP

#include "KokkosArray_View.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"

#include <sstream>
#include <fstream>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 */
template< typename ValueType , class DeviceType >
class LexicographicBlockSparse3Tensor {
public:

  typedef DeviceType                       device_type;
  typedef typename device_type::size_type  size_type;
  typedef ValueType                        value_type;

private:

  typedef KokkosArray::View< size_type[][6] , device_type > coord_array_type;
  typedef KokkosArray::View< value_type[], device_type >    value_array_type;

  coord_array_type  m_coord;
  value_array_type  m_value;
  size_type         m_dimension;
  size_type         m_flops;

public:

  inline
  ~LexicographicBlockSparse3Tensor() {}

  inline
  LexicographicBlockSparse3Tensor() :
    m_coord(),
    m_value(),
    m_dimension(),
    m_flops(0){}

  inline
  LexicographicBlockSparse3Tensor(
    const LexicographicBlockSparse3Tensor & rhs) :
    m_coord(rhs.m_coord),
    m_value(rhs.m_value),
    m_dimension(rhs.m_dimension),
    m_flops(rhs.m_flops) {}

  inline
  LexicographicBlockSparse3Tensor &operator=(
    const LexicographicBlockSparse3Tensor & rhs)
  {
    m_coord = rhs.m_coord;
    m_value = rhs.m_value;
    m_dimension = rhs.m_dimension;
    m_flops = rhs.m_flops;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_dimension; }

  /** \brief  Number of coordinates. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_coord() const { return m_coord.dimension_0(); }

  /** \brief  Number of values. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_value() const { return m_value.dimension_0(); }

  /** \brief   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type get_i_begin(const size_type entry) const {
    return m_coord(entry,0);
  }

  /** \brief   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type get_j_begin(const size_type entry) const {
    return m_coord(entry,1);
  }

  /** \brief   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type get_k_begin(const size_type entry) const {
    return m_coord(entry,2);
  }

  /** \brief   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type get_i_size(const size_type entry) const {
    return m_coord(entry,3);
  }

  /** \brief   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type get_j_size(const size_type entry) const {
    return m_coord(entry,4);
  }

  /** \brief   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type get_k_size(const size_type entry) const {
    return m_coord(entry,5);
  }

  /** \brief Block information for entry 'entry' */
   KOKKOSARRAY_INLINE_FUNCTION
   void block(const size_type entry,
              size_type& i_begin, size_type& j_begin, size_type& k_begin,
              size_type& i_size,  size_type& j_size,  size_type& k_size) const {
     i_begin = m_coord(entry,0);
     j_begin = m_coord(entry,1);
     k_begin = m_coord(entry,2);
     i_size  = m_coord(entry,3);
     j_size  = m_coord(entry,4);
     k_size  = m_coord(entry,5);
   }

  /** \brief  Cijk for entry 'entry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type& value(const size_type entry) const
  { return m_value(entry); }

  /** \brief Number of non-zero's */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_non_zeros() const { return m_value.dimension_0(); }

  /** \brief Number flop's per multiply-add */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_flops() const { return m_flops; }

  template <typename OrdinalType>
  static LexicographicBlockSparse3Tensor
  create(const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
         const Stokhos::LTBSparse3Tensor<OrdinalType,ValueType>& Cijk)
  {
    using Teuchos::Array;
    using Teuchos::RCP;

    // Allocate tensor data
    LexicographicBlockSparse3Tensor tensor ;
    tensor.m_dimension = basis.size();
    tensor.m_coord = coord_array_type( "coord" , Cijk.num_leafs() );
    tensor.m_value = value_array_type( "value" , Cijk.num_entries() );

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror host_coord =
      KokkosArray::create_mirror_view( tensor.m_coord );
    typename value_array_type::HostMirror host_value =
      KokkosArray::create_mirror_view( tensor.m_value );

    // Fill flat 3 tensor
    typedef Stokhos::LTBSparse3Tensor<OrdinalType,ValueType> Cijk_type;
    typedef typename Cijk_type::CijkNode node_type;
    Array< RCP<const node_type> > node_stack;
    Array< OrdinalType > index_stack;
    node_stack.push_back(Cijk.getHeadNode());
    index_stack.push_back(0);
    RCP<const node_type> node;
    OrdinalType child_index;
    OrdinalType coord_index = 0;
    OrdinalType value_index = 0;
    tensor.m_flops = 0;
    while (node_stack.size() > 0) {
      node = node_stack.back();
      child_index = index_stack.back();

      // Leaf
      if (node->is_leaf) {
        host_coord(coord_index, 0) = node->i_begin;
        host_coord(coord_index, 1) = node->j_begin;
        host_coord(coord_index, 2) = node->k_begin;
        host_coord(coord_index, 3) = node->i_size;
        host_coord(coord_index, 4) = node->j_size;
        host_coord(coord_index, 5) = node->k_size;
        ++coord_index;
        for (OrdinalType i=0; i<node->my_num_entries; ++i)
          host_value(value_index++) = node->values[i];
        tensor.m_flops += 3*node->my_num_entries + node->i_size;
        node_stack.pop_back();
        index_stack.pop_back();
      }

      // More children to process -- process them first
      else if (child_index < node->children.size()) {
        ++index_stack.back();
        node = node->children[child_index];
        node_stack.push_back(node);
        index_stack.push_back(0);
      }

      // No more children
      else {
        node_stack.pop_back();
        index_stack.pop_back();
      }

    }
    TEUCHOS_ASSERT(coord_index == Cijk.num_leafs());
    TEUCHOS_ASSERT(value_index == Cijk.num_entries());

    /*
    // Save block volumes to file
    static int index = 0;
    std::stringstream file_name;
    file_name << "cijk_vol_" << index << ".txt";
    std::ofstream file(file_name.str().c_str());
    for (size_type i=0; i<coord_index; ++i) {
      size_type vol = host_coord(i,3) * host_coord(i,4) * host_coord(i,5);
      file << vol << std::endl;
    }
    file.close();
    index++;
    */

    // Copy data to device if necessary
    KokkosArray::deep_copy( tensor.m_coord , host_coord );
    KokkosArray::deep_copy( tensor.m_value , host_value );

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType >
LexicographicBlockSparse3Tensor<ValueType, Device>
create_lexicographic_block_sparse_3_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::LTBSparse3Tensor<OrdinalType,ValueType>& Cijk )
{
  return LexicographicBlockSparse3Tensor<ValueType, Device>::create(Cijk );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP */
