// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP
#define STOKHOS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"
#include "Teuchos_ParameterList.hpp"

#include <sstream>
#include <fstream>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 */
template< typename ValueType , class ExecutionSpace >
class LexicographicBlockSparse3Tensor {
public:

  typedef ExecutionSpace                       execution_space;
  typedef typename execution_space::size_type  size_type;
  typedef ValueType                        value_type;

private:

  typedef Kokkos::View< int[][7] , Kokkos::LayoutRight, execution_space >       coord_array_type;
  typedef Kokkos::View< value_type[], execution_space >    value_array_type;

  coord_array_type  m_coord;
  value_array_type  m_value;
  size_type         m_dimension;
  size_type         m_flops;
  bool              m_symmetric;

public:

  inline
  ~LexicographicBlockSparse3Tensor() {}

  inline
  LexicographicBlockSparse3Tensor() :
    m_coord(),
    m_value(),
    m_dimension(),
    m_flops(0),
    m_symmetric(false) {}

  inline
  LexicographicBlockSparse3Tensor(
    const LexicographicBlockSparse3Tensor & rhs) :
    m_coord(rhs.m_coord),
    m_value(rhs.m_value),
    m_dimension(rhs.m_dimension),
    m_flops(rhs.m_flops),
    m_symmetric(rhs.m_symmetric) {}

  inline
  LexicographicBlockSparse3Tensor &operator=(
    const LexicographicBlockSparse3Tensor & rhs)
  {
    m_coord = rhs.m_coord;
    m_value = rhs.m_value;
    m_dimension = rhs.m_dimension;
    m_flops = rhs.m_flops;
    m_symmetric = rhs.m_symmetric;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_dimension; }

  /** \brief  Number of coordinates. */
  KOKKOS_INLINE_FUNCTION
  size_type num_coord() const { return m_coord.extent(0); }

  /** \brief  Number of values. */
  KOKKOS_INLINE_FUNCTION
  size_type num_value() const { return m_value.extent(0); }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_i_begin(const size_type entry) const {
    return m_coord(entry,0);
  }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_j_begin(const size_type entry) const {
    return m_coord(entry,1);
  }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_k_begin(const size_type entry) const {
    return m_coord(entry,2);
  }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_p_i(const size_type entry) const {
    return m_coord(entry,3);
  }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_p_j(const size_type entry) const {
    return m_coord(entry,4);
  }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_p_k(const size_type entry) const {
    return m_coord(entry,5);
  }

  /** \brief   */
  KOKKOS_INLINE_FUNCTION
  int get_j_eq_k(const size_type entry) const {
    return m_coord(entry,6);
  }

  /** \brief  Cijk for entry 'entry' */
  KOKKOS_INLINE_FUNCTION
  const value_type& value(const size_type entry) const
  { return m_value(entry); }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const { return m_value.extent(0); }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const { return m_flops; }

  /** \brief Is PDF symmetric */
  KOKKOS_INLINE_FUNCTION
  bool symmetric() const { return m_symmetric; }

  template <typename OrdinalType>
  static LexicographicBlockSparse3Tensor
  create(const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
         const Stokhos::LTBSparse3Tensor<OrdinalType,ValueType>& Cijk,
         const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    using Teuchos::Array;
    using Teuchos::RCP;

    // Allocate tensor data
    LexicographicBlockSparse3Tensor tensor ;
    tensor.m_dimension = basis.size();
    tensor.m_symmetric = Cijk.symmetric();
    tensor.m_coord = coord_array_type( "coord" , Cijk.num_leafs() );
    tensor.m_value = value_array_type( "value" , Cijk.num_entries() );

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror host_coord =
      Kokkos::create_mirror_view( tensor.m_coord );
    typename value_array_type::HostMirror host_value =
      Kokkos::create_mirror_view( tensor.m_value );

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
        host_coord(coord_index, 3) = node->p_i;
        host_coord(coord_index, 4) = node->p_j;
        host_coord(coord_index, 5) = node->p_k;
        host_coord(coord_index, 6) = node->parent_j_equals_k;
        ++coord_index;
        for (OrdinalType i=0; i<node->my_num_entries; ++i)
          host_value(value_index++) = node->values[i];
        tensor.m_flops += 5*node->my_num_entries + node->i_size;
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
      int vol = host_coord(i,3) * host_coord(i,4) * host_coord(i,5);
      file << vol << std::endl;
    }
    file.close();
    index++;
    */

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_coord , host_coord );
    Kokkos::deep_copy( tensor.m_value , host_value );

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType >
LexicographicBlockSparse3Tensor<ValueType, Device>
create_lexicographic_block_sparse_3_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::LTBSparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList())
{
  return LexicographicBlockSparse3Tensor<ValueType, Device>::create(
    basis, Cijk, params);
}

  template < typename ValueType, typename Device >
class BlockMultiply< LexicographicBlockSparse3Tensor< ValueType , Device > >
{
public:

  typedef typename Device::size_type size_type ;
  typedef LexicographicBlockSparse3Tensor< ValueType , Device > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nBlock = tensor.num_coord();

    if (tensor.symmetric()) {

      // Loop over coordinate blocks
      size_type value_entry = 0;
      for ( size_type block = 0; block < nBlock; ++block) {
        const int i_begin = tensor.get_i_begin(block);
        const int j_begin = tensor.get_j_begin(block);
        const int k_begin = tensor.get_k_begin(block);
        const int p_i = tensor.get_p_i(block);
        const int p_j = tensor.get_p_j(block);
        const int p_k = tensor.get_p_k(block);
        VectorValue * const y_block = y + i_begin;
        const MatrixValue * const a_j_block = a + j_begin;
        const VectorValue * const x_k_block = x + k_begin;
        const MatrixValue * const a_k_block = a + k_begin;
        const VectorValue * const x_j_block = x + j_begin;

        // for (int i=0; i<=p_i; ++i) {
        //   VectorValue ytmp = 0;
        //   for (int j=0; j<=p_j; ++j) {
        //     int k0 = j_eq_k != 0 ? j : 0;
        //     if (symmetric && (k0 % 2 != (i+j) % 2)) ++k0;
        //     for (int k=k0; k<=p_k; k+=k_inc) {
        //       ytmp += tensor.value(value_entry++) *
        //         ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
        //     }
        //   }
        //   y_block[i] += ytmp ;
        // }

        if (tensor.get_j_eq_k(block) != 0) {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              int k0 = j%2 != (i+j)%2 ? j+1 : j;
              for (int k=k0; k<=p_k; k+=2) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
        else {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              for (int k=(i+j)%2; k<=p_k; k+=2) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
      }

    }

    else {

      // Loop over coordinate blocks
      size_type value_entry = 0;
      for ( size_type block = 0; block < nBlock; ++block) {
        const int i_begin = tensor.get_i_begin(block);
        const int j_begin = tensor.get_j_begin(block);
        const int k_begin = tensor.get_k_begin(block);
        const int p_i = tensor.get_p_i(block);
        const int p_j = tensor.get_p_j(block);
        const int p_k = tensor.get_p_k(block);
        VectorValue * const y_block = y + i_begin;
        const MatrixValue * const a_j_block = a + j_begin;
        const VectorValue * const x_k_block = x + k_begin;
        const MatrixValue * const a_k_block = a + k_begin;
        const VectorValue * const x_j_block = x + j_begin;

        // for (int i=0; i<=p_i; ++i) {
        //   VectorValue ytmp = 0;
        //   for (int j=0; j<=p_j; ++j) {
        //     int k0 = j_eq_k != 0 ? j : 0;
        //     if (symmetric && (k0 % 2 != (i+j) % 2)) ++k0;
        //     for (int k=k0; k<=p_k; k+=k_inc) {
        //       ytmp += tensor.value(value_entry++) *
        //         ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
        //     }
        //   }
        //   y_block[i] += ytmp ;
        // }

        if (tensor.get_j_eq_k(block) != 0) {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              for (int k=j; k<=p_k; ++k) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
        else {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              for (int k=0; k<=p_k; ++k) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
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

#endif /* #ifndef STOKHOS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP */
