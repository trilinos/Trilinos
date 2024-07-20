// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_COO_PRODUCT_TENSOR_HPP
#define STOKHOS_COO_PRODUCT_TENSOR_HPP

#include <type_traits>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_ParameterList.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor using 'COO'-like storage format.
 *
 * Entries are sorted by stochastic 'row', requiring segmented-reduction
 * based on row index.
 *
 * Third template argument determines whether to store (i,j,k) triple in a
 * single value.
 */
template< typename ValueType, class ExecutionSpace, bool PackIndex >
class CooProductTensor {};

/** \brief  Specialization of CooProductTensor for packed (i,j,k)
 */
template< typename ValueType, class ExecutionSpace >
class CooProductTensor<ValueType,ExecutionSpace,true> {
public:

  typedef ExecutionSpace                       execution_space;
  typedef typename execution_space::size_type  size_type;
  typedef ValueType                        value_type;

private:

  // Number of bits available for each index
  static const size_type bits = (sizeof(size_type)*8) / 3;

  // Mask for packing index
  static const size_type mask = (1 << bits)-1;

  typedef Kokkos::View< value_type[], execution_space >  vec_type;
  typedef Kokkos::View< size_type[], execution_space > coord_array_type;
  typedef Kokkos::View< value_type[], execution_space > value_array_type;

  coord_array_type   m_coord;
  value_array_type   m_value;
  size_type          m_dim;
  size_type          m_flops;

  // Pack (i,j,k) into a single integer
  static size_type
  pack( const size_type i, const size_type j, const size_type k ) {
    const size_type ii = i & mask;
    const size_type jj = j & mask;
    const size_type kk = k & mask;
    size_type ijk = ii | (jj << bits) | (kk << 2*bits);
    return ijk;
  }

  KOKKOS_INLINE_FUNCTION
  void unpack( size_type ijk, size_type& i, size_type& j, size_type& k ) const {
    i = ijk & mask; ijk >>= bits;
    j = ijk & mask;
    k = ijk >> bits;
  }

public:

  //! Maximum index storable by packed approach
  static const size_type max_index = 1 << bits;

  inline
  ~CooProductTensor() {}

  inline
  CooProductTensor() :
    m_coord(),
    m_value(),
    m_dim(0),
    m_flops(0) {}

  inline
  CooProductTensor( const CooProductTensor & rhs ) :
    m_coord( rhs.m_coord ),
    m_value( rhs.m_value ),
    m_dim( rhs.m_dim ),
    m_flops( rhs.m_flops ) {}

  inline
  CooProductTensor & operator = ( const CooProductTensor & rhs )
  {
    m_coord = rhs.m_coord;
    m_value = rhs.m_value;
    m_dim = rhs.m_dim;
    m_flops = rhs.m_flops;
    return *this;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_dim; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const { return m_coord.extent(0); }

  /** \brief  Get (i,j,k) coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  void coord( const size_type entry,
              size_type& i, size_type& j, size_type& k ) const {
    unpack(m_coord(entry), i, j, k);
  }

  /** \brief  Value of an entry */
  KOKKOS_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value(entry); }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const { return m_coord.extent(0); }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const { return m_flops; }

  template <typename OrdinalType>
  static CooProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    typedef typename Cijk_type::i_iterator i_iterator;
    typedef typename Cijk_type::ik_iterator k_iterator;
    typedef typename Cijk_type::ikj_iterator j_iterator;

    // Compute entry count
    size_type entry_count = 0;
    for (i_iterator i_it = Cijk.i_begin(); i_it!=Cijk.i_end(); ++i_it) {
      for (k_iterator k_it = Cijk.k_begin(i_it); k_it != Cijk.k_end(i_it);
           ++k_it) {
        OrdinalType k = index(k_it);
        for (j_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it);
             ++j_it) {
          OrdinalType j = index(j_it);
          if (j >= k) {
            ++entry_count;
          }
        }
      }
    }

    // Align entry_count
#if defined( KOKKOS_ENABLE_CUDA )
    enum { Align = std::is_same<ExecutionSpace,Kokkos::Cuda>::value ? 32 : 1 };
#else
    enum { Align = 1 };
#endif

    entry_count = (entry_count+Align-1) & ~(Align-1);
    TEUCHOS_ASSERT(entry_count % Align == 0);

    // Allocate tensor data
    CooProductTensor tensor;
    tensor.m_coord = coord_array_type( "tensor_coord", entry_count );
    tensor.m_value = value_array_type( "tensor_value", entry_count );
    tensor.m_dim = basis.size();
    tensor.m_flops = 5*entry_count + tensor.m_dim;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( tensor.m_coord );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );

    size_type n = 0;
    OrdinalType i=0, j=0, k=0;
    for (i_iterator i_it = Cijk.i_begin(); i_it!=Cijk.i_end(); ++i_it) {
      i = index(i_it);
      for (k_iterator k_it = Cijk.k_begin(i_it); k_it != Cijk.k_end(i_it);
           ++k_it) {
        k = index(k_it);
        for (j_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it);
             ++j_it) {
          j = index(j_it);
          ValueType c = Stokhos::value(j_it);
          if (j >= k) {
            host_value(n) = (j != k) ? c : 0.5*c;
            host_coord(n) = pack(i,j,k);
            ++n;
          }
        }
      }
    }
    for (; n < entry_count; ++n) {
      host_value(n) = 0.0;
      host_coord(n) = pack(i,j,k);
    }

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_coord, host_coord );
    Kokkos::deep_copy( tensor.m_value, host_value );

    return tensor;
  }

  void print(std::ostream& os) const {
    size_type num_entry = entry_count();
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( m_coord );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( m_value );
    Kokkos::deep_copy( host_coord, m_coord );
    Kokkos::deep_copy( host_value, m_value );

    os << "CooProductTensor:  dim = "
       << dimension() << ", entry_count = "
       << num_entry << std::endl
       << "Entries: i j k : cijk" << std::endl;
    for (size_type l=0; l<num_entry; ++l) {
      size_type i,j,k;
      unpack(host_coord(l), i, j, k);
      ValueType cijk = host_value(l);
      os << "\t " << l << " : " << i << " " << j << " " << k << " = " << cijk
         << std::endl;
      if (l > 0 && l % 32 == 0) os << std::endl;
    }
  }
};

/** \brief  Specialization of CooProductTensor for unpacked (i,j,k)
 */
template< typename ValueType, class ExecutionSpace>
class CooProductTensor<ValueType,ExecutionSpace,false> {
public:

  typedef ExecutionSpace                       execution_space;
  typedef typename execution_space::size_type  size_type;
  typedef ValueType                        value_type;

private:

  typedef Kokkos::View< value_type[], execution_space >  vec_type;
  typedef Kokkos::View< size_type[][3], execution_space > coord_array_type;
  typedef Kokkos::View< value_type[], execution_space > value_array_type;

  coord_array_type   m_coord;
  value_array_type   m_value;
  size_type          m_dim;
  size_type          m_flops;

public:

  inline
  ~CooProductTensor() {}

  inline
  CooProductTensor() :
    m_coord(),
    m_value(),
    m_dim(0),
    m_flops(0) {}

  inline
  CooProductTensor( const CooProductTensor & rhs ) :
    m_coord( rhs.m_coord ),
    m_value( rhs.m_value ),
    m_dim( rhs.m_dim ),
    m_flops( rhs.m_flops ) {}

  inline
  CooProductTensor & operator = ( const CooProductTensor & rhs )
  {
    m_coord = rhs.m_coord;
    m_value = rhs.m_value;
    m_dim = rhs.m_dim;
    m_flops = rhs.m_flops;
    return *this;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_dim; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const { return m_coord.extent(0); }

  /** \brief  Get (i,j,k) coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  void coord( const size_type entry,
              size_type& i, size_type& j, size_type& k ) const {
    i = m_coord(entry,0);
    j = m_coord(entry,1);
    k = m_coord(entry,2);
  }

  /** \brief  Value of an entry */
  KOKKOS_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const { return m_coord.extent(0); }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const { return m_flops; }

  template <typename OrdinalType>
  static CooProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    typedef typename Cijk_type::i_iterator i_iterator;
    typedef typename Cijk_type::ik_iterator k_iterator;
    typedef typename Cijk_type::ikj_iterator j_iterator;

    // Compute entry count
    size_type entry_count = 0;
    for (i_iterator i_it = Cijk.i_begin(); i_it!=Cijk.i_end(); ++i_it) {
      for (k_iterator k_it = Cijk.k_begin(i_it); k_it != Cijk.k_end(i_it);
           ++k_it) {
        OrdinalType k = index(k_it);
        for (j_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it);
             ++j_it) {
          OrdinalType j = index(j_it);
          if (j >= k) {
            ++entry_count;
          }
        }
      }
    }

    // Align entry_count
#if defined( KOKKOS_ENABLE_CUDA )
    enum { Align = std::is_same<ExecutionSpace,Kokkos::Cuda>::value ? 32 : 1 };
#else
    enum { Align = 1 };
#endif

    entry_count = (entry_count+Align-1) & ~(Align-1);
    TEUCHOS_ASSERT(entry_count % Align == 0);

    // Allocate tensor data
    CooProductTensor tensor;
    tensor.m_coord = coord_array_type( "tensor_coord", entry_count );
    tensor.m_value = value_array_type( "tensor_value", entry_count );
    tensor.m_dim = basis.size();
    tensor.m_flops = 5*entry_count + tensor.m_dim;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( tensor.m_coord );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );

    // Set entries
    size_type n = 0;
    OrdinalType i=0, j=0, k=0;
    for (i_iterator i_it = Cijk.i_begin(); i_it!=Cijk.i_end(); ++i_it) {
      i = index(i_it);
      for (k_iterator k_it = Cijk.k_begin(i_it); k_it != Cijk.k_end(i_it);
           ++k_it) {
        k = index(k_it);
        for (j_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it);
             ++j_it) {
          j = index(j_it);
          ValueType c = Stokhos::value(j_it);
          if (j >= k) {
            host_value(n) = (j != k) ? c : 0.5*c;
            host_coord(n,0) = i;
            host_coord(n,1) = j;
            host_coord(n,2) = k;
            ++n;
          }
        }
      }
    }
    for (; n < entry_count; ++n) {
      host_value(n) = 0.0;
      host_coord(n,0) = i;
      host_coord(n,1) = j;
      host_coord(n,2) = k;
    }

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_coord, host_coord );
    Kokkos::deep_copy( tensor.m_value, host_value );

    return tensor;
  }

  void print(std::ostream& os) const {
    size_type num_entry = entry_count();
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( m_coord );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( m_value );
    Kokkos::deep_copy( host_coord, m_coord );
    Kokkos::deep_copy( host_value, m_value );

    os << "CooProductTensor:  dim = "
       << dimension() << ", entry_count = "
       << num_entry << std::endl
       << "Entries: i j k : cijk" << std::endl;
    for (size_type l=0; l<num_entry; ++l) {
      size_type i = host_coord(l,0);
      size_type j = host_coord(l,1);
      size_type k = host_coord(l,2);
      ValueType cijk = host_value(l);
      os << "\t " << l << " : " << i << " " << j << " " << k << " = " << cijk
         << std::endl;
      if (l > 0 && l % 32 == 0) os << std::endl;
    }
  }
};

template< class Device, bool Pack, typename OrdinalType, typename ValueType >
CooProductTensor<ValueType, Device, Pack>
create_coo_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList())
{
  return CooProductTensor<ValueType, Device, Pack>::create(
    basis, Cijk, params );
}

template <typename ValueType, typename Device, bool Pack>
std::ostream& operator << (
  std::ostream& os, const CooProductTensor<ValueType,Device,Pack>& tensor)
{
  tensor.print(os);
  return os;
}

template< typename ValueType, typename Device, bool Pack >
class BlockMultiply< CooProductTensor< ValueType, Device, Pack > >
{
public:

  typedef typename Device::size_type size_type;
  typedef CooProductTensor< ValueType, Device, Pack > tensor_type;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor,
                     const MatrixValue * const a,
                     const VectorValue * const x,
                           VectorValue * const y )
  {
    const size_type nEntry = tensor.entry_count();
    size_type i = 0, j = 0, k = 0, i_prev = -1;
    VectorValue val = 0.0, carry_val = 0.0;
    for ( size_type entry = 0 ; entry < nEntry ; ++entry ) {
      tensor.coord(entry, i, j, k);
      val = tensor.value(entry) * ( a[j] * x[k] + a[k] * x[j] );
      if (i == i_prev)
        carry_val += val;
      else {
        y[i_prev] += carry_val;
        carry_val = val;
      }
      i_prev = i;
    }
    y[i] += carry_val;
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

#endif /* #ifndef STOKHOS_COO_PRODUCT_TENSOR_HPP */
