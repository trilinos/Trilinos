// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_FLAT_SPARSE_3_TENSOR_HPP
#define STOKHOS_FLAT_SPARSE_3_TENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_ParameterList.hpp"


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 */
template< typename ValueType , class ExecutionSpace >
class FlatSparse3Tensor {
public:

  typedef ExecutionSpace                       execution_space ;
  typedef typename execution_space::size_type  size_type ;
  typedef ValueType                        value_type ;

private:

  typedef Kokkos::View< size_type[] , execution_space > coord_array_type ;
  typedef Kokkos::View< value_type[], execution_space > value_array_type ;
  typedef Kokkos::View< size_type[], execution_space > entry_array_type ;
  typedef Kokkos::View< size_type[], execution_space > row_map_array_type ;

  coord_array_type   m_k_coord ;
  coord_array_type   m_j_coord ;
  value_array_type   m_value ;
  entry_array_type   m_num_k ;
  entry_array_type   m_num_j ;
  row_map_array_type m_k_row_map ;
  row_map_array_type m_j_row_map ;
  size_type          m_nnz ;
  size_type          m_flops ;

public:

  inline
  ~FlatSparse3Tensor() {}

  inline
  FlatSparse3Tensor() :
    m_k_coord() ,
    m_j_coord() ,
    m_value() ,
    m_num_k() ,
    m_num_j() ,
    m_k_row_map() ,
    m_j_row_map() ,
    m_nnz(0) ,
    m_flops(0) {}

  inline
  FlatSparse3Tensor( const FlatSparse3Tensor & rhs ) :
    m_k_coord( rhs.m_k_coord ) ,
    m_j_coord( rhs.m_j_coord ) ,
    m_value( rhs.m_value ) ,
    m_num_k( rhs.m_num_k ) ,
    m_num_j( rhs.m_num_j ) ,
    m_k_row_map( rhs.m_k_row_map ) ,
    m_j_row_map( rhs.m_j_row_map ) ,
    m_nnz( rhs.m_nnz ) ,
    m_flops( rhs.m_flops ) {}

  inline
  FlatSparse3Tensor & operator = ( const FlatSparse3Tensor & rhs )
  {
    m_k_coord = rhs.m_k_coord ;
    m_j_coord = rhs.m_j_coord ;
    m_value = rhs.m_value ;
    m_num_k = rhs.m_num_k ;
    m_num_j = rhs.m_num_j ;
    m_k_row_map = rhs.m_k_row_map ;
    m_j_row_map = rhs.m_j_row_map ;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_k_row_map.extent(0) - 1 ; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const
  { return m_j_coord.extent(0); }

  /** \brief  Begin k entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type k_begin( size_type i ) const
  { return m_k_row_map[i]; }

  /** \brief  End k entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type k_end( size_type i ) const
  { return m_k_row_map[i] + m_num_k(i); }

  /** \brief  Number of k entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type num_k( size_type i ) const
  { return m_num_k(i); }

  /** \brief  k coordinate for k entry 'kEntry' */
  KOKKOS_INLINE_FUNCTION
  const size_type& k_coord( const size_type kEntry ) const
  { return m_k_coord( kEntry ); }

  /** \brief  Begin j entries with a k entry 'kEntry' */
  KOKKOS_INLINE_FUNCTION
  size_type j_begin( size_type kEntry ) const
  { return m_j_row_map[kEntry]; }

  /** \brief  End j entries with a k entry 'kEntry' */
  KOKKOS_INLINE_FUNCTION
  size_type j_end( size_type kEntry ) const
  { return m_j_row_map[kEntry] + m_num_j(kEntry); }

  /** \brief  Number of j entries with a k entry 'kEntry' */
  KOKKOS_INLINE_FUNCTION
  size_type num_j( size_type kEntry ) const
  { return m_num_j(kEntry); }

  /** \brief  j coordinate for j entry 'jEntry' */
  KOKKOS_INLINE_FUNCTION
  const size_type& j_coord( const size_type jEntry ) const
  { return m_j_coord( jEntry ); }

  /** \brief  Value for j entry 'jEntry' */
  KOKKOS_INLINE_FUNCTION
  const value_type & value( const size_type jEntry ) const
  { return m_value( jEntry ); }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const
  { return m_flops; }

  template <typename OrdinalType>
  static FlatSparse3Tensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;

    // Compute number of k's for each i
    const size_type dimension = basis.size();
    std::vector< size_t > k_coord_work( dimension , (size_t) 0 );
    size_type k_entry_count = 0 ;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin();
         i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      k_coord_work[i] = Cijk.num_k(i_it);
      k_entry_count += Cijk.num_k(i_it);
    }

    // Compute number of j's for each i and k
    std::vector< size_t > j_coord_work( k_entry_count , (size_t) 0 );
    size_type j_entry_count = 0 ;
    size_type k_entry = 0 ;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin();
         i_it!=Cijk.i_end(); ++i_it) {
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it);
           k_it != Cijk.k_end(i_it); ++k_it, ++k_entry) {
        OrdinalType k = index(k_it);
        for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it);
             j_it != Cijk.j_end(k_it); ++j_it) {
          OrdinalType j = index(j_it);
          if (j >= k) {
            ++j_coord_work[k_entry];
            ++j_entry_count;
          }
        }
      }
    }

    /*
    // Pad each row to have size divisible by alignment size
    enum { Align = std::is_same<ExecutionSpace,Kokkos::Cuda>::value ? 32 : 2 };
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      const size_t rem = coord_work[i] % Align;
      if (rem > 0) {
        const size_t pad = Align - rem;
        coord_work[i] += pad;
        entry_count += pad;
      }
    }
    */

    // Allocate tensor data
    FlatSparse3Tensor tensor ;
    tensor.m_k_coord = coord_array_type( "k_coord" , k_entry_count );
    tensor.m_j_coord = coord_array_type( "j_coord" , j_entry_count );
    tensor.m_value = value_array_type( "value" , j_entry_count );
    tensor.m_num_k = entry_array_type( "num_k" , dimension );
    tensor.m_num_j = entry_array_type( "num_j" , k_entry_count );
    tensor.m_k_row_map = row_map_array_type( "k_row_map" , dimension+1 );
    tensor.m_j_row_map = row_map_array_type( "j_row_map" , k_entry_count+1 );

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_k_coord = Kokkos::create_mirror_view( tensor.m_k_coord );
    typename coord_array_type::HostMirror
      host_j_coord = Kokkos::create_mirror_view( tensor.m_j_coord );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_k = Kokkos::create_mirror_view( tensor.m_num_k );
    typename entry_array_type::HostMirror
      host_num_j = Kokkos::create_mirror_view( tensor.m_num_j );
    typename entry_array_type::HostMirror
      host_k_row_map = Kokkos::create_mirror_view( tensor.m_k_row_map );
    typename entry_array_type::HostMirror
      host_j_row_map = Kokkos::create_mirror_view( tensor.m_j_row_map );

    // Compute k row map
    size_type sum = 0;
    host_k_row_map(0) = 0;
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      sum += k_coord_work[i];
      host_k_row_map(i+1) = sum;
      host_num_k(i) = 0;
    }

    // Compute j row map
    sum = 0;
    host_j_row_map(0) = 0;
    for ( size_type i = 0 ; i < k_entry_count ; ++i ) {
      sum += j_coord_work[i];
      host_j_row_map(i+1) = sum;
      host_num_j(i) = 0;
    }

    for ( size_type i = 0 ; i < dimension ; ++i ) {
      k_coord_work[i] = host_k_row_map[i];
    }
    for ( size_type i = 0 ; i < k_entry_count ; ++i ) {
      j_coord_work[i] = host_j_row_map[i];
    }

    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin();
         i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it);
           k_it != Cijk.k_end(i_it); ++k_it) {
        OrdinalType k = index(k_it);
        const size_type kEntry = k_coord_work[i];
        ++k_coord_work[i];
        host_k_coord(kEntry) = k ;
        ++host_num_k(i);
        for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it);
             j_it != Cijk.j_end(k_it); ++j_it) {
          OrdinalType j = index(j_it);
          ValueType c = Stokhos::value(j_it);
          if (j >= k) {
            const size_type jEntry = j_coord_work[kEntry];
            ++j_coord_work[kEntry];
            host_value(jEntry) = (j != k) ? c : 0.5*c;
            host_j_coord(jEntry) = j ;
            ++host_num_j(kEntry);
            ++tensor.m_nnz;
          }
        }
      }
    }

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_k_coord , host_k_coord );
    Kokkos::deep_copy( tensor.m_j_coord , host_j_coord );
    Kokkos::deep_copy( tensor.m_value , host_value );
    Kokkos::deep_copy( tensor.m_num_k , host_num_k );
    Kokkos::deep_copy( tensor.m_num_j , host_num_j );
    Kokkos::deep_copy( tensor.m_k_row_map , host_k_row_map );
    Kokkos::deep_copy( tensor.m_j_row_map , host_j_row_map );

    tensor.m_flops = 5*tensor.m_nnz + dimension;

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType >
FlatSparse3Tensor<ValueType, Device>
create_flat_sparse_3_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList() )
{
  return FlatSparse3Tensor<ValueType, Device>::create( basis, Cijk, params );
}

template <typename ValueType, typename Device>
class BlockMultiply< FlatSparse3Tensor< ValueType , Device > >
{
public:

  typedef typename Device::size_type size_type ;
  typedef FlatSparse3Tensor< ValueType , Device > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {

    const size_type nDim = tensor.dimension();

    // Loop over i
    for ( size_type i = 0; i < nDim; ++i) {
      VectorValue ytmp = 0;

      // Loop over k for this i
      const size_type nk = tensor.num_k(i);
      const size_type kBeg = tensor.k_begin(i);
      const size_type kEnd = kBeg + nk;
      for (size_type kEntry = kBeg; kEntry < kEnd; ++kEntry) {
        const size_type k = tensor.k_coord(kEntry);
        const MatrixValue ak = a[k];
        const VectorValue xk = x[k];

        // Loop over j for this i,k
        const size_type nj = tensor.num_j(kEntry);
        const size_type jBeg = tensor.j_begin(kEntry);
        const size_type jEnd = jBeg + nj;
        for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
          const size_type j = tensor.j_coord(jEntry);
          ytmp += tensor.value(jEntry) * ( a[j] * xk + ak * x[j] );
        }
      }

      y[i] += ytmp ;
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

#endif /* #ifndef STOKHOS_FLAT_SPARSE_3_TENSOR_HPP */
