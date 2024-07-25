// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CRSPRODUCTTENSOR_HPP
#define STOKHOS_CRSPRODUCTTENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"
#include "Stokhos_TinyVec.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 *
 *  This allows product tensor multiplication to be partitioned
 *  on a given coordinate values.
 *
 *  for ( size_type i = 0 ; i < p.dimension() ; ++i ) {
 *    y[i] = 0 ;
 *    for ( size_type e = p.entry_begin(i) ;
 *                    e < p.entry_end(i) ; ++e ) {
 *      const size_type j = p.coord(e,0);
 *      const size_type k = p.coord(e,1);
 *      Scalar tmp = a[j] * x[k] ; if ( j != k ) tmp += a[k] * x[j] ;
 *      y[i] += p.value(e) * tmp ;
 *    }
 *  }
 */
template< typename ValueType, class ExecutionSpace, class Memory = void >
class CrsProductTensor {
public:

  typedef ExecutionSpace  execution_space;
  typedef int         size_type;
  typedef ValueType   value_type;
  typedef Memory      memory_type;

  typedef typename Kokkos::ViewTraits< size_type*, execution_space,void,void >::host_mirror_space host_mirror_space ;
  typedef CrsProductTensor<value_type, host_mirror_space> HostMirror;

// Vectorsize used in multiply algorithm
#if defined(__AVX__)
  static const size_type host_vectorsize = 32/sizeof(value_type);
  static const bool use_intrinsics = true;
  static const size_type num_entry_align = 1;
#elif defined(__MIC__)
  static const size_type host_vectorsize = 16;
  static const bool use_intrinsics = true;
  static const size_type num_entry_align = 8; // avoid use of mask instructions
#else
  static const size_type host_vectorsize = 2;
  static const bool use_intrinsics = false;
  static const size_type num_entry_align = 1;
#endif
  static const size_type cuda_vectorsize = 32;
  static const bool is_cuda =
#if defined( KOKKOS_ENABLE_CUDA )
    std::is_same<ExecutionSpace,Kokkos::Cuda>::value;
#else
    false ;
#endif
  static const size_type vectorsize = is_cuda ? cuda_vectorsize : host_vectorsize;

  // Alignment in terms of number of entries of CRS rows
  static const size_type tensor_align = vectorsize;

private:

  template <class, class, class> friend class CrsProductTensor;

  typedef Kokkos::View< value_type*, Kokkos::LayoutLeft, execution_space, memory_type >  vec_type;
  typedef Kokkos::View< size_type*, Kokkos::LayoutLeft, execution_space, memory_type > coord_array_type;
  typedef Kokkos::View< size_type*[2], Kokkos::LayoutLeft, execution_space, memory_type > coord2_array_type;
  typedef Kokkos::View< value_type*, Kokkos::LayoutLeft, execution_space, memory_type > value_array_type;
  typedef Kokkos::View< size_type*, Kokkos::LayoutLeft, execution_space, memory_type > entry_array_type;
  typedef Kokkos::View< size_type*, Kokkos::LayoutLeft, execution_space, memory_type > row_map_array_type;

  coord_array_type   m_coord;
  coord2_array_type  m_coord2;
  value_array_type   m_value;
  entry_array_type   m_num_entry;
  row_map_array_type m_row_map;
  size_type          m_dim;
  size_type          m_entry_max;
  size_type          m_nnz;
  size_type          m_flops;
  size_type          m_avg_entries_per_row;

  struct CijkRowCount {
    unsigned count;
    unsigned basis;

    CijkRowCount()
      : count(0)
      , basis(0)
      {}
  };

  struct CompareCijkRowCount {
    bool operator() (const CijkRowCount& a, const CijkRowCount& b) const {
      return a.count < b.count;
    }
  };

public:

  KOKKOS_INLINE_FUNCTION
  ~CrsProductTensor() {}

  KOKKOS_INLINE_FUNCTION
  CrsProductTensor() :
    m_coord(),
    m_coord2(),
    m_value(),
    m_num_entry(),
    m_row_map(),
    m_dim(0),
    m_entry_max(0),
    m_nnz(0),
    m_flops(0),
    m_avg_entries_per_row(0) {}

  template <class M>
  KOKKOS_INLINE_FUNCTION
  CrsProductTensor( const CrsProductTensor<value_type,execution_space,M> & rhs ) :
    m_coord( rhs.m_coord ),
    m_coord2( rhs.m_coord2 ),
    m_value( rhs.m_value ),
    m_num_entry( rhs.m_num_entry ),
    m_row_map( rhs.m_row_map ),
    m_dim( rhs.m_dim ),
    m_entry_max( rhs.m_entry_max ),
    m_nnz( rhs.m_nnz ),
    m_flops( rhs.m_flops ),
    m_avg_entries_per_row( rhs.m_avg_entries_per_row ) {}

  template <class M>
  KOKKOS_INLINE_FUNCTION
  CrsProductTensor &
  operator = ( const CrsProductTensor<value_type,execution_space,M> & rhs )
  {
    m_coord = rhs.m_coord;
    m_coord2 = rhs.m_coord2;
    m_value = rhs.m_value;
    m_num_entry = rhs.m_num_entry;
    m_row_map = rhs.m_row_map;
    m_dim = rhs.m_dim;
    m_entry_max = rhs.m_entry_max;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    m_avg_entries_per_row = rhs.m_avg_entries_per_row;
    return *this;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_dim; }

  /** \brief  Is the tensor empty. */
  KOKKOS_INLINE_FUNCTION
  bool is_empty() const { return m_dim == 0; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const
  { return m_coord.extent(0); }

  /** \brief  Maximum sparse entries for any coordinate */
  KOKKOS_INLINE_FUNCTION
  size_type entry_maximum() const
  { return m_entry_max; }

  /** \brief  Begin entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type entry_begin( size_type i ) const
  { return m_row_map[i]; }

  /** \brief  End entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type entry_end( size_type i ) const
  { return m_row_map[i] + m_num_entry(i); }

  /** \brief  Number of entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type num_entry( size_type i ) const
  { return m_num_entry(i); }

  /** \brief  Coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  const size_type& coord( const size_type entry, const size_type c ) const
  { return m_coord2( entry, c ); }

  /** \brief  Coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  const size_type& coord( const size_type entry ) const
  { return m_coord( entry ); }

  /** \brief  Value of an entry */
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

  /** \brief Number average number of entries per row */
  KOKKOS_INLINE_FUNCTION
  size_type avg_entries_per_row() const
  { return m_avg_entries_per_row; }

  template <typename OrdinalType>
  static CrsProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;

    // Note (etp 01/08/15  Commenting out the sorting as it causes a really
    // weird compiler error when compiling with NVCC.  It seems to think the
    // < in CompareCijkRowCount() is part of a template parameter.  We don't
    // seem to use this option, so I am just commenting it out.

    // bool sort_nnz = false;
    // if (params.isParameter("Sort Nonzeros"))
    //   sort_nnz = params.get<bool>("Sort Nonzeros");

    // Compute number of non-zeros for each i
    const size_type dimension = basis.size();
    std::vector< size_t > coord_work( dimension, (size_t) 0 );
    size_type entry_count = 0;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin();
         i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it);
           k_it != Cijk.k_end(i_it); ++k_it) {
        OrdinalType k = index(k_it);
        for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it);
             j_it != Cijk.j_end(k_it); ++j_it) {
          OrdinalType j = index(j_it);
          if (j >= k) {
            ++coord_work[i];
            ++entry_count;
          }
        }
      }
    }

    // Compute average nonzeros per row (must be before padding)
    size_type avg_entries_per_row = entry_count / dimension;

    // Pad each row to have size divisible by alignment size
    for ( size_type i = 0; i < dimension; ++i ) {
      const size_t rem = coord_work[i] % tensor_align;
      if (rem > 0) {
        const size_t pad = tensor_align - rem;
        coord_work[i] += pad;
        entry_count += pad;
      }
    }

    // Sort based on number of non-zeros
    std::vector< CijkRowCount > row_count( dimension );
    for ( size_type i = 0; i < dimension; ++i ) {
      row_count[i].count = coord_work[i];
      row_count[i].basis = i;
    }

    // Note (etp 01/08/15  See above.

    // if (sort_nnz)
    //   std::sort( row_count.begin(), row_count.end(), CompareCijkRowCount() );
    std::vector<size_type> sorted_row_map( dimension );
    for ( size_type i = 0; i < dimension; ++i ) {
      coord_work[i] = row_count[i].count;
      sorted_row_map[ row_count[i].basis ] = i;
    }

    // Allocate tensor data
    // coord and coord2 are initialized to zero because otherwise we get
    // seg faults in the MIC algorithm when processing the tail of each
    // tensor row.  Not quite sure why as the coord loads are padded to
    // length 16 and are masked for the remainder (unless it does load x[j]
    // anyway and masks off the result, so j needs to be valid).
    CrsProductTensor tensor;
    tensor.m_coord = coord_array_type("tensor_coord", entry_count );
    tensor.m_coord2 = coord2_array_type( "tensor_coord2", entry_count );
    tensor.m_value = value_array_type( Kokkos::ViewAllocateWithoutInitializing("tensor_value"), entry_count );
    tensor.m_num_entry = entry_array_type( Kokkos::ViewAllocateWithoutInitializing("tensor_num_entry"), dimension );
    tensor.m_row_map = row_map_array_type( Kokkos::ViewAllocateWithoutInitializing("tensor_row_map"), dimension+1 );
    tensor.m_dim = dimension;
    tensor.m_entry_max = 0;
    tensor.m_avg_entries_per_row = avg_entries_per_row;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( tensor.m_coord );
    typename coord2_array_type::HostMirror
      host_coord2 = Kokkos::create_mirror_view( tensor.m_coord2 );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_entry = Kokkos::create_mirror_view( tensor.m_num_entry );
    typename entry_array_type::HostMirror
      host_row_map = Kokkos::create_mirror_view( tensor.m_row_map );

    // Compute row map
    size_type sum = 0;
    host_row_map(0) = 0;
    for ( size_type i = 0; i < dimension; ++i ) {
      sum += coord_work[i];
      host_row_map(i+1) = sum;
      host_num_entry(i) = 0;
    }

    for ( size_type iCoord = 0; iCoord < dimension; ++iCoord ) {
      coord_work[iCoord] = host_row_map[iCoord];
    }

    // Initialize values and coordinates to zero since we will have extra
    // ones for alignment
    Kokkos::deep_copy( host_value, 0.0 );
    Kokkos::deep_copy( host_coord, 0 );
    Kokkos::deep_copy( host_coord2, 0 );

    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin();
         i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      const size_type row = sorted_row_map[i];
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it);
           k_it != Cijk.k_end(i_it); ++k_it) {
        OrdinalType k = index(k_it);
        for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it);
             j_it != Cijk.j_end(k_it); ++j_it) {
          OrdinalType j = index(j_it);
          ValueType c = Stokhos::value(j_it);
          if (j >= k) {
            const size_type n = coord_work[row]; ++coord_work[row];
            host_value(n) = (j != k) ? c : 0.5*c;
            host_coord2(n,0) = j;
            host_coord2(n,1) = k;
            host_coord(n) = ( k << 16 ) | j;
            ++host_num_entry(row);
            ++tensor.m_nnz;
          }
        }
      }
      // Align num_entry
      host_num_entry(row) =
        (host_num_entry(row) + num_entry_align-1) & ~(num_entry_align-1);
    }

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_coord, host_coord );
    Kokkos::deep_copy( tensor.m_coord2, host_coord2 );
    Kokkos::deep_copy( tensor.m_value, host_value );
    Kokkos::deep_copy( tensor.m_num_entry, host_num_entry );
    Kokkos::deep_copy( tensor.m_row_map, host_row_map );

    for ( size_type i = 0; i < dimension; ++i ) {
      tensor.m_entry_max = std::max( tensor.m_entry_max, host_num_entry(i) );
    }

    tensor.m_flops = 5*tensor.m_nnz + dimension;

    return tensor;
  }

  static CrsProductTensor createMeanBased()
  {
    const size_type dimension = 1;
    const size_type entry_count = 1;

    // Allocate tensor data
    // coord and coord2 are initialized to zero because otherwise we get
    // seg faults in the MIC algorithm when processing the tail of each
    // tensor row.  Not quite sure why as the coord loads are padded to
    // length 16 and are masked for the remainder (unless it does load x[j]
    // anyway and masks off the result, so j needs to be valid).
    CrsProductTensor tensor;
    tensor.m_coord = coord_array_type("tensor_coord", entry_count );
    tensor.m_coord2 = coord2_array_type( "tensor_coord2", entry_count );
    tensor.m_value = value_array_type( Kokkos::ViewAllocateWithoutInitializing("tensor_value"), entry_count );
    tensor.m_num_entry = entry_array_type( Kokkos::ViewAllocateWithoutInitializing("tensor_num_entry"), dimension );
    tensor.m_row_map = row_map_array_type( Kokkos::ViewAllocateWithoutInitializing("tensor_row_map"), dimension+1 );
    tensor.m_dim = dimension;
    tensor.m_entry_max = 1;
    tensor.m_avg_entries_per_row = 1;
    tensor.m_nnz = 1;
    tensor.m_flops = 5*tensor.m_nnz + dimension;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( tensor.m_coord );
    typename coord2_array_type::HostMirror
      host_coord2 = Kokkos::create_mirror_view( tensor.m_coord2 );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_entry = Kokkos::create_mirror_view( tensor.m_num_entry );
    typename entry_array_type::HostMirror
      host_row_map = Kokkos::create_mirror_view( tensor.m_row_map );

    // Compute row map
    host_row_map(0) = 0;
    host_row_map(1) = 1;
    host_num_entry(0) = 1;

    // Compute tensor values
    host_value(0) = 0.5;
    host_coord2(0,0) = 0;
    host_coord2(0,1) = 0;
    host_coord(0) = 0;

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_coord, host_coord );
    Kokkos::deep_copy( tensor.m_coord2, host_coord2 );
    Kokkos::deep_copy( tensor.m_value, host_value );
    Kokkos::deep_copy( tensor.m_num_entry, host_num_entry );
    Kokkos::deep_copy( tensor.m_row_map, host_row_map );

    return tensor;
  }

  static HostMirror
  create_mirror_view( const CrsProductTensor& tensor ) {
    HostMirror host_tensor;

    host_tensor.m_coord     = Kokkos::create_mirror_view( tensor.m_coord );
    host_tensor.m_coord2    = Kokkos::create_mirror_view( tensor.m_coord2 );
    host_tensor.m_value     = Kokkos::create_mirror_view( tensor.m_value );
    host_tensor.m_num_entry = Kokkos::create_mirror_view( tensor.m_num_entry );
    host_tensor.m_row_map   = Kokkos::create_mirror_view( tensor.m_row_map );

    host_tensor.m_dim                 = tensor.m_dim;
    host_tensor.m_entry_max           = tensor.m_entry_max;
    host_tensor.m_avg_entries_per_row = tensor.m_avg_entries_per_row;
    host_tensor.m_nnz                 = tensor.m_nnz;
    host_tensor.m_flops               = tensor.m_flops;

    return host_tensor;
  }

  template < class DstDevice, class DstMemory >
  static void
  deep_copy( const CrsProductTensor<ValueType,DstDevice,DstMemory>& dst ,
             const CrsProductTensor& src ) {
    Kokkos::deep_copy( dst.m_coord,     src.m_coord );
    Kokkos::deep_copy( dst.m_coord2,    src.m_coord2 );
    Kokkos::deep_copy( dst.m_value,     src.m_value );
    Kokkos::deep_copy( dst.m_num_entry, src.m_num_entry );
    Kokkos::deep_copy( dst.m_row_map,   src.m_row_map );
  }

};

template< class Device, typename OrdinalType, typename ValueType>
CrsProductTensor<ValueType, Device>
create_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList())
{
  return CrsProductTensor<ValueType, Device>::create(basis, Cijk, params );
}

template< class Device, typename OrdinalType, typename ValueType,
          class Memory >
CrsProductTensor<ValueType, Device, Memory>
create_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList())
{
  return CrsProductTensor<ValueType, Device, Memory>::create(
    basis, Cijk, params );
}

template< class Device, typename OrdinalType, typename ValueType>
CrsProductTensor<ValueType, Device>
create_mean_based_product_tensor()
{
  return CrsProductTensor<ValueType, Device>::createMeanBased();
}

template< class Device, typename OrdinalType, typename ValueType,
          class Memory >
CrsProductTensor<ValueType, Device, Memory>
create_mean_based_product_tensor()
{
  return CrsProductTensor<ValueType, Device, Memory>::createMeanBased();
}

template < class ValueType, class Device, class Memory >
inline
typename CrsProductTensor<ValueType,Device,Memory>::HostMirror
create_mirror_view( const CrsProductTensor<ValueType,Device,Memory> & src )
{
  return CrsProductTensor<ValueType,Device,Memory>::create_mirror_view( src );
}

  template < class ValueType,
             class DstDevice, class DstMemory,
             class SrcDevice, class SrcMemory >
void
deep_copy( const CrsProductTensor<ValueType,DstDevice,DstMemory> & dst ,
           const CrsProductTensor<ValueType,SrcDevice,SrcMemory> & src )
{
  return CrsProductTensor<ValueType,SrcDevice,SrcMemory>::deep_copy( dst, src );
}

template < typename ValueType, typename Device >
class BlockMultiply< CrsProductTensor< ValueType , Device > >
{
public:

  typedef Device execution_space;
  typedef CrsProductTensor< ValueType , execution_space > tensor_type ;
  typedef typename tensor_type::size_type size_type ;

// Whether to use manual or auto-vectorization
#ifdef __MIC__
#define USE_AUTO_VECTORIZATION 1
#else
#define USE_AUTO_VECTORIZATION 0
#endif

#if defined(__INTEL_COMPILER) && USE_AUTO_VECTORIZATION

  // Version leveraging intel vectorization
  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y ,
                     const VectorValue & alpha = VectorValue(1) )
  {
    // The intel compiler doesn't seem to be able to vectorize through
    // the coord() calls, so extract pointers
    const size_type * cj = &tensor.coord(0,0);
    const size_type * ck = &tensor.coord(0,1);
    const size_type nDim = tensor.dimension();

    for ( size_type iy = 0 ; iy < nDim ; ++iy ) {
      const size_type nEntry = tensor.num_entry(iy);
      const size_type iEntryBeg = tensor.entry_begin(iy);
      const size_type iEntryEnd = iEntryBeg + nEntry;
      VectorValue ytmp = 0;

#pragma simd vectorlength(tensor_type::vectorsize)
#pragma ivdep
#pragma vector aligned
      for (size_type iEntry = iEntryBeg; iEntry<iEntryEnd; ++iEntry) {
        const size_type j    = cj[iEntry]; //tensor.coord(iEntry,0);
        const size_type k    = ck[iEntry]; //tensor.coord(iEntry,1);
        ytmp += tensor.value(iEntry) * ( a[j] * x[k] + a[k] * x[j] );
      }

      y[iy] += alpha * ytmp ;
    }
  }

#elif defined(__MIC__)

  // Version specific to MIC architecture using manual vectorization
  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y ,
                     const VectorValue & alpha = VectorValue(1) )
  {
    const size_type nDim = tensor.dimension();
    for ( size_type iy = 0 ; iy < nDim ; ++iy ) {

      const size_type nEntry = tensor.num_entry(iy);
      const size_type iEntryBeg = tensor.entry_begin(iy);
      const size_type iEntryEnd = iEntryBeg + nEntry;
            size_type iEntry    = iEntryBeg;

      VectorValue ytmp = 0 ;

      const size_type nBlock = nEntry / tensor_type::vectorsize;
      const size_type nEntryB = nBlock * tensor_type::vectorsize;
      const size_type iEnd = iEntryBeg + nEntryB;

      typedef TinyVec<ValueType,tensor_type::vectorsize,tensor_type::use_intrinsics> TV;
      TV vy;
      vy.zero();
      for (size_type block=0; block<nBlock; ++block, iEntry+=tensor_type::vectorsize) {
        const size_type *j = &tensor.coord(iEntry,0);
        const size_type *k = &tensor.coord(iEntry,1);
        TV aj(a, j), ak(a, k), xj(x, j), xk(x, k),
          c(&(tensor.value(iEntry)));

        // vy += c * ( aj * xk + ak * xj)
        aj.times_equal(xk);
        aj.multiply_add(ak, xj);
        vy.multiply_add(c, aj);

      }
      ytmp += vy.sum();

      // The number of nonzeros is always constrained to be a multiple of 8

      const size_type rem = iEntryEnd-iEntry;
      if (rem >= 8) {
        typedef TinyVec<ValueType,8,tensor_type::use_intrinsics> TV2;
        const size_type *j = &tensor.coord(iEntry,0);
        const size_type *k = &tensor.coord(iEntry,1);
        TV2 aj(a, j), ak(a, k), xj(x, j), xk(x, k),
          c(&(tensor.value(iEntry)));

        // vy += c * ( aj * xk + ak * xj)
        aj.times_equal(xk);
        aj.multiply_add(ak, xj);
        aj.times_equal(c);
        ytmp += aj.sum();
      }

      y[iy] += alpha * ytmp ;
    }
  }

#else

  // General version
  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y ,
                     const VectorValue & alpha = VectorValue(1) )
  {
    const size_type nDim = tensor.dimension();
    for ( size_type iy = 0 ; iy < nDim ; ++iy ) {

      const size_type nEntry = tensor.num_entry(iy);
      const size_type iEntryBeg = tensor.entry_begin(iy);
      const size_type iEntryEnd = iEntryBeg + nEntry;
            size_type iEntry    = iEntryBeg;

      VectorValue ytmp = 0 ;

      // Do entries with a blocked loop of size vectorsize
      if (tensor_type::vectorsize > 1 && nEntry >= tensor_type::vectorsize) {
        const size_type nBlock = nEntry / tensor_type::vectorsize;
        const size_type nEntryB = nBlock * tensor_type::vectorsize;
        const size_type iEnd = iEntryBeg + nEntryB;

        typedef TinyVec<ValueType,tensor_type::vectorsize,tensor_type::use_intrinsics> TV;
        TV vy;
        vy.zero();
        for (; iEntry<iEnd; iEntry+=tensor_type::vectorsize) {
          const size_type *j = &tensor.coord(iEntry,0);
          const size_type *k = &tensor.coord(iEntry,1);
          TV aj(a, j), ak(a, k), xj(x, j), xk(x, k), c(&(tensor.value(iEntry)));

          // vy += c * ( aj * xk + ak * xj)
          aj.times_equal(xk);
          aj.multiply_add(ak, xj);
          vy.multiply_add(c, aj);
        }
        ytmp += vy.sum();
      }

      // Do remaining entries with a scalar loop
      for ( ; iEntry<iEntryEnd; ++iEntry) {
        const size_type j = tensor.coord(iEntry,0);
        const size_type k = tensor.coord(iEntry,1);

        ytmp += tensor.value(iEntry) * ( a[j] * x[k] + a[k] * x[j] );
      }

      y[iy] += alpha * ytmp ;
    }
  }
#endif

  KOKKOS_INLINE_FUNCTION
  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  KOKKOS_INLINE_FUNCTION
  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

// Specialization of Multiply< BlockCrsMatrix< BlockSpec, ... > > > for
// CrsProductTensor, which provides a version that processes blocks of FEM
// columns together to reduce the number of global reads of the sparse 3 tensor

// Even though this isn't specific to Threads, templating on Device creates a
// duplicate specialization error for Cuda.  Need to see if we can fix this,
// or put the implementation in another class easily specialized for Threads,
// OpenMP, ...
template< typename ValueType , typename MatrixValue , typename VectorValue ,
          typename Device >
class MultiplyImpl {
public:

  typedef Device execution_space ;
  typedef CrsProductTensor< ValueType , execution_space > tensor_type;
  typedef StochasticProductTensor< ValueType, tensor_type, execution_space > BlockSpec;
  typedef typename BlockSpec::size_type size_type ;
  typedef Kokkos::View< VectorValue** , Kokkos::LayoutLeft , execution_space > block_vector_type ;
  typedef BlockCrsMatrix< BlockSpec , MatrixValue , execution_space >  matrix_type ;

  const matrix_type  m_A ;
  const block_vector_type  m_x ;
  const block_vector_type  m_y ;

  MultiplyImpl( const matrix_type & A ,
                const block_vector_type & x ,
                const block_vector_type & y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  {}

  //--------------------------------------------------------------------------
  //  A( storage_size( m_A.block.size() ) , m_A.graph.row_map.size() );
  //  x( m_A.block.dimension() , m_A.graph.row_map.first_count() );
  //  y( m_A.block.dimension() , m_A.graph.row_map.first_count() );
  //

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iBlockRow ) const
  {
    // Prefer that y[ m_A.block.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    VectorValue * const y = & m_y(0,iBlockRow);

    const size_type iEntryBegin = m_A.graph.row_map[ iBlockRow ];
    const size_type iEntryEnd   = m_A.graph.row_map[ iBlockRow + 1 ];

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = 0 ; j < m_A.block.dimension() ; ++j ) { y[j] = 0 ; }

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      const VectorValue * const x = & m_x( 0 , m_A.graph.entries(iEntry) );
      const MatrixValue * const a = & m_A.values( 0 , iEntry );

      BlockMultiply< BlockSpec >::apply( m_A.block , a , x , y );
    }

  }

  /*
   * Compute work range = (begin, end) such that adjacent threads write to
   * separate cache lines
   */
  KOKKOS_INLINE_FUNCTION
  std::pair< size_type , size_type >
  compute_work_range( const size_type work_count ,
                      const size_type thread_count ,
                      const size_type thread_rank ) const
  {
    enum { work_align = 64 / sizeof(VectorValue) };
    enum { work_shift = Stokhos::power_of_two< work_align >::value };
    enum { work_mask  = work_align - 1 };

    const size_type work_per_thread =
      ( ( ( ( work_count + work_mask ) >> work_shift ) + thread_count - 1 ) /
        thread_count ) << work_shift ;

    const size_type work_begin =
      std::min( thread_rank * work_per_thread , work_count );
    const size_type work_end   =
      std::min( work_begin + work_per_thread , work_count );

    return std::make_pair( work_begin , work_end );
  }

#if defined(__MIC__)

  // A MIC-specific version of the block-multiply algorithm, where block here
  // means processing multiple FEM columns at a time
  KOKKOS_INLINE_FUNCTION
  void operator()( const typename Kokkos::TeamPolicy< execution_space >::member_type & device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A.graph.row_map.extent(0)-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    std::pair<size_type,size_type> work_range =
      compute_work_range(m_A.block.dimension(), num_thread, thread_idx);

    // Prefer that y[ m_A.block.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    VectorValue * const y = & m_y(0,iBlockRow);

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = work_range.first ; j < work_range.second ; ++j )
      y[j] = 0 ;

    const tensor_type& tensor = m_A.block.tensor();

    const size_type iBlockEntryBeg = m_A.graph.row_map[ iBlockRow ];
    const size_type iBlockEntryEnd = m_A.graph.row_map[ iBlockRow + 1 ];
    const size_type BlockSize = 9;
    const size_type numBlock =
      (iBlockEntryEnd-iBlockEntryBeg+BlockSize-1) / BlockSize;

    const MatrixValue* sh_A[BlockSize];
    const VectorValue* sh_x[BlockSize];

    size_type iBlockEntry = iBlockEntryBeg;
    for (size_type block = 0; block<numBlock; ++block, iBlockEntry+=BlockSize) {
      const size_type block_size =
        block == numBlock-1 ? iBlockEntryEnd-iBlockEntry : BlockSize;

      for ( size_type col = 0; col < block_size; ++col ) {
        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry + col );
        sh_x[col] = & m_x( 0 , iBlockColumn );
        sh_A[col] = & m_A.values( 0 , iBlockEntry + col );
      }

      for ( size_type iy = work_range.first ; iy < work_range.second ; ++iy ) {

        const size_type nEntry = tensor.num_entry(iy);
        const size_type iEntryBeg = tensor.entry_begin(iy);
        const size_type iEntryEnd = iEntryBeg + nEntry;
              size_type iEntry    = iEntryBeg;

        VectorValue ytmp = 0 ;

        // Do entries with a blocked loop of size blocksize
        const size_type nBlock = nEntry / tensor_type::vectorsize;
        const size_type nEntryB = nBlock * tensor_type::vectorsize;
        const size_type iEnd = iEntryBeg + nEntryB;

        typedef TinyVec<ValueType,tensor_type::vectorsize,tensor_type::use_intrinsics> ValTV;
        typedef TinyVec<MatrixValue,tensor_type::vectorsize,tensor_type::use_intrinsics> MatTV;
        typedef TinyVec<VectorValue,tensor_type::vectorsize,tensor_type::use_intrinsics> VecTV;
        VecTV vy;
        vy.zero();
        for (size_type block=0; block<nBlock; ++block, iEntry+=tensor_type::vectorsize) {
          const size_type *j = &tensor.coord(iEntry,0);
          const size_type *k = &tensor.coord(iEntry,1);
          ValTV c(&(tensor.value(iEntry)));

          for ( size_type col = 0; col < block_size; ++col ) {
            MatTV aj(sh_A[col], j), ak(sh_A[col], k);
            VecTV xj(sh_x[col], j), xk(sh_x[col], k);

            // vy += c * ( aj * xk + ak * xj)
            aj.times_equal(xk);
            aj.multiply_add(ak, xj);
            vy.multiply_add(c, aj);
          }
        }
        ytmp += vy.sum();

        // The number of nonzeros is always constrained to be a multiple of 8

        const size_type rem = iEntryEnd-iEntry;
        if (rem >= 8) {
          typedef TinyVec<ValueType,8,tensor_type::use_intrinsics> ValTV2;
          typedef TinyVec<MatrixValue,8,tensor_type::use_intrinsics> MatTV2;
          typedef TinyVec<VectorValue,8,tensor_type::use_intrinsics> VecTV2;
          const size_type *j = &tensor.coord(iEntry,0);
          const size_type *k = &tensor.coord(iEntry,1);
          ValTV2 c(&(tensor.value(iEntry)));

          for ( size_type col = 0; col < block_size; ++col ) {
            MatTV2 aj(sh_A[col], j), ak(sh_A[col], k);
            VecTV2 xj(sh_x[col], j), xk(sh_x[col], k);

            // vy += c * ( aj * xk + ak * xj)
            aj.times_equal(xk);
            aj.multiply_add(ak, xj);
            aj.times_equal(c);
            ytmp += aj.sum();
          }
        }

        y[iy] += ytmp ;
      }

      // Add a team barrier to keep the thread team in-sync before going on
      // to the next block
      device.team_barrier();
    }

  }

#else

  // A general hand-vectorized version of the block multiply algorithm, where
  // block here means processing multiple FEM columns at a time.  Note that
  // auto-vectorization of a block algorithm doesn't work, because the
  // stochastic loop is not the inner-most loop.
  KOKKOS_INLINE_FUNCTION
  void operator()( const typename Kokkos::TeamPolicy< execution_space >::member_type & device ) const
  {
    const size_type iBlockRow = device.league_rank();

    // Check for valid row
    const size_type row_count = m_A.graph.row_map.extent(0)-1;
    if (iBlockRow >= row_count)
      return;

    const size_type num_thread = device.team_size();
    const size_type thread_idx = device.team_rank();
    std::pair<size_type,size_type> work_range =
      compute_work_range(m_A.block.dimension(), num_thread, thread_idx);

    // Prefer that y[ m_A.block.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    VectorValue * const y = & m_y(0,iBlockRow);

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = work_range.first ; j < work_range.second ; ++j )
      y[j] = 0 ;

    const tensor_type& tensor = m_A.block.tensor();

    const size_type iBlockEntryBeg = m_A.graph.row_map[ iBlockRow ];
    const size_type iBlockEntryEnd = m_A.graph.row_map[ iBlockRow + 1 ];
    const size_type BlockSize = 14;
    const size_type numBlock =
      (iBlockEntryEnd-iBlockEntryBeg+BlockSize-1) / BlockSize;

    const MatrixValue* sh_A[BlockSize];
    const VectorValue* sh_x[BlockSize];

    size_type iBlockEntry = iBlockEntryBeg;
    for (size_type block = 0; block<numBlock; ++block, iBlockEntry+=BlockSize) {
      const size_type block_size =
        block == numBlock-1 ? iBlockEntryEnd-iBlockEntry : BlockSize;

      for ( size_type col = 0; col < block_size; ++col ) {
        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry + col );
        sh_x[col] = & m_x( 0 , iBlockColumn );
        sh_A[col] = & m_A.values( 0 , iBlockEntry + col );
      }

      for ( size_type iy = work_range.first ; iy < work_range.second ; ++iy ) {

        const size_type nEntry = tensor.num_entry(iy);
        const size_type iEntryBeg = tensor.entry_begin(iy);
        const size_type iEntryEnd = iEntryBeg + nEntry;
              size_type iEntry    = iEntryBeg;

        VectorValue ytmp = 0 ;

        // Do entries with a blocked loop of size blocksize
        if (tensor_type::vectorsize > 1 && nEntry >= tensor_type::vectorsize) {
          const size_type nBlock = nEntry / tensor_type::vectorsize;
          const size_type nEntryB = nBlock * tensor_type::vectorsize;
          const size_type iEnd = iEntryBeg + nEntryB;

          typedef TinyVec<ValueType,tensor_type::vectorsize,tensor_type::use_intrinsics> ValTV;
          typedef TinyVec<MatrixValue,tensor_type::vectorsize,tensor_type::use_intrinsics> MatTV;
          typedef TinyVec<VectorValue,tensor_type::vectorsize,tensor_type::use_intrinsics> VecTV;
          VecTV vy;
          vy.zero();
          for (; iEntry<iEnd; iEntry+=tensor_type::vectorsize) {
            const size_type *j = &tensor.coord(iEntry,0);
            const size_type *k = &tensor.coord(iEntry,1);
            ValTV c(&(tensor.value(iEntry)));

            for ( size_type col = 0; col < block_size; ++col ) {
              MatTV aj(sh_A[col], j), ak(sh_A[col], k);
              VecTV xj(sh_x[col], j), xk(sh_x[col], k);

              // vy += c * ( aj * xk + ak * xj)
              aj.times_equal(xk);
              aj.multiply_add(ak, xj);
              vy.multiply_add(c, aj);
            }
          }
          ytmp += vy.sum();
        }

        // Do remaining entries with a scalar loop
        for ( ; iEntry<iEntryEnd; ++iEntry) {
          const size_type j = tensor.coord(iEntry,0);
          const size_type k = tensor.coord(iEntry,1);
          ValueType cijk = tensor.value(iEntry);

          for ( size_type col = 0; col < block_size; ++col ) {
            ytmp += cijk * ( sh_A[col][j] * sh_x[col][k] +
                             sh_A[col][k] * sh_x[col][j] );
          }

        }

        y[iy] += ytmp ;
      }

      // Add a team barrier to keep the thread team in-sync before going on
      // to the next block
      device.team_barrier();
    }

  }

#endif

  static void apply( const matrix_type & A ,
                     const block_vector_type & x ,
                     const block_vector_type & y )
  {
    // Generally the block algorithm seems to perform better on the MIC,
    // as long as the stochastic size isn't too big, but doesn't perform
    // any better on the CPU (probably because the CPU has a fat L3 cache
    // to store the sparse 3 tensor).
#ifdef __MIC__
    const bool use_block_algorithm = true;
#else
    const bool use_block_algorithm = false;
#endif

    const size_t row_count = A.graph.row_map.extent(0) - 1 ;
    if (use_block_algorithm) {
#ifdef __MIC__
      const size_t team_size = 4;  // 4 hyperthreads for MIC
#else
      const size_t team_size = 2;  // 2 for everything else
#endif
      const size_t league_size = row_count;
      Kokkos::TeamPolicy< execution_space > config(league_size, team_size);
      Kokkos::parallel_for( config , MultiplyImpl(A,x,y) );
    }
    else {
      Kokkos::parallel_for( row_count , MultiplyImpl(A,x,y) );
    }
  }
};

//----------------------------------------------------------------------------

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// Inject some functions into the Kokkos namespace
namespace Kokkos {

  using Stokhos::create_mirror_view;
  using Stokhos::deep_copy;

} // namespace Kokkos

#endif /* #ifndef STOKHOS_CRSPRODUCTTENSOR_HPP */
