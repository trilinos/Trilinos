// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_Sparse3TensorPartition.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_TinyVec.hpp"


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

template< typename ValueType, class ExecutionSpace >
class SimpleTiledCrsProductTensor {
public:

  typedef ExecutionSpace execution_space;
  typedef int size_type;
  typedef ValueType value_type;

// Vectorsize used in multiply algorithm
#if defined(__AVX__)
  static const size_type host_vectorsize = 32/sizeof(value_type);
  static const bool use_intrinsics = true;
#elif defined(__MIC__)
  static const size_type host_vectorsize = 16;
  static const bool use_intrinsics = true;
#else
  static const size_type host_vectorsize = 2;
  static const bool use_intrinsics = false;
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

  typedef Kokkos::View< value_type[], execution_space >  vec_type;
  typedef Kokkos::View< value_type[], execution_space > value_array_type;
  typedef Kokkos::View< size_type[], execution_space > coord_array_type;
  typedef Kokkos::View< size_type[][2], Kokkos::LayoutLeft, execution_space > coord2_array_type;
  typedef Kokkos::View< size_type*, execution_space > i_begin_type;
  typedef Kokkos::View< size_type*, execution_space > i_size_type;
  typedef Kokkos::View< size_type*, execution_space > num_j_type;
  typedef Kokkos::View< size_type**, execution_space > j_begin_type;
  typedef Kokkos::View< size_type**, execution_space > j_size_type;
  typedef Kokkos::View< size_type**, execution_space > num_k_type;
  typedef Kokkos::View< size_type***, execution_space > k_begin_type;
  typedef Kokkos::View< size_type***, execution_space > k_size_type;

  typedef Kokkos::View< size_type****,  Kokkos::LayoutRight, execution_space > row_map_type;
  typedef Kokkos::View< size_type****,  Kokkos::LayoutRight, execution_space > num_entry_type;

  value_array_type   m_value;
  coord_array_type   m_coord;
  coord2_array_type  m_coord2;
  i_begin_type       m_i_begin;
  i_size_type        m_i_size;
  num_j_type         m_num_j;
  j_begin_type       m_j_begin;
  j_size_type        m_j_size;
  num_k_type         m_num_k;
  k_begin_type       m_k_begin;
  k_size_type        m_k_size;
  row_map_type       m_row_map;
  num_entry_type     m_num_entry;
  size_type          m_dimension;
  size_type          m_max_i_tile_size;
  size_type          m_max_jk_tile_size;
  size_type          m_num_i;
  size_type          m_nnz;
  size_type          m_flops;

  struct Coord {
    size_type i, j, k;
    value_type cijk;
  };

  template <typename coord_t>
  struct Tile {
    size_type lower, upper;
    Teuchos::Array<coord_t> parts;
  };

  typedef Tile<Coord> KTile;
  typedef Tile<KTile> JTile;
  typedef Tile<JTile> ITile;

public:

  inline
  ~SimpleTiledCrsProductTensor() {}

  inline
  SimpleTiledCrsProductTensor() :
    m_value(),
    m_coord(),
    m_coord2(),
    m_i_begin(),
    m_i_size(),
    m_num_j(),
    m_j_begin(),
    m_j_size(),
    m_num_k(),
    m_k_begin(),
    m_k_size(),
    m_row_map(),
    m_num_entry(),
    m_dimension(0),
    m_max_i_tile_size(0),
    m_max_jk_tile_size(0),
    m_num_i(0),
    m_nnz(0),
    m_flops(0) {}

  inline
  SimpleTiledCrsProductTensor(const SimpleTiledCrsProductTensor & rhs) :
    m_value(rhs.m_value),
    m_coord(rhs.m_coord),
    m_coord2(rhs.m_coord2),
    m_i_begin(rhs.m_i_begin),
    m_i_size(rhs.m_i_size),
    m_num_j(rhs.m_num_j),
    m_j_begin(rhs.m_j_begin),
    m_j_size(rhs.m_j_size),
    m_num_k(rhs.m_num_k),
    m_k_begin(rhs.m_k_begin),
    m_k_size(rhs.m_k_size),
    m_row_map(rhs.m_row_map),
    m_num_entry(rhs.m_num_entry),
    m_dimension(rhs.m_dimension),
    m_max_i_tile_size(rhs.m_max_i_tile_size),
    m_max_jk_tile_size(rhs.m_max_jk_tile_size),
    m_num_i(rhs.m_num_i),
    m_nnz(rhs.m_nnz),
    m_flops(rhs.m_flops) {}

  inline
  SimpleTiledCrsProductTensor& operator=(
    const SimpleTiledCrsProductTensor & rhs)
  {
    m_value = rhs.m_value;
    m_coord = rhs.m_coord;
    m_coord2 = rhs.m_coord2;
    m_i_begin = rhs.m_i_begin;
    m_i_size = rhs.m_i_size;
    m_num_j = rhs.m_num_j;
    m_j_begin = rhs.m_j_begin;
    m_j_size = rhs.m_j_size;
    m_num_k = rhs.m_num_k;
    m_k_begin = rhs.m_k_begin;
    m_k_size = rhs.m_k_size;
    m_row_map = rhs.m_row_map;
    m_num_entry = rhs.m_num_entry;
    m_dimension = rhs.m_dimension;
    m_max_i_tile_size = rhs.m_max_i_tile_size;
    m_max_jk_tile_size = rhs.m_max_jk_tile_size;
    m_num_i = rhs.m_num_i;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    return *this;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_dimension; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const { return m_coord.extent(0); }

  /** \brief Number i-tiles */
  KOKKOS_INLINE_FUNCTION
  size_type num_i_tiles() const { return m_num_i; }

  /** \brief  Begin entries with for tile 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type i_begin(const size_type i) const { return m_i_begin(i); }

  /** \brief  Number of entries with for tile 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type i_size(const size_type i) const { return m_i_size(i); }

  /** \brief Number j-tiles */
  KOKKOS_INLINE_FUNCTION
  size_type num_j_tiles(const size_type i) const { return m_num_j(i); }

  /** \brief  Begin entries with for tile 'i,j' */
  KOKKOS_INLINE_FUNCTION
  size_type j_begin(const size_type i, const size_type j) const {
    return m_j_begin(i,j);
  }

  /** \brief  Number of entries with for tile 'i,j' */
  KOKKOS_INLINE_FUNCTION
  size_type j_size(const size_type i, const size_type j) const {
    return m_j_size(i,j);
  }

  /** \brief Number k-tiles */
  KOKKOS_INLINE_FUNCTION
  size_type num_k_tiles(const size_type i, const size_type j) const {
    return m_num_k(i,j); }

  /** \brief  Begin entries with for tile 'i,j,k' */
  KOKKOS_INLINE_FUNCTION
  size_type k_begin(const size_type i, const size_type j,
                    const size_type k) const {
    return m_k_begin(i,j,k);
  }

  /** \brief  Number of entries with for tile 'i,j' */
  KOKKOS_INLINE_FUNCTION
  size_type k_size(const size_type i, const size_type j,
                   const size_type k) const {
    return m_k_size(i,j,k);
  }

  /** \brief  Number of entries for tile (i,j,k) and row r */
  KOKKOS_INLINE_FUNCTION
  size_type num_entry(const size_type i, const size_type j,
                      const size_type k, const size_type r) const {
    return m_num_entry(i,j,k,r);
  }

  /** \brief  Begin entries for tile (i,j,k) and row r */
  KOKKOS_INLINE_FUNCTION
  size_type entry_begin(const size_type i, const size_type j,
                        const size_type k, const size_type r) const {
    return m_row_map(i,j,k,r);
  }

  /** \brief  End entries for tile (i,j,k) and row r */
  KOKKOS_INLINE_FUNCTION
  size_type entry_end(const size_type i, const size_type j,
                      const size_type k, const size_type r) const {
    return m_row_map(i,j,k,r) + m_num_entry(i,j,k,r);
  }

  /** \brief  Coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  const size_type& coord(const size_type entry, const size_type c) const {
    return m_coord2(entry, c);
  }

  /** \brief  Coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  const size_type& coord(const size_type entry) const {
    return m_coord(entry);
  }

  /** \brief  Value of an entry */
  KOKKOS_INLINE_FUNCTION
  const value_type & value(const size_type entry) const {
    return m_value(entry);
  }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const
  { return m_flops; }

  /** \brief Max size of any i tile */
  KOKKOS_INLINE_FUNCTION
  size_type max_i_tile_size() const { return m_max_i_tile_size; }

  /** \brief Max size of any j/k tile */
  KOKKOS_INLINE_FUNCTION
  size_type max_jk_tile_size() const { return m_max_jk_tile_size; }

  template <typename OrdinalType>
  static SimpleTiledCrsProductTensor
  create(const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params)
  {
    using Teuchos::rcp;
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::Array;

    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    typedef typename Cijk_type::i_iterator i_iterator;
    typedef typename Cijk_type::ik_iterator ik_iterator;
    typedef typename Cijk_type::ikj_iterator ikj_iterator;

    const size_type i_tile_size = params.get<OrdinalType>("Tile Size");

    // Build 2-way symmetric Cijk tensor
    Cijk_type Cijk_sym;
    i_iterator i_begin = Cijk.i_begin();
    i_iterator i_end = Cijk.i_end();
    for (i_iterator i_it=i_begin; i_it!=i_end; ++i_it) {
      OrdinalType i = index(i_it);
      ik_iterator k_begin = Cijk.k_begin(i_it);
      ik_iterator k_end = Cijk.k_end(i_it);
      for (ik_iterator k_it = k_begin; k_it != k_end; ++k_it) {
        OrdinalType k = index(k_it);
        ikj_iterator j_begin = Cijk.j_begin(k_it);
        ikj_iterator j_end = Cijk.j_end(k_it);
        for (ikj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          OrdinalType j = index(j_it);
          if (k <= j) {
            ValueType c = Stokhos::value(j_it);
            Cijk_sym.add_term(i, j, k, c);
          }
        }
      }
    }
    Cijk_sym.fillComplete();

    // First partition based on i
    size_type j_tile_size = i_tile_size / 2;
    size_type basis_size = basis.size();
    size_type num_i_parts = (basis_size + i_tile_size-1) / i_tile_size;
    //size_type its = basis_size / num_i_parts;
    size_type its = i_tile_size;
    Array<ITile> i_tiles(num_i_parts);
    for (size_type i=0; i<num_i_parts; ++i) {
      i_tiles[i].lower = i*its;
      i_tiles[i].upper = std::min(basis_size, (i+1)*its);
      i_tiles[i].parts.resize(1);
      i_tiles[i].parts[0].lower = basis_size;
      i_tiles[i].parts[0].upper = 0;
    }

    // Next partition j
    size_type max_jk_tile_size = 0;
    for (i_iterator i_it=Cijk_sym.i_begin(); i_it!=Cijk_sym.i_end(); ++i_it) {
      OrdinalType i = index(i_it);

      // Find which part i belongs to
      size_type idx = 0;
      while (idx < num_i_parts && i >= i_tiles[idx].lower) ++idx;
      --idx;
      TEUCHOS_ASSERT(idx >= 0 && idx < num_i_parts);

      ik_iterator k_begin = Cijk_sym.k_begin(i_it);
      ik_iterator k_end = Cijk_sym.k_end(i_it);
      for (ik_iterator k_it = k_begin; k_it != k_end; ++k_it) {
        OrdinalType j = index(k_it);  // using symmetry to interchange j and k

        if (j < i_tiles[idx].parts[0].lower)
          i_tiles[idx].parts[0].lower = j;
        if (j > i_tiles[idx].parts[0].upper)
          i_tiles[idx].parts[0].upper = j;
      }
    }
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type lower = i_tiles[idx].parts[0].lower;
      size_type upper = i_tiles[idx].parts[0].upper;
      size_type range = upper - lower + 1;
      size_type num_j_parts = (range + j_tile_size-1) / j_tile_size;
      //size_type jts = range / num_j_parts;
      size_type jts = j_tile_size;
      max_jk_tile_size = std::max(max_jk_tile_size, jts);
      Array<JTile> j_tiles(num_j_parts);
      for (size_type j=0; j<num_j_parts; ++j) {
        j_tiles[j].lower = lower + j*jts;
        j_tiles[j].upper = std::min(upper+1, lower + (j+1)*jts);
        j_tiles[j].parts.resize(1);
        j_tiles[j].parts[0].lower = basis_size;
        j_tiles[j].parts[0].upper = 0;
      }
      i_tiles[idx].parts.swap(j_tiles);
    }

    // Now partition k
    for (i_iterator i_it=Cijk_sym.i_begin(); i_it!=Cijk_sym.i_end(); ++i_it) {
      OrdinalType i = index(i_it);

      // Find which part i belongs to
      size_type idx = 0;
      while (idx < num_i_parts && i >= i_tiles[idx].lower) ++idx;
      --idx;
      TEUCHOS_ASSERT(idx >= 0 && idx < num_i_parts);

      ik_iterator k_begin = Cijk_sym.k_begin(i_it);
      ik_iterator k_end = Cijk_sym.k_end(i_it);
      for (ik_iterator k_it = k_begin; k_it != k_end; ++k_it) {
        OrdinalType j = index(k_it);  // using symmetry to interchange j and k

        // Find which part j belongs to
        size_type num_j_parts = i_tiles[idx].parts.size();
        size_type jdx = 0;
        while (jdx < num_j_parts && j >= i_tiles[idx].parts[jdx].lower) ++jdx;
        --jdx;
        TEUCHOS_ASSERT(jdx >= 0 && jdx < num_j_parts);

        ikj_iterator j_begin = Cijk_sym.j_begin(k_it);
        ikj_iterator j_end = Cijk_sym.j_end(k_it);
        for (ikj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          OrdinalType k = index(j_it);  // using symmetry to interchange j and k
          ValueType cijk = Stokhos::value(j_it);
          if (k >= j) {
            Coord coord;
            coord.i = i; coord.j = j; coord.k = k; coord.cijk = cijk;
            i_tiles[idx].parts[jdx].parts[0].parts.push_back(coord);
            if (k < i_tiles[idx].parts[jdx].parts[0].lower)
              i_tiles[idx].parts[jdx].parts[0].lower = k;
            if (k > i_tiles[idx].parts[jdx].parts[0].upper)
              i_tiles[idx].parts[jdx].parts[0].upper = k;
          }
        }
      }
    }

    // Now need to divide up k-parts based on lower/upper bounds
    size_type num_coord = 0;
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type num_j_parts = i_tiles[idx].parts.size();
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        size_type lower = i_tiles[idx].parts[jdx].parts[0].lower;
        size_type upper = i_tiles[idx].parts[jdx].parts[0].upper;
        size_type range = upper - lower + 1;
        size_type num_k_parts = (range + j_tile_size-1) / j_tile_size;
        //size_type kts = range / num_k_parts;
        size_type kts = j_tile_size;
        max_jk_tile_size = std::max(max_jk_tile_size, kts);
        Array<KTile> k_tiles(num_k_parts);
        for (size_type k=0; k<num_k_parts; ++k) {
          k_tiles[k].lower = lower + k*kts;
          k_tiles[k].upper = std::min(upper+1, lower + (k+1)*kts);
        }
        size_type num_k = i_tiles[idx].parts[jdx].parts[0].parts.size();
        for (size_type l=0; l<num_k; ++l) {
          size_type i = i_tiles[idx].parts[jdx].parts[0].parts[l].i;
          size_type j = i_tiles[idx].parts[jdx].parts[0].parts[l].j;
          size_type k = i_tiles[idx].parts[jdx].parts[0].parts[l].k;
          value_type cijk = i_tiles[idx].parts[jdx].parts[0].parts[l].cijk;

          // Find which part k belongs to
          size_type kdx = 0;
          while (kdx < num_k_parts && k >= k_tiles[kdx].lower) ++kdx;
          --kdx;
          TEUCHOS_ASSERT(kdx >= 0 && kdx < num_k_parts);

          Coord coord;
          coord.i = i; coord.j = j; coord.k = k; coord.cijk = cijk;
          k_tiles[kdx].parts.push_back(coord);
          ++num_coord;
          if (j != k) ++num_coord;
        }

        // Eliminate parts with zero size
        Array<KTile> k_tiles2;
        for (size_type k=0; k<num_k_parts; ++k) {
          if (k_tiles[k].parts.size() > 0)
            k_tiles2.push_back(k_tiles[k]);
        }
        i_tiles[idx].parts[jdx].parts.swap(k_tiles2);
      }
    }
    TEUCHOS_ASSERT(num_coord == Cijk.num_entries());

    // Compute number of non-zeros for each row in each part
    size_type total_num_rows = 0, max_num_rows = 0, entry_count = 0;
    size_type max_num_j_parts = 0, max_num_k_parts = 0;
    Array< Array< Array< Array<size_type> > > > coord_work(num_i_parts);
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type num_j_parts = i_tiles[idx].parts.size();
      max_num_j_parts = std::max(max_num_j_parts, num_j_parts);
      coord_work[idx].resize(num_j_parts);
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        size_type num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        max_num_k_parts = std::max(max_num_k_parts, num_k_parts);
        coord_work[idx][jdx].resize(num_k_parts);
        for (size_type kdx=0; kdx<num_k_parts; ++kdx) {
          size_type num_rows = i_tiles[idx].upper - i_tiles[idx].lower + 1;
          total_num_rows += num_rows;
          max_num_rows = std::max(max_num_rows, num_rows);
          coord_work[idx][jdx][kdx].resize(num_rows, 0);

          size_type nc = i_tiles[idx].parts[jdx].parts[kdx].parts.size();
          for (size_type c=0; c<nc; ++c) {
            size_type i = i_tiles[idx].parts[jdx].parts[kdx].parts[c].i;
            size_type i_begin = i_tiles[idx].lower;
            ++(coord_work[idx][jdx][kdx][i-i_begin]);
            ++entry_count;
          }
        }
      }
    }

    // Pad each row to have size divisible by alignment size
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type num_j_parts = i_tiles[idx].parts.size();
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        size_type num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        for (size_type kdx=0; kdx<num_k_parts; ++kdx) {
          size_type sz = coord_work[idx][jdx][kdx].size();
          for (size_type i = 0; i < sz; ++i) {
            const size_t rem = coord_work[idx][jdx][kdx][i] % tensor_align;
            if (rem > 0) {
              const size_t pad = tensor_align - rem;
              coord_work[idx][jdx][kdx][i] += pad;
              entry_count += pad;
            }
          }
        }
      }
    }

    // Allocate tensor data
    SimpleTiledCrsProductTensor tensor;
    tensor.m_value = value_array_type("value", entry_count);
    tensor.m_coord = coord_array_type("coord", entry_count);
    tensor.m_coord2 = coord2_array_type("coord2", entry_count);
    tensor.m_i_begin = i_begin_type("i_begin", num_i_parts);
    tensor.m_i_size = i_size_type("i_size", num_i_parts);
    tensor.m_num_j = num_j_type("num_j", num_i_parts);
    tensor.m_j_begin = j_begin_type("j_begin", num_i_parts, max_num_j_parts);
    tensor.m_j_size = j_size_type("j_size", num_i_parts, max_num_j_parts);
    tensor.m_num_k = num_k_type("num_k", num_i_parts, max_num_j_parts);
    tensor.m_k_begin = k_begin_type("k_begin", num_i_parts, max_num_j_parts,
                                    max_num_k_parts);
    tensor.m_k_size = k_size_type("k_size", num_i_parts, max_num_j_parts,
                                  max_num_k_parts);
    tensor.m_row_map = row_map_type("row_map", num_i_parts,
                                    max_num_j_parts, max_num_k_parts,
                                    max_num_rows+1);
    tensor.m_num_entry = num_entry_type("num_entry", num_i_parts,
                                        max_num_j_parts, max_num_k_parts,
                                        max_num_rows);
    tensor.m_dimension = basis.size();
    tensor.m_max_i_tile_size = i_tile_size;
    tensor.m_max_jk_tile_size = max_jk_tile_size;
    tensor.m_num_i = num_i_parts;

    // Create mirror, is a view if is host memory
    typename value_array_type::HostMirror host_value =
      Kokkos::create_mirror_view(tensor.m_value);
    typename coord_array_type::HostMirror host_coord =
      Kokkos::create_mirror_view(tensor.m_coord);
    typename coord2_array_type::HostMirror host_coord2 =
      Kokkos::create_mirror_view(tensor.m_coord2);
    typename i_begin_type::HostMirror host_i_begin =
      Kokkos::create_mirror_view(tensor.m_i_begin);
    typename i_size_type::HostMirror host_i_size =
      Kokkos::create_mirror_view(tensor.m_i_size);
    typename num_j_type::HostMirror host_num_j =
      Kokkos::create_mirror_view(tensor.m_num_j);
    typename j_begin_type::HostMirror host_j_begin =
      Kokkos::create_mirror_view(tensor.m_j_begin);
    typename j_size_type::HostMirror host_j_size =
      Kokkos::create_mirror_view(tensor.m_j_size);
    typename num_k_type::HostMirror host_num_k =
      Kokkos::create_mirror_view(tensor.m_num_k);
    typename k_begin_type::HostMirror host_k_begin =
      Kokkos::create_mirror_view(tensor.m_k_begin);
    typename k_size_type::HostMirror host_k_size =
      Kokkos::create_mirror_view(tensor.m_k_size);
    typename row_map_type::HostMirror host_row_map =
      Kokkos::create_mirror_view(tensor.m_row_map);
    typename num_entry_type::HostMirror host_num_entry =
      Kokkos::create_mirror_view(tensor.m_num_entry);

    // Compute row map
    size_type sum = 0;
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type num_j_parts = i_tiles[idx].parts.size();
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        size_type num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        for (size_type kdx=0; kdx<num_k_parts; ++kdx) {
          size_type nc = coord_work[idx][jdx][kdx].size();
          host_row_map(idx,jdx,kdx,0) = sum;
          for (size_type t=0; t<nc; ++t) {
            sum += coord_work[idx][jdx][kdx][t];
            host_row_map(idx,jdx,kdx,t+1) = sum;
            host_num_entry(idx,jdx,kdx,t) = 0;
          }
        }
      }
    }

    // Copy per part row offsets back into coord_work
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type num_j_parts = i_tiles[idx].parts.size();
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        size_type num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        for (size_type kdx=0; kdx<num_k_parts; ++kdx) {
          size_type nc = coord_work[idx][jdx][kdx].size();
          for (size_type t=0; t<nc; ++t) {
            coord_work[idx][jdx][kdx][t] = host_row_map(idx,jdx,kdx,t);
          }
        }
      }
    }

    // Fill in coordinate and value arrays
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      host_i_begin(idx) = i_tiles[idx].lower;
      host_i_size(idx) = i_tiles[idx].upper - i_tiles[idx].lower;
      TEUCHOS_ASSERT(host_i_size(idx) <= i_tile_size);
      size_type num_j_parts = i_tiles[idx].parts.size();
      host_num_j(idx) = num_j_parts;
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        host_j_begin(idx,jdx) = i_tiles[idx].parts[jdx].lower;
        host_j_size(idx,jdx) = i_tiles[idx].parts[jdx].upper -
          i_tiles[idx].parts[jdx].lower;
        TEUCHOS_ASSERT(host_j_size(idx,jdx) <= max_jk_tile_size);
        size_type num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        host_num_k(idx,jdx) = num_k_parts;
        for (size_type kdx=0; kdx<num_k_parts; ++kdx) {
          host_k_begin(idx,jdx,kdx) = i_tiles[idx].parts[jdx].parts[kdx].lower;
          host_k_size(idx,jdx,kdx) = i_tiles[idx].parts[jdx].parts[kdx].upper -
            i_tiles[idx].parts[jdx].parts[kdx].lower;
          TEUCHOS_ASSERT(host_k_size(idx,jdx,kdx) <= max_jk_tile_size);

          size_type nc = i_tiles[idx].parts[jdx].parts[kdx].parts.size();
          for (size_type t=0; t<nc; ++t) {
            Coord s = i_tiles[idx].parts[jdx].parts[kdx].parts[t];
            const size_type i = s.i;
            const size_type j = s.j;
            const size_type k = s.k;
            const value_type c = s.cijk;

            const size_type row = i - host_i_begin(idx);
            const size_type n = coord_work[idx][jdx][kdx][row];
            ++coord_work[idx][jdx][kdx][row];

            host_value(n) = (j != k) ? c : 0.5*c;
            host_coord2(n,0) = j - host_j_begin(idx,jdx);
            host_coord2(n,1) = k - host_k_begin(idx,jdx,kdx);
            host_coord(n) = (host_coord2(n,1) << 16) | host_coord2(n,0);

            ++host_num_entry(idx,jdx,kdx,row);
            ++tensor.m_nnz;
          }
        }
      }
    }

    // Copy data to device if necessary
    Kokkos::deep_copy(tensor.m_value, host_value);
    Kokkos::deep_copy(tensor.m_coord, host_coord);
    Kokkos::deep_copy(tensor.m_coord2, host_coord2);
    Kokkos::deep_copy(tensor.m_i_begin, host_i_begin);
    Kokkos::deep_copy(tensor.m_i_size, host_i_size);
    Kokkos::deep_copy(tensor.m_num_j, host_num_j);
    Kokkos::deep_copy(tensor.m_j_begin, host_j_begin);
    Kokkos::deep_copy(tensor.m_j_size, host_j_size);
    Kokkos::deep_copy(tensor.m_num_k, host_num_k);
    Kokkos::deep_copy(tensor.m_k_begin, host_k_begin);
    Kokkos::deep_copy(tensor.m_k_size, host_k_size);
    Kokkos::deep_copy(tensor.m_row_map, host_row_map);
    Kokkos::deep_copy(tensor.m_num_entry, host_num_entry);

    tensor.m_flops = 0;
    for (size_type idx=0; idx<num_i_parts; ++idx) {
      size_type num_j_parts = i_tiles[idx].parts.size();
      for (size_type jdx=0; jdx<num_j_parts; ++jdx) {
        size_type num_k_parts = i_tiles[idx].parts[jdx].parts.size();
        for (size_type kdx=0; kdx<num_k_parts; ++kdx) {
          for (size_type i = 0; i < host_i_size(idx); ++i) {
            tensor.m_flops += 5*host_num_entry(idx,jdx,kdx,i) + 1;
          }
        }
      }
    }

    return tensor;
  }
};

template< class Device, typename OrdinalType, typename ValueType >
SimpleTiledCrsProductTensor<ValueType, Device>
create_simple_tiled_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params)
{
  return SimpleTiledCrsProductTensor<ValueType, Device>::create(
    basis, Cijk, params);
}

template < typename ValueType, typename Device >
class BlockMultiply< SimpleTiledCrsProductTensor< ValueType , Device > >
{
public:

  typedef typename Device::size_type size_type ;
  typedef SimpleTiledCrsProductTensor< ValueType , Device > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_INLINE_FUNCTION
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type block_size = 2;
    typedef TinyVec<ValueType,block_size,false> TV;

    const size_type n_i_tile = tensor.num_i_tiles();
    for (size_type i_tile = 0; i_tile<n_i_tile; ++i_tile) {
      const size_type i_begin = tensor.i_begin(i_tile);
      const size_type i_size  = tensor.i_size(i_tile);

      const size_type n_j_tile = tensor.num_j_tiles(i_tile);
      for (size_type j_tile = 0; j_tile<n_j_tile; ++j_tile) {
        const size_type j_begin = tensor.j_begin(i_tile, j_tile);
        //const size_type j_size  = tensor.j_size(i_tile, j_tile);

        const size_type n_k_tile = tensor.num_k_tiles(i_tile, j_tile);
        for (size_type k_tile = 0; k_tile<n_k_tile; ++k_tile) {
          const size_type k_begin = tensor.k_begin(i_tile, j_tile, k_tile);
          //const size_type k_size  = tensor.k_size(i_tile, j_tile, k_tile);

          for (size_type i=0; i<i_size; ++i) {

            const size_type nEntry =
              tensor.num_entry(i_tile,j_tile,k_tile,i);
            const size_type iEntryBeg =
              tensor.entry_begin(i_tile,j_tile,k_tile,i);
            const size_type iEntryEnd = iEntryBeg + nEntry;
            size_type iEntry = iEntryBeg;

            VectorValue ytmp = 0 ;

            // Do entries with a blocked loop of size block_size
            if (block_size > 1) {
              const size_type nBlock = nEntry / block_size;
              const size_type nEntryB = nBlock * block_size;
              const size_type iEnd = iEntryBeg + nEntryB;

              TV vy;
              vy.zero();
              int j[block_size], k[block_size];

              for ( ; iEntry < iEnd ; iEntry += block_size ) {

                for (size_type ii=0; ii<block_size; ++ii) {
                  j[ii] = tensor.coord(iEntry+ii,0) + j_begin;
                  k[ii] = tensor.coord(iEntry+ii,1) + k_begin;
                }
                TV aj(a, j), ak(a, k), xj(x, j), xk(x, k),
                  c(&(tensor.value(iEntry)));

                // vy += c * ( aj * xk + ak * xj)
                aj.times_equal(xk);
                ak.times_equal(xj);
                aj.plus_equal(ak);
                c.times_equal(aj);
                vy.plus_equal(c);

              }

              ytmp += vy.sum();
            }

            // Do remaining entries with a scalar loop
            for ( ; iEntry<iEntryEnd; ++iEntry) {
              const size_type j = tensor.coord(iEntry,0) + j_begin;
              const size_type k = tensor.coord(iEntry,1) + k_begin;

              ytmp += tensor.value(iEntry) * ( a[j] * x[k] + a[k] * x[j] );
            }

            y[i+i_begin] += ytmp;

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

#endif /* #ifndef STOKHOS_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP */
