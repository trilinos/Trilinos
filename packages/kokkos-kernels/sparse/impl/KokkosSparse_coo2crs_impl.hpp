//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef KOKKOSSPARSE_COO2CRS_IMPL_HPP
#define KOKKOSSPARSE_COO2CRS_IMPL_HPP

#include <Kokkos_StdAlgorithms.hpp>
#include "Kokkos_UnorderedMap.hpp"
#include "KokkosKernels_Utils.hpp"

namespace KokkosSparse {
namespace Impl {
template <class DimType, class RowViewType, class ColViewType, class DataViewType, bool InsertMode = false>
class Coo2Crs {
 private:
  using RowViewScalarType  = typename RowViewType::value_type;
  using ColViewScalarType  = typename ColViewType::value_type;
  using DataViewScalarType = typename DataViewType::value_type;
  using CrsST              = DataViewScalarType;
  using CrsOT              = RowViewScalarType;
  using CrsET              = typename DataViewType::execution_space;
  using CrsMT              = void;
  using CrsSzT             = ColViewScalarType;
  using CrsType            = CrsMatrix<CrsST, CrsOT, CrsET, CrsMT, CrsSzT>;
  using CrsValsViewType    = typename CrsType::values_type;
  using CrsRowMapViewType  = typename CrsType::row_map_type::non_const_type;
  using CrsColIdViewType   = typename CrsType::index_type;

  using UmapValueViewType = Kokkos::View<CrsST *, CrsET>;
  using UmapOpTypes       = Kokkos::UnorderedMapInsertOpTypes<UmapValueViewType, CrsOT>;
  using UmapOpType        = typename UmapOpTypes::AtomicAdd;

  // Make public for Kokkos::View
 public:
  using UmapHasherType  = typename Kokkos::pod_hash<CrsOT>;
  using UmapEqualToType = typename Kokkos::pod_equal_to<CrsOT>;
  using UmapType        = Kokkos::UnorderedMap<CrsOT, CrsST, CrsET, UmapHasherType, UmapEqualToType>;
  using UmapMemorySpace = typename UmapType::device_type::memory_space;

  // Public for kokkos policies
  struct coo2crsRp1 {};
  struct rowmapRp1 {};
  struct copyTp1 {};
  struct copyRp1 {};

  using copyTp1Pt         = Kokkos::TeamPolicy<copyTp1, CrsET>;
  using copyTp1MemberType = typename copyTp1Pt::member_type;

 private:
  using CrsRowMapView       = Kokkos::View<CrsOT *, CrsET>;
  using CrsRowMapAtomicView = Kokkos::View<CrsOT *, CrsET, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  using CrsValuesView       = Kokkos::View<CrsST *, CrsET>;
  using CrsColIdsView       = Kokkos::View<CrsSzT *, CrsET>;

  // Needed since Kokkos::Bitset cannot be accessed on the host
  using BmapViewType = Kokkos::View<bool *, CrsET, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  using Bitset       = Kokkos::Bitset<CrsET>;

  CrsRowMapView m_crs_row_map;
  CrsRowMapAtomicView m_crs_row_map_tmp;
  CrsValuesView m_crs_vals;
  CrsColIdsView m_crs_col_ids;
  UmapType *m_umaps;
  BmapViewType m_capacity_bmap;
  Bitset m_tuple_bmap;
  UmapOpType m_insert_op;
  CrsOT m_nrows;
  CrsOT m_ncols;
  RowViewType m_row;
  ColViewType m_col;
  DataViewType m_data;
  CrsSzT m_nnz;

  int m_n_tuples;

 public:
  KOKKOS_INLINE_FUNCTION
  void operator()(const coo2crsRp1 &, const int &idx) const {
    auto i           = m_row(idx);
    auto j           = m_col(idx);
    auto is_inserted = m_tuple_bmap.test(idx);

    if (i >= m_nrows || j >= m_ncols) {
      Kokkos::abort("tuple is out of bounds");
    } else if (!is_inserted && i >= 0 && j >= 0) {
      if (m_umaps[i].insert(j, m_data(idx), m_insert_op).failed()) {
        m_capacity_bmap(i) = true;  // hmap at index i reached capacity
      } else {
        m_tuple_bmap.set(idx);  // checklist of inserted tuples
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const copyRp1 &, const int &i) const {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif  // KOKKOS_ENABLE_PRAGMA_UN
    for (int j = 0; j < m_ncols; j++) {
      if (m_umaps[i].exists(j)) {
        auto umap_idx         = m_umaps[i].find(j);
        auto offset           = m_crs_row_map_tmp(i)++;
        m_crs_vals(offset)    = m_umaps[i].value_at(umap_idx);
        m_crs_col_ids(offset) = m_umaps[i].key_at(umap_idx);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const copyTp1 &, const copyTp1MemberType &member) const {
    auto row_idx = member.league_rank();
    auto cpy_beg = m_crs_row_map(row_idx);
    auto cpy_end = m_crs_row_map(row_idx + 1);
    auto cpy_len = cpy_end - cpy_beg;

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, cpy_len), [&](const CrsOT &i) {
      auto offset           = i + cpy_beg;
      m_crs_vals(offset)    = m_umaps[i].value_at(i);
      m_crs_col_ids(offset) = m_umaps[i].key_at(i);
    });
  }

  Coo2Crs(DimType m, DimType n, RowViewType row, ColViewType col, DataViewType data) {
    m_n_tuples = data.extent(0);
    m_nrows    = m;
    m_ncols    = n;
    m_row      = row;
    m_col      = col;
    m_data     = data;

    typename UmapType::size_type arg_capacity_hint = m_nrows > 0 ? (m_n_tuples / m_nrows / 4) : 16;
    typename UmapType::hasher_type arg_hasher;
    typename UmapType::equal_to_type arg_equal_to;
    arg_capacity_hint = arg_capacity_hint < 16 ? 16 : arg_capacity_hint;

    // Record of whether capacity was reached in any unordered map
    m_capacity_bmap                                          = BmapViewType("m_capacity_bmap", m_nrows);
    typename BmapViewType::HostMirror m_capacity_bmap_mirror = Kokkos::create_mirror_view(m_capacity_bmap);

    // Track which tuples have been processed
    m_tuple_bmap = Bitset(m_n_tuples);

    m_crs_row_map = CrsRowMapView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_crs_row_map"), m_nrows + 1);

    // Memory management notes for `umap_ptrs` and `m_umaps`:
    //   `umap_ptrs` is a two dimensional array. The first dimension contains
    //   pointers to mixed-memory (host and device memory). The second
    //   dimension is the array of UnorderedMap objects. Some of the object
    //   methods are callable from only the device (device-callable), others
    //   are callable from only the host. Some of the host-callable methods,
    //   such as rehash are intended to be observable on the device.
    //   See Kokkos::UnorderedMap for details.
    //
    //   `m_umaps` is a single dimension array of device memory. This array
    //   contains a shallow copy of all the UnorderedMap members that are
    //   allocated manually below.
    //
    //   Any time a host-callable method with device observable results is
    //   invoked, we must shallow-copy the given `umap_ptrs` member back to
    //   the device.
    //
    //   However, since we are using shallow copies of objects of type
    //   UnorderedMap, we do not need to copy the device memory back to
    //   the host before using a host-callable method.

    // Setup a nrows length array of Unordered Maps
    m_umaps =
        reinterpret_cast<UmapType *>(Kokkos::kokkos_malloc<UmapMemorySpace>("m_umaps", m_nrows * sizeof(UmapType)));

    auto shallow_copy_to_device = [](UmapType *dst, UmapType const *src, std::size_t cnt) {
      std::size_t nn = cnt / sizeof(UmapType);
      Kokkos::deep_copy(Kokkos::View<UmapType *, UmapMemorySpace>(dst, nn),
                        Kokkos::View<UmapType const *, Kokkos::HostSpace>(src, nn));
    };

    UmapType **umap_ptrs = new UmapType *[m_nrows];
    // TODO: use host-level parallel_for with tag rowmapRp1
    for (int i = 0; i < m_nrows; i++) {
      umap_ptrs[i] = new UmapType(arg_capacity_hint, arg_hasher, arg_equal_to);
      shallow_copy_to_device(m_umaps + i, umap_ptrs[i], sizeof(UmapType));
    }

    using coo2crsRp1Pt = Kokkos::RangePolicy<coo2crsRp1, CrsET>;
    bool rehashed      = true;
    while (rehashed) {
      Kokkos::parallel_for("coo2crsRp1", coo2crsRp1Pt(0, m_n_tuples), *this);

      CrsET().fence();  // Wait for bitmap writes to land
      Kokkos::deep_copy(m_capacity_bmap_mirror, m_capacity_bmap);
      CrsET().fence();

      rehashed = false;
      // TODO: covert to host-level parallel for.
      for (int i = 0; i < m_nrows; i++) {
        if (m_capacity_bmap_mirror(i)) {
          umap_ptrs[i]->rehash(umap_ptrs[i]->capacity() * 2);
          rehashed                  = true;
          m_capacity_bmap_mirror(i) = false;
          shallow_copy_to_device(m_umaps + i, umap_ptrs[i], sizeof(UmapType));
        }
      }
      Kokkos::deep_copy(m_capacity_bmap, m_capacity_bmap_mirror);
      CrsET().fence();
    }

    typename CrsRowMapView::HostMirror m_crs_row_map_h = Kokkos::create_mirror_view(m_crs_row_map);

    // TODO: convert to host-level parallel_for / prefix sum
    m_crs_row_map_h(0) = 0;
    for (int i = 1; i < m_nrows + 1; i++) {
      auto adj_i         = i - 1;
      auto sz            = umap_ptrs[adj_i]->size();
      m_crs_row_map_h(i) = m_crs_row_map_h(adj_i) + sz;
    }

    m_crs_row_map_tmp =
        CrsRowMapAtomicView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_crs_row_map_tmp"), m_nrows + 1);
    Kokkos::deep_copy(m_crs_row_map, m_crs_row_map_h);
    Kokkos::deep_copy(m_crs_row_map_tmp, m_crs_row_map_h);
    CrsET().fence();

    m_nnz = m_crs_row_map_h(m_nrows);

    m_crs_vals    = CrsValuesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_crs_vals"), m_nnz);
    m_crs_col_ids = CrsColIdsView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_crs_col_ids"), m_nnz);

    using copyRp1Pt = Kokkos::RangePolicy<copyRp1, CrsET>;
    Kokkos::parallel_for("copyRp1", copyRp1Pt(0, m_nrows), *this);
    CrsET().fence();

    // Cleanup
    for (int i = 0; i < m_nrows; i++) {
      delete umap_ptrs[i];
    }
    delete[] umap_ptrs;
    Kokkos::kokkos_free<UmapMemorySpace>(m_umaps);
  }

  CrsType get_crsMat() { return CrsType("coo2crs", m_nrows, m_ncols, m_nnz, m_crs_vals, m_crs_row_map, m_crs_col_ids); }
};
}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_COO2CRS_IMPL_HPP
