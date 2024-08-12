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

#ifndef KOKKOSSPARSE_IMPL_SPILUK_SYMBOLIC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_SYMBOLIC_HPP_

/// \file KokkosSparse_spiluk_symbolic_impl.hpp
/// \brief Implementation of the symbolic phase of sparse ILU(k).

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_spiluk_handle.hpp>
#include <KokkosSparse_SortCrs.hpp>
#include <KokkosKernels_Error.hpp>

// #define SYMBOLIC_OUTPUT_INFO

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

template <class IlukHandle, class RowMapType, class EntriesType, class LevelType1, class LevelType2, class LevelType3,
          class size_type>
void level_sched(IlukHandle& thandle, const RowMapType row_map, const EntriesType entries, LevelType1& level_list,
                 LevelType2& level_ptr, LevelType3& level_idx, size_type& nlevels) {
  // Scheduling currently compute on host

  using nnz_lno_t = typename IlukHandle::nnz_lno_t;

  size_type nrows = thandle.get_nrows();

  nlevels      = 0;
  level_ptr(0) = 0;

  for (size_type i = 0; i < nrows; ++i) {
    size_type l        = 0;
    size_type rowstart = row_map(i);
    size_type rowend   = row_map(i + 1);
    for (size_type j = rowstart; j < rowend; ++j) {
      nnz_lno_t col = entries(j);
      l             = std::max(l, level_list(col));
    }
    level_list(i) = l + 1;
    level_ptr(l + 1) += 1;
    nlevels = std::max(nlevels, l + 1);
  }

  for (size_type i = 1; i <= nlevels; ++i) {
    level_ptr(i) += level_ptr(i - 1);
  }

  for (size_type i = 0; i < nrows; i++) {
    level_idx(level_ptr(level_list(i) - 1)) = i;
    level_ptr(level_list(i) - 1) += 1;
  }

  if (nlevels > 0) {  // note: to avoid wrapping around to the max of size_t
                      // when nlevels = 0.
    for (size_type i = nlevels - 1; i > 0; --i) {
      level_ptr(i) = level_ptr(i - 1);
    }
  }

  level_ptr(0) = 0;

  // Find the maximum number of rows of levels
  size_type maxrows = 0;
  for (size_type i = 0; i < nlevels; ++i) {
    size_type lnrows = level_ptr(i + 1) - level_ptr(i);
    if (maxrows < lnrows) {
      maxrows = lnrows;
    }
  }

  thandle.set_num_levels(nlevels);
  thandle.set_level_maxrows(maxrows);
}

// SEQLVLSCHD_TP1 algorithm (chunks)
template <class IlukHandle, class RowMapType, class EntriesType, class LevelType1, class LevelType2, class LevelType3,
          class size_type>
void level_sched_tp(IlukHandle& thandle, const RowMapType row_map, const EntriesType entries, LevelType1& level_list,
                    LevelType2& level_ptr, LevelType3& level_idx, size_type& nlevels, int nstreams = 1) {
  // Scheduling currently compute on host

  using nnz_lno_t           = typename IlukHandle::nnz_lno_t;
  using nnz_lno_view_host_t = typename IlukHandle::nnz_lno_view_host_t;

  size_type nrows = thandle.get_nrows();

  nlevels      = 0;
  level_ptr(0) = 0;

  for (size_type i = 0; i < nrows; ++i) {
    size_type l        = 0;
    size_type rowstart = row_map(i);
    size_type rowend   = row_map(i + 1);
    for (size_type j = rowstart; j < rowend; ++j) {
      nnz_lno_t col = entries(j);
      l             = std::max(l, level_list(col));
    }
    level_list(i) = l + 1;
    level_ptr(l + 1) += 1;
    nlevels = std::max(nlevels, l + 1);
  }

  for (size_type i = 1; i <= nlevels; ++i) {
    level_ptr(i) += level_ptr(i - 1);
  }

  for (size_type i = 0; i < nrows; i++) {
    level_idx(level_ptr(level_list(i) - 1)) = i;
    level_ptr(level_list(i) - 1) += 1;
  }

  if (nlevels > 0) {  // note: to avoid wrapping around to the max of size_t
                      // when nlevels = 0.
    for (size_type i = nlevels - 1; i > 0; --i) {
      level_ptr(i) = level_ptr(i - 1);
    }
  }

  level_ptr(0) = 0;

  // Find max rows, number of chunks, max rows of chunks across levels
  thandle.alloc_level_nchunks(nlevels);
  thandle.alloc_level_nrowsperchunk(nlevels);
  nnz_lno_view_host_t lnchunks       = thandle.get_level_nchunks();
  nnz_lno_view_host_t lnrowsperchunk = thandle.get_level_nrowsperchunk();

#ifdef KOKKOS_ENABLE_CUDA
  using memory_space = typename IlukHandle::memory_space;
  size_t avail_byte  = 0;
  if (std::is_same<memory_space, Kokkos::CudaSpace>::value) {
    size_t free_byte, total_byte;
    KokkosKernels::Impl::kk_get_free_total_memory<memory_space>(free_byte, total_byte);
    avail_byte = static_cast<size_t>(0.85 * static_cast<double>(free_byte) / static_cast<double>(nstreams));
  }
#endif

  size_type maxrows         = 0;
  size_type maxrowsperchunk = 0;
  for (size_type i = 0; i < nlevels; ++i) {
    size_type lnrows = level_ptr(i + 1) - level_ptr(i);
    if (maxrows < lnrows) {
      maxrows = lnrows;
    }
#ifdef KOKKOS_ENABLE_CUDA
    size_t required_size = static_cast<size_t>(lnrows) * nrows * sizeof(nnz_lno_t);
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value) {
      lnchunks(i)       = required_size / avail_byte + 1;
      lnrowsperchunk(i) = (lnrows % lnchunks(i) == 0) ? (lnrows / lnchunks(i)) : (lnrows / lnchunks(i) + 1);
    } else
#endif
    {
      // Workaround to fix unused-parameter nstreams error
      lnchunks(i)       = static_cast<nnz_lno_t>(nstreams / nstreams);
      lnrowsperchunk(i) = lnrows;
    }
    if (maxrowsperchunk < static_cast<size_type>(lnrowsperchunk(i))) {
      maxrowsperchunk = lnrowsperchunk(i);
    }
  }

  thandle.set_num_levels(nlevels);
  thandle.set_level_maxrows(maxrows);
  thandle.set_level_maxrowsperchunk(maxrowsperchunk);
}

// Linear Search for the smallest row index
template <class size_type, class nnz_lno_t, class ViewType>
size_type search_col_index(nnz_lno_t j, size_type lenl, ViewType h_iL, ViewType h_llev, ViewType h_iw) {
  nnz_lno_t irow = h_iL(j);
  nnz_lno_t ipos = j;

  // Find the smallest col index
  for (size_type k = j + 1; k < lenl; ++k) {
    if (h_iL(k) < irow) {
      irow = h_iL(k);
      ipos = k;
    }
  }

  if (ipos != j) {  // Swap entries
    nnz_lno_t row = h_iL(j);
    h_iL(j)       = h_iL(ipos);
    h_iL(ipos)    = row;

    nnz_lno_t t  = h_llev(j);
    h_llev(j)    = h_llev(ipos);
    h_llev(ipos) = t;

    h_iw(irow) = j;
    h_iw(row)  = ipos;
  }
  return ((size_type)irow);
}

template <class IlukHandle, class ARowMapType, class AEntriesType, class LRowMapType, class LEntriesType,
          class URowMapType, class UEntriesType>
void iluk_symbolic(IlukHandle& thandle, const typename IlukHandle::const_nnz_lno_t& fill_lev,
                   const ARowMapType& A_row_map_d, const AEntriesType& A_entries_d, LRowMapType& L_row_map_d,
                   LEntriesType& L_entries_d, URowMapType& U_row_map_d, UEntriesType& U_entries_d, int nstreams = 1) {
  if (thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1)
  /*   || thandle.get_algorithm() ==
     KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHED_TP2 )*/
  {
    // Scheduling and symbolic phase currently compute on host - need host copy
    // of all views

    using size_type = typename IlukHandle::size_type;
    using nnz_lno_t = typename IlukHandle::nnz_lno_t;

    using HandleDeviceEntriesType = typename IlukHandle::nnz_lno_view_t;
    using HandleDeviceRowMapType  = typename IlukHandle::nnz_row_view_t;

    // typedef typename IlukHandle::signed_integral_t signed_integral_t;

    size_type nrows = thandle.get_nrows();

    auto A_row_map = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_row_map_d);
    auto A_entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_entries_d);
    auto L_row_map = Kokkos::create_mirror_view(Kokkos::HostSpace(), L_row_map_d);
    auto L_entries = Kokkos::create_mirror_view(Kokkos::HostSpace(), L_entries_d);
    auto U_row_map = Kokkos::create_mirror_view(Kokkos::HostSpace(), U_row_map_d);
    auto U_entries = Kokkos::create_mirror_view(Kokkos::HostSpace(), U_entries_d);

    HandleDeviceRowMapType dlevel_list = thandle.get_level_list();
    auto level_list                    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dlevel_list);

    HandleDeviceEntriesType dlevel_ptr = thandle.get_level_ptr();
    auto level_ptr                     = thandle.get_host_level_ptr();

    HandleDeviceEntriesType dlevel_idx = thandle.get_level_idx();
    auto level_idx                     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dlevel_idx);

    size_type nlev = 0;

    // Level scheduling on A???
    // level_sched<IlukHandle, AHostRowMapType, AHostEntriesType,
    // HandleHostRowMapType, HandleHostEntriesType, size_type >
    //                                    (thandle, A_row_map, A_entries, nrows,
    //                                    level_list, level_ptr, level_idx,
    //                                    nlev);
    // level_sched (thandle, A_row_map, A_entries, nrows, level_list, level_ptr,
    // level_idx, nlev);

    // Symbolic phase
    // Kokkos::resize(L_row_map_d, nrows-3);// error: static assertion failed:
    // Can only resize managed views Kokkos::resize(L_entries_d,
    // L_entries_d.extent(0)-3); thandle.set_nnzL(L_entries_d.extent(0)+5);

    using HostTmpViewType = Kokkos::View<nnz_lno_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;

    HostTmpViewType h_lev("h_lev", thandle.get_nnzU());
    HostTmpViewType h_iw("h_iw", nrows);
    HostTmpViewType h_iL("h_iL", nrows);
    HostTmpViewType h_llev("h_llev", nrows);

    size_type cntL = 0;
    size_type cntU = 0;
    size_type iU, ulev, lenu, lenl;

    L_row_map(0) = 0;
    U_row_map(0) = 0;

    Kokkos::deep_copy(h_iw, nnz_lno_t(-1));

    // Main loop
    for (size_type i = 0; i < nrows; ++i) {
      iU   = i;
      ulev = i;
      lenl = lenu = 0;

      // Unpack the ith row
      size_type k1 = A_row_map(i);
      size_type k2 = A_row_map(i + 1);

      for (size_type k = k1; k < k2; ++k) {
        size_type col = static_cast<size_type>(A_entries(k));
        if (col < nrows) {  // Ignore column elements that are not in the square
                            // matrix
          if (col > i) {    // U part
            h_iw(col)           = lenu;
            h_iL(iU + lenu)     = col;
            h_llev(ulev + lenu) = 0;
            lenu++;
          } else if (col < i) {  // L part
            h_iw(col)    = lenl;
            h_iL(lenl)   = col;
            h_llev(lenl) = 0;
            lenl++;
          }
        }
      }

      // Eliminate rows
      nnz_lno_t j = -1;
      while (static_cast<size_type>(++j) < lenl) {
        size_type row  = search_col_index(j, lenl, h_iL, h_llev, h_iw);
        nnz_lno_t jlev = h_llev(j);
        k1             = U_row_map(row) + 1;
        k2             = U_row_map(row + 1);
        for (size_type k = k1; k < k2; ++k) {
          size_type col  = static_cast<size_type>(U_entries(k));
          nnz_lno_t lev1 = jlev + h_lev(k) + 1;
          if (lev1 > fill_lev) continue;
          nnz_lno_t ipos = h_iw(col);
          if (ipos == -1) {  // Fill-in
            if (col > i) {   // U part
              h_iw(col)           = lenu;
              h_iL(iU + lenu)     = col;
              h_llev(ulev + lenu) = lev1;
              lenu++;
            } else if (col < i) {  // L part
              h_iw(col)    = lenl;
              h_iL(lenl)   = col;
              h_llev(lenl) = lev1;
              lenl++;
            }
          } else {  // Not a fill-in
            if (col > i)
              h_llev(ulev + ipos) = std::min(h_llev(ulev + ipos), lev1);
            else if (col < i)
              h_llev(ipos) = std::min(h_llev(ipos), lev1);
          }
        }
      }

      // Reset iw
      for (size_type k = 0; k < lenl; ++k) h_iw(h_iL(k)) = -1;
      for (size_type k = 0; k < lenu; ++k) h_iw(h_iL(iU + k)) = -1;

      // Copy U part+diag and levels
      if (cntU + lenu + 1 > static_cast<size_type>(U_entries_d.extent(0))) {
        // size_type newsize = (size_type)(U_entries_d.extent(0)*EXPAND_FACT);
        // Kokkos::resize(h_lev, newsize);
        // Kokkos::resize(U_entries, newsize);
        // Kokkos::resize(U_entries_d, newsize);
        std::ostringstream os;
        os << "KokkosSparse::Experimental::spiluk_symbolic: U_entries's extent "
              "must be larger than "
           << U_entries_d.extent(0) << ", must be at least " << cntU + lenu + 1;
        KokkosKernels::Impl::throw_runtime_exception(os.str());
      }
      // U diag entry
      U_entries(cntU) = i;
      cntU++;
      // U part
      for (size_type k = 0; k < lenu; ++k) {
        U_entries(cntU) = h_iL(iU + k);
        h_lev(cntU)     = h_llev(ulev + k);
        cntU++;
      }
      U_row_map(i + 1) = cntU;

      // Copy L part
      if (cntL + lenl + 1 > static_cast<size_type>(L_entries_d.extent(0))) {
        // size_type newsize = (size_type) (L_entries_d.extent(0)*EXPAND_FACT);
        // Kokkos::resize(L_entries, newsize);
        // Kokkos::resize(L_entries_d, newsize);
        std::ostringstream os;
        os << "KokkosSparse::Experimental::spiluk_symbolic: L_entries's extent "
              "must be larger than "
           << L_entries_d.extent(0) << ", must be at least " << cntL + lenl + 1;
        KokkosKernels::Impl::throw_runtime_exception(os.str());
      }
      for (size_type k = 0; k < lenl; ++k) {
        L_entries(cntL) = h_iL(k);
        cntL++;
      }
      // L diag entry
      L_entries(cntL) = i;
      cntL++;
      L_row_map(i + 1) = cntL;
    }  // End main loop i

    thandle.set_nnzL(cntL);
    thandle.set_nnzU(cntU);

    // Sort
    KokkosSparse::sort_crs_graph<Kokkos::DefaultHostExecutionSpace, decltype(L_row_map), decltype(L_entries)>(
        L_row_map, L_entries);
    KokkosSparse::sort_crs_graph<Kokkos::DefaultHostExecutionSpace, decltype(U_row_map), decltype(U_entries)>(
        U_row_map, U_entries);

    // Level scheduling on L
    if (thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1) {
      level_sched_tp(thandle, L_row_map, L_entries, level_list, level_ptr, level_idx, nlev, nstreams);
      thandle.alloc_iw(thandle.get_level_maxrowsperchunk(), nrows);
    } else {
      level_sched(thandle, L_row_map, L_entries, level_list, level_ptr, level_idx, nlev);
      thandle.alloc_iw(thandle.get_level_maxrows(), nrows);
    }

    Kokkos::deep_copy(dlevel_ptr, level_ptr);
    Kokkos::deep_copy(dlevel_idx, level_idx);
    Kokkos::deep_copy(dlevel_list, level_list);

    Kokkos::deep_copy(L_row_map_d, L_row_map);
    Kokkos::deep_copy(L_entries_d, L_entries);
    Kokkos::deep_copy(U_row_map_d, U_row_map);
    Kokkos::deep_copy(U_entries_d, U_entries);

    thandle.set_symbolic_complete();

    // Output check
#ifdef SYMBOLIC_OUTPUT_INFO
    std::cout << "  ILU(k) fill_level: " << fill_lev << std::endl;
    std::cout << "  symbolic complete: " << thandle.is_symbolic_complete() << std::endl;
    std::cout << "  num levels: " << thandle.get_num_levels() << std::endl;
    std::cout << "  max num rows among levels: " << thandle.get_level_maxrows() << std::endl;
    std::cout << "  max num rows among chunks among levels: " << thandle.get_level_maxrowsperchunk() << std::endl;

    std::cout << "  iluk_symbolic result: " << std::endl;

    std::cout << "  level_list = ";
    for (size_type i = 0; i < nrows; ++i) {
      std::cout << level_list(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  level_ptr = ";
    for (size_type i = 0; i < nlev + 1; ++i) {
      std::cout << level_ptr(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  level_idx = ";
    for (size_type i = 0; i < nrows; ++i) {
      std::cout << level_idx(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  nnzL: " << thandle.get_nnzL() << std::endl;
    std::cout << "  L_row_map = ";
    for (size_type i = 0; i < nrows + 1; ++i) {
      std::cout << L_row_map(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  L_entries = ";
    for (size_type i = 0; i < thandle.get_nnzL(); ++i) {
      std::cout << L_entries(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  nnzU: " << thandle.get_nnzU() << std::endl;
    std::cout << "  U_row_map = ";
    for (size_type i = 0; i < nrows + 1; ++i) {
      std::cout << U_row_map(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  U_entries = ";
    for (size_type i = 0; i < thandle.get_nnzU(); ++i) {
      std::cout << U_entries(i) << " ";
    }
    std::cout << std::endl;
#endif
  }
}  // end iluk_symbolic

}  // namespace Experimental
}  // namespace Impl
}  // namespace KokkosSparse

#endif
