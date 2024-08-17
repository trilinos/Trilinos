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

#ifndef _KOKKOS_SPADD_NUMERIC_IMPL_HPP
#define _KOKKOS_SPADD_NUMERIC_IMPL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_Sorting.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosSparse {
namespace Impl {

template <typename size_type, typename ordinal_type, typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
          typename AcolindsT, typename BcolindsT, typename CcolindsT, typename AvaluesT, typename BvaluesT,
          typename CvaluesT, typename AscalarT, typename BscalarT>
struct SortedNumericSumFunctor {
  using CscalarT = typename CvaluesT::non_const_value_type;

  SortedNumericSumFunctor(const ArowptrsT& Arowptrs_, const BrowptrsT& Browptrs_, const CrowptrsT& Crowptrs_,
                          const AcolindsT& Acolinds_, const BcolindsT& Bcolinds_, const CcolindsT& Ccolinds_,
                          const AvaluesT& Avalues_, const BvaluesT& Bvalues_, const CvaluesT& Cvalues_,
                          const AscalarT alpha_, const BscalarT beta_)
      : Arowptrs(Arowptrs_),
        Browptrs(Browptrs_),
        Crowptrs(Crowptrs_),
        Acolinds(Acolinds_),
        Bcolinds(Bcolinds_),
        Ccolinds(Ccolinds_),
        Avalues(Avalues_),
        Bvalues(Bvalues_),
        Cvalues(Cvalues_),
        alpha(alpha_),
        beta(beta_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
    const ordinal_type ORDINAL_MAX = Kokkos::ArithTraits<ordinal_type>::max();

    // count the union of nonzeros in Arow and Brow
    size_type ai        = 0;
    size_type bi        = 0;
    size_type Arowstart = Arowptrs(i);
    size_type Arowlen   = Arowptrs(i + 1) - Arowstart;
    size_type Browstart = Browptrs(i);
    size_type Browlen   = Browptrs(i + 1) - Browstart;
    ordinal_type Acol   = (Arowlen == 0) ? ORDINAL_MAX : Acolinds(Arowstart);
    ordinal_type Bcol   = (Browlen == 0) ? ORDINAL_MAX : Bcolinds(Browstart);
    size_type Coffset   = Crowptrs(i);
    while (Acol != ORDINAL_MAX || Bcol != ORDINAL_MAX) {
      ordinal_type Ccol = (Acol < Bcol) ? Acol : Bcol;
      // Eat all entries in both A and B which have this column
      // This also results in Acol/Bcol being updated to following entries for
      // next loop iter
      CscalarT accum = Kokkos::ArithTraits<CscalarT>::zero();
      while (Acol == Ccol) {
        accum += static_cast<CscalarT>(alpha * Avalues(Arowstart + ai));
        ai++;
        if (ai == Arowlen)
          Acol = ORDINAL_MAX;
        else
          Acol = Acolinds(Arowstart + ai);
      }
      while (Bcol == Ccol) {
        accum += static_cast<CscalarT>(beta * Bvalues(Browstart + bi));
        bi++;
        if (bi == Browlen)
          Bcol = ORDINAL_MAX;
        else
          Bcol = Bcolinds(Browstart + bi);
      }
      Ccolinds(Coffset) = Ccol;
      Cvalues(Coffset)  = accum;
      Coffset++;
    }
  }

  const ArowptrsT Arowptrs;
  const BrowptrsT Browptrs;
  const CrowptrsT Crowptrs;
  const AcolindsT Acolinds;
  const BcolindsT Bcolinds;
  CcolindsT Ccolinds;
  const AvaluesT Avalues;
  const BvaluesT Bvalues;
  CvaluesT Cvalues;
  const AscalarT alpha;
  const BscalarT beta;
};

template <typename size_type, typename ordinal_type, typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
          typename AcolindsT, typename BcolindsT, typename CcolindsT, typename AvaluesT, typename BvaluesT,
          typename CvaluesT, typename AscalarT, typename BscalarT>
struct UnsortedNumericSumFunctor {
  using CscalarT = typename CvaluesT::non_const_value_type;

  UnsortedNumericSumFunctor(const ArowptrsT Arowptrs_, const BrowptrsT Browptrs_, const CrowptrsT Crowptrs_,
                            const AcolindsT Acolinds_, const BcolindsT Bcolinds_, CcolindsT Ccolinds_,
                            const AvaluesT Avalues_, const BvaluesT Bvalues_, CvaluesT Cvalues_, const AscalarT alpha_,
                            const BscalarT beta_, const CcolindsT Apos_, const CcolindsT Bpos_)
      : Arowptrs(Arowptrs_),
        Browptrs(Browptrs_),
        Crowptrs(Crowptrs_),
        Acolinds(Acolinds_),
        Bcolinds(Bcolinds_),
        Ccolinds(Ccolinds_),
        Avalues(Avalues_),
        Bvalues(Bvalues_),
        Cvalues(Cvalues_),
        alpha(alpha_),
        beta(beta_),
        Apos(Apos_),
        Bpos(Bpos_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
    size_type CrowStart = Crowptrs(i);
    size_type CrowEnd   = Crowptrs(i + 1);
    size_type ArowStart = Arowptrs(i);
    size_type ArowEnd   = Arowptrs(i + 1);
    size_type BrowStart = Browptrs(i);
    size_type BrowEnd   = Browptrs(i + 1);
    for (size_type j = CrowStart; j < CrowEnd; j++) Cvalues(j) = Kokkos::ArithTraits<CscalarT>::zero();
    // add in A entries, while setting C colinds
    for (size_type j = ArowStart; j < ArowEnd; j++) {
      Cvalues(CrowStart + Apos(j)) += alpha * Avalues(j);
      Ccolinds(CrowStart + Apos(j)) = Acolinds(j);
    }
    // add in B entries, while setting C colinds
    for (size_type j = BrowStart; j < BrowEnd; j++) {
      Cvalues(CrowStart + Bpos(j)) += beta * Bvalues(j);
      Ccolinds(CrowStart + Bpos(j)) = Bcolinds(j);
    }
  }
  const ArowptrsT Arowptrs;
  const BrowptrsT Browptrs;
  const CrowptrsT Crowptrs;
  const AcolindsT Acolinds;
  const BcolindsT Bcolinds;
  CcolindsT Ccolinds;
  const AvaluesT Avalues;
  const BvaluesT Bvalues;
  CvaluesT Cvalues;
  const AscalarT alpha;
  const BscalarT beta;
  const CcolindsT Apos;
  const CcolindsT Bpos;
};

// Two types are the same (ignoring const)
template <typename T, typename U>
constexpr bool spadd_numeric_same_type =
    std::is_same_v<typename std::remove_const_t<T>, typename std::remove_const_t<U>>;

template <typename execution_space, typename KernelHandle, typename alno_row_view_t, typename alno_nnz_view_t,
          typename ascalar_t, typename ascalar_nnz_view_t, typename blno_row_view_t, typename blno_nnz_view_t,
          typename bscalar_t, typename bscalar_nnz_view_t, typename clno_row_view_t, typename clno_nnz_view_t,
          typename cscalar_nnz_view_t>
void spadd_numeric_impl(const execution_space& exec, KernelHandle* kernel_handle, const alno_row_view_t a_rowmap,
                        const alno_nnz_view_t a_entries, const ascalar_nnz_view_t a_values, const ascalar_t alpha,
                        const blno_row_view_t b_rowmap, const blno_nnz_view_t b_entries,
                        const bscalar_nnz_view_t b_values, const bscalar_t beta, const clno_row_view_t c_rowmap,
                        clno_nnz_view_t c_entries, cscalar_nnz_view_t c_values) {
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_type;
  // Check that A/B/C data types match KernelHandle types, and that C data types
  // are nonconst (doesn't matter if A/B types are const)
  static_assert(spadd_numeric_same_type<ascalar_t, scalar_type>, "A scalar type must match handle scalar type");
  static_assert(spadd_numeric_same_type<bscalar_t, scalar_type>, "B scalar type must match handle scalar type");
  static_assert(spadd_numeric_same_type<typename alno_row_view_t::value_type, size_type>,
                "add_symbolic: A size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(spadd_numeric_same_type<typename blno_row_view_t::value_type, size_type>,
                "add_symbolic: B size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(spadd_numeric_same_type<typename clno_row_view_t::non_const_value_type, size_type>,
                "add_symbolic: C size_type must match KernelHandle size_type)");
  static_assert(spadd_numeric_same_type<typename alno_nnz_view_t::value_type, ordinal_type>,
                "add_symbolic: A entry type must match KernelHandle entry type "
                "(aka nnz_lno_t, and const doesn't matter)");
  static_assert(spadd_numeric_same_type<typename blno_nnz_view_t::value_type, ordinal_type>,
                "add_symbolic: B entry type must match KernelHandle entry type "
                "(aka nnz_lno_t, and const doesn't matter)");
  static_assert(spadd_numeric_same_type<typename clno_nnz_view_t::value_type, ordinal_type>,
                "add_symbolic: C entry type must match KernelHandle entry type "
                "(aka nnz_lno_t)");
  static_assert(std::is_same_v<typename clno_nnz_view_t::non_const_value_type, typename clno_nnz_view_t::value_type>,
                "add_symbolic: C entry type must not be const");
  static_assert(spadd_numeric_same_type<typename ascalar_nnz_view_t::value_type, scalar_type>,
                "add_symbolic: A scalar type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");
  static_assert(spadd_numeric_same_type<typename bscalar_nnz_view_t::value_type, scalar_type>,
                "add_symbolic: B scalar type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");
  static_assert(spadd_numeric_same_type<typename cscalar_nnz_view_t::value_type, scalar_type>,
                "add_symbolic: C scalar type must match KernelHandle entry type (aka "
                "nnz_lno_t)");
  static_assert(
      std::is_same_v<typename cscalar_nnz_view_t::non_const_value_type, typename cscalar_nnz_view_t::value_type>,
      "add_symbolic: C scalar type must not be const");
  typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
  auto addHandle = kernel_handle->get_spadd_handle();
  // rowmap length can be 0 or 1 if #rows is 0.
  // Otherwise, it's always #rows+1.
  if (a_rowmap.extent(0) == 0 || a_rowmap.extent(0) == 1) {
    addHandle->set_call_numeric();
    return;
  }
  ordinal_type nrows = a_rowmap.extent(0) - 1;
  if (addHandle->is_input_sorted()) {
    SortedNumericSumFunctor<size_type, ordinal_type, alno_row_view_t, blno_row_view_t, clno_row_view_t, alno_nnz_view_t,
                            blno_nnz_view_t, clno_nnz_view_t, ascalar_nnz_view_t, bscalar_nnz_view_t,
                            cscalar_nnz_view_t, ascalar_t, bscalar_t>
        sortedNumeric(a_rowmap, b_rowmap, c_rowmap, a_entries, b_entries, c_entries, a_values, b_values, c_values,
                      alpha, beta);
    Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputSorted", range_type(exec, 0, nrows), sortedNumeric);
  } else {
    // use a_pos and b_pos (set in the handle by symbolic) to quickly compute C
    // entries and values
    UnsortedNumericSumFunctor<size_type, ordinal_type, alno_row_view_t, blno_row_view_t, clno_row_view_t,
                              alno_nnz_view_t, blno_nnz_view_t, clno_nnz_view_t, ascalar_nnz_view_t, bscalar_nnz_view_t,
                              cscalar_nnz_view_t, ascalar_t, bscalar_t>
        unsortedNumeric(a_rowmap, b_rowmap, c_rowmap, a_entries, b_entries, c_entries, a_values, b_values, c_values,
                        alpha, beta, addHandle->get_a_pos(), addHandle->get_b_pos());
    Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputNotSorted", range_type(exec, 0, nrows), unsortedNumeric);
  }
  addHandle->set_call_numeric();
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif
