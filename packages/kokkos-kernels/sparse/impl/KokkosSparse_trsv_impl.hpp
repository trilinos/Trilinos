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

#ifndef KOKKOSSPARSE_TRSV_IMPL_HPP_
#define KOKKOSSPARSE_TRSV_IMPL_HPP_

/// \file KokkosSparse_trsv_impl.hpp
/// \brief Implementation(s) of sequential sparse triangular solve.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_Gesv.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "KokkosBlas1_set.hpp"

namespace KokkosSparse {
namespace Impl {
namespace Sequential {

template <class CrsMatrixType, class DomainMultiVectorType, class RangeMultiVectorType>
struct TrsvWrap {
  using offset_type = typename CrsMatrixType::row_map_type::non_const_value_type;
  using lno_t       = typename CrsMatrixType::index_type::non_const_value_type;
  using scalar_t    = typename CrsMatrixType::values_type::non_const_value_type;
  using device_t    = typename CrsMatrixType::device_type;
  using sview_1d    = typename Kokkos::View<scalar_t*, device_t>;
  using STS         = Kokkos::ArithTraits<scalar_t>;

  static inline void manual_copy(RangeMultiVectorType X, DomainMultiVectorType Y) {
    auto numRows = X.extent(0);
    auto numVecs = X.extent(1);
    for (decltype(numRows) i = 0; i < numRows; ++i) {
      for (decltype(numVecs) j = 0; j < numVecs; ++j) {
        X(i, j) = Y(i, j);
      }
    }
  }

  struct CommonUnblocked {
    CommonUnblocked(const lno_t block_size) {
      KK_REQUIRE_MSG(block_size == 1, "Tried to use block_size>1 for non-block-enabled Common");
    }

    scalar_t zero() { return STS::zero(); }

    template <typename ValuesView>
    scalar_t get(const ValuesView& vals, const offset_type i) {
      return vals(i);
    }

    void pluseq(scalar_t& lhs, const scalar_t& rhs) { lhs += rhs; }

    void gemv(RangeMultiVectorType X, const scalar_t& A, const lno_t r, const lno_t c, const lno_t j,
              const char = 'N') {
      X(r, j) -= A * X(c, j);
    }

    template <bool IsLower, bool Transpose = false>
    void divide(RangeMultiVectorType X, const scalar_t& A, const lno_t r, const lno_t j) {
      X(r, j) /= A;
    }
  };

  struct CommonBlocked {
    // BSR data is in LayoutRight!
    using Layout = Kokkos::LayoutRight;

    using UBlock = Kokkos::View<scalar_t**, Layout, typename CrsMatrixType::device_type,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

    using Block = Kokkos::View<scalar_t**, Layout, typename CrsMatrixType::device_type,
                               Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

    using Vector =
        Kokkos::View<scalar_t*, typename CrsMatrixType::device_type, Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

    using UVector = Kokkos::View<scalar_t*, typename CrsMatrixType::device_type,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

    lno_t m_block_size;
    lno_t m_block_items;
    Vector m_ones;
    Block m_data;
    Block m_tmp;    // Needed for SerialGesv
    UBlock m_utmp;  // Needed for SerialGesv
    Vector m_vec_data1;
    Vector m_vec_data2;

    CommonBlocked(const lno_t block_size)
        : m_block_size(block_size),
          m_block_items(block_size * block_size),
          m_ones("ones", block_size),
          m_data("m_data", block_size, block_size),
          m_tmp("m_tmp", block_size, block_size + 4),
          m_utmp(m_tmp.data(), block_size, block_size + 4),
          m_vec_data1("m_vec_data1", block_size),
          m_vec_data2("m_vec_data2", block_size) {
      Kokkos::deep_copy(m_ones, 1.0);
    }

    UBlock zero() {
      UBlock block(m_data.data(), m_block_size, m_block_size);
      KokkosBlas::SerialSet::invoke(STS::zero(), block);
      return block;
    }

    template <typename ValuesView>
    UBlock get(const ValuesView& vals, const offset_type i) {
      scalar_t* data = const_cast<scalar_t*>(vals.data());
      UBlock rv(data + (i * m_block_items), m_block_size, m_block_size);
      return rv;
    }

    void pluseq(UBlock& lhs, const UBlock& rhs) { KokkosBatched::SerialAxpy::invoke(m_ones, rhs, lhs); }

    void gemv(RangeMultiVectorType X, const UBlock& A, const lno_t r, const lno_t c, const lno_t j,
              const char transpose = 'N') {
      // Create and populate x and y
      UVector x(m_vec_data1.data(), m_block_size);
      UVector y(m_vec_data2.data(), m_block_size);
      for (lno_t b = 0; b < m_block_size; ++b) {
        x(b) = X(c * m_block_size + b, j);
        y(b) = X(r * m_block_size + b, j);
      }

      KokkosBlas::Experimental::serial_gemv(transpose, -1, A, x, 1, y);

      for (lno_t b = 0; b < m_block_size; ++b) {
        X(r * m_block_size + b, j) = y(b);
      }
    }

    template <bool IsLower, bool Transpose = false>
    void divide(RangeMultiVectorType X, const UBlock& A, const lno_t r, const lno_t j) {
      UVector x(m_vec_data1.data(), m_block_size);
      UVector y(m_vec_data2.data(), m_block_size);
      for (lno_t b = 0; b < m_block_size; ++b) {
        y(b) = X(r * m_block_size + b, j);
      }

      // if StaticPivoting is used, there are compiler errors related to
      // comparing complex and non-complex.
      using Algo = KokkosBatched::Gesv::NoPivoting;

      KokkosBatched::SerialGesv<Algo>::invoke(A, x, y, m_utmp);

      for (lno_t b = 0; b < m_block_size; ++b) {
        X(r * m_block_size + b, j) = x(b);
      }
    }
  };

  using CommonOps = std::conditional_t<KokkosSparse::Experimental::is_bsr_matrix<CrsMatrixType>::value, CommonBlocked,
                                       CommonUnblocked>;

  static void lowerTriSolveCsrUnitDiag(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    for (lno_t r = 0; r < numRows; ++r) {
      const offset_type beg = ptr(r);
      const offset_type end = ptr(r + 1);
      for (offset_type k = beg; k < end; ++k) {
        const scalar_t A_rc = val(k);
        const lno_t c       = ind(k);
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current row r
    }    // for each row r
  }

  static void lowerTriSolveCsr(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    CommonOps co(block_size);

    manual_copy(X, Y);

    for (lno_t r = 0; r < numRows; ++r) {
      auto A_rr             = co.zero();
      const offset_type beg = ptr(r);
      const offset_type end = ptr(r + 1);

      for (offset_type k = beg; k < end; ++k) {
        const auto A_rc = co.get(val, k);
        const lno_t c   = ind(k);
        // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
        // has equal local row and column indices.  That may not
        // necessarily hold, depending on the row and column Maps.  The
        // way to fix this would be for Tpetra::CrsMatrix to remember
        // the local column index of the diagonal entry (if there is
        // one) in each row, and pass that along to this function.
        if (r == c) {
          co.pluseq(A_rr, A_rc);
        } else {
          for (lno_t j = 0; j < numVecs; ++j) {
            co.gemv(X, A_rc, r, c, j);
          }
        }
      }  // for each entry A_rc in the current row r
      for (lno_t j = 0; j < numVecs; ++j) {
        co.template divide<true>(X, A_rr, r, j);
      }
    }  // for each row r
  }

  static void upperTriSolveCsrUnitDiag(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    // If lno_t is unsigned and numRows is 0, the loop
    // below will have entirely the wrong number of iterations.
    if (numRows == 0) {
      return;
    }

    // Don't use r >= 0 as the test, because that fails if
    // lno_t is unsigned.  We do r == 0 (last
    // iteration) below.
    for (lno_t r = numRows - 1; r != 0; --r) {
      const offset_type beg = ptr(r);
      const offset_type end = ptr(r + 1);
      for (offset_type k = beg; k < end; ++k) {
        const scalar_t A_rc = val(k);
        const lno_t c       = ind(k);
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current row r
    }    // for each row r

    // Last iteration: r = 0.
    {
      const lno_t r         = 0;
      const offset_type beg = ptr(r);
      const offset_type end = ptr(r + 1);
      for (offset_type k = beg; k < end; ++k) {
        const scalar_t A_rc = val(k);
        const lno_t c       = ind(k);
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current row r
    }    // last iteration: r = 0
  }

  static void upperTriSolveCsr(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    CommonOps co(block_size);

    manual_copy(X, Y);

    // If lno_t is unsigned and numRows is 0, the loop
    // below will have entirely the wrong number of iterations.
    if (numRows == 0) {
      return;
    }

    // Don't use r >= 0 as the test, because that fails if
    // lno_t is unsigned.  We do r == 0 (last
    // iteration) below.
    for (lno_t r = numRows - 1; r != 0; --r) {
      const offset_type beg = ptr(r);
      const offset_type end = ptr(r + 1);
      auto A_rr             = co.zero();
      for (offset_type k = beg; k < end; ++k) {
        const auto A_rc = co.get(val, k);
        const lno_t c   = ind(k);
        if (r == c) {
          co.pluseq(A_rr, A_rc);
        } else {
          for (lno_t j = 0; j < numVecs; ++j) {
            co.gemv(X, A_rc, r, c, j);
          }
        }
      }  // for each entry A_rc in the current row r
      for (lno_t j = 0; j < numVecs; ++j) {
        co.template divide<false>(X, A_rr, r, j);
      }
    }  // for each row r

    // Last iteration: r = 0.
    {
      const lno_t r         = 0;
      const offset_type beg = ptr(r);
      const offset_type end = ptr(r + 1);
      auto A_rr             = co.zero();
      for (offset_type k = beg; k < end; ++k) {
        const auto A_rc = co.get(val, k);
        const lno_t c   = ind(k);
        if (r == c) {
          co.pluseq(A_rr, A_rc);
        } else {
          for (lno_t j = 0; j < numVecs; ++j) {
            co.gemv(X, A_rc, r, c, j);
          }
        }
      }  // for each entry A_rc in the current row r
      for (lno_t j = 0; j < numVecs; ++j) {
        co.template divide<false>(X, A_rr, r, j);
      }
    }  // last iteration: r = 0
  }

  static void upperTriSolveCscUnitDiag(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    // If lno_t is unsigned and numCols is 0, the loop
    // below will have entirely the wrong number of iterations.
    if (numCols == 0) {
      return;
    }

    // Don't use c >= 0 as the test, because that fails if
    // lno_t is unsigned.  We do c == 0 (last
    // iteration) below.
    for (lno_t c = numCols - 1; c != 0; --c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const scalar_t A_rc = val(k);
        const lno_t r       = ind(k);
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c

    // Last iteration: c = 0.
    {
      const lno_t c         = 0;
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const scalar_t A_rc = val(k);
        const lno_t r       = ind(k);
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current column c
    }
  }

  static void upperTriSolveCsc(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    manual_copy(X, Y);

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    // If lno_t is unsigned and numCols is 0, the loop
    // below will have entirely the wrong number of iterations.
    if (numCols == 0) {
      return;
    }

    // Don't use c >= 0 as the test, because that fails if
    // lno_t is unsigned.  We do c == 0 (last
    // iteration) below.
    for (lno_t c = numCols - 1; c != 0; --c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = end - 1; k >= beg; --k) {
        const lno_t r   = ind(k);
        const auto A_rc = val(k);
        /*(vqd 20 Jul 2020) This assumes that the diagonal entry
          has equal local row and column indices.  That may not
          necessarily hold, depending on the row and column Maps.  See
          note above.*/
        for (lno_t j = 0; j < numVecs; ++j) {
          if (r == c) {
            X(c, j) = X(c, j) / A_rc;
          } else {
            X(r, j) -= A_rc * X(c, j);
          }
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c

    // Last iteration: c = 0.
    {
      const offset_type beg = ptr(0);
      const auto A_rc       = val(beg);
      /*(vqd 20 Jul 2020) This assumes that the diagonal entry
        has equal local row and column indices.  That may not
        necessarily hold, depending on the row and column Maps.  See
        note above.*/
      for (lno_t j = 0; j < numVecs; ++j) {
        X(0, j) = X(0, j) / A_rc;
      }
    }
  }

  static void lowerTriSolveCscUnitDiag(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    for (lno_t c = 0; c < numCols; ++c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = val(k);
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c
  }

  static void upperTriSolveCscUnitDiagConj(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    // If lno_t is unsigned and numCols is 0, the loop
    // below will have entirely the wrong number of iterations.
    if (numCols == 0) {
      return;
    }

    // Don't use c >= 0 as the test, because that fails if
    // lno_t is unsigned.  We do c == 0 (last
    // iteration) below.
    for (lno_t c = numCols - 1; c != 0; --c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = STS::conj(val(k));
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c

    // Last iteration: c = 0.
    {
      const lno_t c         = 0;
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = STS::conj(val(k));
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current column c
    }
  }

  static void upperTriSolveCscConj(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    // If lno_t is unsigned and numCols is 0, the loop
    // below will have entirely the wrong number of iterations.
    if (numCols == 0) {
      return;
    }

    // Don't use c >= 0 as the test, because that fails if
    // lno_t is unsigned.  We do c == 0 (last
    // iteration) below.
    for (lno_t c = numCols - 1; c != 0; --c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = end - 1; k >= beg; --k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = STS::conj(val(k));
        /*(vqd 20 Jul 2020) This assumes that the diagonal entry
          has equal local row and column indices.  That may not
          necessarily hold, depending on the row and column Maps.  See
          note above.*/
        for (lno_t j = 0; j < numVecs; ++j) {
          if (r == c) {
            X(c, j) = X(c, j) / A_rc;
          } else {
            X(r, j) -= A_rc * X(c, j);
          }
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c

    // Last iteration: c = 0.
    {
      const offset_type beg = ptr(0);
      const scalar_t A_rc   = STS::conj(val(beg));
      /*(vqd 20 Jul 2020) This assumes that the diagonal entry
        has equal local row and column indices.  That may not
        necessarily hold, depending on the row and column Maps.  See
        note above.*/
      for (lno_t j = 0; j < numVecs; ++j) {
        X(0, j) = X(0, j) / A_rc;
      }
    }
  }

  static void lowerTriSolveCsc(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    for (lno_t c = 0; c < numCols; ++c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = val(k);
        /*(vqd 20 Jul 2020) This assumes that the diagonal entry
          has equal local row and column indices.  That may not
          necessarily hold, depending on the row and column Maps.  See
          note above.*/
        for (lno_t j = 0; j < numVecs; ++j) {
          if (r == c) {
            X(c, j) = X(c, j) / A_rc;
          } else {
            X(r, j) -= A_rc * X(c, j);
          }
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c
  }

  static void lowerTriSolveCscUnitDiagConj(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    for (lno_t c = 0; c < numCols; ++c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = STS::conj(val(k));
        for (lno_t j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c
  }

  static void lowerTriSolveCscConj(RangeMultiVectorType X, const CrsMatrixType& A, DomainMultiVectorType Y) {
    const lno_t numRows = A.numRows();
    if (numRows == 0) return;

    const lno_t numCols                      = A.numCols();
    const lno_t numPointRows                 = A.numPointRows();
    const lno_t block_size                   = numPointRows / numRows;
    const lno_t numVecs                      = X.extent(1);
    typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
    typename CrsMatrixType::index_type ind   = A.graph.entries;
    typename CrsMatrixType::values_type val  = A.values;

    KK_REQUIRE_MSG(block_size == 1, "BSRs not support for this function yet");

    manual_copy(X, Y);

    for (lno_t c = 0; c < numCols; ++c) {
      const offset_type beg = ptr(c);
      const offset_type end = ptr(c + 1);
      for (offset_type k = beg; k < end; ++k) {
        const lno_t r       = ind(k);
        const scalar_t A_rc = STS::conj(val(k));
        /*(vqd 20 Jul 2020) This assumes that the diagonal entry
          has equal local row and column indices.  That may not
          necessarily hold, depending on the row and column Maps.  See
          note above.*/
        for (lno_t j = 0; j < numVecs; ++j) {
          if (r == c) {
            X(c, j) = X(c, j) / A_rc;
          } else {
            X(r, j) -= A_rc * X(c, j);
          }
        }
      }  // for each entry A_rc in the current column c
    }    // for each column c
  }
};

}  // namespace Sequential
}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_TRSV_IMPL_HPP_
