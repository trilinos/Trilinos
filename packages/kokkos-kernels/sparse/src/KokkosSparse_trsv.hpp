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

/// \file KokkosSparse_trsv.hpp
/// \brief Local sparse triangular solve
///
/// This file provides KokkosSparse::trsv.  This function performs a
/// local (no MPI) sparse triangular solve on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_TRSV_HPP_
#define KOKKOSSPARSE_TRSV_HPP_

#include <type_traits>

#include "KokkosSparse_trsv_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosSparse {

/// \brief Solve the triangular sparse linear system Op(A) x = b.
///
/// \tparam AMatrix KokkosSparse::CrsMatrix specialization.
/// \tparam BMV The type of the input (right-hand side) (multi)vector.
/// \tparam XMV The type of the output (left-hand side) (multi)vector.
///
/// \param uplo [in] "U" (for upper triangular) or "L" (for lower
///   triangular).
/// \param trans [in] "C" (for conjugate transpose), "T" (for
///   transpose), or "N" (for no transpose).
/// \param diag [in] "U" (for implicit unit diagonal) or "N" (for
///   not).
/// \param A [in] The input matrix A; must be upper triangular or
///   lower triangular.
/// \param b [in] The input (right-hand side) (multi)vector.
/// \param x [in] The output (left-hand side) (multi)vector.
template <class AMatrix, class BMV, class XMV>
void trsv(const char uplo[], const char trans[], const char diag[], const AMatrix& A, const BMV& b, const XMV& x) {
  // FIXME (mfh 23 Apr 2015) Need to implement rank-1 version of this function.
  static_assert(BMV::rank == 2,
                "KokkosBlas::trsv: Rank-1 version of this "
                "function has not yet been implemented.");

  static_assert(Kokkos::is_view<BMV>::value, "KokkosBlas::trsv: b is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value, "KokkosBlas::trsv: x is not a Kokkos::View.");
  static_assert((int)BMV::rank == (int)XMV::rank, "KokkosBlas::trsv: The ranks of b and x do not match.");
  static_assert(BMV::rank == 1 || BMV::rank == 2, "KokkosBlas::trsv: b and x must both either have rank 1, or rank 2.");
  static_assert(std::is_same<typename XMV::value_type, typename XMV::non_const_value_type>::value,
                "KokkosBlas::trsv: The output x must be nonconst.");

  static_assert(
      KokkosSparse::is_crs_matrix<AMatrix>::value || KokkosSparse::Experimental::is_bsr_matrix<AMatrix>::value,
      "KokkosBlas::trsv: A is not a CRS or BSR matrix.");

  // The following three code lines have been moved up by Massimiliano Lupo
  // Pasini
  typedef typename BMV::size_type size_type;
  const size_type numRows = static_cast<size_type>(A.numPointRows());
  const size_type numCols = static_cast<size_type>(A.numPointCols());
  const size_type zero    = static_cast<size_type>(0);

  if (zero != numRows && uplo[0] != 'U' && uplo[0] != 'u' && uplo[0] != 'L' && uplo[0] != 'l') {
    std::ostringstream os;
    os << "Invalid uplo[0] = \'" << uplo << "\'";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (zero != numRows && trans[0] != 'C' && trans[0] != 'c' && trans[0] != 'T' && trans[0] != 't' && trans[0] != 'N' &&
      trans[0] != 'n') {
    std::ostringstream os;
    os << "Invalid trans[0] = \'" << trans << "\'";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (zero != numRows && diag[0] != 'U' && diag[0] != 'u' && diag[0] != 'N' && diag[0] != 'n') {
    std::ostringstream os;
    os << "Invalid diag[0] = \'" << diag << "\'";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  /*  typedef typename BMV::size_type size_type;
    const size_type numRows = static_cast<size_type> (A.numRows ());
    const size_type numCols = static_cast<size_type> (A.numCols ());*/

  const bool transpose = trans[0] != 'N' && trans[0] != 'n';
  if (!transpose && (numCols != x.extent(0) || numRows != b.extent(0))) {
    std::ostringstream os;
    os << "Dimensions do not match (non-transpose case).  "
       << "A is " << numRows << " x " << numCols << ", x is " << x.extent(0) << " x " << x.extent(1) << ", and b is "
       << b.extent(0) << " x " << b.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (transpose && (numRows != x.extent(0) || numCols != b.extent(0))) {
    std::ostringstream os;
    os << "Dimensions do not match (transpose or conjugate transpose case).  "
       << "A is " << numRows << " x " << numCols << ", x is " << x.extent(0) << " x " << x.extent(1) << ", and b is "
       << b.extent(0) << " x " << b.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using AMatrix_Bsr_Internal =
      KokkosSparse::Experimental::BsrMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                                            typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                                            typename AMatrix::const_size_type>;

  using AMatrix_Internal = std::conditional_t<
      KokkosSparse::is_crs_matrix<AMatrix>::value,
      KokkosSparse::CrsMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                              typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                              typename AMatrix::const_size_type>,
      AMatrix_Bsr_Internal>;

  AMatrix_Internal A_i(A);

  typedef Kokkos::View<typename BMV::const_value_type**, typename BMV::array_layout, typename BMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      BMV_Internal;

  typedef Kokkos::View<typename XMV::non_const_value_type**, typename XMV::array_layout, typename XMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XMV_Internal;

  BMV_Internal b_i = b;
  XMV_Internal x_i = x;

  KokkosSparse::Impl::TRSV<AMatrix_Internal, BMV_Internal, XMV_Internal>::trsv(uplo, trans, diag, A_i, b_i, x_i);
}

}  // namespace KokkosSparse

#endif  // KOKKOS_SPARSE_TRSV_HPP_
