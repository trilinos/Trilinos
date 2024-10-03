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

#ifndef KOKKOS_SPARSE_IMPL_GETDIAGCOPYWITHOFFSETS_HPP_
#define KOKKOS_SPARSE_IMPL_GETDIAGCOPYWITHOFFSETS_HPP_

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_OrdinalTraits.hpp"
#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosSparse {
namespace Impl {

/// \brief Functor that implements much of the two-argument overload
///   of Tpetra::CrsMatrix::getLocalDiagCopy, for the case where the
///   matrix is fill complete.
///
/// The two-argument version of Tpetra::CrsMatrix::getLocalDiagCopy is
/// optimized by already having the offsets.  The getLocalDiagOffsets
/// method computes the offsets.
///
/// \tparam CrsMatrixType Specialization of KokkosSparse::CrsMatrix
/// \tparam OffsetType The type of offsets (the entries of \c offsets)
template <class DiagType, class OffsetsType, class CrsMatrixType>
struct CrsMatrixGetDiagCopyWithOffsetsFunctor {
  typedef typename CrsMatrixType::ordinal_type LO;  // local ordinal type
  typedef typename CrsMatrixType::device_type device_type;
  typedef typename CrsMatrixType::value_type scalar_type;
  typedef typename OffsetsType::non_const_value_type offset_type;

  /// \brief Constructor
  ///
  /// \param D [out] 1-D Kokkos::View into which to store the matrix's
  ///   diagonal.
  /// \param offsets [in] Offsets, precomputed using
  ///   Tpetra::CrsMatrix::getLocalDiagOffsets.
  /// \param A [in] The sparse matrix from which to get the diagonal.
  CrsMatrixGetDiagCopyWithOffsetsFunctor(const DiagType& D, const OffsetsType& offsets, const CrsMatrixType& A)
      : D_(D), offsets_(offsets), A_(A) {
    static_assert(Kokkos::is_view<DiagType>::value, "The DiagType template parameter must be a Kokkos::View.");
    static_assert(static_cast<int>(DiagType::rank) == 1, "The DiagType template parameter must be a 1-D Kokkos::View.");
    static_assert(std::is_same<typename DiagType::value_type, typename DiagType::non_const_value_type>::value,
                  "The DiagType template parameter must be a nonconst Kokkos::View.");
    static_assert(Kokkos::is_view<OffsetsType>::value, "The OffsetsType template parameter must be a Kokkos::View.");
    static_assert(static_cast<int>(OffsetsType::rank) == 1,
                  "The OffsetsType template parameter must be a 1-D Kokkos::View.");
  }

  /// \brief Operator for Kokkos::parallel_for.
  ///
  /// \param lclRow [in] The current (local) row of the sparse matrix.
  KOKKOS_INLINE_FUNCTION void operator()(const LO& lclRow) const {
    const offset_type INV  = KokkosSparse::OrdinalTraits<offset_type>::invalid();
    const scalar_type ZERO = Kokkos::ArithTraits<scalar_type>::zero();

    // If the row lacks a stored diagonal entry, then its value is zero.
    D_(lclRow)               = ZERO;
    const offset_type offset = offsets_(lclRow);
    if (offset != INV) {
      auto curRow = A_.rowConst(lclRow);
      D_(lclRow)  = curRow.value(offset);
    }
  }

 private:
  //! 1-D Kokkos::View into which to store the matrix's diagonal.
  DiagType D_;
  //! Offsets, precomputed using Tpetra::CrsMatrix::getLocalDiagOffsets.
  OffsetsType offsets_;
  //! The sparse matrix from which to get the diagonal.
  CrsMatrixType A_;
};

template <class ScalarType, class OrdinalType, class DeviceType, class OffsetType>
struct CrsMatrixGetDiagCopyWithOffsets {
  typedef ScalarType scalar_type;
  typedef OrdinalType ordinal_type;
  typedef DeviceType device_type;
  typedef OffsetType offset_type;
  typedef ::KokkosSparse::CrsMatrix<scalar_type, ordinal_type, device_type, void, offset_type> crs_matrix_type;
  typedef Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type, Kokkos::MemoryUnmanaged> diag_type;
  typedef Kokkos::View<const size_t*, device_type, Kokkos::MemoryUnmanaged> offsets_type;
  static void getDiagCopy(const diag_type& D, const offsets_type& offsets, const crs_matrix_type& A) {
    typedef typename device_type::execution_space execution_space;
    const ordinal_type numRows = static_cast<ordinal_type>(D.extent(0));
    CrsMatrixGetDiagCopyWithOffsetsFunctor<diag_type, offsets_type, crs_matrix_type> functor(D, offsets, A);
    typedef Kokkos::RangePolicy<execution_space, ordinal_type> policy_type;
    Kokkos::parallel_for(policy_type(0, numRows), functor);
  }
};

//
// Macro for declarations of full specialization of
// KokkosSparse::Impl::CrsMatrixGetDiagCopyWithOffsets.  This is NOT for
// users!!!
//
#define KOKKOSSPARSE_IMPL_GETDIAGCOPYWITHOFFSETS_DECL(SCALAR, ORDINAL, EXEC_SPACE, MEM_SPACE, OFFSET)            \
  extern template struct CrsMatrixGetDiagCopyWithOffsets<SCALAR, ORDINAL, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                                         OFFSET>;

//
// Macro for definitions of full specialization of
// KokkosSparse::Impl::CrsMatrixGetDiagCopyWithOffsets.  This is NOT for
// users!!!
//
#define KOKKOSSPARSE_IMPL_GETDIAGCOPYWITHOFFSETS_DEF(SCALAR, ORDINAL, EXEC_SPACE, MEM_SPACE, OFFSET) \
  template struct CrsMatrixGetDiagCopyWithOffsets<SCALAR, ORDINAL, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, OFFSET>;

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOS_SPARSE_IMPL_GETDIAGCOPYWITHOFFSETS_HPP_
