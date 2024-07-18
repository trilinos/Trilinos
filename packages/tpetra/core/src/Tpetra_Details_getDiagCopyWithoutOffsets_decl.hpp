// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_DECL_HPP
#define TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_DECL_HPP

/// \file Tpetra_Details_getDiagCopyWithoutOffsets_decl.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::getDiagCopyWithoutOffsets, and declaration
///   (only) of
///   Tpetra::Details::getDiagCopyWithoutOffsetsNotFillcomplete.
///
/// This header file, and the functions declared in it, are
/// implementation details of Tpetra::CrsMatrix.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_RowMatrix_decl.hpp"
#include "Tpetra_Vector_decl.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {

/// \brief Functor that implements much of the one-argument overload
///   of Tpetra::CrsMatrix::getLocalDiagCopy, for the case where the
///   matrix is fill complete.
///
/// \warning This is an implementation detail of
///   getDiagCopyWithoutOffsets.  See that function's documentation
///   first, before considering direct use of this functor.
///
/// \tparam DiagType 1-D nonconst Kokkos::View
/// \tparam CrsMatrixType Specialization of KokkosSparse::CrsMatrix
/// \tparam LocalMapType Specialization of Tpetra::Details::LocalMap;
///   type of the "local" part of a Tpetra::Map
template<class DiagType,
         class LocalMapType,
         class CrsMatrixType>
struct CrsMatrixGetDiagCopyFunctor {
  typedef typename LocalMapType::local_ordinal_type LO; // local ordinal type
  typedef typename LocalMapType::global_ordinal_type GO; // global ordinal type
  typedef typename CrsMatrixType::device_type device_type;
  typedef typename CrsMatrixType::value_type scalar_type;
  typedef typename CrsMatrixType::size_type offset_type;

  //! The result of the reduction; number of errors.
  typedef LO value_type;

  /// \brief Constructor
  ///
  /// \param D [out] 1-D Kokkos::View into which to store the matrix's
  ///   diagonal.
  /// \param rowMap [in] Local part of the Tpetra row Map.
  /// \param colMap [in] Local part of the Tpetra column Map.
  /// \param A [in] The sparse matrix from which to get the diagonal.
  CrsMatrixGetDiagCopyFunctor (const DiagType& D,
                               const LocalMapType& rowMap,
                               const LocalMapType& colMap,
                               const CrsMatrixType& A) :
    D_ (D), rowMap_ (rowMap), colMap_ (colMap), A_ (A)
  {}

  /// \brief Operator for Kokkos::parallel_for.
  ///
  /// \param lclRowInd [in] Index of current (local) row of the sparse matrix.
  KOKKOS_FUNCTION void
  operator () (const LO& lclRowInd, value_type& errCount) const
  {
    const LO INV = Tpetra::Details::OrdinalTraits<LO>::invalid ();
    const scalar_type ZERO =
      Kokkos::ArithTraits<scalar_type>::zero ();

    // If the row lacks a stored diagonal entry, then its value is zero.
    D_(lclRowInd) = ZERO;
    const GO gblInd = rowMap_.getGlobalElement (lclRowInd);
    const LO lclColInd = colMap_.getLocalElement (gblInd);

    if (lclColInd != INV) {
      auto curRow = A_.rowConst (lclRowInd);

      // FIXME (mfh 12 May 2016) Use binary search when the row is
      // long enough.  findRelOffset currently lives in KokkosKernels
      // (in tpetra/kernels/src/Kokkos_Sparse_findRelOffset.hpp).
      LO offset = 0;
      const LO numEnt = curRow.length;
      for ( ; offset < numEnt; ++offset) {
        if (curRow.colidx(offset) == lclColInd) {
          break;
        }
      }

      if (offset == numEnt) {
        ++errCount;
      }
      else {
        D_(lclRowInd) = curRow.value(offset);
      }
    }
  }

private:
  //! 1-D Kokkos::View into which to store the matrix's diagonal.
  DiagType D_;
  //! Local part of the Tpetra row Map.
  LocalMapType rowMap_;
  //! Local part of the Tpetra column Map.
  LocalMapType colMap_;
  //! The sparse matrix from which to get the diagonal.
  CrsMatrixType A_;
};


/// \brief Given a locally indexed, local sparse matrix, and
///   corresponding local row and column Maps, extract the matrix's
///   diagonal entries into a 1-D Kokkos::View.
///
/// \warning This is an implementation detail of Tpetra::CrsMatrix.
///   This function may disappear or change its interface at any time.
///
/// This function implements much of the one-argument overload of
/// Tpetra::CrsMatrix::getLocalDiagCopy, for the case where the matrix
/// is fill complete.  The function computes offsets of diagonal
/// entries inline, and does not store them.  If you want to store the
/// offsets, call computeOffsets() instead.
///
/// \tparam DiagType 1-D nonconst Kokkos::View
/// \tparam CrsMatrixType Specialization of KokkosSparse::CrsMatrix
/// \tparam LocalMapType Specialization of Tpetra::Details::LocalMap;
///   type of the "local" part of a Tpetra::Map
///
/// \param D [out] 1-D Kokkos::View to which to write the diagonal entries.
/// \param rowMap [in] "Local" part of the sparse matrix's row Map.
/// \param colMap [in] "Local" part of the sparse matrix's column Map.
/// \param A [in] The sparse matrix.
template<class DiagType,
         class LocalMapType,
         class CrsMatrixType>
static typename LocalMapType::local_ordinal_type
getDiagCopyWithoutOffsets (const DiagType& D,
                           const LocalMapType& rowMap,
                           const LocalMapType& colMap,
                           const CrsMatrixType& A)
{
  static_assert (Kokkos::is_view<DiagType>::value,
                 "DiagType must be a Kokkos::View.");
  static_assert (static_cast<int> (DiagType::rank) == 1,
                 "DiagType must be a 1-D Kokkos::View.");
  static_assert (std::is_same<typename DiagType::value_type, typename DiagType::non_const_value_type>::value,
                 "DiagType must be a nonconst Kokkos::View.");
  typedef typename LocalMapType::local_ordinal_type LO;
  typedef typename CrsMatrixType::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> policy_type;

  typedef Kokkos::View<typename DiagType::non_const_value_type*,
    typename DiagType::array_layout,
    typename DiagType::device_type,
    Kokkos::MemoryUnmanaged> diag_type;
  diag_type D_out = D;
  CrsMatrixGetDiagCopyFunctor<diag_type, LocalMapType, CrsMatrixType>
    functor (D_out, rowMap, colMap, A);
  const LO numRows = static_cast<LO> (D.extent (0));
  LO errCount = 0;
  Kokkos::parallel_reduce (policy_type (0, numRows), functor, errCount);
  return errCount;
}

/// \brief Given a locally indexed, global sparse matrix, extract the
///   matrix's diagonal entries into a Tpetra::Vector.
///
/// \warning This is an implementation detail of Tpetra::CrsMatrix.
///   This function may disappear or change its interface at any time.
///
/// This function is a work-around for Github Issue #499.  It
/// implements one-argument Tpetra::CrsMatrix::getLocalDiagCopy for
/// the case where the matrix is not fill complete.  The function
/// computes offsets of diagonal entries inline, and does not store
/// them.  If you want to store the offsets, call computeOffsets()
/// instead.
///
/// \tparam SC Same as first template parameter (Scalar) of
///   Tpetra::CrsMatrix and Tpetra::Vector.
/// \tparam LO Same as second template parameter (LocalOrdinal) of
///   Tpetra::CrsMatrix and Tpetra::Vector.
/// \tparam GO Same as third template parameter (GlobalOrdinal) of
///   Tpetra::CrsMatrix and Tpetra::Vector.
/// \tparam NT Same as fourth template parameter (Node) of
///   Tpetra::CrsMatrix and Tpetra::Vector.
///
/// \param diag [out] Tpetra::Vector to which to write the diagonal
///   entries.  Its Map must be the same (in the sense of
///   Tpetra::Map::isSameAs()) as the row Map of \c A.
/// \param A [in] The sparse matrix.  Must be a Tpetra::RowMatrix (the
///   base class of Tpetra::CrsMatrix), must be locally indexed, and
///   must have row views.
/// \param debug [in] Whether to do extra run-time checks.  This costs
///   MPI communication.  The default is false in a release build, and
///   true in a debug build.
///
/// We pass in the sparse matrix as a Tpetra::RowMatrix because the
/// implementation of Tpetra::CrsMatrix uses this function, and we
/// want to avoid a circular header dependency.  On the other hand,
/// the implementation does not actually depend on Tpetra::CrsMatrix.
template<class SC, class LO, class GO, class NT>
LO
getLocalDiagCopyWithoutOffsetsNotFillComplete ( ::Tpetra::Vector<SC, LO, GO, NT>& diag,
                                                const ::Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                                                const bool debug =
#ifdef HAVE_TPETRA_DEBUG
                                                true);
#else // ! HAVE_TPETRA_DEBUG
                                                false);
#endif // HAVE_TPETRA_DEBUG

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETDIAGCOPYWITHOUTOFFSETS_DECL_HPP
