/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_UnitTestEquilibration.cpp
/// \brief Unit test for Ifpack2's equilibration

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_copyConvert.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <memory>
#include <type_traits>

namespace { // (anonymous)

/// \brief Struct storing results of computeRowAndColumnOneNorms.
///
/// Compute this as the return value of computeRowAndColumnOneNorms.
/// Use it either as input to leftAndOrRightScaleCrsMatrix (e.g., for
/// equilibration, or balancing in the symmetric / Hermitian positive
/// definite case), or to analyze the matrix (e.g., whether it is
/// diagonally dominant).
///
/// If this is the return value of computeRowAndColumnOneNorms,
/// results are always global.  computeLocalRowAndColumnOneNorms only
/// computes local results, where "local" means "to an MPI process."
/// globalizeRowOneNorms and globalizeColumnOneNorms "globalize" the
/// results, so they refer to the entire matrix, globally distributed
/// over an MPI communicator.
template<class ScalarType, class DeviceType>
struct EquilibrationInfo {
  using val_type = typename Kokkos::ArithTraits<ScalarType>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using device_type = typename DeviceType::device_type;
  using host_device_type = typename Kokkos::View<mag_type*, device_type>::HostMirror::device_type;
  using HostMirror = EquilibrationInfo<val_type, host_device_type>;

  EquilibrationInfo (const std::size_t lclNumRows,
                     const std::size_t lclNumCols,
                     const bool assumeSymmetric_) :
    rowNorms (Kokkos::View<mag_type*, device_type> ("rowNorms", lclNumRows)),
    rowDiagonalEntries (Kokkos::View<val_type*, device_type> ("rowDiagonalEntries", lclNumRows)),
    colNorms (Kokkos::View<mag_type*, device_type> ("colNorms", lclNumCols)),
    colDiagonalEntries (Kokkos::View<val_type*, device_type> ("colDiagonalEntries",
                                                              assumeSymmetric_ ?
                                                              std::size_t (0) :
                                                              lclNumCols)),
    rowScaledColNorms (Kokkos::View<mag_type*, device_type> ("rowScaledColNorms",
                                                             assumeSymmetric_ ?
                                                             std::size_t (0) :
                                                             lclNumCols)),
    assumeSymmetric (assumeSymmetric_)
  {}

  EquilibrationInfo (const Kokkos::View<mag_type*, device_type>& rowNorms_,
                     const Kokkos::View<val_type*, device_type>& rowDiagonalEntries_,
                     const Kokkos::View<mag_type*, device_type>& colNorms_,
                     const Kokkos::View<val_type*, device_type>& colDiagonalEntries_,
                     const Kokkos::View<mag_type*, device_type>& rowScaledColNorms_,
                     const bool assumeSymmetric_) :
    rowNorms (rowNorms_),
    rowDiagonalEntries (rowDiagonalEntries_),
    colNorms (colNorms_),
    colDiagonalEntries (colDiagonalEntries_),
    rowScaledColNorms (rowScaledColNorms_),
    assumeSymmetric (assumeSymmetric_)
  {}

  //! Deep-copy src into *this.
  template<class SrcDeviceType>
  void
  assign (const EquilibrationInfo<ScalarType, SrcDeviceType>& src)
  {
    Kokkos::deep_copy (rowNorms, src.rowNorms);
    Kokkos::deep_copy (rowDiagonalEntries, src.rowDiagonalEntries);
    Kokkos::deep_copy (colNorms, src.colNorms);
    if (src.colDiagonalEntries.extent (0) == 0) {
      colDiagonalEntries =
        Kokkos::View<val_type*, device_type> ("colDiagonalEntries", 0);
    }
    else {
      Kokkos::deep_copy (colDiagonalEntries, src.colDiagonalEntries);
    }
    if (src.rowScaledColNorms.extent (0) == 0) {
      rowScaledColNorms =
        Kokkos::View<mag_type*, device_type> ("rowScaledColNorms", 0);
    }
    else {
      Kokkos::deep_copy (rowScaledColNorms, src.rowScaledColNorms);
    }
  }

  typename EquilibrationInfo<val_type, device_type>::HostMirror
  createMirrorView ()
  {
    auto rowNorms_h = Kokkos::create_mirror_view (rowNorms);
    auto rowDiagonalEntries_h = Kokkos::create_mirror_view (rowDiagonalEntries);
    auto colNorms_h = Kokkos::create_mirror_view (colNorms);
    auto colDiagonalEntries_h = Kokkos::create_mirror_view (colDiagonalEntries);
    auto rowScaledColNorms_h = Kokkos::create_mirror_view (rowScaledColNorms);

    return HostMirror {rowNorms_h, rowDiagonalEntries_h, colNorms_h,
        colDiagonalEntries_h, rowScaledColNorms_h, assumeSymmetric};
  }

  // We call a row a "diagonally dominant row" if the absolute value
  // of the diagonal entry is >= the sum of the absolute values of the
  // off-diagonal entries.  The row norm is the sum of those two
  // things, so this means diagAbsVal >= rowNorm - diagAbsVal.  Ditto
  // for a column.

  //! One-norms of the matrix's rows, distributed via the row Map.
  Kokkos::View<mag_type*, device_type> rowNorms;

  //! Diagonal entries of the matrix, distributed via the row Map.
  Kokkos::View<val_type*, device_type> rowDiagonalEntries;

  /// \brief One-norms of the matrix's columns, distributed via the
  ///   column Map.
  ///
  /// If assumeSymmetric is true, this is just a redistributed version
  /// of rowNorms.
  Kokkos::View<mag_type*, device_type> colNorms;

  /// \brief Diagonal entries of the matrix, distributed via the column Map.
  ///
  /// Only use this if assumeSymmetric is false.
  Kokkos::View<val_type*, device_type> colDiagonalEntries;

  /// \brief One-norms of the matrix's columns, after the matrix's
  ///   rows have been scaled by rowNorms.
  ///
  /// This is only valid if assumeSymmetric is false.
  ///
  /// For the nonsymmetric case, we imitate LAPACK's DGEEQU, in doing
  /// the row scaling first, then the column scaling.  Thus, the
  /// column norms are "scaled" by the row norms (above).  We still
  /// keep the unscaled column norms (colNorms; see above) in that
  /// case, because they are diagnostic.
  Kokkos::View<mag_type*, device_type> rowScaledColNorms;

  /// \brief Whether to assume that the matrix is (globally) symmetric.
  ///
  /// This affects whether colDiagonalEntries and rowScaledColNorms
  /// are valid.
  bool assumeSymmetric;
};

template<class SC, class LO, class GO, class NT>
std::size_t
lclMaxNumEntriesRowMatrix (const Tpetra::RowMatrix<SC, LO, GO, NT>& A)
{
  const auto& rowMap = * (A.getRowMap ());
  const LO lclNumRows = static_cast<LO> (rowMap.getNodeNumElements ());

  std::size_t maxNumEnt {0};
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const std::size_t numEnt = A.getNumEntriesInLocalRow (lclRow);
    maxNumEnt = numEnt > maxNumEnt ? numEnt : maxNumEnt;
  }
  return maxNumEnt;
}

template<class SC, class LO, class GO, class NT>
void
forEachLocalRowMatrixRow (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                          const LO lclNumRows,
                          const std::size_t maxNumEnt,
                          std::function<void (const LO lclRow,
                                              const Teuchos::ArrayView<LO>& /* ind */,
                                              const Teuchos::ArrayView<SC>& /* val */,
                                              std::size_t /* numEnt */ )> doForEachRow)
{
  Teuchos::Array<LO> indBuf (maxNumEnt);
  Teuchos::Array<SC> valBuf (maxNumEnt);

  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    std::size_t numEnt = A.getNumEntriesInLocalRow (lclRow);
    Teuchos::ArrayView<LO> ind = indBuf.view (0, numEnt);
    Teuchos::ArrayView<SC> val = valBuf.view (0, numEnt);
    A.getLocalRowCopy (lclRow, ind, val, numEnt);
    doForEachRow (lclRow, ind, val, numEnt);
  }
}

template<class SC, class LO, class GO, class NT>
void
forEachLocalRowMatrixRow (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                          std::function<void (const LO lclRow,
                                              const Teuchos::ArrayView<LO>& /* ind */,
                                              const Teuchos::ArrayView<SC>& /* val */,
                                              std::size_t /* numEnt */ )> doForEachRow)
{
  const auto& rowMap = * (A.getRowMap ());
  const LO lclNumRows = static_cast<LO> (rowMap.getNodeNumElements ());
  const std::size_t maxNumEnt = lclMaxNumEntriesRowMatrix (A);

  forEachLocalRowMatrixRow<SC, LO, GO, NT> (A, lclNumRows, maxNumEnt, doForEachRow);
}

/// \brief For a given Tpetra::RowMatrix that is not a
///   Tpetra::CrsMatrix, assume that result.rowNorms has been computed
///   (and globalized), and compute result.rowScaledColNorms.
template<class SC, class LO, class GO, class NT>
void
computeLocalRowScaledColumnNorms_RowMatrix (EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                                                              typename NT::device_type>& result,
                                            const Tpetra::RowMatrix<SC, LO, GO, NT>& A)
{
  using KAT = Kokkos::ArithTraits<SC>;
  using mag_type = typename KAT::mag_type;

  auto rowNorms_h = Kokkos::create_mirror_view (result.rowNorms);
  Kokkos::deep_copy (rowNorms_h, result.rowNorms);
  auto rowScaledColNorms_h = Kokkos::create_mirror_view (result.rowScaledColNorms);

  forEachLocalRowMatrixRow<SC, LO, GO, NT> (A,
    [&] (const LO lclRow,
         const Teuchos::ArrayView<LO>& ind,
         const Teuchos::ArrayView<SC>& val,
         std::size_t numEnt) {
      const mag_type rowNorm = rowNorms_h[lclRow];
      for (std::size_t k = 0; k < numEnt; ++k) {
        const mag_type matrixAbsVal = KAT::abs (val[k]);
        const LO lclCol = ind[k];

        rowScaledColNorms_h[lclCol] += matrixAbsVal / rowNorm;
      }
    });
  Kokkos::deep_copy (result.rowScaledColNorms, rowScaledColNorms_h);
}

/// \brief Implementation of computeLocalRowAndColumnOneNorms for a
///   Tpetra::RowMatrix that is NOT a Tpetra::CrsMatrix.
template<class SC, class LO, class GO, class NT>
EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type, typename NT::device_type>
computeLocalRowAndColumnOneNorms_RowMatrix (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                                            const bool assumeSymmetric)
{
  using KAT = Kokkos::ArithTraits<SC>;
  using val_type = typename KAT::val_type;
  using mag_type = typename KAT::mag_type;
  using device_type = typename NT::device_type;

  const auto& rowMap = * (A.getRowMap ());
  const auto& colMap = * (A.getColMap ());
  const LO lclNumRows = static_cast<LO> (rowMap.getNodeNumElements ());
  const LO lclNumCols = static_cast<LO> (colMap.getNodeNumElements ());

  EquilibrationInfo<val_type, device_type> result
    (lclNumRows, lclNumCols, assumeSymmetric);
  auto result_h = result.createMirrorView ();

  forEachLocalRowMatrixRow<SC, LO, GO, NT> (A,
    [&] (const LO lclRow,
         const Teuchos::ArrayView<LO>& ind,
         const Teuchos::ArrayView<SC>& val,
         std::size_t numEnt) {
      mag_type rowNorm {0.0};
      val_type diagVal {0.0};
      const GO gblRow = rowMap.getGlobalElement (lclRow);
      // OK if invalid(); then we simply won't find the diagonal entry.
      const GO lclDiagColInd = colMap.getLocalElement (gblRow);
      for (std::size_t k = 0; k < numEnt; ++k) {
        const mag_type matrixAbsVal = KAT::abs (val[k]);
        rowNorm += matrixAbsVal;
        const LO lclCol = ind[k];
        if (lclCol == lclDiagColInd) {
          diagVal += val[k]; // repeats count additively
        }
        if (! assumeSymmetric) {
          result_h.colNorms[lclCol] += matrixAbsVal;
        }
      }

      // NOTE (mfh 24 May 2018) We could actually compute local
      // rowScaledColNorms in situ at this point, if ! assumeSymmetric
      // and row Map is the same as range Map (so that the local row
      // norms are the same as the global row norms).

      result_h.rowDiagonalEntries[lclRow] += diagVal;
      result_h.rowNorms[lclRow] = rowNorm;
      if (! assumeSymmetric &&
          lclDiagColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
        result_h.colDiagonalEntries[lclDiagColInd] += diagVal;
      }
    });

  result.assign (result_h);
  return result;
}


template<class SC, class LO, class GO, class NT>
class ComputeLocalRowScaledColumnNorms {
public:
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using device_type = typename crs_matrix_type::device_type;

  ComputeLocalRowScaledColumnNorms (const Kokkos::View<mag_type*, device_type>& rowScaledColNorms,
                                    const Kokkos::View<const mag_type*, device_type>& rowNorms,
                                    const crs_matrix_type& A) :
    rowScaledColNorms_ (rowScaledColNorms),
    rowNorms_ (rowNorms),
    A_lcl_ (A.getLocalMatrix ())
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const LO lclRow) const {
    using KAT = Kokkos::ArithTraits<val_type>;

    const auto curRow = A_lcl_.rowConst (lclRow);
    const mag_type rowNorm = rowNorms_[lclRow];
    const LO numEnt = curRow.length;
    for (LO k = 0; k < numEnt; ++k) {
      const mag_type matrixAbsVal = KAT::abs (curRow.value(k));
      const LO lclCol = curRow.colidx(k);

      Kokkos::atomic_add (&rowScaledColNorms_[lclCol], matrixAbsVal / rowNorm);
    }
  }

  static void
  run (const Kokkos::View<mag_type*, device_type>& rowScaledColNorms,
       const Kokkos::View<const mag_type*, device_type>& rowNorms,
       const crs_matrix_type& A)
  {
    using execution_space = typename device_type::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    using functor_type = ComputeLocalRowScaledColumnNorms<SC, LO, GO, NT>;

    functor_type functor (rowScaledColNorms, rowNorms, A);
    const LO lclNumRows =
      static_cast<LO> (A.getRowMap ()->getNodeNumElements ());
    Kokkos::parallel_for ("computeLocalRowScaledColumnNorms",
                          range_type (0, lclNumRows), functor);
  }

private:
  Kokkos::View<mag_type*, device_type> rowScaledColNorms_;
  Kokkos::View<const mag_type*, device_type> rowNorms_;

  using local_matrix_type = typename crs_matrix_type::local_matrix_type;
  local_matrix_type A_lcl_;
};

template<class SC, class LO, class GO, class NT>
void
computeLocalRowScaledColumnNorms_CrsMatrix (EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                                                              typename NT::device_type>& result,
                                            const Tpetra::CrsMatrix<SC, LO, GO, NT>& A)
{
  using impl_type = ComputeLocalRowScaledColumnNorms<SC, LO, GO, NT>;
  impl_type::run (result.rowScaledColNorms, result.rowNorms, A);
}

template<class SC, class LO, class GO, class NT>
void
computeLocalRowScaledColumnNorms (EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                                                    typename NT::device_type>& result,
                                  const Tpetra::RowMatrix<SC, LO, GO, NT>& A)
{
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using device_type = typename NT::device_type;

  auto colMapPtr = A.getColMap ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (colMapPtr.get () == nullptr, std::invalid_argument,
     "computeLocalRowScaledColumnNorms: "
     "Input matrix A must have a nonnull column Map.");
  const LO lclNumCols = static_cast<LO> (colMapPtr->getNodeNumElements ());
  if (static_cast<std::size_t> (result.rowScaledColNorms.extent (0)) !=
      static_cast<std::size_t> (lclNumCols)) {
    result.rowScaledColNorms =
      Kokkos::View<mag_type*, device_type> ("rowScaledColNorms", lclNumCols);
  }

  const crs_matrix_type* A_crs = dynamic_cast<const crs_matrix_type*> (&A);
  if (A_crs == nullptr) {
    computeLocalRowScaledColumnNorms_RowMatrix (result, A);
  }
  else {
    computeLocalRowScaledColumnNorms_CrsMatrix (result, *A_crs);
  }
}


template<class SC, class LO, class GO, class NT>
class ComputeLocalRowAndColumnOneNorms {
public:
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using result_type = EquilibrationInfo<val_type, typename NT::device_type>;

public:
  ComputeLocalRowAndColumnOneNorms (const crs_matrix_type& A,
                                    const bool assumeSymmetric) :
    A_lcl_ (A.getLocalMatrix ()),
    rowMap_ (A.getRowMap ()->getLocalMap ()),
    colMap_ (A.getColMap ()->getLocalMap ()),
    result_ (result_type (static_cast<LO> (rowMap_.getNodeNumElements ()),
                          static_cast<LO> (colMap_.getNodeNumElements ()),
                          assumeSymmetric))
  {}

  result_type getResult () const
  {
    return result_;
  }

  KOKKOS_INLINE_FUNCTION void operator () (const LO lclRow) const
  {
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;

    const GO gblRow = rowMap_.getGlobalElement (lclRow);
    // OK if invalid(); then we simply won't find the diagonal entry.
    const GO lclDiagColInd = colMap_.getLocalElement (gblRow);

    const auto curRow = A_lcl_.rowConst (lclRow);
    const LO numEnt = curRow.length;
    const bool assumeSymmetric = result_.assumeSymmetric;

    mag_type rowNorm {0.0};
    val_type diagVal {0.0};
    for (LO k = 0; k < numEnt; ++k) {
      const mag_type matrixAbsVal = KAT::abs (curRow.value (k));
      rowNorm += matrixAbsVal;
      const LO lclCol = curRow.colidx (k);
      if (lclCol == lclDiagColInd) {
        diagVal = curRow.value (k); // assume no repeats
      }
      if (! assumeSymmetric) {
        Kokkos::atomic_add (&(result_.colNorms[lclCol]), matrixAbsVal);
      }
    }

    // NOTE (mfh 24 May 2018) We could actually compute local
    // rowScaledColNorms in situ at this point, if ! assumeSymmetric
    // and row Map is the same as range Map (so that the local row
    // norms are the same as the global row norms).

    result_.rowDiagonalEntries[lclRow] = diagVal;
    result_.rowNorms[lclRow] = rowNorm;
    if (! assumeSymmetric &&
        lclDiagColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
      // Don't need an atomic update here, since this lclDiagColInd is
      // a one-to-one function of lclRow.
      result_.colDiagonalEntries[lclDiagColInd] += diagVal;
    }
  }

private:
  using local_matrix_type = typename crs_matrix_type::local_matrix_type;
  local_matrix_type A_lcl_;

  using local_map_type = typename ::Tpetra::Map<LO, GO, NT>::local_map_type;
  local_map_type rowMap_;
  local_map_type colMap_;

  result_type result_;
};


/// \brief Implementation of computeLocalRowAndColumnOneNorms for a
///   Tpetra::CrsMatrix.
template<class SC, class LO, class GO, class NT>
EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type, typename NT::device_type>
computeLocalRowAndColumnOneNorms_CrsMatrix (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                                            const bool assumeSymmetric)
{
  using execution_space = typename NT::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  using functor_type = ComputeLocalRowAndColumnOneNorms<SC, LO, GO, NT>;

  functor_type functor (A, assumeSymmetric);
  const LO lclNumRows = static_cast<LO> (A.getRowMap ()->getNodeNumElements ());
  Kokkos::parallel_for ("computeLocalRowAndColumnOneNorms",
                        range_type (0, lclNumRows), functor);
  return functor.getResult ();
}

/// \brief Compute LOCAL row and column one-norms ("row sums" etc.) of
///   the input sparse matrix A.  Optionally, also compute row-scaled
///   column norms (in the manner of LAPACK's DGEEQU routine).
///
/// \param A [in] The input sparse matrix A.
///
/// \param assumeSymmetric [in] Whether to assume that the matrix A is
///   (globally) symmetric.  If so, don't compute row-scaled column
///   norms separately from row norms.
///
/// This function will only compute (local) row-scaled column norms in
/// the same pass as row norms, if and only if BOTH of the following
/// conditions hold:
/// <ol>
/// <li> When A's row Map and range Map are the same (so that local
///      row norms == global row norms, so that it's correct to scale
///      by row norms in the same pass over the local matrix as
///      computing the row norms) </li>
/// <li> When the matrix is nonsymmetric (otherwise the row norms
///      suffice) </li>
/// </ol>
template<class SC, class LO, class GO, class NT>
EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type, typename NT::device_type>
computeLocalRowAndColumnOneNorms (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                                  const bool assumeSymmetric)
{
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  const crs_matrix_type* A_crs = dynamic_cast<const crs_matrix_type*> (&A);

  if (A_crs == nullptr) {
    return computeLocalRowAndColumnOneNorms_RowMatrix (A, assumeSymmetric);
  }
  else {
    return computeLocalRowAndColumnOneNorms_CrsMatrix (*A_crs, assumeSymmetric);
  }
}

template<class SC, class LO, class GO, class NT>
typename Tpetra::MultiVector<SC, LO, GO, NT>::dual_view_type::t_dev
getLocalView_2d (const Tpetra::MultiVector<SC, LO, GO, NT>& X)
{
  return X.template getLocalView<typename NT::device_type::memory_space> ();
}

template<class SC, class LO, class GO, class NT>
auto getLocalView_1d (const Tpetra::MultiVector<SC, LO, GO, NT>& X,
                      const LO whichColumn)
  -> decltype (Kokkos::subview (getLocalView_2d (X), Kokkos::ALL (), whichColumn))
{
  if (X.isConstantStride ()) {
    return Kokkos::subview (getLocalView_2d (X), Kokkos::ALL (), whichColumn);
  }
  else {
    auto X_whichColumn = X.getVector (whichColumn);
    return Kokkos::subview (getLocalView_2d (*X_whichColumn), Kokkos::ALL (), 0);
  }
}

template<class SC, class LO, class GO, class NT, class ViewValueType>
void
copy1DViewIntoMultiVectorColumn (Tpetra::MultiVector<SC, LO, GO, NT>& X,
                                 const LO whichColumn,
                                 const Kokkos::View<ViewValueType*, typename NT::device_type>& view)
{
  using dev_memory_space = typename NT::device_type::memory_space;
  // MultiVector always starts sync'd to device.
  X.template modify<dev_memory_space> ();
  auto X_lcl = getLocalView_1d (X, whichColumn);
  Tpetra::Details::copyConvert (X_lcl, view);
}

template<class SC, class LO, class GO, class NT, class ViewValueType>
void
copyMultiVectorColumnInto1DView (const Kokkos::View<ViewValueType*, typename NT::device_type>& view,
                                 Tpetra::MultiVector<SC, LO, GO, NT>& X,
                                 const LO whichColumn)
{
  using dev_memory_space = typename NT::device_type::memory_space;
  X.template sync<dev_memory_space> ();
  auto X_lcl = getLocalView_1d (X, whichColumn);
  Tpetra::Details::copyConvert (view, X_lcl);
}

template<class SC, class LO, class GO, class NT>
void
globalizeRowOneNorms (EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                                        typename NT::device_type>& equib,
                      const Tpetra::RowMatrix<SC, LO, GO, NT>& A)
{
  using mv_type = Tpetra::MultiVector<SC, LO, GO, NT>;

  auto G = A.getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.get () == nullptr, std::invalid_argument,
     "globalizeRowOneNorms: Input RowMatrix A must have a nonnull graph "
     "(that is, getGraph() must return nonnull).");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! G->isFillComplete (), std::invalid_argument,
     "globalizeRowOneNorms: Input CrsGraph G must be fillComplete.");

  auto exp = G->getExporter ();
  if (! exp.is_null ()) {
    // If the matrix has an overlapping row Map, first Export the
    // local row norms with ADD CombineMode to a range Map Vector to
    // get the global row norms, then reverse them back with REPLACE
    // CombineMode to the row Map Vector.  Ditto for the local row
    // diagonal entries.  Use SC instead of mag_type, so we can
    // communicate both row norms and row diagonal entries at once.

    // FIXME (mfh 16 May 2018) Clever DualView tricks could possibly
    // avoid the local copy here.
    mv_type rowMapMV (G->getRowMap (), 2, false);

    copy1DViewIntoMultiVectorColumn (rowMapMV, 0, equib.rowNorms);
    copy1DViewIntoMultiVectorColumn (rowMapMV, 1, equib.rowDiagonalEntries);
    {
      mv_type rangeMapMV (G->getRangeMap (), 2, true);
      rangeMapMV.doExport (rowMapMV, *exp, Tpetra::ADD); // forward mode
      rowMapMV.doImport (rangeMapMV, *exp, Tpetra::REPLACE); // reverse mode
    }
    copyMultiVectorColumnInto1DView (equib.rowNorms, rowMapMV, 0);
    copyMultiVectorColumnInto1DView (equib.rowDiagonalEntries, rowMapMV, 1);
  }
}

template<class SC, class LO, class GO, class NT>
void
globalizeColumnOneNorms (EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                                           typename NT::device_type>& equib,
                         const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                         const bool assumeSymmetric) // if so, use row norms
{
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using mv_type = Tpetra::MultiVector<mag_type, LO, GO, NT>;
  using device_type = typename NT::device_type;

  auto G = A.getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.get () == nullptr, std::invalid_argument,
     "globalizeRowOneNorms: Input RowMatrix A must have a nonnull graph "
     "(that is, getGraph() must return nonnull).");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! G->isFillComplete (), std::invalid_argument,
     "globalizeRowOneNorms: Input CrsGraph G must be fillComplete.");

  auto imp = G->getImporter ();
  if (assumeSymmetric) {
    const LO numCols = 2;
    // Redistribute local row info to global column info.

    // Get the data into a MultiVector on the domain Map.
    mv_type rowNorms_domMap (G->getDomainMap (), numCols, false);
    const bool rowMapSameAsDomainMap = G->getRowMap ()->isSameAs (* (G->getDomainMap ()));
    if (rowMapSameAsDomainMap) {
      copy1DViewIntoMultiVectorColumn (rowNorms_domMap, 0, equib.rowNorms);
      copy1DViewIntoMultiVectorColumn (rowNorms_domMap, 1, equib.rowDiagonalEntries);
    }
    else {
      // This is not a common case; it would normally arise when the
      // matrix has an overlapping row Map.
      Tpetra::Export<LO, GO, NT> rowToDom (G->getRowMap (), G->getDomainMap ());
      mv_type rowNorms_rowMap (G->getRowMap (), numCols, true);
      copy1DViewIntoMultiVectorColumn (rowNorms_rowMap, 0, equib.rowNorms);
      copy1DViewIntoMultiVectorColumn (rowNorms_rowMap, 1, equib.rowDiagonalEntries);
      rowNorms_domMap.doExport (rowNorms_rowMap, rowToDom, Tpetra::REPLACE);
    }

    // Use the existing Import to redistribute the row norms from the
    // domain Map to the column Map.
    std::unique_ptr<mv_type> rowNorms_colMap;
    if (imp.is_null ()) {
      // Shallow copy of rowNorms_domMap, since the two must have the
      // same Maps.
      rowNorms_colMap =
        std::unique_ptr<mv_type> (new mv_type (G->getColMap (),
                                               rowNorms_domMap.getDualView ()));
    }
    else {
      rowNorms_colMap =
        std::unique_ptr<mv_type> (new mv_type (G->getColMap (), numCols, true));
      rowNorms_colMap->doImport (rowNorms_domMap, *imp, Tpetra::REPLACE);
    }

    // Make sure the result has allocations of the right size.
    const LO lclNumCols =
      static_cast<LO> (G->getColMap ()->getNodeNumElements ());
    if (static_cast<LO> (equib.colNorms.extent (0)) != lclNumCols) {
      equib.colNorms =
        Kokkos::View<mag_type*, device_type> ("colNorms", lclNumCols);
    }
    if (static_cast<LO> (equib.colDiagonalEntries.extent (0)) != lclNumCols) {
      equib.colDiagonalEntries =
        Kokkos::View<val_type*, device_type> ("colDiagonalEntries", lclNumCols);
    }

    // Copy row norms and diagonal entries, appropriately
    // redistributed, into column norms resp. diagonal entries.
    copyMultiVectorColumnInto1DView (equib.colNorms, *rowNorms_colMap, 0);
    copyMultiVectorColumnInto1DView (equib.colDiagonalEntries, *rowNorms_colMap, 1);
  }
  else {
    if (! imp.is_null ()) {
      const LO numCols = 3;
      // If the matrix has an overlapping column Map (this is usually
      // the case), first Export (reverse-mode Import) the local info
      // to a domain Map Vector to get the global info, then Import
      // them back with REPLACE CombineMode to the column Map Vector.
      // Ditto for the row-scaled column norms.

      // FIXME (mfh 16 May 2018) Clever DualView tricks could possibly
      // avoid the local copy here.
      mv_type colMapMV (G->getColMap (), numCols, false);

      copy1DViewIntoMultiVectorColumn (colMapMV, 0, equib.colNorms);
      copy1DViewIntoMultiVectorColumn (colMapMV, 1, equib.colDiagonalEntries);
      copy1DViewIntoMultiVectorColumn (colMapMV, 2, equib.rowScaledColNorms);
      {
        mv_type domainMapMV (G->getDomainMap (), numCols, true);
        domainMapMV.doExport (colMapMV, *imp, Tpetra::ADD); // reverse mode
        colMapMV.doImport (domainMapMV, *imp, Tpetra::REPLACE); // forward mode
      }
      copyMultiVectorColumnInto1DView (equib.colNorms, colMapMV, 0);
      copyMultiVectorColumnInto1DView (equib.colDiagonalEntries, colMapMV, 1);
      copyMultiVectorColumnInto1DView (equib.rowScaledColNorms, colMapMV, 2);
    }
  }
}

/// \brief Compute global row and column one-norms ("row sums" and
///   "column sums") of the input sparse matrix A, in a way suitable
///   for equilibration.
///
/// \note USERS: This is a function you want.
///
/// \note For AztecOO users: If you set assumeSymmetric=true, this
///   function should behave like setting the <tt>AZ_scaling</tt>
///   option to <tt>AZ_sym_row_sum</tt>.
///
/// \note This function is collective over A's communicator, and may
///   need to communicate, depending on A's Maps.
///
/// For the nonsymmetric case, this function works like a sparse
/// version of LAPACK's DGEEQU routine, except that it uses one norms
/// (sums of absolute values) instead of infinity norms (maximum
/// absolute value).  The resulting row and column scaling is NOT
/// symmetric.  For the symmetric case, this function computes the row
/// norms and uses those for the column norms.  The resulting scaling
/// is symmetric IF you take square roots.
///
/// \param A [in] The input sparse matrix A.
///
/// \param assumeSymmetric [in] Whether to assume that the matrix A is
///   (globally) symmetric.  If so, don't compute row-scaled column
///   norms separately from row norms.
///
/// \return Input to leftAndOrRightScaleCrsMatrix (which see).
template<class SC, class LO, class GO, class NT>
EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                  typename NT::device_type>
computeRowAndColumnOneNorms (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                             const bool assumeSymmetric)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::invalid_argument,
     "computeRowAndColumnOneNorms: Input matrix A must be fillComplete.");
  auto result = computeLocalRowAndColumnOneNorms (A, assumeSymmetric);

  globalizeRowOneNorms (result, A);
  if (! assumeSymmetric) {
    // Row-norm-scaled column norms are trivial if the matrix is
    // symmetric, since the row norms and column norms are the same in
    // that case.
    computeLocalRowScaledColumnNorms (result, A);
  }
  globalizeColumnOneNorms (result, A, assumeSymmetric);
  return result;
}


/// \brief Kokkos::parallel_for functor that left-scales a KokkosSparse::CrsMatrix.
///
/// \tparam LocalSparseMatrixType KokkosSparse::CrsMatrix specialization.
template<class LocalSparseMatrixType>
class LeftScaleLocalCrsMatrix {
public:
  using local_matrix_type = LocalSparseMatrixType;
  using val_type = typename local_matrix_type::value_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using device_type = typename local_matrix_type::device_type;

  LeftScaleLocalCrsMatrix (const local_matrix_type& A_lcl,
                           const EquilibrationInfo<val_type, device_type>& result) :
    A_lcl_ (A_lcl),
    rowNorms_ (result.rowNorms),
    assumeSymmetric_ (result.assumeSymmetric)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const typename local_matrix_type::ordinal_type lclRow) const
  {
    using LO = typename local_matrix_type::ordinal_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    const mag_type curRowNorm = rowNorms_(lclRow);
    // Dividing by zero gives unpleasant results, so don't do it.
    if (curRowNorm > KAM::zero ()) {
      const mag_type scalingFactor = assumeSymmetric_ ?
        KAM::sqrt (curRowNorm) : curRowNorm;
      auto curRow = A_lcl_.row (lclRow);
      const LO numEnt = curRow.length;
      for (LO k = 0; k < numEnt; ++k) {
        curRow.value (k) = curRow.value(k) / scalingFactor;
      }
    }
  }

  /// \brief leftAndOrRightScaleCrsMatrix should call this.
  ///
  /// Create a functor instance and invoke it in a Kokkos::parallel_for.
  static void
  run (const local_matrix_type& A_lcl,
       const EquilibrationInfo<val_type, device_type>& result)
  {
    using execution_space = typename device_type::execution_space;
    using LO = typename local_matrix_type::ordinal_type;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;

    LeftScaleLocalCrsMatrix<local_matrix_type> functor (A_lcl, result);
    const LO lclNumRows = A_lcl.numRows ();
    Kokkos::parallel_for ("leftScaleLocalCrsMatrix",
                          range_type (0, lclNumRows), functor);
  }

private:
  local_matrix_type A_lcl_;
  Kokkos::View<const mag_type*, device_type> rowNorms_;
  bool assumeSymmetric_;
};

/// \brief Kokkos::parallel_for functor that right-scales a KokkosSparse::CrsMatrix.
///
/// \tparam LocalSparseMatrixType KokkosSparse::CrsMatrix specialization.
template<class LocalSparseMatrixType>
class RightScaleLocalCrsMatrix {
public:
  using local_matrix_type = LocalSparseMatrixType;
  using val_type = typename local_matrix_type::value_type;
  using device_type = typename local_matrix_type::device_type;

  RightScaleLocalCrsMatrix (const local_matrix_type& A_lcl,
                            const EquilibrationInfo<val_type, device_type>& result) :
    A_lcl_ (A_lcl),
    scalingFactors_ (result.assumeSymmetric ? result.colNorms : result.rowScaledColNorms),
    assumeSymmetric_ (result.assumeSymmetric)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const typename local_matrix_type::ordinal_type lclRow) const
  {
    using LO = typename local_matrix_type::ordinal_type;
    using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    auto curRow = A_lcl_.row (lclRow);
    const LO numEnt = curRow.length;
    for (LO k = 0; k < numEnt; ++k) {
      const LO lclColInd = curRow.colidx(k);
      const mag_type curColNorm = scalingFactors_(lclColInd);
      if (curColNorm != KAM::zero ()) {
        const mag_type scalingFactor = assumeSymmetric_ ?
          KAM::sqrt (curColNorm) : curColNorm;
        curRow.value(k) = curRow.value(k) / scalingFactor;
      }
    }
  }

  /// \brief leftAndOrRightScaleCrsMatrix should call this.
  ///
  /// Create a functor instance and invoke it in a Kokkos::parallel_for.
  static void
  run (const local_matrix_type& A_lcl,
       const EquilibrationInfo<val_type, device_type>& result)
  {
    using execution_space = typename device_type::execution_space;
    using LO = typename local_matrix_type::ordinal_type;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;

    RightScaleLocalCrsMatrix<local_matrix_type> functor (A_lcl, result);
    Kokkos::parallel_for ("rightScaleLocalCrsMatrix",
                          range_type (0, A_lcl.numRows ()), functor);
  }

private:
  local_matrix_type A_lcl_;

  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  Kokkos::View<const mag_type*, device_type> scalingFactors_;
  bool assumeSymmetric_;
};


/// \brief Left-scale and/or right-scale (in that order) the entries
///   of the input Tpetra::CrsMatrix A.
///
/// \note USERS: This is a function you want.
///
/// \param A [in/out] The sparse matrix A to scale.  It must have a
///   valid KokkosSparse::CrsMatrix.  This is true if fillComplete has
///   been called on it at least once, or if the matrix was created
///   with a local sparse matrix.
///
/// \param equib [in] Return value of computeRowAndColumnNorms (which
///   see), called on the input matrix A.
///
/// \param leftScale [in] Whether to left-scale A.  Left scaling
///   happens first.
///
/// \param rightScale [in] Whether to right-scale A.  Right scaling
///   happens last.
template<class SC, class LO, class GO, class NT>
void
leftAndOrRightScaleCrsMatrix (Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                              const EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type, typename NT::device_type>& equib,
                              const bool leftScale,
                              const bool rightScale)
{
  if (! leftScale && ! rightScale) {
    return;
  }

  const bool A_fillComplete_on_input = A.isFillComplete ();
  if (! A_fillComplete_on_input) {
    // Make sure that A has a valid local matrix.  It might not if it
    // was not created with a local matrix, and if fillComplete has
    // never been called on it before.  It's invalid if the local
    // number of rows is zero, but shouldn't be.
    auto A_lcl = A.getLocalMatrix ();
    const LO lclNumRows =
      static_cast<LO> (A.getRowMap ()->getNodeNumElements ());
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_lcl.numRows () != lclNumRows, std::invalid_argument,
       "leftAndOrRightScaleCrsMatrix: Local matrix is not valid.  "
       "This means that A was not created with a local matrix, "
       "and that fillComplete has never yet been called on A before.  "
       "Please call fillComplete on A at least once first "
       "before calling this method.");
  }
  else {
    A.resumeFill ();
  }

  if (leftScale) {
    using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using impl_type = LeftScaleLocalCrsMatrix<local_matrix_type>;
    impl_type::run (A.getLocalMatrix (), equib);
  }
  if (rightScale) {
    using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using impl_type = RightScaleLocalCrsMatrix<local_matrix_type>;
    impl_type::run (A.getLocalMatrix (), equib);
  }

  if (A_fillComplete_on_input) { // put A back how we found it
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList ();
    params->set ("No Nonlocal Changes", true);
    A.fillComplete (A.getDomainMap (), A.getRangeMap (), params);
  }
}

// NOTE (mfh 21 May 2018) This function is only for the test!
template<class ValueType>
bool
near (const ValueType& x,
      const ValueType& y,
      const typename Kokkos::ArithTraits<ValueType>::mag_type& factor)
{
  const auto eps = Kokkos::ArithTraits<ValueType>::eps ();
  const auto absDiff = Kokkos::ArithTraits<ValueType>::abs (x - y);
  return absDiff <= factor * eps;
}

// NOTE (mfh 21 May 2018) This function is only for the test!
template<class SC, class LO, class GO, class NT>
Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
deepCopyFillCompleteCrsMatrix (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A)
{
  using Teuchos::RCP;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;

  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::invalid_argument,
     "deepCopyFillCompleteCrsMatrix: Input matrix A must be fillComplete.");
  RCP<crs_matrix_type> A_copy (new crs_matrix_type (A.getCrsGraph ()));
  auto A_copy_lcl = A_copy->getLocalMatrix ();
  auto A_lcl = A.getLocalMatrix ();
  Kokkos::deep_copy (A_copy_lcl.values, A_lcl.values);
  A_copy->fillComplete (A.getDomainMap (), A.getRangeMap ());
  return A_copy;
}

// NOTE (mfh 24 May 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
void
testCrsMatrixEquality (bool& success,
                       Teuchos::FancyOStream& out,
                       const Tpetra::CrsMatrix<SC, LO, GO, NT>& A_expected,
                       const Tpetra::CrsMatrix<SC, LO, GO, NT>& A_actual)
{
  using std::endl;
  using mag_type = typename Kokkos::ArithTraits<SC>::mag_type;
  const mag_type toleranceFactor = 10.0; // factor of eps

  auto A_expected_lcl = A_expected.getLocalMatrix ();
  auto ptr_h = Kokkos::create_mirror_view (A_expected_lcl.graph.row_map);
  Kokkos::deep_copy (ptr_h, A_expected_lcl.graph.row_map);

  auto expected_val_h = Kokkos::create_mirror_view (A_expected_lcl.values);
  Kokkos::deep_copy (expected_val_h, A_expected_lcl.values);

  auto A_actual_lcl = A_actual.getLocalMatrix ();
  auto actual_val_h = Kokkos::create_mirror_view (A_actual_lcl.values);
  Kokkos::deep_copy (actual_val_h, A_actual_lcl.values);

  using size_type = typename decltype (A_actual_lcl.graph)::size_type;

  if (A_expected_lcl.numRows () != A_actual_lcl.numRows ()) {
    TEST_EQUALITY( A_expected_lcl.numRows (), A_actual_lcl.numRows () );
    return;
  }
  const LO lclNumRows = A_expected_lcl.numRows ();

  bool ok = true;
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
      if (! near (actual_val_h[k],
                  expected_val_h[k],
                  toleranceFactor)) {
        ok = false;
        break;
      }
    }
  }
  if (! ok) {
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      out << "lclRow: " << lclRow << endl;
      Teuchos::OSTab tab2 (out);

      using size_type = typename decltype (A_actual_lcl.graph)::size_type;
      out << "Expected values: [";
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        out << expected_val_h[k];
        if (k + size_type (1) < ptr_h[lclRow+1]) {
          out << ", ";
        }
      }
      out << "]" << endl
          << "Actual values: [";
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        out << actual_val_h[k];
        if (k + size_type (1) < ptr_h[lclRow+1]) {
          out << ", ";
        }
      }
      out << "]" << endl;
    }
  }
}


// NOTE (mfh 23 May 2018) This struct is only for the test.
template<class SC, class LO, class GO, class NT>
struct EquilibrationTest {
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;

  Teuchos::RCP<crs_matrix_type> A;
  std::vector<mag_type> lclRowNorms;
  std::vector<val_type> lclRowDiagonalEntries;
  std::vector<mag_type> lclColNorms;
  std::vector<val_type> lclColDiagonalEntries;
  std::vector<mag_type> lclRowScaledColNorms;
  std::vector<mag_type> gblRowNorms;
  std::vector<mag_type> gblColNorms;
  std::vector<mag_type> gblRowScaledColNorms;
  Teuchos::RCP<crs_matrix_type> A_leftScaled;
  Teuchos::RCP<crs_matrix_type> A_rightScaled;
};

// NOTE (mfh 23 May 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
void
testEquilibration (Teuchos::FancyOStream& out,
                   bool& success,
                   const EquilibrationTest<SC, LO, GO, NT>& test,
                   const bool assumeSymmetric)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using row_matrix_type = Tpetra::RowMatrix<SC, LO, GO, NT>;
  using mag_type = typename Kokkos::ArithTraits<SC>::mag_type;

  const mag_type toleranceFactor = 10.0; // factor of eps
  int lclSuccess = success ? 1 : 0; // for reduceAll (see below)
  int gblSuccess = 0; // output argument for reduceAll (see below)
  auto comm = test.A->getMap ()->getComm ();

  out << "Test equilibration: assumeSymmetric="
      << (assumeSymmetric ? "true" : "false") << endl;
  Teuchos::OSTab tab1 (out);

  const LO lclNumRows =
    static_cast<LO> (test.A->getRowMap ()->getNodeNumElements ());
  RCP<const map_type> colMap = test.A->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());

  // Test computeLocalRowAndColumnOneNorms (CrsMatrix)
  {
    out << "Test computeLocalRowAndColumnOneNorms (CrsMatrix)" << endl;
    Teuchos::OSTab tab2 (out);
    auto result0 = computeLocalRowAndColumnOneNorms (* (test.A), assumeSymmetric);

    {
      out << "Test local row norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowNorms_h = Kokkos::create_mirror_view (result0.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result0.rowNorms);

      std::ostringstream os;
      os << "Expected local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.lclRowNorms[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowNorms_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.lclRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row norms

    {
      out << "Test local row diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowDiagonalEntries_h =
        Kokkos::create_mirror_view (result0.rowDiagonalEntries);
      Kokkos::deep_copy (rowDiagonalEntries_h, result0.rowDiagonalEntries);

      std::ostringstream os;
      os << "Expected local rowDiagonalEntries: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.lclRowDiagonalEntries[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowDiagonalEntries: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowDiagonalEntries_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowDiagonalEntries_h[lclRow],
                    test.lclRowDiagonalEntries[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row diagonal entries

    if (! assumeSymmetric) {
      out << "Test local column norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto colNorms_h = Kokkos::create_mirror_view (result0.colNorms);
      Kokkos::deep_copy (colNorms_h, result0.colNorms);

      std::ostringstream os;
      os << "Expected local colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.lclColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << colNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.lclColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local column norms

    if (! assumeSymmetric) {
      out << "Test local column diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);

      auto colDiagonalEntries_h =
        Kokkos::create_mirror_view (result0.colDiagonalEntries);
      Kokkos::deep_copy (colDiagonalEntries_h, result0.colDiagonalEntries);

      std::ostringstream os;
      os << "Expected local colDiagonalEntries: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.lclColDiagonalEntries[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local colDiagonalEntries: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << colDiagonalEntries_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colDiagonalEntries_h[lclCol],
                    test.lclColDiagonalEntries[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local column diagonal entries

    out << "Test globalizing row one norms" << endl;
    globalizeRowOneNorms (result0, * (test.A));

    {
      out << "Test global row norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowNorms_h = Kokkos::create_mirror_view (result0.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result0.rowNorms);

      std::ostringstream os;
      os << "Expected global rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.gblRowNorms[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowNorms_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.gblRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row norms

    if (! assumeSymmetric) {
      out << "Test computeLocalRowScaledColumnNorms" << endl;
      computeLocalRowScaledColumnNorms (result0, * (test.A));

      auto rowScaledColNorms_h =
        Kokkos::create_mirror_view (result0.rowScaledColNorms);
      Kokkos::deep_copy (rowScaledColNorms_h, result0.rowScaledColNorms);

      std::ostringstream os;
      os << "Expected local rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.lclRowScaledColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << rowScaledColNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (rowScaledColNorms_h[lclCol],
                    test.lclRowScaledColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test computeLocalRowScaledColumnNorms

    out << "Test globalizing column one norms" << endl;
    globalizeColumnOneNorms (result0, * (test.A), assumeSymmetric);

    {
      out << "Test global column norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto colNorms_h = Kokkos::create_mirror_view (result0.colNorms);
      Kokkos::deep_copy (colNorms_h, result0.colNorms);

      std::ostringstream os;
      os << "Expected global colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.gblColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << colNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.gblColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global column norms

    if (! assumeSymmetric) {
      out << "Test global row-scaled column norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowScaledColNorms_h =
        Kokkos::create_mirror_view (result0.rowScaledColNorms);
      Kokkos::deep_copy (rowScaledColNorms_h, result0.rowScaledColNorms);

      std::ostringstream os;
      os << "Expected global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.gblRowScaledColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << rowScaledColNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (rowScaledColNorms_h[lclCol],
                    test.gblRowScaledColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row-scaled column norms

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED on some process!" << endl;
      return;
    }
  } // test computeLocalRowAndColumnNorms (CrsMatrix)

  // Test computeLocalRowAndColumnNorms (RowMatrix)
  {
    out << "Test computeLocalRowAndColumnOneNorms (RowMatrix)" << endl;
    Teuchos::OSTab tab2 (out);
    auto result1 =
      computeLocalRowAndColumnOneNorms_RowMatrix (* (test.A), assumeSymmetric);

    {
      out << "Test local row norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowNorms_h = Kokkos::create_mirror_view (result1.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result1.rowNorms);

      std::ostringstream os;
      os << "Expected local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.lclRowNorms[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowNorms_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.lclRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row norms

    {
      out << "Test local row diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowDiagonalEntries_h =
        Kokkos::create_mirror_view (result1.rowDiagonalEntries);
      Kokkos::deep_copy (rowDiagonalEntries_h, result1.rowDiagonalEntries);

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowDiagonalEntries_h[lclRow],
                    test.lclRowDiagonalEntries[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row diagonal entries

    if (! assumeSymmetric) {
      out << "Test local column norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto colNorms_h = Kokkos::create_mirror_view (result1.colNorms);
      Kokkos::deep_copy (colNorms_h, result1.colNorms);

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.lclColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // ! assumeSymmetric

    if (! assumeSymmetric) {
      out << "Test local column diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);
      auto colDiagonalEntries_h =
        Kokkos::create_mirror_view (result1.colDiagonalEntries);
      Kokkos::deep_copy (colDiagonalEntries_h, result1.colDiagonalEntries);

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colDiagonalEntries_h[lclCol],
                    test.lclColDiagonalEntries[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // ! assumeSymmetric

    // We've already tested globalize{Row,Column}OneNorms above;
    // neither depends on whether the matrix is a CrsMatrix.

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED on some process!" << endl;
      return;
    }
  } // test computeLocalRowAndColumnNorms (RowMatrix)

  // Test computeRowAndColumnOneNorms
  {
    out << "Test computeRowAndColumnNorms" << endl;
    Teuchos::OSTab tab2 (out);
    // FIXME (mfh 24 May 2018) Why the cast?
    auto result2 = computeRowAndColumnOneNorms (static_cast<row_matrix_type&> (* (test.A)),
                                                assumeSymmetric);
    {
      out << "Test whether assumeSymmetric got communicated" << endl;
      Teuchos::OSTab tab3 (out);
      TEST_EQUALITY( assumeSymmetric, result2.assumeSymmetric );
    }

    {
      out << "Test global row norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowNorms_h = Kokkos::create_mirror_view (result2.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result2.rowNorms);

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.gblRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row norms

    {
      out << "Test global column norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto colNorms_h = Kokkos::create_mirror_view (result2.colNorms);
      Kokkos::deep_copy (colNorms_h, result2.colNorms);

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.gblColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global column norms

    if (! assumeSymmetric) {
      out << "Test global row-scaled column norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowScaledColNorms_h =
        Kokkos::create_mirror_view (result2.rowScaledColNorms);
      Kokkos::deep_copy (rowScaledColNorms_h, result2.rowScaledColNorms);

      std::ostringstream os;
      os << "Expected global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.gblRowScaledColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << rowScaledColNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (rowScaledColNorms_h[lclCol],
                    test.gblRowScaledColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row-scaled column norms

    {
      out << "Test deepCopyFillCompleteCrsMatrix" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      testCrsMatrixEquality (success, out, *(test.A), *A_copy);
    }

    {
      out << "Test left-scaling CrsMatrix" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      leftAndOrRightScaleCrsMatrix (*A_copy, result2, true, false);
      testCrsMatrixEquality (success, out, *(test.A_leftScaled), *A_copy);
    }

    {
      out << "Test right-scaling CrsMatrix" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      leftAndOrRightScaleCrsMatrix (*A_copy, result2, false, true);
      testCrsMatrixEquality (success, out, *(test.A_rightScaled), *A_copy);
    }
  } // test computeRowAndColumnOneNorms

  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Test FAILED on some process!" << endl;
  }
}

// NOTE (mfh 23 May 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
EquilibrationTest<SC, LO, GO, NT>
makeSymmetricPositiveDefiniteTridiagonalMatrixTest (Teuchos::FancyOStream& out,
                                                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                                    const bool assumeSymmetric)
{
  using Teuchos::RCP;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KAT = Kokkos::ArithTraits<SC>;
  using val_type = typename KAT::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using KAM = Kokkos::ArithTraits<mag_type>;

  const LO lclNumRows = 5;
  const GO gblNumRows =
    static_cast<GO> (comm->getSize ()) * static_cast<GO> (lclNumRows);
  const GO indexBase = 0;

  out << "Create symmetric positive definite tridiagonal matrix problem"
      << endl;
  Teuchos::OSTab tab0 (out);

  out << "Create Maps" << endl;
  RCP<const map_type> rowMap =
    rcp (new map_type (static_cast<GST> (gblNumRows),
                       static_cast<size_t> (lclNumRows),
                       indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Create CrsGraph" << endl;
  const size_t maxNumEntPerRow = 3;
  RCP<crs_graph_type> G =
    rcp (new crs_graph_type (rowMap, maxNumEntPerRow, Tpetra::StaticProfile));
  std::vector<GO> globalIndices (maxNumEntPerRow);

  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);

    if (gblRow == 0) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow;
      globalIndices[1] = gblRow+1;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
    else if (gblRow == gblNumRows - GO (1)) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
    else {
      const LO numEnt = 3;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      globalIndices[2] = gblRow+1;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
  }
  G->fillComplete (domMap, ranMap);

  const SC diagVal {2.0};
  const SC offDiagVal {-1.0};
  std::vector<SC> firstRowValues {diagVal, offDiagVal};
  std::vector<SC> middleRowValues {offDiagVal, diagVal, offDiagVal};
  std::vector<SC> lastRowValues {offDiagVal, diagVal};

  out << "Create test CrsMatrix" << endl;
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (G));
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);

    if (gblRow == 0) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow;
      globalIndices[1] = gblRow+1;
      A->replaceGlobalValues (gblRow, numEnt, firstRowValues.data (), globalIndices.data ());
    }
    else if (gblRow == gblNumRows - GO (1)) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      A->replaceGlobalValues (gblRow, numEnt, lastRowValues.data (), globalIndices.data ());
    }
    else {
      const LO numEnt = 3;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      globalIndices[2] = gblRow+1;
      A->replaceGlobalValues (gblRow, numEnt, middleRowValues.data (), globalIndices.data ());
    }
  }
  A->fillComplete (domMap, ranMap);

  RCP<const map_type> colMap = G->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());
  const GO gblNumCols = static_cast<GO> (G->getDomainMap ()->getGlobalNumElements ());

  const mag_type diagAbsVal = KAT::abs (diagVal);
  const mag_type offDiagAbsVal = KAT::abs (offDiagVal);

  out << "Compute local row norms and diagonal entries" << endl;
  std::vector<val_type> lclRowDiagonalEntries (lclNumRows);
  std::vector<mag_type> lclRowNorms (lclNumRows);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);
    const mag_type expectedRowNorm =
      (gblRow == 0 || gblRow == gblNumRows - GO (1)) ?
      offDiagAbsVal + diagAbsVal :
      offDiagAbsVal + diagAbsVal + offDiagAbsVal;

    lclRowDiagonalEntries[lclRow] = diagVal;
    lclRowNorms[lclRow] = expectedRowNorm;
  }

  // For this matrix, the global row norms are the same as the local
  // row norms, since the matrix's row Map and range Map are the same.
  // This may not be the case for matrices with an overlapping or
  // permuted row Map.
  out << "Compute global row norms" << endl;
  std::vector<mag_type> gblRowNorms (lclRowNorms.begin (), lclRowNorms.end ());

  Teuchos::RCP<crs_matrix_type> A_leftScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);

      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const GO gblCol = colMap->getGlobalElement (lclCol);
        const val_type expectedUnscaledVal =
          (gblRow == gblCol) ? diagVal : offDiagVal;
        const mag_type rowNorm = gblRowNorms[lclRow];
        const mag_type scalingFactor = assumeSymmetric ?
          KAM::sqrt (rowNorm) : rowNorm;
        val_h[k] = expectedUnscaledVal / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  if (assumeSymmetric) {
    out << "assumeSymmetric=true: Skip local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }
  else {
    out << "assumeSymmetric=false: Compute local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }

  std::vector<mag_type> lclColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<val_type> lclColDiagonalEntries
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<mag_type> lclRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);

  if (! assumeSymmetric) {
    // Columns are a little more complicated, since in the usual
    // distributed case, the local column norms may not be the same as
    // the global column norms.
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      const GO gblCol = colMap->getGlobalElement (lclCol);
      const GO gblRow = gblCol;

      mag_type expectedColNorm {0.0};
      val_type expectedColDiagonalVal {0.0};

      // These are local column norms, not global column norms.
      if (rowMap->isNodeGlobalElement (gblRow-1)) {
        expectedColNorm += offDiagAbsVal;
      }
      if (rowMap->isNodeGlobalElement (gblRow)) {
        expectedColNorm += diagAbsVal;
        expectedColDiagonalVal += diagVal;
      }
      if (rowMap->isNodeGlobalElement (gblRow+1)) {
        expectedColNorm += offDiagAbsVal;
      }

      lclColNorms[lclCol] = expectedColNorm;
      lclColDiagonalEntries[lclCol] = expectedColDiagonalVal;
    } // for each local column

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      mag_type lclRowScaledColNorm {0.0};

      if (gblRow == rowMap->getMinAllGlobalIndex ()) {
        const LO diagLclColInd = colMap->getLocalElement (gblCol);
        const LO rightLclColInd = colMap->getLocalElement (gblCol + GO (1));

        lclRowScaledColNorms[diagLclColInd] += diagAbsVal / gblRowNorms[lclRow];
        lclRowScaledColNorms[rightLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
      }
      else if (gblRow == rowMap->getMaxAllGlobalIndex ()) {
        const LO leftLclColInd = colMap->getLocalElement (gblCol - GO (1));
        const LO diagLclColInd = colMap->getLocalElement (gblCol);

        lclRowScaledColNorms[leftLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
        lclRowScaledColNorms[diagLclColInd] += diagAbsVal / gblRowNorms[lclRow];
      }
      else {
        const LO leftLclColInd = colMap->getLocalElement (gblCol - GO (1));
        const LO diagLclColInd = colMap->getLocalElement (gblCol);
        const LO rightLclColInd = colMap->getLocalElement (gblCol + GO (1));

        lclRowScaledColNorms[leftLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
        lclRowScaledColNorms[diagLclColInd] += diagAbsVal / gblRowNorms[lclRow];
        lclRowScaledColNorms[rightLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
      }
    }
  } // ! assumeSymmetric

  out << "Compute global column norms" << endl;
  std::vector<mag_type> gblColNorms (lclNumCols);
  // The matrix is symmetric, so this holds regardless of assumeSymmetric.
  for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
    const GO gblCol = colMap->getGlobalElement (lclCol);

    mag_type expectedColNorm {0.0};
    if (gblCol == 0 || gblCol == gblNumCols - GO (1)) {
      expectedColNorm = diagAbsVal + offDiagAbsVal;
    }
    else {
      expectedColNorm = offDiagAbsVal + diagAbsVal + offDiagAbsVal;
    }
    gblColNorms[lclCol] = expectedColNorm;
  }

  if (assumeSymmetric) {
    out << "Skip global row-scaled column norms" << endl;
  }
  else {
    out << "Compute global row-scaled column norms" << endl;
  }
  std::vector<mag_type> gblRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  if (! assumeSymmetric && colMap->getGlobalNumElements () != 0) {
    const GO gblMinGblCol = domMap->getMinAllGlobalIndex ();
    const GO gblMaxGblCol = domMap->getMaxAllGlobalIndex ();

    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      const GO gblCol = colMap->getGlobalElement (lclCol);

      if (gblCol == gblMinGblCol || gblCol == gblMaxGblCol) {
        gblRowScaledColNorms[lclCol] = 11.0 / 12.0; // 2/3 + 1/4
      }
      else if (gblCol == gblMinGblCol + GO (1) ||
               gblCol == gblMaxGblCol - GO (1)) {
        gblRowScaledColNorms[lclCol] = 13.0 / 12.0; // 1/3 + 1/2 + 1/4
      }
      else {
        gblRowScaledColNorms[lclCol] = 1.0; // 1/4 + 1/2 + 1/4
      }
    }
  }

  Teuchos::RCP<crs_matrix_type> A_rightScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const mag_type colNorm = assumeSymmetric ?
          gblColNorms[lclCol] : gblRowScaledColNorms[lclCol];
        // const mag_type rowScalingFactor = assumeSymmetric ?
        //   KAM::sqrt (gblRowNorms[lclRow]) : gblRowNorms[lclRow];
        const mag_type colScalingFactor = assumeSymmetric ?
          KAM::sqrt (colNorm) : colNorm;
        // const mag_type scalingFactor = rowScalingFactor * colScalingFactor;
        const mag_type scalingFactor = colScalingFactor;
        val_h[k] = val_h[k] / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  return EquilibrationTest<SC, LO, GO, NT> {
    A,
    lclRowNorms,
    lclRowDiagonalEntries,
    lclColNorms,
    lclColDiagonalEntries,
    lclRowScaledColNorms,
    gblRowNorms,
    gblColNorms,
    gblRowScaledColNorms,
    A_leftScaled,
    A_rightScaled
  };
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Equilibration, Test0, SC, LO, GO, NT)
{
  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)
  using std::endl;

  //const bool debugMode = true;
  const bool debugMode = false;
  Teuchos::RCP<Teuchos::FancyOStream> errStrmPtr = debugMode ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) : Teuchos::null;
  Teuchos::FancyOStream& testOut = debugMode ? *errStrmPtr : out;

  testOut << "Test equilibration" << endl;
  Teuchos::OSTab tab0 (testOut);
  auto comm = Tpetra::getDefaultComm ();

  for (bool assumeSymmetric : {false, true}) {
    // The default FancyOStream 'out' only prints to Process 0 by
    // default.  Thus, we gather up output from each process into a
    // single string, and print it on Process 0 at the end.
    Teuchos::RCP<std::ostringstream> osPtr (new std::ostringstream);
    Teuchos::RCP<Teuchos::FancyOStream> curOutPtr = Teuchos::getFancyOStream (osPtr);
    Teuchos::FancyOStream& curOut = *curOutPtr;

    curOut << endl << ">>> Process " << comm->getRank () << ":" << endl << endl
           << "assumeSymmetric=" << (assumeSymmetric ? "true" : "false")
           << endl;
    Teuchos::OSTab tab1 (curOut);
    auto test0 =
      makeSymmetricPositiveDefiniteTridiagonalMatrixTest<SC, LO, GO, NT> (curOut, comm,
                                                                          assumeSymmetric);
    bool curSuccess = true;
    testEquilibration (curOut, curSuccess, test0, assumeSymmetric);
    success = success && curSuccess;
    Tpetra::Details::gathervPrint (testOut, osPtr->str (), *comm);

    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED on some process!" << endl;
      return;
    }
  }
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Equilibration, Test0, SC, LO, GO, NT )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type;
// using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type;
// using default_node_type = Tpetra::Map<>::node_type;

IFPACK2_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

// UNIT_TEST_GROUP( double, default_local_ordinal_type, default_global_ordinal_type, default_node_type )

} // namespace (anonymous)

