// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_COMPUTEROWANDCOLUMNONENORMS_DEF_HPP
#define TPETRA_COMPUTEROWANDCOLUMNONENORMS_DEF_HPP

/// \file Tpetra_computeRowAndColumnOneNorms_def.hpp
/// \brief Definition of Tpetra::computeRowAndColumnOneNorms.
///
/// For the declaration of this function and its public Doxygen
/// documentation, please see
/// Tpetra_computeRowAndColumnOneNorms_decl.hpp in this directory.

#include "Tpetra_Details_copyConvert.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <memory>

namespace Tpetra {
namespace Details {

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
  using KAM = Kokkos::ArithTraits<mag_type>;
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
        const val_type matrixVal = val[k];
        if (KAT::isInf (matrixVal)) {
          result_h.foundInf = true;
        }
        if (KAT::isNan (matrixVal)) {
          result_h.foundNan = true;
        }
        const mag_type matrixAbsVal = KAT::abs (matrixVal);
        rowNorm += matrixAbsVal;
        const LO lclCol = ind[k];
        if (lclCol == lclDiagColInd) {
          diagVal += val[k]; // repeats count additively
        }
        if (! assumeSymmetric) {
          result_h.colNorms[lclCol] += matrixAbsVal;
        }
      } // for each entry in row

      // This is a local result.  If the matrix has an overlapping
      // row Map, then the global result might differ.
      if (diagVal == KAT::zero ()) {
        result_h.foundZeroDiag = true;
      }
      if (rowNorm == KAM::zero ()) {
        result_h.foundZeroRowNorm = true;
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

// Kokkos::parallel_reduce functor that is part of the implementation
// of computeLocalRowAndColumnOneNorms_CrsMatrix.
template<class SC, class LO, class GO, class NT>
class ComputeLocalRowAndColumnOneNorms {
public:
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using equib_info_type = EquilibrationInfo<val_type, typename NT::device_type>;
  using local_matrix_type = typename ::Tpetra::CrsMatrix<SC, LO, GO, NT>::local_matrix_type;
  using local_map_type = typename ::Tpetra::Map<LO, GO, NT>::local_map_type;

public:
  ComputeLocalRowAndColumnOneNorms (const equib_info_type& equib,   // in/out
                                    const local_matrix_type& A_lcl, // in
                                    const local_map_type& rowMap,   // in
                                    const local_map_type& colMap) : // in
    equib_ (equib),
    A_lcl_ (A_lcl),
    rowMap_ (rowMap),
    colMap_ (colMap)
  {}

  // (result & 1) != 0 means "found Inf."
  // (result & 2) != 0 means "found NaN."
  // (result & 4) != 0 means "found zero diag."
  // (result & 8) != 0 means "found zero row norm."
  // Pack into a single int so the reduction is cheaper,
  // esp. on GPU.
  using value_type = int;

  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    dst = 0;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    dst |= src;
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const LO lclRow, value_type& dst) const
  {
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    const GO gblRow = rowMap_.getGlobalElement (lclRow);
    // OK if invalid(); then we simply won't find the diagonal entry.
    const GO lclDiagColInd = colMap_.getLocalElement (gblRow);

    const auto curRow = A_lcl_.rowConst (lclRow);
    const LO numEnt = curRow.length;
    const bool assumeSymmetric = equib_.assumeSymmetric;

    mag_type rowNorm {0.0};
    val_type diagVal {0.0};

    for (LO k = 0; k < numEnt; ++k) {
      const val_type matrixVal = curRow.value (k);
      if (KAT::isInf (matrixVal)) {
        dst |= 1;
      }
      if (KAT::isNan (matrixVal)) {
        dst |= 2;
      }
      const mag_type matrixAbsVal = KAT::abs (matrixVal);
      rowNorm += matrixAbsVal;
      const LO lclCol = curRow.colidx (k);
      if (lclCol == lclDiagColInd) {
        diagVal = curRow.value (k); // assume no repeats
      }
      if (! assumeSymmetric) {
        Kokkos::atomic_add (&(equib_.colNorms[lclCol]), matrixAbsVal);
      }
    } // for each entry in row

    // This is a local result.  If the matrix has an overlapping
    // row Map, then the global result might differ.
    if (diagVal == KAT::zero ()) {
      dst |= 4;
    }
    if (rowNorm == KAM::zero ()) {
      dst |= 8;
    }
    // NOTE (mfh 24 May 2018) We could actually compute local
    // rowScaledColNorms in situ at this point, if ! assumeSymmetric
    // and row Map is the same as range Map (so that the local row
    // norms are the same as the global row norms).
    equib_.rowDiagonalEntries[lclRow] = diagVal;
    equib_.rowNorms[lclRow] = rowNorm;
    if (! assumeSymmetric &&
        lclDiagColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
      // Don't need an atomic update here, since this lclDiagColInd is
      // a one-to-one function of lclRow.
      equib_.colDiagonalEntries[lclDiagColInd] += diagVal;
    }
  }

private:
  equib_info_type equib_;
  local_matrix_type A_lcl_;
  local_map_type rowMap_;
  local_map_type colMap_;
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
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using device_type = typename NT::device_type;
  using equib_info_type = EquilibrationInfo<val_type, device_type>;

  const LO lclNumRows = static_cast<LO> (A.getRowMap ()->getNodeNumElements ());
  const LO lclNumCols = static_cast<LO> (A.getColMap ()->getNodeNumElements ());
  equib_info_type equib (lclNumRows, lclNumCols, assumeSymmetric);

  functor_type functor (equib, A.getLocalMatrix (),
                        A.getRowMap ()->getLocalMap (),
                        A.getColMap ()->getLocalMap ());
  int result = 0;
  Kokkos::parallel_reduce ("computeLocalRowAndColumnOneNorms",
                           range_type (0, lclNumRows), functor,
                           result);
  equib.foundInf = (result & 1) != 0;
  equib.foundNan = (result & 2) != 0;
  equib.foundZeroDiag = (result & 4) != 0;
  equib.foundZeroRowNorm = (result & 8) != 0;
  return equib;
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

template<class OneDViewType, class IndexType>
class FindZero {
public:
  static_assert (OneDViewType::Rank == 1,
                 "OneDViewType must be a rank-1 Kokkos::View.");
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  FindZero (const OneDViewType& x) : x_ (x) {}
  // Kokkos historically didn't like bool reduction results on CUDA,
  // so we use int as the reduction result type.
  KOKKOS_INLINE_FUNCTION void
  operator () (const IndexType i, int& result) const {
    using val_type = typename OneDViewType::non_const_value_type;
    result = (x_(i) == Kokkos::ArithTraits<val_type>::zero ()) ? 1 : result;
  }
private:
  OneDViewType x_;
};

template<class OneDViewType>
bool findZero (const OneDViewType& x)
{
  using view_type = typename OneDViewType::const_type;
  using execution_space = typename view_type::execution_space;
  using size_type = typename view_type::size_type;
  using functor_type = FindZero<view_type, size_type>;

  Kokkos::RangePolicy<execution_space, size_type> range (0, x.extent (0));
  range.set (Kokkos::ChunkSize (500)); // adjust as needed

  int foundZero = 0;
  Kokkos::parallel_reduce ("findZero", range, functor_type (x), foundZero);
  return foundZero == 1;
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

    // It's not common for users to solve linear systems with a
    // nontrival Export, so it's OK for this to cost an additional
    // pass over rowDiagonalEntries.
    equib.foundZeroDiag = findZero (equib.rowDiagonalEntries);
    equib.foundZeroRowNorm = findZero (equib.rowNorms);
  }

  constexpr int allReduceCount = 4;
  int lclNaughtyMatrix[allReduceCount];
  lclNaughtyMatrix[0] = equib.foundInf ? 1 : 0;
  lclNaughtyMatrix[1] = equib.foundNan ? 1 : 0;
  lclNaughtyMatrix[2] = equib.foundZeroDiag ? 1 : 0;
  lclNaughtyMatrix[3] = equib.foundZeroRowNorm ? 1 : 0;

  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;
  auto comm = G->getComm ();
  int gblNaughtyMatrix[allReduceCount];
  reduceAll<int, int> (*comm, REDUCE_MAX, allReduceCount,
                       lclNaughtyMatrix, gblNaughtyMatrix);

  equib.foundInf = gblNaughtyMatrix[0] == 1;
  equib.foundNan = gblNaughtyMatrix[1] == 1;
  equib.foundZeroDiag = gblNaughtyMatrix[2] == 1;
  equib.foundZeroRowNorm = gblNaughtyMatrix[3] == 1;
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
     "globalizeColumnOneNorms: Input RowMatrix A must have a nonnull graph "
     "(that is, getGraph() must return nonnull).");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! G->isFillComplete (), std::invalid_argument,
     "globalizeColumnOneNorms: Input CrsGraph G must be fillComplete.");

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

} // namespace Details

template<class SC, class LO, class GO, class NT>
Details::EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                           typename NT::device_type>
computeRowAndColumnOneNorms (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                             const bool assumeSymmetric)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::invalid_argument,
     "computeRowAndColumnOneNorms: Input matrix A must be fillComplete.");
  auto result = Details::computeLocalRowAndColumnOneNorms (A, assumeSymmetric);

  Details::globalizeRowOneNorms (result, A);
  if (! assumeSymmetric) {
    // Row-norm-scaled column norms are trivial if the matrix is
    // symmetric, since the row norms and column norms are the same in
    // that case.
    Details::computeLocalRowScaledColumnNorms (result, A);
  }
  Details::globalizeColumnOneNorms (result, A, assumeSymmetric);
  return result;
}

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_COMPUTEROWANDCOLUMNONENORMS_INSTANT(SC,LO,GO,NT) \
  template Details::EquilibrationInfo<Kokkos::ArithTraits<SC>::val_type, NT::device_type> \
  computeRowAndColumnOneNorms (const Tpetra::RowMatrix<SC, LO, GO, NT>& A, \
                               const bool assumeSymmetric);

#endif // TPETRA_COMPUTEROWANDCOLUMNONENORMS_DEF_HPP
