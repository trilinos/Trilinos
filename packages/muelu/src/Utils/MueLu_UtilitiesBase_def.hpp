// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_UTILITIESBASE_DEF_HPP
#define MUELU_UTILITIESBASE_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_UtilitiesBase_decl.hpp"

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_getDiagCopy.hpp>

#include <Xpetra_BlockedVector.hpp>
#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_ExportFactory.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_StridedMap.hpp>

#include "MueLu_Exceptions.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include <KokkosKernels_Handle.hpp>
#include <KokkosGraph_RCM.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Crs2Op(RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  if (Op.is_null())
    return Teuchos::null;
  return rcp(new CrsMatrixWrap(Op));
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
removeSmallEntries(Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                   const typename Teuchos::ScalarTraits<Scalar>::magnitudeType threshold,
                   const bool keepDiagonal) {
  using crs_matrix      = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using row_ptr_type    = typename crs_matrix::local_graph_type::row_map_type::non_const_type;
  using col_idx_type    = typename crs_matrix::local_graph_type::entries_type::non_const_type;
  using vals_type       = typename crs_matrix::local_matrix_type::values_type;
  using execution_space = typename crs_matrix::local_matrix_type::execution_space;

  using ATS      = Kokkos::ArithTraits<Scalar>;
  using impl_SC  = typename ATS::val_type;
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;

  auto lclA = A->getLocalMatrixDevice();

  auto rowptr = row_ptr_type("rowptr", lclA.numRows() + 1);
  col_idx_type idx;
  vals_type vals;
  LocalOrdinal nnz;

  if (keepDiagonal) {
    auto lclRowMap = A->getRowMap()->getLocalMap();
    auto lclColMap = A->getColMap()->getLocalMap();
    Kokkos::parallel_scan(
        "removeSmallEntries::rowptr",
        Kokkos::RangePolicy<LocalOrdinal, execution_space>(0, lclA.numRows()),
        KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {
          auto row      = lclA.row(rlid);
          auto rowInCol = lclColMap.getLocalElement(lclRowMap.getGlobalElement(rlid));
          for (LocalOrdinal k = 0; k < row.length; ++k) {
            if ((impl_ATS::magnitude(row.value(k)) > threshold) || (row.colidx(k) == rowInCol)) {
              partial_nnz += 1;
            }
          }
          if (is_final)
            rowptr(rlid + 1) = partial_nnz;
        },
        nnz);

    idx  = col_idx_type("idx", nnz);
    vals = vals_type("vals", nnz);

    Kokkos::parallel_for(
        "removeSmallEntries::indicesValues",
        Kokkos::RangePolicy<LocalOrdinal, execution_space>(0, lclA.numRows()),
        KOKKOS_LAMBDA(const LocalOrdinal rlid) {
          auto row      = lclA.row(rlid);
          auto rowInCol = lclColMap.getLocalElement(lclRowMap.getGlobalElement(rlid));
          auto I        = rowptr(rlid);
          for (LocalOrdinal k = 0; k < row.length; ++k) {
            if ((impl_ATS::magnitude(row.value(k)) > threshold) || (row.colidx(k) == rowInCol)) {
              idx(I)  = row.colidx(k);
              vals(I) = row.value(k);
              I += 1;
            }
          }
        });

    Kokkos::fence();
  } else {
    Kokkos::parallel_scan(
        "removeSmallEntries::rowptr",
        Kokkos::RangePolicy<LocalOrdinal, execution_space>(0, lclA.numRows()),
        KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {
          auto row = lclA.row(rlid);
          for (LocalOrdinal k = 0; k < row.length; ++k) {
            if (impl_ATS::magnitude(row.value(k)) > threshold) {
              partial_nnz += 1;
            }
          }
          if (is_final)
            rowptr(rlid + 1) = partial_nnz;
        },
        nnz);

    idx  = col_idx_type("idx", nnz);
    vals = vals_type("vals", nnz);

    Kokkos::parallel_for(
        "removeSmallEntries::indicesValues",
        Kokkos::RangePolicy<LocalOrdinal, execution_space>(0, lclA.numRows()),
        KOKKOS_LAMBDA(const LocalOrdinal rlid) {
          auto row = lclA.row(rlid);
          auto I   = rowptr(rlid);
          for (LocalOrdinal k = 0; k < row.length; ++k) {
            if (impl_ATS::magnitude(row.value(k)) > threshold) {
              idx(I)  = row.colidx(k);
              vals(I) = row.value(k);
              I += 1;
            }
          }
        });

    Kokkos::fence();
  }

  auto lclNewA = typename crs_matrix::local_matrix_type("thresholdedMatrix", lclA.numRows(), lclA.numCols(), nnz, vals, rowptr, idx);
  auto newA    = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclNewA, A->getRowMap(), A->getColMap(), A->getDomainMap(), A->getRangeMap());

  return newA;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetThresholdedMatrix(const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ain, const typename Teuchos::ScalarTraits<Scalar>::magnitudeType threshold, const bool keepDiagonal, const GlobalOrdinal expectedNNZperRow) {
  auto crsWrap = rcp_dynamic_cast<CrsMatrixWrap>(Ain);
  if (!crsWrap.is_null()) {
    auto crsMat      = crsWrap->getCrsMatrix();
    auto filteredMat = removeSmallEntries(crsMat, threshold, keepDiagonal);
    return rcp_static_cast<CrsMatrixWrap>(filteredMat);
  }

  RCP<const Map> rowmap   = Ain->getRowMap();
  RCP<const Map> colmap   = Ain->getColMap();
  RCP<CrsMatrixWrap> Aout = rcp(new CrsMatrixWrap(rowmap, expectedNNZperRow <= 0 ? Ain->getGlobalMaxNumRowEntries() : expectedNNZperRow));
  // loop over local rows
  for (size_t row = 0; row < Ain->getLocalNumRows(); row++) {
    size_t nnz = Ain->getNumEntriesInLocalRow(row);

    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    Ain->getLocalRowView(row, indices, vals);

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");

    Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(), Teuchos::ScalarTraits<GlobalOrdinal>::zero());
    Teuchos::ArrayRCP<Scalar> valout(indices.size(), Teuchos::ScalarTraits<Scalar>::zero());
    size_t nNonzeros = 0;
    if (keepDiagonal) {
      GlobalOrdinal glbRow   = rowmap->getGlobalElement(row);
      LocalOrdinal lclColIdx = colmap->getLocalElement(glbRow);
      for (size_t i = 0; i < (size_t)indices.size(); i++) {
        if (Teuchos::ScalarTraits<Scalar>::magnitude(vals[i]) > Teuchos::ScalarTraits<Scalar>::magnitude(threshold) || indices[i] == lclColIdx) {
          indout[nNonzeros] = colmap->getGlobalElement(indices[i]);  // LID -> GID (column)
          valout[nNonzeros] = vals[i];
          nNonzeros++;
        }
      }
    } else
      for (size_t i = 0; i < (size_t)indices.size(); i++) {
        if (Teuchos::ScalarTraits<Scalar>::magnitude(vals[i]) > Teuchos::ScalarTraits<Scalar>::magnitude(threshold)) {
          indout[nNonzeros] = colmap->getGlobalElement(indices[i]);  // LID -> GID (column)
          valout[nNonzeros] = vals[i];
          nNonzeros++;
        }
      }

    indout.resize(nNonzeros);
    valout.resize(nNonzeros);

    Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0, indout.size()), valout.view(0, valout.size()));
  }
  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

  return Aout;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetThresholdedGraph(const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A, const Magnitude threshold, const GlobalOrdinal expectedNNZperRow) {
  using STS                     = Teuchos::ScalarTraits<Scalar>;
  RCP<CrsGraph> sparsityPattern = CrsGraphFactory::Build(A->getRowMap(), expectedNNZperRow <= 0 ? A->getGlobalMaxNumRowEntries() : expectedNNZperRow);

  RCP<Vector> diag         = GetMatrixOverlappedDiagonal(*A);
  ArrayRCP<const Scalar> D = diag->getData(0);

  for (size_t row = 0; row < A->getLocalNumRows(); row++) {
    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    A->getLocalRowView(row, indices, vals);

    GlobalOrdinal globalRow = A->getRowMap()->getGlobalElement(row);
    LocalOrdinal col        = A->getColMap()->getLocalElement(globalRow);

    const Scalar Dk = STS::magnitude(D[col]) > 0.0 ? STS::magnitude(D[col]) : 1.0;
    Array<GlobalOrdinal> indicesNew;

    for (size_t i = 0; i < size_t(indices.size()); i++)
      // keep diagonal per default
      if (col == indices[i] || STS::magnitude(STS::squareroot(Dk) * vals[i] * STS::squareroot(Dk)) > STS::magnitude(threshold))
        indicesNew.append(A->getColMap()->getGlobalElement(indices[i]));

    sparsityPattern->insertGlobalIndices(globalRow, ArrayView<const GlobalOrdinal>(indicesNew.data(), indicesNew.length()));
  }
  sparsityPattern->fillComplete();

  return sparsityPattern;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixDiagonal_arcp(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  size_t numRows = A.getRowMap()->getLocalNumElements();
  Teuchos::ArrayRCP<Scalar> diag(numRows);
  Teuchos::ArrayView<const LocalOrdinal> cols;
  Teuchos::ArrayView<const Scalar> vals;
  for (size_t i = 0; i < numRows; ++i) {
    A.getLocalRowView(i, cols, vals);
    LocalOrdinal j = 0;
    for (; j < cols.size(); ++j) {
      if (Teuchos::as<size_t>(cols[j]) == i) {
        diag[i] = vals[j];
        break;
      }
    }
    if (j == cols.size()) {
      // Diagonal entry is absent
      diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }
  }
  return diag;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixDiagonal(const Matrix& A) {
  const auto rowMap = A.getRowMap();
  auto diag         = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, true);

  const CrsMatrixWrap* crsOp = dynamic_cast<const CrsMatrixWrap*>(&A);
  if ((crsOp != NULL) && (rowMap->lib() == Xpetra::UseTpetra)) {
    using device_type = typename CrsGraph::device_type;
    Kokkos::View<size_t*, device_type> offsets("offsets", rowMap->getLocalNumElements());
    crsOp->getCrsGraph()->getLocalDiagOffsets(offsets);
    crsOp->getCrsMatrix()->getLocalDiagCopy(*diag, offsets);
  } else {
    A.getLocalDiagCopy(*diag);
  }

  return diag;
}

// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
// RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
// UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
// GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & A, Magnitude tol, Scalar valReplacement) {
//   Teuchos::TimeMonitor MM = *Teuchos::TimeMonitor::getNewTimer("UtilitiesBase::GetMatrixDiagonalInverse");

//   RCP<const Map> rowMap = A.getRowMap();
//   RCP<Vector> diag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap,true);

//   A.getLocalDiagCopy(*diag);

//   RCP<Vector> inv = MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetInverse(diag, tol, valReplacement);

//   return inv;
// }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                             typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol,
                             Scalar valReplacement,
                             const bool doLumped) {
  Teuchos::TimeMonitor MM = *Teuchos::TimeMonitor::getNewTimer("Utilities::GetMatrixDiagonalInverse");

  RCP<const BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpFromRef(A));
  if (!bA.is_null()) {
    RCP<const Map> rowMap = A.getRowMap();
    RCP<Vector> diag      = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, true);
    A.getLocalDiagCopy(*diag);
    RCP<Vector> inv = MueLu::UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetInverse(diag, tol, valReplacement);
    return inv;
  }

  // Some useful type definitions
  using local_matrix_type = typename Matrix::local_matrix_type;
  // using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using value_type      = typename local_matrix_type::value_type;
  using values_type     = typename local_matrix_type::values_type;
  using scalar_type     = typename values_type::non_const_value_type;
  using ordinal_type    = typename local_matrix_type::ordinal_type;
  using execution_space = typename local_matrix_type::execution_space;
  // using memory_space      = typename local_matrix_type::memory_space;
  // Be careful with this one, if using Kokkos::ArithTraits<Scalar>
  // you are likely to run into errors when handling std::complex<>
  // a good way to work around that is to use the following:
  // using KAT = Kokkos::ArithTraits<Kokkos::ArithTraits<Scalar>::val_type> >
  // here we have: value_type = Kokkos::ArithTraits<Scalar>::val_type
  using KAT = Kokkos::ArithTraits<value_type>;

  // Get/Create distributed objects
  RCP<const Map> rowMap = A.getRowMap();
  RCP<Vector> diag      = VectorFactory::Build(rowMap, false);

  // Now generate local objects
  local_matrix_type localMatrix = A.getLocalMatrixDevice();
  auto diagVals                 = diag->getDeviceLocalView(Xpetra::Access::ReadWrite);

  ordinal_type numRows = localMatrix.graph.numRows();

  scalar_type valReplacement_dev = valReplacement;

  // Note: 2019-11-21, LBV
  // This could be implemented with a TeamPolicy over the rows
  // and a TeamVectorRange over the entries in a row if performance
  // becomes more important here.
  if (!doLumped)
    Kokkos::parallel_for(
        "Utilities::GetMatrixDiagonalInverse",
        Kokkos::RangePolicy<ordinal_type, execution_space>(0, numRows),
        KOKKOS_LAMBDA(const ordinal_type rowIdx) {
          bool foundDiagEntry = false;
          auto myRow          = localMatrix.rowConst(rowIdx);
          for (ordinal_type entryIdx = 0; entryIdx < myRow.length; ++entryIdx) {
            if (myRow.colidx(entryIdx) == rowIdx) {
              foundDiagEntry = true;
              if (KAT::magnitude(myRow.value(entryIdx)) > KAT::magnitude(tol)) {
                diagVals(rowIdx, 0) = KAT::one() / myRow.value(entryIdx);
              } else {
                diagVals(rowIdx, 0) = valReplacement_dev;
              }
              break;
            }
          }

          if (!foundDiagEntry) {
            diagVals(rowIdx, 0) = KAT::zero();
          }
        });
  else
    Kokkos::parallel_for(
        "Utilities::GetMatrixDiagonalInverse",
        Kokkos::RangePolicy<ordinal_type, execution_space>(0, numRows),
        KOKKOS_LAMBDA(const ordinal_type rowIdx) {
          auto myRow = localMatrix.rowConst(rowIdx);
          for (ordinal_type entryIdx = 0; entryIdx < myRow.length; ++entryIdx) {
            diagVals(rowIdx, 0) += KAT::magnitude(myRow.value(entryIdx));
          }
          if (KAT::magnitude(diagVals(rowIdx, 0)) > KAT::magnitude(tol))
            diagVals(rowIdx, 0) = KAT::one() / diagVals(rowIdx, 0);
          else
            diagVals(rowIdx, 0) = valReplacement_dev;
        });

  return diag;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetLumpedMatrixDiagonal(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> const& A, const bool doReciprocal,
                            Magnitude tol,
                            Scalar valReplacement,
                            const bool replaceSingleEntryRowWithZero,
                            const bool useAverageAbsDiagVal) {
  typedef Teuchos::ScalarTraits<Scalar> TST;

  RCP<Vector> diag  = Teuchos::null;
  const Scalar zero = TST::zero();
  const Scalar one  = TST::one();
  const Scalar two  = one + one;

  Teuchos::RCP<const Matrix> rcpA = Teuchos::rcpFromRef(A);

  RCP<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bA =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(rcpA);
  if (bA == Teuchos::null) {
    RCP<const Map> rowMap = rcpA->getRowMap();
    diag                  = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, true);

    if (rowMap->lib() == Xpetra::UnderlyingLib::UseTpetra) {
      Teuchos::TimeMonitor MM = *Teuchos::TimeMonitor::getNewTimer("UtilitiesBase::GetLumpedMatrixDiagonal (Kokkos implementation)");
      // Implement using Kokkos
      using local_vector_type = typename Vector::dual_view_type::t_dev_um;
      using local_matrix_type = typename Matrix::local_matrix_type;
      using execution_space   = typename local_vector_type::execution_space;
      // using rowmap_type       = typename local_matrix_type::row_map_type;
      // using entries_type      = typename local_matrix_type::index_type;
      using values_type = typename local_matrix_type::values_type;
      using scalar_type = typename values_type::non_const_value_type;
      using mag_type    = typename Kokkos::ArithTraits<scalar_type>::mag_type;
      using KAT_S       = typename Kokkos::ArithTraits<scalar_type>;
      using KAT_M       = typename Kokkos::ArithTraits<mag_type>;
      using size_type   = typename local_matrix_type::non_const_size_type;

      local_vector_type diag_dev      = diag->getDeviceLocalView(Xpetra::Access::OverwriteAll);
      local_matrix_type local_mat_dev = rcpA->getLocalMatrixDevice();
      Kokkos::RangePolicy<execution_space, int> my_policy(0, static_cast<int>(diag_dev.extent(0)));
      scalar_type valReplacement_dev = valReplacement;

      if (doReciprocal) {
        Kokkos::View<int*, execution_space> nnzPerRow("nnz per rows", diag_dev.extent(0));
        Kokkos::View<scalar_type*, execution_space> regSum("regSum", diag_dev.extent(0));
        Kokkos::View<mag_type, execution_space> avgAbsDiagVal_dev("avgAbsDiagVal");
        Kokkos::View<int, execution_space> numDiagsEqualToOne_dev("numDiagsEqualToOne");

        {
          Teuchos::TimeMonitor MMM = *Teuchos::TimeMonitor::getNewTimer("GetLumpedMatrixDiagonal: parallel_for (doReciprocal)");
          Kokkos::parallel_for(
              "GetLumpedMatrixDiagonal", my_policy,
              KOKKOS_LAMBDA(const int rowIdx) {
                diag_dev(rowIdx, 0) = KAT_S::zero();
                for (size_type entryIdx = local_mat_dev.graph.row_map(rowIdx);
                     entryIdx < local_mat_dev.graph.row_map(rowIdx + 1);
                     ++entryIdx) {
                  regSum(rowIdx) += local_mat_dev.values(entryIdx);
                  if (KAT_M::zero() < KAT_S::abs(local_mat_dev.values(entryIdx))) {
                    ++nnzPerRow(rowIdx);
                  }
                  diag_dev(rowIdx, 0) += KAT_S::abs(local_mat_dev.values(entryIdx));
                  if (rowIdx == local_mat_dev.graph.entries(entryIdx)) {
                    Kokkos::atomic_add(&avgAbsDiagVal_dev(), KAT_S::abs(local_mat_dev.values(entryIdx)));
                  }
                }

                if (nnzPerRow(rowIdx) == 1 && KAT_S::magnitude(diag_dev(rowIdx, 0)) == KAT_M::one()) {
                  Kokkos::atomic_add(&numDiagsEqualToOne_dev(), 1);
                }
              });
        }
        if (useAverageAbsDiagVal) {
          Teuchos::TimeMonitor MMM                                                   = *Teuchos::TimeMonitor::getNewTimer("GetLumpedMatrixDiagonal: useAverageAbsDiagVal");
          typename Kokkos::View<mag_type, execution_space>::HostMirror avgAbsDiagVal = Kokkos::create_mirror_view(avgAbsDiagVal_dev);
          Kokkos::deep_copy(avgAbsDiagVal, avgAbsDiagVal_dev);
          int numDiagsEqualToOne;
          Kokkos::deep_copy(numDiagsEqualToOne, numDiagsEqualToOne_dev);

          tol = TST::magnitude(100 * Teuchos::ScalarTraits<Scalar>::eps()) * (avgAbsDiagVal() - numDiagsEqualToOne) / (rowMap->getLocalNumElements() - numDiagsEqualToOne);
        }

        {
          Teuchos::TimeMonitor MMM = *Teuchos::TimeMonitor::getNewTimer("ComputeLumpedDiagonalInverse: parallel_for (doReciprocal)");
          Kokkos::parallel_for(
              "ComputeLumpedDiagonalInverse", my_policy,
              KOKKOS_LAMBDA(const int rowIdx) {
                if (replaceSingleEntryRowWithZero && nnzPerRow(rowIdx) <= 1) {
                  diag_dev(rowIdx, 0) = KAT_S::zero();
                } else if ((diag_dev(rowIdx, 0) != KAT_S::zero()) && (KAT_S::magnitude(diag_dev(rowIdx, 0)) < KAT_S::magnitude(2 * regSum(rowIdx)))) {
                  diag_dev(rowIdx, 0) = KAT_S::one() / KAT_S::magnitude(2 * regSum(rowIdx));
                } else {
                  if (KAT_S::magnitude(diag_dev(rowIdx, 0)) > tol) {
                    diag_dev(rowIdx, 0) = KAT_S::one() / diag_dev(rowIdx, 0);
                  } else {
                    diag_dev(rowIdx, 0) = valReplacement_dev;
                  }
                }
              });
        }

      } else {
        Teuchos::TimeMonitor MMM = *Teuchos::TimeMonitor::getNewTimer("GetLumpedMatrixDiagonal: parallel_for");
        Kokkos::parallel_for(
            "GetLumpedMatrixDiagonal", my_policy,
            KOKKOS_LAMBDA(const int rowIdx) {
              diag_dev(rowIdx, 0) = KAT_S::zero();
              for (size_type entryIdx = local_mat_dev.graph.row_map(rowIdx);
                   entryIdx < local_mat_dev.graph.row_map(rowIdx + 1);
                   ++entryIdx) {
                diag_dev(rowIdx, 0) += KAT_S::magnitude(local_mat_dev.values(entryIdx));
              }
            });
      }
    } else {
      // Implement using Teuchos
      Teuchos::TimeMonitor MMM  = *Teuchos::TimeMonitor::getNewTimer("UtilitiesBase: GetLumpedMatrixDiagonal: (Teuchos implementation)");
      ArrayRCP<Scalar> diagVals = diag->getDataNonConst(0);
      Teuchos::Array<Scalar> regSum(diag->getLocalLength());
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;

      std::vector<int> nnzPerRow(rowMap->getLocalNumElements());

      // FIXME 2021-10-22 JHU   If this is called with doReciprocal=false, what should the correct behavior be?  Currently,
      // FIXME 2021-10-22 JHU   the diagonal entry is set to be the sum of the absolute values of the row entries.

      const Magnitude zeroMagn = TST::magnitude(zero);
      Magnitude avgAbsDiagVal  = TST::magnitude(zero);
      int numDiagsEqualToOne   = 0;
      for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
        nnzPerRow[i] = 0;
        rcpA->getLocalRowView(i, cols, vals);
        diagVals[i] = zero;
        for (LocalOrdinal j = 0; j < cols.size(); ++j) {
          regSum[i] += vals[j];
          const Magnitude rowEntryMagn = TST::magnitude(vals[j]);
          if (rowEntryMagn > zeroMagn)
            nnzPerRow[i]++;
          diagVals[i] += rowEntryMagn;
          if (static_cast<size_t>(cols[j]) == i)
            avgAbsDiagVal += rowEntryMagn;
        }
        if (nnzPerRow[i] == 1 && TST::magnitude(diagVals[i]) == 1.)
          numDiagsEqualToOne++;
      }
      if (useAverageAbsDiagVal)
        tol = TST::magnitude(100 * Teuchos::ScalarTraits<Scalar>::eps()) * (avgAbsDiagVal - numDiagsEqualToOne) / (rowMap->getLocalNumElements() - numDiagsEqualToOne);
      if (doReciprocal) {
        for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
          if (replaceSingleEntryRowWithZero && nnzPerRow[i] <= static_cast<int>(1))
            diagVals[i] = zero;
          else if ((diagVals[i] != zero) && (TST::magnitude(diagVals[i]) < TST::magnitude(two * regSum[i])))
            diagVals[i] = one / TST::magnitude((two * regSum[i]));
          else {
            if (TST::magnitude(diagVals[i]) > tol)
              diagVals[i] = one / diagVals[i];
            else {
              diagVals[i] = valReplacement;
            }
          }
        }
      }
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(doReciprocal, Xpetra::Exceptions::RuntimeError,
                               "UtilitiesBase::GetLumpedMatrixDiagonal(): extracting reciprocal of diagonal of a blocked matrix is not supported");
    diag = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(bA->getRangeMapExtractor()->getFullMap(), true);

    for (size_t row = 0; row < bA->Rows(); ++row) {
      for (size_t col = 0; col < bA->Cols(); ++col) {
        if (!bA->getMatrix(row, col).is_null()) {
          // if we are in Thyra mode, but the block (row,row) is again a blocked operator, we have to use (pseudo) Xpetra-style GIDs with offset!
          bool bThyraMode      = bA->getRangeMapExtractor()->getThyraMode() && (Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(bA->getMatrix(row, col)) == Teuchos::null);
          RCP<Vector> ddtemp   = bA->getRangeMapExtractor()->ExtractVector(diag, row, bThyraMode);
          RCP<const Vector> dd = GetLumpedMatrixDiagonal(*(bA->getMatrix(row, col)));
          ddtemp->update(Teuchos::as<Scalar>(1.0), *dd, Teuchos::as<Scalar>(1.0));
          bA->getRangeMapExtractor()->InsertVector(ddtemp, row, diag, bThyraMode);
        }
      }
    }
  }

  return diag;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixMaxMinusOffDiagonal(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  // Get/Create distributed objects
  RCP<const Map> rowMap = A.getRowMap();
  auto diag             = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, false);

  // Implement using Kokkos
  using local_vector_type = typename Vector::dual_view_type::t_dev_um;
  using local_matrix_type = typename Matrix::local_matrix_type;
  using execution_space   = typename local_vector_type::execution_space;
  using values_type       = typename local_matrix_type::values_type;
  using scalar_type       = typename values_type::non_const_value_type;
  using KAT_S             = typename Kokkos::ArithTraits<scalar_type>;

  auto diag_dev      = diag->getDeviceLocalView(Xpetra::Access::OverwriteAll);
  auto local_mat_dev = A.getLocalMatrixDevice();
  Kokkos::RangePolicy<execution_space, int> my_policy(0, static_cast<int>(diag_dev.extent(0)));

  Kokkos::parallel_for(
      "GetMatrixMaxMinusOffDiagonal", my_policy,
      KOKKOS_LAMBDA(const LocalOrdinal rowIdx) {
        auto mymax = KAT_S::zero();
        auto row   = local_mat_dev.rowConst(rowIdx);
        for (LocalOrdinal entryIdx = 0; entryIdx < row.length; ++entryIdx) {
          if (rowIdx != row.colidx(entryIdx)) {
            if (KAT_S::real(mymax) < -KAT_S::real(row.value(entryIdx)))
              mymax = -KAT_S::real(row.value(entryIdx));
          }
        }
        diag_dev(rowIdx, 0) = mymax;
      });

  return diag;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixMaxMinusOffDiagonal(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber) {
  TEUCHOS_TEST_FOR_EXCEPTION(!A.getColMap()->isSameAs(*BlockNumber.getMap()), std::runtime_error, "GetMatrixMaxMinusOffDiagonal: BlockNumber must match's A's column map.");

  // Get/Create distributed objects
  RCP<const Map> rowMap = A.getRowMap();
  auto diag             = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, false);

  // Implement using Kokkos
  using local_vector_type = typename Vector::dual_view_type::t_dev_um;
  using local_matrix_type = typename Matrix::local_matrix_type;
  using execution_space   = typename local_vector_type::execution_space;
  using values_type       = typename local_matrix_type::values_type;
  using scalar_type       = typename values_type::non_const_value_type;
  using KAT_S             = typename Kokkos::ArithTraits<scalar_type>;

  auto diag_dev        = diag->getDeviceLocalView(Xpetra::Access::OverwriteAll);
  auto local_mat_dev   = A.getLocalMatrixDevice();
  auto local_block_dev = BlockNumber.getDeviceLocalView(Xpetra::Access::ReadOnly);
  Kokkos::RangePolicy<execution_space, int> my_policy(0, static_cast<int>(diag_dev.extent(0)));

  Kokkos::parallel_for(
      "GetMatrixMaxMinusOffDiagonal", my_policy,
      KOKKOS_LAMBDA(const LocalOrdinal rowIdx) {
        auto mymax = KAT_S::zero();
        auto row   = local_mat_dev.row(rowIdx);
        for (LocalOrdinal entryIdx = 0; entryIdx < row.length; ++entryIdx) {
          if ((rowIdx != row.colidx(entryIdx)) && (local_block_dev(rowIdx, 0) == local_block_dev(row.colidx(entryIdx), 0))) {
            if (KAT_S::real(mymax) < -KAT_S::real(row.value(entryIdx)))
              mymax = -KAT_S::real(row.value(entryIdx));
          }
        }
        diag_dev(rowIdx, 0) = mymax;
      });

  return diag;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetInverse(Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol, Scalar valReplacement) {
  RCP<Vector> ret = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(v->getMap(), true);

  // check whether input vector "v" is a BlockedVector
  RCP<const BlockedVector> bv = Teuchos::rcp_dynamic_cast<const BlockedVector>(v);
  if (bv.is_null() == false) {
    RCP<BlockedVector> bret = Teuchos::rcp_dynamic_cast<BlockedVector>(ret);
    TEUCHOS_TEST_FOR_EXCEPTION(bret.is_null() == true, MueLu::Exceptions::RuntimeError, "MueLu::UtilitiesBase::GetInverse: return vector should be of type BlockedVector");
    RCP<const BlockedMap> bmap = bv->getBlockedMap();
    for (size_t r = 0; r < bmap->getNumMaps(); ++r) {
      RCP<const MultiVector> submvec = bv->getMultiVector(r, bmap->getThyraMode());
      RCP<const Vector> subvec       = submvec->getVector(0);
      RCP<Vector> subvecinf          = MueLu::UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetInverse(subvec, tol, valReplacement);
      bret->setMultiVector(r, subvecinf, bmap->getThyraMode());
    }
    return ret;
  }

  // v is an {Epetra,Tpetra}Vector: work with the underlying raw data
  ArrayRCP<Scalar> retVals         = ret->getDataNonConst(0);
  ArrayRCP<const Scalar> inputVals = v->getData(0);
  for (size_t i = 0; i < v->getMap()->getLocalNumElements(); ++i) {
    if (Teuchos::ScalarTraits<Scalar>::magnitude(inputVals[i]) > tol)
      retVals[i] = Teuchos::ScalarTraits<Scalar>::one() / inputVals[i];
    else
      retVals[i] = valReplacement;
  }
  return ret;
}

// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
// RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
// UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
// GetMatrixOverlappedDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & A) {
//   RCP<const Map> rowMap = A.getRowMap(), colMap = A.getColMap();

//   // Undo block map (if we have one)
//   RCP<const BlockedMap> browMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rowMap);
//   if(!browMap.is_null()) rowMap = browMap->getMap();

//   RCP<Vector> localDiag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap);
//   try {
//     const CrsMatrixWrap* crsOp = dynamic_cast<const CrsMatrixWrap*>(&A);
//     if (crsOp == NULL) {
//       throw Exceptions::RuntimeError("cast to CrsMatrixWrap failed");
//     }
//     Teuchos::ArrayRCP<size_t> offsets;
//     crsOp->getLocalDiagOffsets(offsets);
//     crsOp->getLocalDiagCopy(*localDiag,offsets());
//   }
//   catch (...) {
//     ArrayRCP<Scalar>   localDiagVals = localDiag->getDataNonConst(0);
//     Teuchos::ArrayRCP<Scalar> diagVals = GetMatrixDiagonal(A);
//     for (LocalOrdinal i = 0; i < localDiagVals.size(); i++)
//       localDiagVals[i] = diagVals[i];
//     localDiagVals = diagVals = null;
//   }

//   RCP<Vector> diagonal = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colMap);
//   RCP< const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer;
//   importer = A.getCrsGraph()->getImporter();
//   if (importer == Teuchos::null) {
//     importer = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap, colMap);
//   }
//   diagonal->doImport(*localDiag, *(importer), Xpetra::INSERT);
//   return diagonal;
// }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixOverlappedDiagonal(const Matrix& A) {
  RCP<const Map> rowMap = A.getRowMap(), colMap = A.getColMap();
  RCP<Vector> localDiag      = GetMatrixDiagonal(A);
  RCP<Vector> diagonal       = VectorFactory::Build(colMap);
  RCP<const Import> importer = A.getCrsGraph()->getImporter();
  if (importer == Teuchos::null) {
    importer = ImportFactory::Build(rowMap, colMap);
  }
  diagonal->doImport(*localDiag, *(importer), Xpetra::INSERT);

  return diagonal;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixOverlappedDeletedRowsum(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  using STS = typename Teuchos::ScalarTraits<SC>;

  // Undo block map (if we have one)
  RCP<const Map> rowMap = A.getRowMap(), colMap = A.getColMap();
  RCP<const BlockedMap> browMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rowMap);
  if (!browMap.is_null()) rowMap = browMap->getMap();

  RCP<Vector> local      = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(rowMap);
  RCP<Vector> ghosted    = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(colMap, true);
  ArrayRCP<SC> localVals = local->getDataNonConst(0);

  for (LO row = 0; row < static_cast<LO>(A.getRowMap()->getLocalNumElements()); ++row) {
    size_t nnz = A.getNumEntriesInLocalRow(row);
    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    A.getLocalRowView(row, indices, vals);

    SC si = STS::zero();

    for (LO colID = 0; colID < static_cast<LO>(nnz); colID++) {
      if (indices[colID] != row) {
        si += vals[colID];
      }
    }
    localVals[row] = si;
  }

  RCP<const Xpetra::Import<LO, GO, Node>> importer;
  importer = A.getCrsGraph()->getImporter();
  if (importer == Teuchos::null) {
    importer = Xpetra::ImportFactory<LO, GO, Node>::Build(rowMap, colMap);
  }
  ghosted->doImport(*local, *(importer), Xpetra::INSERT);
  return ghosted;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixOverlappedAbsDeletedRowsum(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  RCP<const Map> rowMap = A.getRowMap(), colMap = A.getColMap();
  using STS              = typename Teuchos::ScalarTraits<Scalar>;
  using MTS              = typename Teuchos::ScalarTraits<Magnitude>;
  using MT               = Magnitude;
  using RealValuedVector = Xpetra::Vector<MT, LO, GO, Node>;

  // Undo block map (if we have one)
  RCP<const BlockedMap> browMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rowMap);
  if (!browMap.is_null()) rowMap = browMap->getMap();

  RCP<RealValuedVector> local   = Xpetra::VectorFactory<MT, LO, GO, Node>::Build(rowMap);
  RCP<RealValuedVector> ghosted = Xpetra::VectorFactory<MT, LO, GO, Node>::Build(colMap, true);
  ArrayRCP<MT> localVals        = local->getDataNonConst(0);

  for (LO rowIdx = 0; rowIdx < static_cast<LO>(A.getRowMap()->getLocalNumElements()); ++rowIdx) {
    size_t nnz = A.getNumEntriesInLocalRow(rowIdx);
    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    A.getLocalRowView(rowIdx, indices, vals);

    MT si = MTS::zero();

    for (LO colID = 0; colID < static_cast<LO>(nnz); ++colID) {
      if (indices[colID] != rowIdx) {
        si += STS::magnitude(vals[colID]);
      }
    }
    localVals[rowIdx] = si;
  }

  RCP<const Xpetra::Import<LO, GO, Node>> importer;
  importer = A.getCrsGraph()->getImporter();
  if (importer == Teuchos::null) {
    importer = Xpetra::ImportFactory<LO, GO, Node>::Build(rowMap, colMap);
  }
  ghosted->doImport(*local, *(importer), Xpetra::INSERT);
  return ghosted;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ResidualNorm(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS) {
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
  const size_t numVecs = X.getNumVectors();
  RCP<MultiVector> RES = Residual(Op, X, RHS);
  Teuchos::Array<Magnitude> norms(numVecs);
  RES->norm2(norms);
  return norms;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ResidualNorm(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS, MultiVector& Resid) {
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
  const size_t numVecs = X.getNumVectors();
  Residual(Op, X, RHS, Resid);
  Teuchos::Array<Magnitude> norms(numVecs);
  Resid.norm2(norms);
  return norms;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Residual(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS) {
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
  const size_t numVecs = X.getNumVectors();
  // TODO Op.getRangeMap should return a BlockedMap if it is a BlockedCrsOperator
  RCP<MultiVector> RES = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(RHS.getMap(), numVecs, false);  // no need to initialize to zero
  Op.residual(X, RHS, *RES);
  return RES;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Residual(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS, MultiVector& Resid) {
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides");
  TEUCHOS_TEST_FOR_EXCEPTION(Resid.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of residual vectors != number of right-hand sides");
  Op.residual(X, RHS, Resid);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    PowerMethod(const Matrix& A, bool scaleByDiag,
                LocalOrdinal niters, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance, bool verbose, unsigned int seed) {
  TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRangeMap()->isSameAs(*(A.getDomainMap()))), Exceptions::Incompatible,
                             "Utils::PowerMethod: operator must have domain and range maps that are equivalent.");

  // power iteration
  RCP<Vector> diagInvVec;
  if (scaleByDiag) {
    diagInvVec = GetMatrixDiagonalInverse(A);
  }

  Scalar lambda = PowerMethod(A, diagInvVec, niters, tolerance, verbose, seed);
  return lambda;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    PowerMethod(const Matrix& A, const RCP<Vector>& diagInvVec,
                LocalOrdinal niters, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance, bool verbose, unsigned int seed) {
  TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRangeMap()->isSameAs(*(A.getDomainMap()))), Exceptions::Incompatible,
                             "Utils::PowerMethod: operator must have domain and range maps that are equivalent.");

  // Create three vectors, fill z with random numbers
  RCP<Vector> q = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A.getDomainMap());
  RCP<Vector> r = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A.getRangeMap());
  RCP<Vector> z = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A.getRangeMap());

  z->setSeed(seed);    // seed random number generator
  z->randomize(true);  // use Xpetra implementation: -> same results for Epetra and Tpetra

  Teuchos::Array<Magnitude> norms(1);

  typedef Teuchos::ScalarTraits<Scalar> STS;

  const Scalar zero = STS::zero(), one = STS::one();

  Scalar lambda      = zero;
  Magnitude residual = STS::magnitude(zero);

  // power iteration
  for (int iter = 0; iter < niters; ++iter) {
    z->norm2(norms);                      // Compute 2-norm of z
    q->update(one / norms[0], *z, zero);  // Set q = z / normz
    A.apply(*q, *z);                      // Compute z = A*q
    if (diagInvVec != Teuchos::null)
      z->elementWiseMultiply(one, *diagInvVec, *z, zero);
    lambda = q->dot(*z);  // Approximate maximum eigenvalue: lamba = dot(q,z)

    if (iter % 100 == 0 || iter + 1 == niters) {
      r->update(1.0, *z, -lambda, *q, zero);  // Compute A*q - lambda*q
      r->norm2(norms);
      residual = STS::magnitude(norms[0] / lambda);
      if (verbose) {
        std::cout << "Iter = " << iter
                  << "  Lambda = " << lambda
                  << "  Residual of A*q - lambda*q = " << residual
                  << std::endl;
      }
    }
    if (residual < tolerance)
      break;
  }
  return lambda;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Teuchos::FancyOStream>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MakeFancy(std::ostream& os) {
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
  return fancy;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Distance2(const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1) {
  const size_t numVectors = v.size();

  Scalar d = Teuchos::ScalarTraits<Scalar>::zero();
  for (size_t j = 0; j < numVectors; j++) {
    d += (v[j][i0] - v[j][i1]) * (v[j][i0] - v[j][i1]);
  }
  return Teuchos::ScalarTraits<Scalar>::magnitude(d);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Distance2(const Teuchos::ArrayView<double>& weight, const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1) {
  const size_t numVectors = v.size();
  using MT                = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  Scalar d = Teuchos::ScalarTraits<Scalar>::zero();
  for (size_t j = 0; j < numVectors; j++) {
    d += Teuchos::as<MT>(weight[j]) * (v[j][i0] - v[j][i1]) * (v[j][i0] - v[j][i1]);
  }
  return Teuchos::ScalarTraits<Scalar>::magnitude(d);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const bool>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DetectDirichletRows(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& tol, bool count_twos_as_dirichlet) {
  LocalOrdinal numRows = A.getLocalNumRows();
  typedef Teuchos::ScalarTraits<Scalar> STS;
  ArrayRCP<bool> boundaryNodes(numRows, true);
  if (count_twos_as_dirichlet) {
    for (LocalOrdinal row = 0; row < numRows; row++) {
      ArrayView<const LocalOrdinal> indices;
      ArrayView<const Scalar> vals;
      A.getLocalRowView(row, indices, vals);
      size_t nnz = A.getNumEntriesInLocalRow(row);
      if (nnz > 2) {
        size_t col;
        for (col = 0; col < nnz; col++)
          if ((indices[col] != row) && STS::magnitude(vals[col]) > tol) {
            if (!boundaryNodes[row])
              break;
            boundaryNodes[row] = false;
          }
        if (col == nnz)
          boundaryNodes[row] = true;
      }
    }
  } else {
    for (LocalOrdinal row = 0; row < numRows; row++) {
      ArrayView<const LocalOrdinal> indices;
      ArrayView<const Scalar> vals;
      A.getLocalRowView(row, indices, vals);
      size_t nnz = A.getNumEntriesInLocalRow(row);
      if (nnz > 1)
        for (size_t col = 0; col < nnz; col++)
          if ((indices[col] != row) && STS::magnitude(vals[col]) > tol) {
            boundaryNodes[row] = false;
            break;
          }
    }
  }
  return boundaryNodes;
}

template <class SC, class LO, class GO, class NO, class memory_space>
Kokkos::View<bool*, memory_space>
DetectDirichletRows_kokkos(const Xpetra::Matrix<SC, LO, GO, NO>& A,
                           const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol,
                           const bool count_twos_as_dirichlet) {
  using impl_scalar_type = typename Kokkos::ArithTraits<SC>::val_type;
  using ATS              = Kokkos::ArithTraits<impl_scalar_type>;
  using range_type       = Kokkos::RangePolicy<LO, typename NO::execution_space>;
  using helpers          = Xpetra::Helpers<SC, LO, GO, NO>;

  Kokkos::View<bool*, typename NO::device_type::memory_space> boundaryNodes;

  if (helpers::isTpetraBlockCrs(A)) {
    const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& Am = helpers::Op2TpetraBlockCrs(A);
    auto b_graph                                     = Am.getCrsGraph().getLocalGraphDevice();
    auto b_rowptr                                    = Am.getCrsGraph().getLocalRowPtrsDevice();
    auto values                                      = Am.getValuesDevice();
    LO numBlockRows                                  = Am.getLocalNumRows();
    const LO stride                                  = Am.getBlockSize() * Am.getBlockSize();

    boundaryNodes = Kokkos::View<bool*, typename NO::device_type::memory_space>(Kokkos::ViewAllocateWithoutInitializing("boundaryNodes"), numBlockRows);

    if (count_twos_as_dirichlet)
      throw Exceptions::RuntimeError("BlockCrs does not support counting twos as Dirichlet");

    Kokkos::parallel_for(
        "MueLu:Utils::DetectDirichletRowsBlockCrs", range_type(0, numBlockRows),
        KOKKOS_LAMBDA(const LO row) {
          auto rowView = b_graph.rowConst(row);
          auto length  = rowView.length;
          LO valstart  = b_rowptr[row] * stride;

          boundaryNodes(row)     = true;
          decltype(length) colID = 0;
          for (; colID < length; colID++) {
            if (rowView.colidx(colID) != row) {
              LO current = valstart + colID * stride;
              for (LO k = 0; k < stride; k++) {
                if (ATS::magnitude(values[current + k]) > tol) {
                  boundaryNodes(row) = false;
                  break;
                }
              }
            }
            if (boundaryNodes(row) == false)
              break;
          }
        });
  } else {
    auto localMatrix = A.getLocalMatrixDevice();
    LO numRows       = A.getLocalNumRows();
    boundaryNodes    = Kokkos::View<bool*, typename NO::device_type::memory_space>(Kokkos::ViewAllocateWithoutInitializing("boundaryNodes"), numRows);

    if (count_twos_as_dirichlet)
      Kokkos::parallel_for(
          "MueLu:Utils::DetectDirichletRows_Twos_As_Dirichlet", range_type(0, numRows),
          KOKKOS_LAMBDA(const LO row) {
            auto rowView = localMatrix.row(row);
            auto length  = rowView.length;

            boundaryNodes(row) = true;
            if (length > 2) {
              decltype(length) colID = 0;
              for (; colID < length; colID++)
                if ((rowView.colidx(colID) != row) &&
                    (ATS::magnitude(rowView.value(colID)) > tol)) {
                  if (!boundaryNodes(row))
                    break;
                  boundaryNodes(row) = false;
                }
              if (colID == length)
                boundaryNodes(row) = true;
            }
          });
    else
      Kokkos::parallel_for(
          "MueLu:Utils::DetectDirichletRows", range_type(0, numRows),
          KOKKOS_LAMBDA(const LO row) {
            auto rowView = localMatrix.row(row);
            auto length  = rowView.length;

            boundaryNodes(row) = true;
            for (decltype(length) colID = 0; colID < length; colID++)
              if ((rowView.colidx(colID) != row) &&
                  (ATS::magnitude(rowView.value(colID)) > tol)) {
                boundaryNodes(row) = false;
                break;
              }
          });
  }
  if constexpr (std::is_same<memory_space, typename NO::device_type::memory_space>::value)
    return boundaryNodes;
  else {
    Kokkos::View<bool*, memory_space> boundaryNodes2(Kokkos::ViewAllocateWithoutInitializing("boundaryNodes"), boundaryNodes.extent(0));
    Kokkos::deep_copy(boundaryNodes2, boundaryNodes);
    return boundaryNodes2;
  }
  // CAG: No idea why this is needed to avoid "warning: missing return statement at end of non-void function"
  Kokkos::View<bool*, memory_space> dummy("dummy", 0);
  return dummy;
}

template <class SC, class LO, class GO, class NO>
Kokkos::View<bool*, typename NO::device_type::memory_space>
UtilitiesBase<SC, LO, GO, NO>::
    DetectDirichletRows_kokkos(const Xpetra::Matrix<SC, LO, GO, NO>& A,
                               const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol,
                               const bool count_twos_as_dirichlet) {
  return MueLu::DetectDirichletRows_kokkos<SC, LO, GO, NO, typename NO::device_type::memory_space>(A, tol, count_twos_as_dirichlet);
}

template <class SC, class LO, class GO, class NO>
Kokkos::View<bool*, typename Kokkos::HostSpace>
UtilitiesBase<SC, LO, GO, NO>::
    DetectDirichletRows_kokkos_host(const Xpetra::Matrix<SC, LO, GO, NO>& A,
                                    const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol,
                                    const bool count_twos_as_dirichlet) {
  return MueLu::DetectDirichletRows_kokkos<SC, LO, GO, NO, typename Kokkos::HostSpace>(A, tol, count_twos_as_dirichlet);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const bool>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DetectDirichletRowsExt(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, bool& bHasZeroDiagonal, const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& tol) {
  // assume that there is no zero diagonal in matrix
  bHasZeroDiagonal = false;

  Teuchos::RCP<Vector> diagVec = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A.getRowMap());
  A.getLocalDiagCopy(*diagVec);
  Teuchos::ArrayRCP<const Scalar> diagVecData = diagVec->getData(0);

  LocalOrdinal numRows = A.getLocalNumRows();
  typedef Teuchos::ScalarTraits<Scalar> STS;
  ArrayRCP<bool> boundaryNodes(numRows, false);
  for (LocalOrdinal row = 0; row < numRows; row++) {
    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    A.getLocalRowView(row, indices, vals);
    size_t nnz    = 0;  // collect nonzeros in row (excluding the diagonal)
    bool bHasDiag = false;
    for (decltype(indices.size()) col = 0; col < indices.size(); col++) {
      if (indices[col] != row) {
        if (STS::magnitude(vals[col] / STS::magnitude(sqrt(STS::magnitude(diagVecData[row]) * STS::magnitude(diagVecData[col])))) > tol) {
          nnz++;
        }
      } else
        bHasDiag = true;  // found a diagonal entry
    }
    if (bHasDiag == false)
      bHasZeroDiagonal = true;  // we found at least one row without a diagonal
    else if (nnz == 0)
      boundaryNodes[row] = true;
  }
  return boundaryNodes;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FindNonZeros(const Teuchos::ArrayRCP<const Scalar> vals,
                 Teuchos::ArrayRCP<bool> nonzeros) {
  TEUCHOS_ASSERT(vals.size() == nonzeros.size());
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  const magnitudeType eps = 2.0 * Teuchos::ScalarTraits<magnitudeType>::eps();
  for (size_t i = 0; i < static_cast<size_t>(vals.size()); i++) {
    nonzeros[i] = (Teuchos::ScalarTraits<Scalar>::magnitude(vals[i]) > eps);
  }
}

// Find Nonzeros in a device view
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FindNonZeros(const typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev_const_um vals,
                 Kokkos::View<bool*, typename Node::device_type> nonzeros) {
  using ATS        = Kokkos::ArithTraits<Scalar>;
  using impl_ATS   = Kokkos::ArithTraits<typename ATS::val_type>;
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;
  TEUCHOS_ASSERT(vals.extent(0) == nonzeros.extent(0));
  const typename ATS::magnitudeType eps = 2.0 * impl_ATS::eps();

  Kokkos::parallel_for(
      "MueLu:Maxwell1::FindNonZeros", range_type(0, vals.extent(0)),
      KOKKOS_LAMBDA(const size_t i) {
        nonzeros(i) = (impl_ATS::magnitude(vals(i, 0)) > eps);
      });
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                  const Teuchos::ArrayRCP<bool>& dirichletRows,
                                  Teuchos::ArrayRCP<bool> dirichletCols,
                                  Teuchos::ArrayRCP<bool> dirichletDomain) {
  const Scalar one                                                 = Teuchos::ScalarTraits<Scalar>::one();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domMap = A.getDomainMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap = A.getRowMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap = A.getColMap();
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletRows.size()) == rowMap->getLocalNumElements());
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletCols.size()) == colMap->getLocalNumElements());
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletDomain.size()) == domMap->getLocalNumElements());
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> myColsToZero = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(colMap, 1, /*zeroOut=*/true);
  // Find all local column indices that are in Dirichlet rows, record in myColsToZero as 1.0
  for (size_t i = 0; i < (size_t)dirichletRows.size(); i++) {
    if (dirichletRows[i]) {
      ArrayView<const LocalOrdinal> indices;
      ArrayView<const Scalar> values;
      A.getLocalRowView(i, indices, values);
      for (size_t j = 0; j < static_cast<size_t>(indices.size()); j++)
        myColsToZero->replaceLocalValue(indices[j], 0, one);
    }
  }

  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> globalColsToZero;
  RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> importer = A.getCrsGraph()->getImporter();
  if (!importer.is_null()) {
    globalColsToZero = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(domMap, 1, /*zeroOut=*/true);
    // export to domain map
    globalColsToZero->doExport(*myColsToZero, *importer, Xpetra::ADD);
    // import to column map
    myColsToZero->doImport(*globalColsToZero, *importer, Xpetra::INSERT);
  } else
    globalColsToZero = myColsToZero;

  FindNonZeros(globalColsToZero->getData(0), dirichletDomain);
  FindNonZeros(myColsToZero->getData(0), dirichletCols);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                  const Kokkos::View<bool*, typename Node::device_type>& dirichletRows,
                                  Kokkos::View<bool*, typename Node::device_type> dirichletCols,
                                  Kokkos::View<bool*, typename Node::device_type> dirichletDomain) {
  using ATS                                                        = Kokkos::ArithTraits<Scalar>;
  using impl_ATS                                                   = Kokkos::ArithTraits<typename ATS::val_type>;
  using range_type                                                 = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domMap = A.getDomainMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap = A.getRowMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap = A.getColMap();
  TEUCHOS_ASSERT(dirichletRows.extent(0) == rowMap->getLocalNumElements());
  TEUCHOS_ASSERT(dirichletCols.extent(0) == colMap->getLocalNumElements());
  TEUCHOS_ASSERT(dirichletDomain.extent(0) == domMap->getLocalNumElements());
  RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> myColsToZero = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(colMap, /*zeroOut=*/true);
  // Find all local column indices that are in Dirichlet rows, record in myColsToZero as 1.0
  auto myColsToZeroView = myColsToZero->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto localMatrix      = A.getLocalMatrixDevice();
  Kokkos::parallel_for(
      "MueLu:Maxwell1::DetectDirichletCols", range_type(0, rowMap->getLocalNumElements()),
      KOKKOS_LAMBDA(const LocalOrdinal row) {
        if (dirichletRows(row)) {
          auto rowView = localMatrix.row(row);
          auto length  = rowView.length;

          for (decltype(length) colID = 0; colID < length; colID++)
            myColsToZeroView(rowView.colidx(colID), 0) = impl_ATS::one();
        }
      });

  RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> globalColsToZero;
  RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> importer = A.getCrsGraph()->getImporter();
  if (!importer.is_null()) {
    globalColsToZero = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(domMap, /*zeroOut=*/true);
    // export to domain map
    globalColsToZero->doExport(*myColsToZero, *importer, Xpetra::ADD);
    // import to column map
    myColsToZero->doImport(*globalColsToZero, *importer, Xpetra::INSERT);
  } else
    globalColsToZero = myColsToZero;
  UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FindNonZeros(globalColsToZero->getDeviceLocalView(Xpetra::Access::ReadOnly), dirichletDomain);
  UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FindNonZeros(myColsToZero->getDeviceLocalView(Xpetra::Access::ReadOnly), dirichletCols);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol, Teuchos::ArrayRCP<bool>& dirichletRows) {
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> MTS;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowmap = A.getRowMap();
  for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(rowmap->getLocalNumElements()); ++row) {
    size_t nnz = A.getNumEntriesInLocalRow(row);
    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    A.getLocalRowView(row, indices, vals);

    Scalar rowsum  = STS::zero();
    Scalar diagval = STS::zero();

    for (LocalOrdinal colID = 0; colID < Teuchos::as<LocalOrdinal>(nnz); colID++) {
      LocalOrdinal col = indices[colID];
      if (row == col)
        diagval = vals[colID];
      rowsum += vals[colID];
    }
    //        printf("A(%d,:) row_sum(point) = %6.4e\n",row,rowsum);
    if (rowSumTol < MTS::one() && STS::magnitude(rowsum) > STS::magnitude(diagval) * rowSumTol) {
      // printf("Row %d triggers rowsum\n",(int)row);
      dirichletRows[row] = true;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber, const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol, Teuchos::ArrayRCP<bool>& dirichletRows) {
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> MTS;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowmap = A.getRowMap();

  TEUCHOS_TEST_FOR_EXCEPTION(!A.getColMap()->isSameAs(*BlockNumber.getMap()), std::runtime_error, "ApplyRowSumCriterion: BlockNumber must match's A's column map.");

  Teuchos::ArrayRCP<const LocalOrdinal> block_id = BlockNumber.getData(0);
  for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(rowmap->getLocalNumElements()); ++row) {
    size_t nnz = A.getNumEntriesInLocalRow(row);
    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    A.getLocalRowView(row, indices, vals);

    Scalar rowsum  = STS::zero();
    Scalar diagval = STS::zero();
    for (LocalOrdinal colID = 0; colID < Teuchos::as<LocalOrdinal>(nnz); colID++) {
      LocalOrdinal col = indices[colID];
      if (row == col)
        diagval = vals[colID];
      if (block_id[row] == block_id[col])
        rowsum += vals[colID];
    }

    //        printf("A(%d,:) row_sum(block) = %6.4e\n",row,rowsum);
    if (rowSumTol < MTS::one() && STS::magnitude(rowsum) > STS::magnitude(diagval) * rowSumTol) {
      // printf("Row %d triggers rowsum\n",(int)row);
      dirichletRows[row] = true;
    }
  }
}

// Applies rowsum criterion
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class memory_space>
void ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                          const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                          Kokkos::View<bool*, memory_space>& dirichletRows) {
  typedef Teuchos::ScalarTraits<Scalar> STS;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowmap = A.getRowMap();

  auto dirichletRowsHost = Kokkos::create_mirror_view(dirichletRows);
  Kokkos::deep_copy(dirichletRowsHost, dirichletRows);

  for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(rowmap->getLocalNumElements()); ++row) {
    size_t nnz = A.getNumEntriesInLocalRow(row);
    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    A.getLocalRowView(row, indices, vals);

    Scalar rowsum  = STS::zero();
    Scalar diagval = STS::zero();
    for (LocalOrdinal colID = 0; colID < Teuchos::as<LocalOrdinal>(nnz); colID++) {
      LocalOrdinal col = indices[colID];
      if (row == col)
        diagval = vals[colID];
      rowsum += vals[colID];
    }
    if (STS::real(rowsum) > STS::magnitude(diagval) * rowSumTol)
      dirichletRowsHost(row) = true;
  }

  Kokkos::deep_copy(dirichletRows, dirichletRowsHost);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                         const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                         Kokkos::View<bool*, typename Node::device_type::memory_space>& dirichletRows) {
  MueLu::ApplyRowSumCriterion<Scalar, LocalOrdinal, GlobalOrdinal, Node, typename Node::device_type::memory_space>(A, rowSumTol, dirichletRows);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyRowSumCriterionHost(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                             const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                             Kokkos::View<bool*, Kokkos::HostSpace>& dirichletRows) {
  MueLu::ApplyRowSumCriterion<Scalar, LocalOrdinal, GlobalOrdinal, Node, Kokkos::HostSpace>(A, rowSumTol, dirichletRows);
}

// Applies rowsum criterion
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class memory_space>
void ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                          const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber,
                          const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                          Kokkos::View<bool*, memory_space>& dirichletRows) {
  typedef Teuchos::ScalarTraits<Scalar> STS;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowmap = A.getRowMap();

  TEUCHOS_TEST_FOR_EXCEPTION(!A.getColMap()->isSameAs(*BlockNumber.getMap()), std::runtime_error, "ApplyRowSumCriterion: BlockNumber must match's A's column map.");

  auto dirichletRowsHost = Kokkos::create_mirror_view(dirichletRows);
  Kokkos::deep_copy(dirichletRowsHost, dirichletRows);

  Teuchos::ArrayRCP<const LocalOrdinal> block_id = BlockNumber.getData(0);
  for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(rowmap->getLocalNumElements()); ++row) {
    size_t nnz = A.getNumEntriesInLocalRow(row);
    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    A.getLocalRowView(row, indices, vals);

    Scalar rowsum  = STS::zero();
    Scalar diagval = STS::zero();
    for (LocalOrdinal colID = 0; colID < Teuchos::as<LocalOrdinal>(nnz); colID++) {
      LocalOrdinal col = indices[colID];
      if (row == col)
        diagval = vals[colID];
      if (block_id[row] == block_id[col])
        rowsum += vals[colID];
    }
    if (STS::real(rowsum) > STS::magnitude(diagval) * rowSumTol)
      dirichletRowsHost(row) = true;
  }

  Kokkos::deep_copy(dirichletRows, dirichletRowsHost);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                         const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber,
                         const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                         Kokkos::View<bool*, typename Node::device_type::memory_space>& dirichletRows) {
  MueLu::ApplyRowSumCriterion<Scalar, LocalOrdinal, GlobalOrdinal, Node, typename Node::device_type::memory_space>(A, BlockNumber, rowSumTol, dirichletRows);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyRowSumCriterionHost(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                             const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber,
                             const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                             Kokkos::View<bool*, Kokkos::HostSpace>& dirichletRows) {
  MueLu::ApplyRowSumCriterion<Scalar, LocalOrdinal, GlobalOrdinal, Node, Kokkos::HostSpace>(A, BlockNumber, rowSumTol, dirichletRows);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const bool>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DetectDirichletCols(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                        const Teuchos::ArrayRCP<const bool>& dirichletRows) {
  Scalar zero                                                                               = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar one                                                                                = Teuchos::ScalarTraits<Scalar>::one();
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domMap                 = A.getDomainMap();
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap                 = A.getColMap();
  Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> myColsToZero = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(colMap, 1);
  myColsToZero->putScalar(zero);
  // Find all local column indices that are in Dirichlet rows, record in myColsToZero as 1.0
  for (size_t i = 0; i < (size_t)dirichletRows.size(); i++) {
    if (dirichletRows[i]) {
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> values;
      A.getLocalRowView(i, indices, values);
      for (size_t j = 0; j < static_cast<size_t>(indices.size()); j++)
        myColsToZero->replaceLocalValue(indices[j], 0, one);
    }
  }

  Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> globalColsToZero = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(domMap, 1);
  globalColsToZero->putScalar(zero);
  Teuchos::RCP<Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node>> exporter = Xpetra::ExportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(colMap, domMap);
  // export to domain map
  globalColsToZero->doExport(*myColsToZero, *exporter, Xpetra::ADD);
  // import to column map
  myColsToZero->doImport(*globalColsToZero, *exporter, Xpetra::INSERT);
  Teuchos::ArrayRCP<const Scalar> myCols = myColsToZero->getData(0);
  Teuchos::ArrayRCP<bool> dirichletCols(colMap->getLocalNumElements(), true);
  Magnitude eps = Teuchos::ScalarTraits<Magnitude>::eps();
  for (size_t i = 0; i < colMap->getLocalNumElements(); i++) {
    dirichletCols[i] = Teuchos::ScalarTraits<Scalar>::magnitude(myCols[i]) > 2.0 * eps;
  }
  return dirichletCols;
}

template <class SC, class LO, class GO, class NO>
Kokkos::View<bool*, typename NO::device_type>
UtilitiesBase<SC, LO, GO, NO>::
    DetectDirichletCols(const Xpetra::Matrix<SC, LO, GO, NO>& A,
                        const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows) {
  using ATS        = Kokkos::ArithTraits<SC>;
  using impl_ATS   = Kokkos::ArithTraits<typename ATS::val_type>;
  using range_type = Kokkos::RangePolicy<LO, typename NO::execution_space>;

  SC zero = ATS::zero();

  auto localMatrix = A.getLocalMatrixDevice();
  LO numRows       = A.getLocalNumRows();

  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> domMap             = A.getDomainMap();
  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> colMap             = A.getColMap();
  Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> myColsToZero = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(colMap, 1);
  myColsToZero->putScalar(zero);
  auto myColsToZeroView = myColsToZero->getDeviceLocalView(Xpetra::Access::ReadWrite);
  // Find all local column indices that are in Dirichlet rows, record in myColsToZero as 1.0
  Kokkos::parallel_for(
      "MueLu:Utils::DetectDirichletCols1", range_type(0, numRows),
      KOKKOS_LAMBDA(const LO row) {
        if (dirichletRows(row)) {
          auto rowView = localMatrix.row(row);
          auto length  = rowView.length;

          for (decltype(length) colID = 0; colID < length; colID++)
            myColsToZeroView(rowView.colidx(colID), 0) = impl_ATS::one();
        }
      });

  Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> globalColsToZero = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(domMap, 1);
  globalColsToZero->putScalar(zero);
  Teuchos::RCP<Xpetra::Export<LO, GO, NO>> exporter = Xpetra::ExportFactory<LO, GO, NO>::Build(colMap, domMap);
  // export to domain map
  globalColsToZero->doExport(*myColsToZero, *exporter, Xpetra::ADD);
  // import to column map
  myColsToZero->doImport(*globalColsToZero, *exporter, Xpetra::INSERT);

  auto myCols          = myColsToZero->getDeviceLocalView(Xpetra::Access::ReadOnly);
  size_t numColEntries = colMap->getLocalNumElements();
  Kokkos::View<bool*, typename NO::device_type> dirichletCols(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), numColEntries);
  const typename ATS::magnitudeType eps = 2.0 * ATS::eps();

  Kokkos::parallel_for(
      "MueLu:Utils::DetectDirichletCols2", range_type(0, numColEntries),
      KOKKOS_LAMBDA(const size_t i) {
        dirichletCols(i) = impl_ATS::magnitude(myCols(i, 0)) > eps;
      });
  return dirichletCols;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Frobenius(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B) {
  // We check only row maps. Column may be different. One would hope that they are the same, as we typically
  // calculate frobenius norm of the specified sparsity pattern with an updated matrix from the previous step,
  // but matrix addition, even when one is submatrix of the other, changes column map (though change may be as
  // simple as couple of elements swapped)
  TEUCHOS_TEST_FOR_EXCEPTION(!A.getRowMap()->isSameAs(*B.getRowMap()), Exceptions::Incompatible, "MueLu::CGSolver::Frobenius: row maps are incompatible");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete() || !B.isFillComplete(), Exceptions::RuntimeError, "Matrices must be fill completed");

  const Map& AColMap = *A.getColMap();
  const Map& BColMap = *B.getColMap();

  Teuchos::ArrayView<const LocalOrdinal> indA, indB;
  Teuchos::ArrayView<const Scalar> valA, valB;
  size_t nnzA = 0, nnzB = 0;

  // We use a simple algorithm
  // for each row we fill valBAll array with the values in the corresponding row of B
  // as such, it serves as both sorted array and as storage, so we don't need to do a
  // tricky problem: "find a value in the row of B corresponding to the specific GID"
  // Once we do that, we translate LID of entries of row of A to LID of B, and multiply
  // corresponding entries.
  // The algorithm should be reasonably cheap, as it does not sort anything, provided
  // that getLocalElement and getGlobalElement functions are reasonably effective. It
  // *is* possible that the costs are hidden in those functions, but if maps are close
  // to linear maps, we should be fine
  Teuchos::Array<Scalar> valBAll(BColMap.getLocalNumElements());

  LocalOrdinal invalid = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero(), f = zero, gf;
  size_t numRows = A.getLocalNumRows();
  for (size_t i = 0; i < numRows; i++) {
    A.getLocalRowView(i, indA, valA);
    B.getLocalRowView(i, indB, valB);
    nnzA = indA.size();
    nnzB = indB.size();

    // Set up array values
    for (size_t j = 0; j < nnzB; j++)
      valBAll[indB[j]] = valB[j];

    for (size_t j = 0; j < nnzA; j++) {
      // The cost of the whole Frobenius dot product function depends on the
      // cost of the getLocalElement and getGlobalElement functions here.
      LocalOrdinal ind = BColMap.getLocalElement(AColMap.getGlobalElement(indA[j]));
      if (ind != invalid)
        f += valBAll[ind] * valA[j];
    }

    // Clean up array values
    for (size_t j = 0; j < nnzB; j++)
      valBAll[indB[j]] = zero;
  }

  MueLu_sumAll(AColMap.getComm(), f, gf);

  return gf;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    SetRandomSeed(const Teuchos::Comm<int>& comm) {
  // Distribute the seeds evenly in [1,maxint-1].  This guarantees nothing
  // about where in random number stream we are, but avoids overflow situations
  // in parallel when multiplying by a PID.  It would be better to use
  // a good parallel random number generator.
  double one = 1.0;
  int maxint = INT_MAX;  //= 2^31-1 = 2147483647 for 32-bit integers
  int mySeed = Teuchos::as<int>((maxint - 1) * (one - (comm.getRank() + 1) / (comm.getSize() + one)));
  if (mySeed < 1 || mySeed == maxint) {
    std::ostringstream errStr;
    errStr << "Error detected with random seed = " << mySeed << ". It should be in the interval [1,2^31-2].";
    throw Exceptions::RuntimeError(errStr.str());
  }
  std::srand(mySeed);
  // For Tpetra, we could use Kokkos' random number generator here.
  Teuchos::ScalarTraits<Scalar>::seedrandom(mySeed);
  // Epetra
  //   MultiVector::Random() -> Epetra_Util::RandomDouble() -> Epetra_Utils::RandomInt()
  // Its own random number generator, based on Seed_. Seed_ is initialized in Epetra_Util constructor with std::rand()
  // So our setting std::srand() affects that too
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FindDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                      std::vector<LocalOrdinal>& dirichletRows, bool count_twos_as_dirichlet) {
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  dirichletRows.resize(0);
  for (size_t i = 0; i < A->getLocalNumRows(); i++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;
    A->getLocalRowView(i, indices, values);
    int nnz = 0;
    for (size_t j = 0; j < (size_t)indices.size(); j++) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(values[j]) > Teuchos::ScalarTraits<MT>::eps()) {
        nnz++;
      }
    }
    if (nnz == 1 || (count_twos_as_dirichlet && nnz == 2)) {
      dirichletRows.push_back(i);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                         const std::vector<LocalOrdinal>& dirichletRows) {
  RCP<const Map> Rmap = A->getRowMap();
  RCP<const Map> Cmap = A->getColMap();
  Scalar one          = Teuchos::ScalarTraits<Scalar>::one();
  Scalar zero         = Teuchos::ScalarTraits<Scalar>::zero();

  for (size_t i = 0; i < dirichletRows.size(); i++) {
    GlobalOrdinal row_gid = Rmap->getGlobalElement(dirichletRows[i]);

    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;
    A->getLocalRowView(dirichletRows[i], indices, values);
    // NOTE: This won't work with fancy node types.
    Scalar* valuesNC = const_cast<Scalar*>(values.getRawPtr());
    for (size_t j = 0; j < (size_t)indices.size(); j++) {
      if (Cmap->getGlobalElement(indices[j]) == row_gid)
        valuesNC[j] = one;
      else
        valuesNC[j] = zero;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                         const Teuchos::ArrayRCP<const bool>& dirichletRows) {
  TEUCHOS_ASSERT(A->isFillComplete());
  RCP<const Map> domMap = A->getDomainMap();
  RCP<const Map> ranMap = A->getRangeMap();
  RCP<const Map> Rmap   = A->getRowMap();
  RCP<const Map> Cmap   = A->getColMap();
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletRows.size()) == Rmap->getLocalNumElements());
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  A->resumeFill();
  for (size_t i = 0; i < (size_t)dirichletRows.size(); i++) {
    if (dirichletRows[i]) {
      GlobalOrdinal row_gid = Rmap->getGlobalElement(i);

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> values;
      A->getLocalRowView(i, indices, values);

      Teuchos::ArrayRCP<Scalar> valuesNC(values.size());
      for (size_t j = 0; j < (size_t)indices.size(); j++) {
        if (Cmap->getGlobalElement(indices[j]) == row_gid)
          valuesNC[j] = one;
        else
          valuesNC[j] = zero;
      }
      A->replaceLocalValues(i, indices, valuesNC());
    }
  }
  A->fillComplete(domMap, ranMap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                         const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows) {
  TEUCHOS_ASSERT(A->isFillComplete());
  using ATS        = Kokkos::ArithTraits<Scalar>;
  using impl_ATS   = Kokkos::ArithTraits<typename ATS::val_type>;
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domMap = A->getDomainMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> ranMap = A->getRangeMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> Rmap   = A->getRowMap();
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> Cmap   = A->getColMap();

  TEUCHOS_ASSERT(static_cast<size_t>(dirichletRows.size()) == Rmap->getLocalNumElements());

  auto localMatrix = A->getLocalMatrixDevice();
  auto localRmap   = Rmap->getLocalMap();
  auto localCmap   = Cmap->getLocalMap();

  Kokkos::parallel_for(
      "MueLu::Utils::ApplyOAZ", range_type(0, dirichletRows.extent(0)),
      KOKKOS_LAMBDA(const LocalOrdinal row) {
        if (dirichletRows(row)) {
          auto rowView = localMatrix.row(row);
          auto length  = rowView.length;
          auto row_gid = localRmap.getGlobalElement(row);
          auto row_lid = localCmap.getLocalElement(row_gid);

          for (decltype(length) colID = 0; colID < length; colID++)
            if (rowView.colidx(colID) == row_lid)
              rowView.value(colID) = impl_ATS::one();
            else
              rowView.value(colID) = impl_ATS::zero();
        }
      });
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                      const std::vector<LocalOrdinal>& dirichletRows,
                      Scalar replaceWith) {
  for (size_t i = 0; i < dirichletRows.size(); i++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;
    A->getLocalRowView(dirichletRows[i], indices, values);
    // NOTE: This won't work with fancy node types.
    Scalar* valuesNC = const_cast<Scalar*>(values.getRawPtr());
    for (size_t j = 0; j < (size_t)indices.size(); j++)
      valuesNC[j] = replaceWith;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                      const Teuchos::ArrayRCP<const bool>& dirichletRows,
                      Scalar replaceWith) {
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletRows.size()) == A->getRowMap()->getLocalNumElements());
  for (size_t i = 0; i < (size_t)dirichletRows.size(); i++) {
    if (dirichletRows[i]) {
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> values;
      A->getLocalRowView(i, indices, values);
      // NOTE: This won't work with fancy node types.
      Scalar* valuesNC = const_cast<Scalar*>(values.getRawPtr());
      for (size_t j = 0; j < (size_t)indices.size(); j++)
        valuesNC[j] = replaceWith;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletRows(Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X,
                      const Teuchos::ArrayRCP<const bool>& dirichletRows,
                      Scalar replaceWith) {
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletRows.size()) == X->getMap()->getLocalNumElements());
  for (size_t i = 0; i < (size_t)dirichletRows.size(); i++) {
    if (dirichletRows[i]) {
      for (size_t j = 0; j < X->getNumVectors(); j++)
        X->replaceLocalValue(i, j, replaceWith);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletRows(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                      const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows,
                      Scalar replaceWith) {
  using ATS        = Kokkos::ArithTraits<Scalar>;
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  typename ATS::val_type impl_replaceWith = replaceWith;

  auto localMatrix     = A->getLocalMatrixDevice();
  LocalOrdinal numRows = A->getLocalNumRows();

  Kokkos::parallel_for(
      "MueLu:Utils::ZeroDirichletRows", range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal row) {
        if (dirichletRows(row)) {
          auto rowView = localMatrix.row(row);
          auto length  = rowView.length;
          for (decltype(length) colID = 0; colID < length; colID++)
            rowView.value(colID) = impl_replaceWith;
        }
      });
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletRows(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X,
                      const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows,
                      Scalar replaceWith) {
  using ATS        = Kokkos::ArithTraits<Scalar>;
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  typename ATS::val_type impl_replaceWith = replaceWith;

  auto myCols    = X->getDeviceLocalView(Xpetra::Access::ReadWrite);
  size_t numVecs = X->getNumVectors();
  Kokkos::parallel_for(
      "MueLu:Utils::ZeroDirichletRows_MV", range_type(0, dirichletRows.size()),
      KOKKOS_LAMBDA(const size_t i) {
        if (dirichletRows(i)) {
          for (size_t j = 0; j < numVecs; j++)
            myCols(i, j) = impl_replaceWith;
        }
      });
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletCols(Teuchos::RCP<Matrix>& A,
                      const Teuchos::ArrayRCP<const bool>& dirichletCols,
                      Scalar replaceWith) {
  TEUCHOS_ASSERT(static_cast<size_t>(dirichletCols.size()) == A->getColMap()->getLocalNumElements());
  for (size_t i = 0; i < A->getLocalNumRows(); i++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;
    A->getLocalRowView(i, indices, values);
    // NOTE: This won't work with fancy node types.
    Scalar* valuesNC = const_cast<Scalar*>(values.getRawPtr());
    for (size_t j = 0; j < static_cast<size_t>(indices.size()); j++)
      if (dirichletCols[indices[j]])
        valuesNC[j] = replaceWith;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ZeroDirichletCols(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                      const Kokkos::View<const bool*, typename Node::device_type>& dirichletCols,
                      Scalar replaceWith) {
  using ATS        = Kokkos::ArithTraits<Scalar>;
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  typename ATS::val_type impl_replaceWith = replaceWith;

  auto localMatrix     = A->getLocalMatrixDevice();
  LocalOrdinal numRows = A->getLocalNumRows();

  Kokkos::parallel_for(
      "MueLu:Utils::ZeroDirichletCols", range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal row) {
        auto rowView = localMatrix.row(row);
        auto length  = rowView.length;
        for (decltype(length) colID = 0; colID < length; colID++)
          if (dirichletCols(rowView.colidx(colID))) {
            rowView.value(colID) = impl_replaceWith;
          }
      });
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FindDirichletRowsAndPropagateToCols(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                        Teuchos::RCP<Xpetra::Vector<int, LocalOrdinal, GlobalOrdinal, Node>>& isDirichletRow,
                                        Teuchos::RCP<Xpetra::Vector<int, LocalOrdinal, GlobalOrdinal, Node>>& isDirichletCol) {
  // Make sure A's RowMap == DomainMap
  if (!A->getRowMap()->isSameAs(*A->getDomainMap())) {
    throw std::runtime_error("UtilitiesBase::FindDirichletRowsAndPropagateToCols row and domain maps must match.");
  }
  RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> importer = A->getCrsGraph()->getImporter();
  bool has_import                                                       = !importer.is_null();

  // Find the Dirichlet rows
  std::vector<LocalOrdinal> dirichletRows;
  FindDirichletRows(A, dirichletRows);

#if 0
    printf("[%d] DirichletRow Ids = ",A->getRowMap()->getComm()->getRank());
    for(size_t i=0; i<(size_t) dirichletRows.size(); i++)
      printf("%d ",dirichletRows[i]);
    printf("\n");
    fflush(stdout);
#endif
  // Allocate all as non-Dirichlet
  isDirichletRow = Xpetra::VectorFactory<int, LocalOrdinal, GlobalOrdinal, Node>::Build(A->getRowMap(), true);
  isDirichletCol = Xpetra::VectorFactory<int, LocalOrdinal, GlobalOrdinal, Node>::Build(A->getColMap(), true);

  {
    Teuchos::ArrayRCP<int> dr_rcp = isDirichletRow->getDataNonConst(0);
    Teuchos::ArrayView<int> dr    = dr_rcp();
    Teuchos::ArrayRCP<int> dc_rcp = isDirichletCol->getDataNonConst(0);
    Teuchos::ArrayView<int> dc    = dc_rcp();
    for (size_t i = 0; i < (size_t)dirichletRows.size(); i++) {
      dr[dirichletRows[i]] = 1;
      if (!has_import) dc[dirichletRows[i]] = 1;
    }
  }

  if (has_import)
    isDirichletCol->doImport(*isDirichletRow, *importer, Xpetra::CombineMode::ADD);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ReplaceNonZerosWithOnes(const RCP<Matrix>& original) {
  using ISC               = typename Kokkos::ArithTraits<Scalar>::val_type;
  using range_type        = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;
  using local_matrix_type = typename CrsMatrix::local_matrix_type;
  using values_type       = typename local_matrix_type::values_type;

  const ISC ONE  = Kokkos::ArithTraits<ISC>::one();
  const ISC ZERO = Kokkos::ArithTraits<ISC>::zero();

  // Copy the values array of the old matrix to a new array, replacing all the non-zeros with one
  auto localMatrix = original->getLocalMatrixDevice();
  TEUCHOS_TEST_FOR_EXCEPTION(!original->hasCrsGraph(), Exceptions::RuntimeError, "ReplaceNonZerosWithOnes: Cannot get CrsGraph");
  values_type new_values("values", localMatrix.nnz());

  Kokkos::parallel_for(
      "ReplaceNonZerosWithOnes", range_type(0, localMatrix.nnz()), KOKKOS_LAMBDA(const size_t i) {
        if (localMatrix.values(i) != ZERO)
          new_values(i) = ONE;
        else
          new_values(i) = ZERO;
      });

  // Build the new matrix
  RCP<Matrix> NewMatrix = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(original->getCrsGraph(), new_values);
  TEUCHOS_TEST_FOR_EXCEPTION(NewMatrix.is_null(), Exceptions::RuntimeError, "ReplaceNonZerosWithOnes: MatrixFactory::Build() did not return matrix");
  NewMatrix->fillComplete(original->getDomainMap(), original->getRangeMap());
  return NewMatrix;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GeneratedBlockedTargetMap(const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>& sourceBlockedMap,
                              const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>& Importer) {
  typedef Xpetra::Vector<int, LocalOrdinal, GlobalOrdinal, Node> IntVector;
  Xpetra::UnderlyingLib lib = sourceBlockedMap.lib();

  // De-stride the map if we have to (might regret this later)
  RCP<const Map> fullMap    = sourceBlockedMap.getMap();
  RCP<const Map> stridedMap = Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>(fullMap);
  if (!stridedMap.is_null()) fullMap = stridedMap->getMap();

  // Initial sanity checking for map compatibil
  const size_t numSubMaps = sourceBlockedMap.getNumMaps();
  if (!Importer.getSourceMap()->isCompatible(*fullMap))
    throw std::runtime_error("GenerateBlockedTargetMap(): Map compatibility error");

  // Build an indicator vector
  RCP<IntVector> block_ids = Xpetra::VectorFactory<int, LocalOrdinal, GlobalOrdinal, Node>::Build(fullMap);

  for (size_t i = 0; i < numSubMaps; i++) {
    RCP<const Map> map = sourceBlockedMap.getMap(i);

    for (size_t j = 0; j < map->getLocalNumElements(); j++) {
      LocalOrdinal jj = fullMap->getLocalElement(map->getGlobalElement(j));
      block_ids->replaceLocalValue(jj, (int)i);
    }
  }

  // Get the block ids for the new map
  RCP<const Map> targetMap     = Importer.getTargetMap();
  RCP<IntVector> new_block_ids = Xpetra::VectorFactory<int, LocalOrdinal, GlobalOrdinal, Node>::Build(targetMap);
  new_block_ids->doImport(*block_ids, Importer, Xpetra::CombineMode::ADD);
  Teuchos::ArrayRCP<const int> dataRCP = new_block_ids->getData(0);
  Teuchos::ArrayView<const int> data   = dataRCP();

  // Get the GIDs for each subblock
  Teuchos::Array<Teuchos::Array<GlobalOrdinal>> elementsInSubMap(numSubMaps);
  for (size_t i = 0; i < targetMap->getLocalNumElements(); i++) {
    elementsInSubMap[data[i]].push_back(targetMap->getGlobalElement(i));
  }

  // Generate the new submaps
  std::vector<RCP<const Map>> subMaps(numSubMaps);
  for (size_t i = 0; i < numSubMaps; i++) {
    subMaps[i] = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), elementsInSubMap[i](), targetMap->getIndexBase(), targetMap->getComm());
  }

  // Build the BlockedMap
  return rcp(new BlockedMap(targetMap, subMaps));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapsAreNested(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& rowMap, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& colMap) {
  ArrayView<const GlobalOrdinal> rowElements = rowMap.getLocalElementList();
  ArrayView<const GlobalOrdinal> colElements = colMap.getLocalElementList();

  const size_t numElements = rowElements.size();

  if (size_t(colElements.size()) < numElements)
    return false;

  bool goodMap = true;
  for (size_t i = 0; i < numElements; i++)
    if (rowElements[i] != colElements[i]) {
      goodMap = false;
      break;
    }

  return goodMap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ReverseCuthillMcKee(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  using local_matrix_type = typename Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using lno_nnz_view_t    = typename local_graph_type::entries_type::non_const_type;
  using device            = typename local_graph_type::device_type;
  using execution_space   = typename local_matrix_type::execution_space;
  using ordinal_type      = typename local_matrix_type::ordinal_type;

  local_graph_type localGraph = Op.getLocalMatrixDevice().graph;

  lno_nnz_view_t rcmOrder = KokkosGraph::Experimental::graph_rcm<device, typename local_graph_type::row_map_type, typename local_graph_type::entries_type, lno_nnz_view_t>(localGraph.row_map, localGraph.entries);

  RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> retval =
      Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(Op.getRowMap());

  // Copy out and reorder data
  auto view1D = Kokkos::subview(retval->getDeviceLocalView(Xpetra::Access::ReadWrite), Kokkos::ALL(), 0);
  Kokkos::parallel_for(
      "Utilities::ReverseCuthillMcKee",
      Kokkos::RangePolicy<ordinal_type, execution_space>(0, localGraph.numRows()),
      KOKKOS_LAMBDA(const ordinal_type rowIdx) {
        view1D(rcmOrder(rowIdx)) = rowIdx;
      });
  return retval;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>
UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CuthillMcKee(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  using local_matrix_type = typename Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using lno_nnz_view_t    = typename local_graph_type::entries_type::non_const_type;
  using device            = typename local_graph_type::device_type;
  using execution_space   = typename local_matrix_type::execution_space;
  using ordinal_type      = typename local_matrix_type::ordinal_type;

  local_graph_type localGraph = Op.getLocalMatrixDevice().graph;
  LocalOrdinal numRows        = localGraph.numRows();

  lno_nnz_view_t rcmOrder = KokkosGraph::Experimental::graph_rcm<device, typename local_graph_type::row_map_type, typename local_graph_type::entries_type, lno_nnz_view_t>(localGraph.row_map, localGraph.entries);

  RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> retval =
      Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(Op.getRowMap());

  // Copy out data
  auto view1D = Kokkos::subview(retval->getDeviceLocalView(Xpetra::Access::ReadWrite), Kokkos::ALL(), 0);
  // Since KokkosKernels produced RCM, also reverse the order of the view to get CM
  Kokkos::parallel_for(
      "Utilities::ReverseCuthillMcKee",
      Kokkos::RangePolicy<ordinal_type, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const ordinal_type rowIdx) {
        view1D(rcmOrder(numRows - 1 - rowIdx)) = rowIdx;
      });
  return retval;
}

}  // namespace MueLu

#define MUELU_UTILITIESBASE_SHORT
#endif  // MUELU_UTILITIESBASE_DEF_HPP

//  LocalWords:  LocalOrdinal
