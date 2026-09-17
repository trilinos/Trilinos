// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DENSECONSTRAINT_DEF_HPP
#define MUELU_DENSECONSTRAINT_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_CrsGraph.hpp>
#include "Teuchos_CommHelpers.hpp"
#include "KokkosKernels_ArithTraits.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_DenseConstraint_decl.hpp"
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MatrixFactory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
DenseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DenseConstraint(const RCP<MultiVector>& B,
                    const RCP<MultiVector>& Bc,
                    RCP<const CrsGraph> Ppattern,
                    const std::string& solverType) {
  this->SetPattern(Ppattern);
  B_  = B;
  Bc_ = Bc;
  this->SetProcRankVerbose(Ppattern->getRowMap()->getComm()->getRank());
  Setup();
  this->PrepareLeastSquaresSolve(solverType, /*detect_singular_blocks=*/false);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DenseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup() {
  using graph_t         = typename CrsGraph::local_graph_type;
  using matrix_t        = typename CrsMatrix::local_matrix_type;
  using lno_view_t      = typename graph_t::row_map_type::non_const_type;
  using lno_nnz_view_t  = typename graph_t::entries_type::non_const_type;
  using scalar_view_t   = typename matrix_t::values_type::non_const_type;
  using execution_space = typename Node::execution_space;
  using range_type      = Kokkos::RangePolicy<LocalOrdinal, execution_space>;

  using ATS            = KokkosKernels::ArithTraits<Scalar>;
  using magnitude_type = typename ATS::magnitudeType;

  auto Ppattern            = this->GetPattern();
  auto B                   = B_;
  auto Bc                  = Bc_;
  const size_t NSDim       = Bc->getNumVectors();
  const size_t numUnknowns = Ppattern->getLocalNumEntries();
  const size_t numRows     = Ppattern->getLocalNumRows();

  TEUCHOS_TEST_FOR_EXCEPTION(!B->getMap()->isSameAs(*Ppattern->getRangeMap()),
                             Exceptions::Incompatible,
                             "Maps are incompatible");
  TEUCHOS_TEST_FOR_EXCEPTION(!Bc->getMap()->isSameAs(*Ppattern->getDomainMap()),
                             Exceptions::Incompatible,
                             "Maps are incompatible");

  auto importer = Ppattern->getImporter();
  RCP<MultiVector> ghostedBc;
  if (!importer.is_null()) {
    ghostedBc = MultiVectorFactory::Build(Ppattern->getColMap(), NSDim);
    ghostedBc->doImport(*Bc, *importer, Xpetra::INSERT);
  } else {
    ghostedBc = Bc;
  }

  auto lib                = Ppattern->getRowMap()->lib();
  auto comm               = Ppattern->getComm();
  GlobalOrdinal indexBase = Ppattern->getRowMap()->getIndexBase();

  auto lclPattern   = Ppattern->getLocalGraphDevice();
  auto lclNullspace = ghostedBc->getLocalViewDevice(Tpetra::Access::ReadOnly);

  const LocalOrdinal numPotentialConstraints =
      static_cast<LocalOrdinal>(numRows * NSDim);

  // Treat only exactly-zero coefficients as structurally inactive.  This
  // skips canonical-nullspace constraints such as y/z translation constraints
  // on an x-displacement row, without discarding small but legitimate
  // coordinate-dependent nullspace entries.
  const magnitude_type zeroTol =
      KokkosKernels::ArithTraits<magnitude_type>::zero();

  // ---------------------------------------------------------------------
  // Build constraints only for nonzero constraint rows.
  //
  // The old implementation created NSDim constraints for every fine scalar
  // row of P, regardless of whether the corresponding constraint row had any
  // nonzero coefficients. For vector elasticity with canonical translations,
  // this creates many trivial 0 = 0 constraints, e.g. y/z translation
  // constraints on x-displacement rows. Those rows make X*X^T singular and
  // cause NaNs in the direct block inverse.
  //
  // A constraint row (i,k) is active iff row i of Ppattern has at least one
  // column j for which Bc(j,k) is nonzero.
  // ---------------------------------------------------------------------

  LocalOrdinal numConstraints = 0;

  Kokkos::parallel_reduce(
      "MueLu::DenseConstraint::countNonzeroConstraintRows",
      range_type(0, numPotentialConstraints),
      KOKKOS_LAMBDA(const LocalOrdinal c, LocalOrdinal& count) {
        const LocalOrdinal i = static_cast<LocalOrdinal>(c / static_cast<LocalOrdinal>(NSDim));
        const LocalOrdinal k = static_cast<LocalOrdinal>(c - i * static_cast<LocalOrdinal>(NSDim));

        bool active = false;

        for (LocalOrdinal jj = static_cast<LocalOrdinal>(lclPattern.row_map(i));
             jj < static_cast<LocalOrdinal>(lclPattern.row_map(i + 1)); ++jj) {
          const auto j = lclPattern.entries(jj);

          if (ATS::magnitude(lclNullspace(j, k)) != zeroTol) {
            active = true;
            break;
          }
        }

        if (active)
          ++count;
      },
      numConstraints);

  Kokkos::fence();

  // Store the original fine row and nullspace-vector index for each retained
  // constraint row.
  lno_nnz_view_t activeRows(
      Kokkos::ViewAllocateWithoutInitializing("DenseConstraint_activeRows"),
      numConstraints);

  lno_nnz_view_t activeModes(
      Kokkos::ViewAllocateWithoutInitializing("DenseConstraint_activeModes"),
      numConstraints);

  LocalOrdinal scanTotal = 0;

  Kokkos::parallel_scan(
      "MueLu::DenseConstraint::fillNonzeroConstraintRows",
      range_type(0, numPotentialConstraints),
      KOKKOS_LAMBDA(const LocalOrdinal c, LocalOrdinal& update, const bool final) {
        const LocalOrdinal i = static_cast<LocalOrdinal>(c / static_cast<LocalOrdinal>(NSDim));
        const LocalOrdinal k = static_cast<LocalOrdinal>(c - i * static_cast<LocalOrdinal>(NSDim));

        bool active = false;

        for (LocalOrdinal jj = static_cast<LocalOrdinal>(lclPattern.row_map(i));
             jj < static_cast<LocalOrdinal>(lclPattern.row_map(i + 1)); ++jj) {
          const auto j = lclPattern.entries(jj);

          if (ATS::magnitude(lclNullspace(j, k)) != zeroTol) {
            active = true;
            break;
          }
        }

        if (active) {
          const LocalOrdinal activeId = update;

          if (final) {
            activeRows(activeId)  = i;
            activeModes(activeId) = k;
          }

          ++update;
        }
      },
      scanTotal);

  Kokkos::fence();

  TEUCHOS_ASSERT_EQUALITY(scanTotal, numConstraints);

  // Count matrix entries in X. Each retained constraint row has one entry for
  // every nonzero in the corresponding Ppattern row.
  LocalOrdinal nnz = 0;

  Kokkos::parallel_reduce(
      "MueLu::DenseConstraint::countConstraintNnz",
      range_type(0, numConstraints),
      KOKKOS_LAMBDA(const LocalOrdinal activeId, LocalOrdinal& count) {
        const LocalOrdinal i = activeRows(activeId);

        const LocalOrdinal rowBegin =
            static_cast<LocalOrdinal>(lclPattern.row_map(i));

        const LocalOrdinal rowEnd =
            static_cast<LocalOrdinal>(lclPattern.row_map(i + 1));

        count += rowEnd - rowBegin;
      },
      nnz);

  Kokkos::fence();

  GlobalOrdinal localConstraintsGO = static_cast<GlobalOrdinal>(numConstraints);
  GlobalOrdinal globalConstraints  = 0;

  Teuchos::reduceAll(
      *comm,
      Teuchos::REDUCE_SUM,
      1,
      &localConstraintsGO,
      &globalConstraints);

  Xpetra::global_size_t global_numConstraints =
      static_cast<Xpetra::global_size_t>(globalConstraints);

  Xpetra::global_size_t global_numUnknowns =
      Ppattern->getGlobalNumEntries();

  auto constraint_rowmap =
      MapFactory::Build(lib,
                        global_numConstraints,
                        static_cast<size_t>(numConstraints),
                        indexBase,
                        comm);

  auto constraint_domainmap =
      MapFactory::Build(lib,
                        global_numUnknowns,
                        numUnknowns,
                        indexBase,
                        comm);

  // The matrix of constraints X has size
  //
  //   global_numConstraints x global_numUnknowns.
  //
  // Empty/trivial constraint rows are omitted.
  RCP<Matrix> X;
  {
    lno_view_t rowptr(
        Kokkos::ViewAllocateWithoutInitializing("constraint_rowptr"),
        static_cast<size_t>(numConstraints) + 1);

    lno_nnz_view_t colind(
        Kokkos::ViewAllocateWithoutInitializing("constraint_indices"),
        nnz);

    scalar_view_t values(
        Kokkos::ViewAllocateWithoutInitializing("constraint_values"),
        nnz);

    LocalOrdinal nnzScanTotal = 0;

    Kokkos::parallel_scan(
        "MueLu::DenseConstraint::buildConstraintRowptr",
        range_type(0, numConstraints),
        KOKKOS_LAMBDA(const LocalOrdinal activeId, LocalOrdinal& update, const bool final) {
          const LocalOrdinal i = activeRows(activeId);

          const LocalOrdinal rowBegin =
              static_cast<LocalOrdinal>(lclPattern.row_map(i));

          const LocalOrdinal rowEnd =
              static_cast<LocalOrdinal>(lclPattern.row_map(i + 1));

          if (final)
            rowptr(activeId) = update;

          update += rowEnd - rowBegin;
        },
        nnzScanTotal);

    Kokkos::parallel_for(
        "MueLu::DenseConstraint::setFinalRowptr",
        Kokkos::RangePolicy<execution_space>(0, 1),
        KOKKOS_LAMBDA(const int) {
          rowptr(numConstraints) = nnzScanTotal;
        });

    Kokkos::fence();

    TEUCHOS_ASSERT_EQUALITY(nnzScanTotal, nnz);

    Kokkos::parallel_for(
        "MueLu::DenseConstraint::buildNullspaceConstraint",
        range_type(0, numConstraints),
        KOKKOS_LAMBDA(const LocalOrdinal activeId) {
          const LocalOrdinal i = activeRows(activeId);
          const LocalOrdinal k = activeModes(activeId);

          const LocalOrdinal rowBegin =
              static_cast<LocalOrdinal>(lclPattern.row_map(i));

          const LocalOrdinal rowEnd =
              static_cast<LocalOrdinal>(lclPattern.row_map(i + 1));

          const LocalOrdinal outBegin = rowptr(activeId);

          LocalOrdinal l = 0;
          for (LocalOrdinal jj = rowBegin; jj < rowEnd; ++jj) {
            const auto j = lclPattern.entries(jj);

            // The unknown vector stores entries of P in local row-storage
            // order.  Thus jj is the vectorized local P-entry ordinal; j is
            // the local column index used to read the ghosted coarse nullspace.
            colind(outBegin + l) = jj;
            values(outBegin + l) = lclNullspace(j, k);

            ++l;
          }
        });

    Kokkos::fence();

    auto lclConstraintGraph = graph_t(colind, rowptr);
    auto lclConstraint      = matrix_t("constraint", numUnknowns, values, lclConstraintGraph);

    X = MatrixFactory::Build(lclConstraint,
                             constraint_rowmap,
                             constraint_domainmap,
                             constraint_domainmap,
                             constraint_rowmap);
  }

  this->SetConstraintsMatrix(X);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type
DenseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FindBlocks(RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>& /*XXt*/) {
  using graph_type      = typename CrsGraph::local_graph_type;
  using execution_space = typename Node::execution_space;
  using range_type      = Kokkos::RangePolicy<LocalOrdinal, execution_space>;

  using ATS            = KokkosKernels::ArithTraits<Scalar>;
  using magnitude_type = typename ATS::magnitudeType;

  auto Ppattern = this->GetPattern();

  const size_t NSDim = Bc_->getNumVectors();

  auto importer = Ppattern->getImporter();
  RCP<MultiVector> ghostedBc;
  if (!importer.is_null()) {
    ghostedBc = MultiVectorFactory::Build(Ppattern->getColMap(), NSDim);
    ghostedBc->doImport(*Bc_, *importer, Xpetra::INSERT);
  } else {
    ghostedBc = Bc_;
  }

  auto lclPattern   = Ppattern->getLocalGraphDevice();
  auto lclNullspace = ghostedBc->getLocalViewDevice(Tpetra::Access::ReadOnly);

  auto X    = this->GetConstraintMatrix();
  auto lclX = X->getLocalMatrixDevice();

  const LocalOrdinal numRows =
      static_cast<LocalOrdinal>(Ppattern->getLocalNumRows());

  const LocalOrdinal numConstraints =
      static_cast<LocalOrdinal>(lclX.numRows());

  const magnitude_type zeroTol =
      KokkosKernels::ArithTraits<magnitude_type>::zero();

  typename graph_type::entries_type::non_const_type activeModeCounts(
      Kokkos::ViewAllocateWithoutInitializing("DenseConstraint_block_activeModeCounts"),
      numRows);

  Kokkos::parallel_for(
      "MueLu::DenseConstraint::countActiveModesPerRow",
      range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal i) {
        LocalOrdinal count = 0;

        for (size_t k = 0; k < NSDim; ++k) {
          bool active = false;

          for (LocalOrdinal jj = static_cast<LocalOrdinal>(lclPattern.row_map(i));
               jj < static_cast<LocalOrdinal>(lclPattern.row_map(i + 1)); ++jj) {
            const auto j = lclPattern.entries(jj);

            if (ATS::magnitude(lclNullspace(j, k)) != zeroTol) {
              active = true;
              break;
            }
          }

          if (active)
            ++count;
        }

        activeModeCounts(i) = count;
      });

  Kokkos::fence();

  LocalOrdinal numBlocks = 0;

  Kokkos::parallel_reduce(
      "MueLu::DenseConstraint::countNonemptyBlocks",
      range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal i, LocalOrdinal& count) {
        if (activeModeCounts(i) > static_cast<LocalOrdinal>(0))
          ++count;
      },
      numBlocks);

  Kokkos::fence();

  typename graph_type::row_map_type::non_const_type fullConstraintOffsets(
      Kokkos::ViewAllocateWithoutInitializing("DenseConstraint_fullConstraintOffsets"),
      static_cast<size_t>(numRows) + 1);

  LocalOrdinal totalConstraints = 0;

  Kokkos::parallel_scan(
      "MueLu::DenseConstraint::scanActiveModeCounts",
      range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal i, LocalOrdinal& update, const bool final) {
        if (final)
          fullConstraintOffsets(i) = update;

        update += activeModeCounts(i);
      },
      totalConstraints);

  Kokkos::parallel_for(
      "MueLu::DenseConstraint::setFinalFullConstraintOffset",
      Kokkos::RangePolicy<execution_space>(0, 1),
      KOKKOS_LAMBDA(const int) {
        fullConstraintOffsets(numRows) = totalConstraints;
      });

  Kokkos::fence();

  TEUCHOS_ASSERT_EQUALITY(totalConstraints, numConstraints);

  typename graph_type::row_map_type::non_const_type rowptr(
      Kokkos::ViewAllocateWithoutInitializing("blocks_rowptr"),
      static_cast<size_t>(numBlocks) + 1);

  typename graph_type::entries_type::non_const_type indices(
      Kokkos::ViewAllocateWithoutInitializing("blocks_indices"),
      numConstraints);

  LocalOrdinal blockScanTotal = 0;

  Kokkos::parallel_scan(
      "MueLu::DenseConstraint::buildVariableDenseBlocks",
      range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal i, LocalOrdinal& update, const bool final) {
        const LocalOrdinal count = activeModeCounts(i);

        if (count > static_cast<LocalOrdinal>(0)) {
          const LocalOrdinal blockId = update;

          if (final) {
            const LocalOrdinal begin =
                static_cast<LocalOrdinal>(fullConstraintOffsets(i));

            const LocalOrdinal end =
                static_cast<LocalOrdinal>(fullConstraintOffsets(i + 1));

            rowptr(blockId) = begin;

            for (LocalOrdinal constraintId = begin; constraintId < end; ++constraintId)
              indices(constraintId) = constraintId;
          }

          ++update;
        }
      },
      blockScanTotal);

  Kokkos::parallel_for(
      "MueLu::DenseConstraint::setFinalBlockRowptr",
      Kokkos::RangePolicy<execution_space>(0, 1),
      KOKKOS_LAMBDA(const int) {
        rowptr(numBlocks) = numConstraints;
      });

  Kokkos::fence();

  TEUCHOS_ASSERT_EQUALITY(blockScanTotal, numBlocks);

  return graph_type(indices, rowptr);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
DenseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ResidualNorm(const RCP<const Matrix> P) const {
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  using ATS            = KokkosKernels::ArithTraits<Scalar>;
  using magnitude_type = typename ATS::magnitudeType;

  const auto one  = Teuchos::ScalarTraits<Scalar>::one();
  const auto zero = Teuchos::ScalarTraits<Scalar>::zero();

  auto residual = MultiVectorFactory::Build(B_->getMap(), B_->getNumVectors());
  P->apply(*Bc_, *residual, Teuchos::NO_TRANS);
  residual->update(one, *B_, -one);

  auto Ppattern = this->GetPattern();

  auto importer = Ppattern->getImporter();
  RCP<MultiVector> ghostedBc;
  if (!importer.is_null()) {
    ghostedBc = MultiVectorFactory::Build(Ppattern->getColMap(), Bc_->getNumVectors());
    ghostedBc->doImport(*Bc_, *importer, Xpetra::INSERT);
  } else {
    ghostedBc = Bc_;
  }

  auto lclGraph     = Ppattern->getLocalGraphDevice();
  auto lclNullspace = ghostedBc->getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto lclRes       = residual->getLocalViewDevice(Tpetra::Access::ReadWrite);

  const size_t numVecs = residual->getNumVectors();
  const LocalOrdinal numRows =
      static_cast<LocalOrdinal>(Ppattern->getLocalNumRows());

  const magnitude_type zeroTol =
      KokkosKernels::ArithTraits<magnitude_type>::zero();

  // Mask residual entries for constraint rows that were skipped in Setup,
  // i.e., rows (i,k) whose coefficients against the Ppattern row are all zero.
  Kokkos::parallel_for(
      "MueLu::DenseConstraint::maskSkippedConstraintResiduals",
      range_type(0, numRows),
      KOKKOS_LAMBDA(const LocalOrdinal i) {
        for (size_t k = 0; k < numVecs; ++k) {
          bool active = false;

          for (LocalOrdinal jj = static_cast<LocalOrdinal>(lclGraph.row_map(i));
               jj < static_cast<LocalOrdinal>(lclGraph.row_map(i + 1)); ++jj) {
            const auto j = lclGraph.entries(jj);

            if (ATS::magnitude(lclNullspace(j, k)) != zeroTol) {
              active = true;
              break;
            }
          }

          if (!active)
            lclRes(i, k) = zero;
        }
      });

  Kokkos::fence();

  Teuchos::Array<MagnitudeType> norms(B_->getNumVectors());
  residual->norm2(norms);

  MagnitudeType residualNorm = Teuchos::ScalarTraits<MagnitudeType>::zero();
  for (size_t k = 0; k < B_->getNumVectors(); ++k) {
    residualNorm += norms[k] * norms[k];
  }

  return Teuchos::ScalarTraits<MagnitudeType>::squareroot(residualNorm);
}

}  // namespace MueLu

#endif  // ifndef MUELU_DENSECONSTRAINT_DEF_HPP
