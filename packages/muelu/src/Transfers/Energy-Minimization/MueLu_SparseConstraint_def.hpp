// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SPARSECONSTRAINT_DEF_HPP
#define MUELU_SPARSECONSTRAINT_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include "Kokkos_Atomic.hpp"
#include "Kokkos_Pair.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Xpetra_MatrixFactory.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_SparseConstraint_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SparseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    SparseConstraint(const RCP<Matrix>& P_nodal,
                     const RCP<Matrix>& D,
                     const RCP<Matrix>& Dc,
                     RCP<const CrsGraph> Ppattern,
                     const std::string& solverType) {
  this->SetPattern(Ppattern);
  P_nodal_ = P_nodal;
  D_       = D;
  Dc_      = Dc;
  Setup();
  this->PrepareLeastSquaresSolve(solverType, /*detect_singular_blocks=*/true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SparseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup() {
  using graph_t        = typename CrsGraph::local_graph_type;
  using matrix_t       = typename CrsMatrix::local_matrix_type;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename matrix_t::values_type::non_const_type;
  using range_type     = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  Monitor m(*this, "Setup");

  auto D        = D_;
  auto Dc       = Dc_;
  auto Ppattern = this->GetPattern();

  // The constraint on Pe (with graph Ppattern) takes the form
  //
  //   Pe * Dc = D * Pn.
  //
  // This means that we have nnz(Pe * Dc) constraints for nnz(Ppattern)
  // unknowns.
  //
  // A single constraint corresponds to an entry (i,j) of Pe * D0c and is
  // written out via the sparse matrix-matrix products between Pe and D0c:
  //
  //  sum_{k} Pe_{i,k} Dc_{k,j} = (D * Pn)_{i,j}
  //
  // We map (i,j) to its offset I in (Pe * D0c) and (i,k) to its offset J in Pe.
  //
  // The constraint matrix X then has the entry
  //   X_{I,J} = Dc_{k,j}.

  auto lib = Ppattern->getRowMap()->lib();

  // If we rebalanced then Dc lives on a smaller communicator than D.
  // Since we need to perform matrix-matrix multiplications with Dc, we construct a version of it that lives on the same communicator.
  auto comm = Ppattern->getRowMap()->getComm();
  if (Dc.is_null() || Dc->getRowMap()->getComm()->getSize() < comm->getSize()) {
    if (Dc.is_null()) {
      Kokkos::View<GlobalOrdinal*, typename Node::memory_space> dummy("", 0);
      auto big_coarse_nodal_map    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, dummy, 0, comm);
      auto big_coarse_edge_map     = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, dummy, 0, comm);
      auto big_coarse_nodal_colmap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, dummy, 0, comm);

      typename Matrix::local_matrix_device_type dummyLocalMatrix;
      big_Dc_ = MatrixFactory::Build(dummyLocalMatrix, big_coarse_edge_map, big_coarse_nodal_colmap, big_coarse_nodal_map, big_coarse_edge_map);

    } else {
      auto big_coarse_nodal_map    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, Dc->getDomainMap()->getMyGlobalIndicesDevice(), 0, comm);
      auto big_coarse_edge_map     = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, Dc->getRangeMap()->getMyGlobalIndicesDevice(), 0, comm);
      auto big_coarse_nodal_colmap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, Dc->getColMap()->getMyGlobalIndicesDevice(), 0, comm);

      big_Dc_ = MatrixFactory::Build(Dc->getLocalMatrixDevice(), big_coarse_edge_map, big_coarse_nodal_colmap, big_coarse_nodal_map, big_coarse_edge_map);
    }
  } else {
    big_Dc_ = Dc;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!D->getRangeMap()->isSameAs(*Ppattern->getRangeMap()),
                             Exceptions::Incompatible,
                             "Maps are incompatible");
  TEUCHOS_TEST_FOR_EXCEPTION(!big_Dc_->getRangeMap()->isSameAs(*Ppattern->getDomainMap()),
                             Exceptions::Incompatible,
                             "Maps are incompatible");

  // Construct auxiliary graph via Ppattern * Dc
  RCP<const CrsGraph> auxGraph;
  {
    const auto one = Teuchos::ScalarTraits<Scalar>::one();
    auto absP      = MatrixFactory::Build(Ppattern);
    absP->setAllToScalar(one);
    absP->fillComplete();

    auto absDc = MatrixFactory::BuildCopy(big_Dc_);
    absDc->setAllToScalar(one);

    auto P_Dc = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absP, false, *absDc, false, this->GetOStream(Statistics2), true, true);
    auxGraph  = P_Dc->getCrsGraph();
  }
  RHS_pattern_ = auxGraph;

  GlobalOrdinal indexBase                     = Ppattern->getRowMap()->getIndexBase();
  const size_t numUnknowns                    = Ppattern->getLocalNumEntries();
  const size_t numRows                        = Ppattern->getLocalNumRows();
  Xpetra::global_size_t global_numConstraints = auxGraph->getGlobalNumEntries();
  Xpetra::global_size_t global_numUnknowns    = Ppattern->getGlobalNumEntries();
  const size_t numConstraints                 = auxGraph->getLocalNumEntries();
  auto constraint_rowmap                      = MapFactory::Build(lib, global_numConstraints, numConstraints, indexBase, comm);
  auto constraint_domainmap                   = MapFactory::Build(lib, global_numUnknowns, numUnknowns, indexBase, comm);

  RCP<Matrix> ghostedDc;
  if (!Ppattern->getImporter().is_null())
    ghostedDc = MatrixFactory::Build(big_Dc_, *Ppattern->getImporter());
  else
    ghostedDc = big_Dc_;

  RCP<Matrix> X;
  {
    auto lclPattern  = Ppattern->getLocalGraphDevice();
    auto lclD0       = ghostedDc->getLocalMatrixDevice();
    auto lclAuxGraph = auxGraph->getLocalGraphDevice();

    // Over-allocate by 1. Makes the logic a bit easier in what follows.
    lno_view_t rowptr("constraint_rowptr", numConstraints + 2);

    Kokkos::parallel_for(
        "sparse_constraint_num_entries_per_row",
        range_type(0, numRows),
        KOKKOS_LAMBDA(const size_t pattern_i) {
          for (size_t pattern_jj = lclPattern.row_map(pattern_i); pattern_jj < lclPattern.row_map(pattern_i + 1); ++pattern_jj) {
            auto pattern_j = lclPattern.entries(pattern_jj);
            // entry (pattern_i, pattern_j) in Ppattern

            for (size_t D0_jj = lclD0.graph.row_map(pattern_j); D0_jj < lclD0.graph.row_map(pattern_j + 1); ++D0_jj) {
              auto D0_j = lclD0.graph.entries(D0_jj);
              // entry (pattern_j, D0_j) in ghosted D0

              // Find entry (pattern_i, D0_j) in tempGraph
              size_t constraint_I;
              for (constraint_I = lclAuxGraph.row_map(pattern_i); constraint_I < lclAuxGraph.row_map(pattern_i + 1); ++constraint_I) {
                if (lclAuxGraph.entries(constraint_I) == D0_j)
                  break;
              }
#ifdef HAVE_MUELU_DEBUG
              if (lclAuxGraph.entries(constraint_I) != D0_j)
                ::Kokkos::abort("Did not find entry in row of tempGraph.");
#endif
              // Need an entry in row constraint_I.
              // We offset by 2 since we do not want to compute the final rowptr just yet.
              // That will happen during fill.
              Kokkos::atomic_add(&rowptr(constraint_I + 2), 1);
            }
          }
        });

    // The usual prefix sum.
    size_t nnz = 0;
    Kokkos::parallel_scan(
        "sparse_constraint_prefix_sum",
        range_type(1, numConstraints + 2),
        KOKKOS_LAMBDA(const size_t constraint_i, size_t& partial_nnz, bool is_final) {
          partial_nnz += rowptr(constraint_i);
          if (is_final)
            rowptr(constraint_i) = partial_nnz;
        },
        nnz);

    // allocate indices and values
    lno_nnz_view_t colind(Kokkos::ViewAllocateWithoutInitializing("constraint_indices"), nnz);
    scalar_view_t values(Kokkos::ViewAllocateWithoutInitializing("constraint_values"), nnz);

    // fill indices and values
    Kokkos::parallel_for(
        "sparse_constraint_fill",
        range_type(0, numRows),
        KOKKOS_LAMBDA(const size_t pattern_i) {
          for (size_t pattern_jj = lclPattern.row_map(pattern_i); pattern_jj < lclPattern.row_map(pattern_i + 1); ++pattern_jj) {
            auto pattern_j = lclPattern.entries(pattern_jj);
            // entry (pattern_i, pattern_j) in Ppattern

            for (size_t D0_jj = lclD0.graph.row_map(pattern_j); D0_jj < lclD0.graph.row_map(pattern_j + 1); ++D0_jj) {
              auto D0_j   = lclD0.graph.entries(D0_jj);
              auto D0_val = lclD0.values(D0_jj);
              // entry (pattern_j, D0_j) in ghosted D0

              // Find entry (pattern_i, D0_j) in tempGraph
              size_t constraint_I;
              for (constraint_I = lclAuxGraph.row_map(pattern_i); constraint_I < lclAuxGraph.row_map(pattern_i + 1); ++constraint_I) {
                if (lclAuxGraph.entries(constraint_I) == D0_j)
                  break;
              }
#ifdef HAVE_MUELU_DEBUG
              if (lclAuxGraph.entries(constraint_I) != D0_j)
                ::Kokkos::abort("Did not find entry in row of tempGraph.");
#endif
              // Enter data into constraint matrix.
              // (constraint_I, pattern_jj) -> D0_val
              // After this the rowptr will be correct.
              // This is why we had to offset the index by 2 earlier on.
              auto constraint_jj    = Kokkos::atomic_fetch_inc(&rowptr(constraint_I + 1));
              colind(constraint_jj) = pattern_jj;
              values(constraint_jj) = D0_val;
            }
          }
        });

    auto lclConstraintGraph = graph_t(colind, Kokkos::subview(rowptr, Kokkos::make_pair(size_t(0), numConstraints + 1)));
    auto lclConstraint      = matrix_t("constraint", numUnknowns, values, lclConstraintGraph);
    X                       = MatrixFactory::Build(lclConstraint, constraint_rowmap, constraint_domainmap, constraint_domainmap, constraint_rowmap);
  }
  this->SetConstraintsMatrix(X);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
SparseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ResidualNorm(const RCP<const Matrix> P) const {
  const auto one = Teuchos::ScalarTraits<Scalar>::one();

  // P*Dc
  RCP<Matrix> temp;
  temp = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P, false,
                                                                                   *big_Dc_, false,
                                                                                   temp,
                                                                                   this->GetOStream(Runtime0), true, true);
  // D*P_nodal
  RCP<Matrix> temp2;
  temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*D_, false,
                                                                                    *P_nodal_, false,
                                                                                    temp2,
                                                                                    this->GetOStream(Runtime0), true, true);

  // D*P_nodal - P*Dc
  RCP<Matrix> residual;
  Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*temp2, false, one,
                                                                                *temp, false, -one,
                                                                                residual,
                                                                                this->GetOStream(Runtime0));
  residual->fillComplete();
  return Teuchos::ScalarTraits<MagnitudeType>::squareroot(Teuchos::ScalarTraits<Scalar>::magnitude(Utilities::Frobenius(*residual, *residual)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SparseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AssignMatrixEntriesToConstraintVector(const Matrix& A,
                                                                                                        MultiVector& vecC) const {
  this->AssignMatrixEntriesToVector(A, RHS_pattern_, vecC);
}

}  // namespace MueLu

#endif  // ifndef MUELU_SPARSECONSTRAINT_DEF_HPP
