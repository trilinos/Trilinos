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
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Xpetra_Access.hpp"
#include "Xpetra_MatrixFactory.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_SparseConstraint_decl.hpp"
#include "MueLu_Utilities.hpp"

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
  this->PrepareLeastSquaresSolve(solverType, /*singular=*/true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SparseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup() {
  using graph_t        = typename CrsGraph::local_graph_type;
  using matrix_t       = typename CrsMatrix::local_matrix_type;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename matrix_t::values_type::non_const_type;
  using range_type     = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  auto D        = D_;
  auto Dc       = Dc_;
  auto Ppattern = this->GetPattern();

  TEUCHOS_TEST_FOR_EXCEPTION(!D->getRangeMap()->isSameAs(*Ppattern->getRangeMap()),
                             Exceptions::Incompatible,
                             "Maps are incompatible");
  TEUCHOS_TEST_FOR_EXCEPTION(!Dc->getRangeMap()->isSameAs(*Ppattern->getDomainMap()),
                             Exceptions::Incompatible,
                             "Maps are incompatible");

  // Construct auxiliary graph via Ppattern * Dc
  RCP<const CrsGraph> auxGraph;
  {
    const auto one = Teuchos::ScalarTraits<Scalar>::one();
    auto absP      = MatrixFactory::Build(Ppattern);
    absP->setAllToScalar(one);
    absP->fillComplete();

    auto absDc = MatrixFactory::BuildCopy(Dc);
    absDc->setAllToScalar(one);

    auto P_Dc = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absP, false, *absDc, false, this->GetOStream(Statistics2), true, true);
    auxGraph  = P_Dc->getCrsGraph();
  }
  RHS_pattern_ = auxGraph;

  auto lib                                    = Ppattern->getRowMap()->lib();
  auto comm                                   = Ppattern->getComm();
  GlobalOrdinal indexBase                     = Ppattern->getRowMap()->getIndexBase();
  const size_t numUnknowns                    = Ppattern->getLocalNumEntries();
  const size_t numRows                        = Ppattern->getLocalNumRows();
  Xpetra::global_size_t global_numConstraints = auxGraph->getGlobalNumEntries();
  Xpetra::global_size_t global_numUnknowns    = Ppattern->getGlobalNumEntries();
  const size_t numConstraints                 = auxGraph->getLocalNumEntries();
  auto constraint_rowmap                      = MapFactory::Build(lib, global_numConstraints, numConstraints, indexBase, comm);
  auto constraint_domainmap                   = MapFactory::Build(lib, global_numUnknowns, numUnknowns, indexBase, comm);

  RCP<Matrix> X;
  {
    auto lclPattern  = Ppattern->getLocalGraphDevice();
    auto lclD0       = Dc->getLocalMatrixDevice();
    auto lclAuxGraph = auxGraph->getLocalGraphDevice();

    lno_view_t rowptr("constraint_rowptr", numConstraints + 1);

    size_t nnz = 0;
    Kokkos::parallel_reduce(
        "sparse_constraint",
        range_type(0, numRows),
        KOKKOS_LAMBDA(const size_t pattern_i, size_t& partial_nnz) {
          for (size_t pattern_jj = lclPattern.row_map(pattern_i); pattern_jj < lclPattern.row_map(pattern_i + 1); ++pattern_jj) {
            auto pattern_j = lclPattern.entries(pattern_jj);
            for (size_t mat_jj = lclD0.graph.row_map(pattern_j); mat_jj < lclD0.graph.row_map(pattern_j + 1); ++mat_jj) {
              auto mat_j = lclD0.graph.entries(mat_jj);
              // Find entry mat_j in row graph_i of tempGraph
              size_t constraint_i;
              for (constraint_i = lclAuxGraph.row_map(pattern_i); constraint_i < lclAuxGraph.row_map(pattern_i + 1); ++constraint_i) {
                if (lclAuxGraph.entries(constraint_i) == mat_j)
                  break;
              }
              partial_nnz += 1;
              if (constraint_i + 2 < numConstraints + 1) {
                Kokkos::atomic_add(&rowptr(constraint_i + 2), 1);
              }
            }
          }
        },
        nnz);

    // todo: merge?
    Kokkos::parallel_scan(
        "sparse_constraint",
        range_type(1, numConstraints + 1),
        KOKKOS_LAMBDA(const size_t constraint_i, size_t& partial_nnz, bool is_final) {
          partial_nnz += rowptr(constraint_i);
          if (is_final)
            rowptr(constraint_i) = partial_nnz;
        });

    lno_nnz_view_t colind(Kokkos::ViewAllocateWithoutInitializing("constraint_indices"), nnz);
    scalar_view_t values(Kokkos::ViewAllocateWithoutInitializing("constraint_values"), nnz);

    Kokkos::parallel_for(
        "sparse_constraint",
        range_type(0, numRows),
        KOKKOS_LAMBDA(const size_t pattern_i) {
          for (size_t pattern_jj = lclPattern.row_map(pattern_i); pattern_jj < lclPattern.row_map(pattern_i + 1); ++pattern_jj) {
            auto pattern_j = lclPattern.entries(pattern_jj);
            for (size_t mat_jj = lclD0.graph.row_map(pattern_j); mat_jj < lclD0.graph.row_map(pattern_j + 1); ++mat_jj) {
              auto mat_j   = lclD0.graph.entries(mat_jj);
              auto mat_val = lclD0.values(mat_jj);
              // Find entry mat_j in row graph_i of tempGraph
              size_t constraint_i;
              for (constraint_i = lclAuxGraph.row_map(pattern_i); constraint_i < lclAuxGraph.row_map(pattern_i + 1); ++constraint_i) {
                if (lclAuxGraph.entries(constraint_i) == mat_j)
                  break;
              }
              auto constraint_jj    = rowptr(constraint_i + 1);
              colind(constraint_jj) = pattern_jj;
              values(constraint_jj) = mat_val;
              ++rowptr(constraint_i + 1);
            }
          }
        });

    auto lclConstraintGraph = graph_t(colind, rowptr);
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
                                                                                   *Dc_, false,
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
