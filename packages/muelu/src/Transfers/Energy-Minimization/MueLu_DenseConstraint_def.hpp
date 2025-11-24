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

#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Atomic.hpp"
#include "Kokkos_Macros.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_DenseConstraint_decl.hpp"
#include "Serial/Kokkos_Serial_Parallel_Range.hpp"
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
  Setup();
  this->PrepareLeastSquaresSolve(solverType, /*detect_singular_blocks=*/false);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DenseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup() {
  using graph_t        = typename CrsGraph::local_graph_type;
  using matrix_t       = typename CrsMatrix::local_matrix_type;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename matrix_t::values_type::non_const_type;
  using range_type     = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  auto Ppattern               = this->GetPattern();
  auto B                      = B_;
  auto Bc                     = Bc_;
  const size_t NSDim          = Bc->getNumVectors();
  const size_t numUnknowns    = Ppattern->getLocalNumEntries();
  const size_t numRows        = Ppattern->getLocalNumRows();
  const size_t numConstraints = numRows * NSDim;

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

  auto lib                                    = Ppattern->getRowMap()->lib();
  auto comm                                   = Ppattern->getComm();
  GlobalOrdinal indexBase                     = Ppattern->getRowMap()->getIndexBase();
  Xpetra::global_size_t global_numConstraints = Ppattern->getGlobalNumRows() * NSDim;
  Xpetra::global_size_t global_numUnknowns    = Ppattern->getGlobalNumEntries();
  auto constraint_rowmap                      = MapFactory::Build(lib, global_numConstraints, numConstraints, indexBase, comm);
  auto constraint_domainmap                   = MapFactory::Build(lib, global_numUnknowns, numUnknowns, indexBase, comm);

  // The matrix of constraints X has size (global_numConstraints, global_numUnknowns).
  RCP<Matrix> X;
  {
    auto lclPattern     = Ppattern->getLocalGraphDevice();
    auto pattern_rowptr = lclPattern.row_map;

    auto lclNullspace = ghostedBc->getLocalViewDevice(Tpetra::Access::ReadOnly);

    lno_view_t rowptr(Kokkos::ViewAllocateWithoutInitializing("constraint_rowptr"), numConstraints + 1);
    LocalOrdinal nnz = NSDim * lclPattern.entries.extent(0);
    lno_nnz_view_t colind(Kokkos::ViewAllocateWithoutInitializing("constraint_indices"), nnz);
    scalar_view_t values(Kokkos::ViewAllocateWithoutInitializing("constraint_values"), nnz);

    Kokkos::parallel_for(
        "nullspace_constraint",
        range_type(0, numRows + 1),
        KOKKOS_LAMBDA(const size_t i) {
          if (i < numRows) {
            for (size_t k = 0; k < NSDim; ++k) {
              rowptr(NSDim * i + k) = NSDim * lclPattern.row_map(i) + k * (lclPattern.row_map(i + 1) - lclPattern.row_map(i));

              size_t l = 0;
              for (size_t jj = lclPattern.row_map(i); jj < lclPattern.row_map(i + 1); ++jj) {
                auto j                            = lclPattern.entries(jj);
                colind(rowptr(NSDim * i + k) + l) = jj;
                values(rowptr(NSDim * i + k) + l) = lclNullspace(j, k);
                ++l;
              }
            }
          } else {
            rowptr(numConstraints) = nnz;
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
DenseConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ResidualNorm(const RCP<const Matrix> P) const {
  const auto one = Teuchos::ScalarTraits<Scalar>::one();

  auto residual = MultiVectorFactory::Build(B_->getMap(), B_->getNumVectors());
  P->apply(*Bc_, *residual, Teuchos::NO_TRANS);
  residual->update(one, *B_, -one);
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
