// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SCALARDROPPINGCLASSICAL_DEF_HPP
#define MUELU_SCALARDROPPINGCLASSICAL_DEF_HPP

#include "MueLu_ScalarDroppingClassical_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure SoC>
void ScalarDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, SoC>::runDroppingFunctors_on_A(matrix_type& A,
                                                                                                       results_view& results,
                                                                                                       rowptr_type& filtered_rowptr,
                                                                                                       LocalOrdinal& nnz_filtered,
                                                                                                       boundary_nodes_type& boundaryNodes,
                                                                                                       const std::string& droppingMethod,
                                                                                                       const magnitudeType threshold,
                                                                                                       const bool aggregationMayCreateDirichlet,
                                                                                                       const bool symmetrizeDroppedGraph,
                                                                                                       const bool useBlocking,
                                                                                                       Level& level,
                                                                                                       const Factory& factory) {
  auto lclA               = A.getLocalMatrixDevice();
  auto preserve_diagonals = Misc::KeepDiagonalFunctor(lclA, results);

  if (droppingMethod == "point-wise") {
    auto dropping = ClassicalDropping::make_drop_functor<SoC>(A, threshold, results);

    if (aggregationMayCreateDirichlet) {
      auto mark_singletons_as_boundary = Misc::MarkSingletonFunctor(lclA, boundaryNodes, results);

      if (symmetrizeDroppedGraph) {
        auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(A, boundaryNodes, results);
        ScalarDroppingClassical::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                     dropping,
                                                     drop_boundaries,
                                                     preserve_diagonals,
                                                     mark_singletons_as_boundary);
      } else {
        auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);
        ScalarDroppingClassical::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                     dropping,
                                                     drop_boundaries,
                                                     preserve_diagonals,
                                                     mark_singletons_as_boundary);
      }
    } else {
      if (symmetrizeDroppedGraph) {
        auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(A, boundaryNodes, results);
        ScalarDroppingClassical::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                     dropping,
                                                     drop_boundaries,
                                                     preserve_diagonals);
      } else {
        auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);
        ScalarDroppingClassical::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                     dropping,
                                                     drop_boundaries,
                                                     preserve_diagonals);
      }
    }
  } else if (droppingMethod == "cut-drop") {
    auto comparison = CutDrop::make_comparison_functor<SoC>(A, results);
    auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

    if (symmetrizeDroppedGraph) {
      auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(A, boundaryNodes, results);
      ScalarDroppingClassical::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                   drop_boundaries,
                                                   preserve_diagonals,
                                                   cut_drop);
    } else {
      auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);
      ScalarDroppingClassical::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                   drop_boundaries,
                                                   preserve_diagonals,
                                                   cut_drop);
    }
  }
}
}  // namespace MueLu

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  MUELU_ETI_SLGN_SoC(MueLu::ScalarDroppingClassical, SC, LO, GO, NO)

#endif
