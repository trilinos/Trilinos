// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SCALARDROPPINGDISTANCELAPLACIAN_DECL_HPP
#define MUELU_SCALARDROPPINGDISTANCELAPLACIAN_DECL_HPP

#include "MueLu_ScalarDroppingBase.hpp"
#include "MueLu_CutDrop.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_DistanceLaplacianDropping.hpp"
#include "MueLu_MatrixConstruction.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure SoC>
class ScalarDroppingDistanceLaplacian : public ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using matrix_type         = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_matrix_type     = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType           = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type   = typename crs_matrix_type::local_matrix_device_type;
  using local_graph_type    = typename GraphType::local_graph_device_type;
  using rowptr_type         = typename local_graph_type::row_map_type::non_const_type;
  using entries_type        = typename local_graph_type::entries_type::non_const_type;
  using values_type         = typename local_matrix_type::values_type::non_const_type;
  using device_type         = typename Node::device_type;
  using memory_space        = typename device_type::memory_space;
  using results_view        = Kokkos::View<DecisionType*, memory_space>;
  using magnitudeType       = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using boundary_nodes_type = typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type;
  using Level               = MueLu::Level;

  template <class DistanceFunctorType>
  static void runDroppingFunctors_on_dlap_inner(matrix_type& A,
                                                results_view& results,
                                                rowptr_type& filtered_rowptr,
                                                LocalOrdinal& nnz_filtered,
                                                boundary_nodes_type& boundaryNodes,
                                                const std::string& droppingMethod,
                                                const magnitudeType threshold,
                                                const bool aggregationMayCreateDirichlet,
                                                const bool symmetrizeDroppedGraph,
                                                const bool useBlocking,
                                                DistanceFunctorType& dist2,
                                                Level& level,
                                                const Factory& factory) {
    auto lclA               = A.getLocalMatrixDevice();
    auto preserve_diagonals = Misc::KeepDiagonalFunctor(lclA, results);

    if (droppingMethod == "point-wise") {
      auto dist_laplacian_dropping = DistanceLaplacian::make_drop_functor<SoC>(A, threshold, dist2, results);

      if (aggregationMayCreateDirichlet) {
        auto mark_singletons_as_boundary = Misc::MarkSingletonFunctor(lclA, boundaryNodes, results);

        if (symmetrizeDroppedGraph) {
          auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(A, boundaryNodes, results);
          ScalarDroppingDistanceLaplacian::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                               dist_laplacian_dropping,
                                                               drop_boundaries,
                                                               preserve_diagonals,
                                                               mark_singletons_as_boundary);
        } else {
          auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);
          ScalarDroppingDistanceLaplacian::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                               dist_laplacian_dropping,
                                                               drop_boundaries,
                                                               preserve_diagonals,
                                                               mark_singletons_as_boundary);
        }
      } else {
        if (symmetrizeDroppedGraph) {
          auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(A, boundaryNodes, results);
          ScalarDroppingDistanceLaplacian::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                               dist_laplacian_dropping,
                                                               drop_boundaries,
                                                               preserve_diagonals);
        } else {
          auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);
          ScalarDroppingDistanceLaplacian::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                               dist_laplacian_dropping,
                                                               drop_boundaries,
                                                               preserve_diagonals);
        }
      }
    } else if (droppingMethod == "cut-drop") {
      auto comparison = CutDrop::make_dlap_comparison_functor<SoC>(A, dist2, results);
      auto cut_drop   = CutDrop::CutDropFunctor(comparison, threshold);

      if (symmetrizeDroppedGraph) {
        auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(A, boundaryNodes, results);
        ScalarDroppingDistanceLaplacian::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                             drop_boundaries,
                                                             preserve_diagonals,
                                                             cut_drop);
      } else {
        auto drop_boundaries = Misc::PointwiseDropBoundaryFunctor(lclA, boundaryNodes, results);
        ScalarDroppingDistanceLaplacian::runDroppingFunctors(A, results, filtered_rowptr, nnz_filtered, useBlocking, level, factory,
                                                             drop_boundaries,
                                                             preserve_diagonals,
                                                             cut_drop);
      }
    }
  }

  static void runDroppingFunctors_on_dlap(matrix_type& A,
                                          results_view& results,
                                          rowptr_type& filtered_rowptr,
                                          LocalOrdinal& nnz_filtered,
                                          boundary_nodes_type& boundaryNodes,
                                          const std::string& droppingMethod,
                                          const magnitudeType threshold,
                                          const bool aggregationMayCreateDirichlet,
                                          const bool symmetrizeDroppedGraph,
                                          const bool useBlocking,
                                          const std::string& distanceLaplacianMetric,
                                          Level& level,
                                          const Factory& factory);
};

}  // namespace MueLu
#endif
