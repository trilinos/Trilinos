// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_VECTORDROPPINGCLASSICAL_DECL_HPP
#define MUELU_VECTORDROPPINGCLASSICAL_DECL_HPP

#include "MueLu_VectorDroppingBase.hpp"
#include "MueLu_BoundaryDetection.hpp"
#include "MueLu_ClassicalDropping.hpp"
#include "MueLu_CutDrop.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_DistanceLaplacianDropping.hpp"
#include "MueLu_MatrixConstruction.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure SoC>
class VectorDroppingClassical : public VectorDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using matrix_type             = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_matrix_type         = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType               = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type       = typename crs_matrix_type::local_matrix_device_type;
  using local_graph_type        = typename GraphType::local_graph_device_type;
  using rowptr_type             = typename local_graph_type::row_map_type::non_const_type;
  using device_type             = typename Node::device_type;
  using memory_space            = typename device_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using magnitudeType           = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using boundary_nodes_type     = typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type;
  using Level                   = MueLu::Level;
  using nnz_count_type          = Kokkos::pair<LocalOrdinal, LocalOrdinal>;
  using block_indices_view_type = typename Kokkos::View<LocalOrdinal*, typename Node::device_type>;

  static void runDroppingFunctors_on_A(matrix_type& A,
                                       matrix_type& mergedA,
                                       LocalOrdinal blkPartSize,
                                       block_indices_view_type& rowTranslation,
                                       block_indices_view_type& colTranslation,
                                       results_view& results,
                                       rowptr_type& filtered_rowptr,
                                       rowptr_type& graph_rowptr,
                                       nnz_count_type& nnz,
                                       boundary_nodes_type& boundaryNodes,
                                       const std::string& droppingMethod,
                                       const magnitudeType threshold,
                                       const bool aggregationMayCreateDirichlet,
                                       const bool symmetrizeDroppedGraph,
                                       const bool useBlocking,
                                       Level& level,
                                       const Factory& factory);
};

}  // namespace MueLu
#endif
