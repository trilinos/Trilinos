// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SCALARDROPPINGCLASSICAL_DECL_HPP
#define MUELU_SCALARDROPPINGCLASSICAL_DECL_HPP

#include "MueLu_ScalarDroppingBase.hpp"
#include "MueLu_ClassicalDropping.hpp"
#include "MueLu_CutDrop.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_MatrixConstruction.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure SoC>
class ScalarDroppingClassical : private ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using matrix_type         = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_matrix_type     = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType           = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type   = typename crs_matrix_type::local_matrix_device_type;
  using local_graph_type    = typename GraphType::local_graph_device_type;
  using rowptr_type         = typename local_graph_type::row_map_type::non_const_type;
  using device_type         = typename Node::device_type;
  using memory_space        = typename device_type::memory_space;
  using results_view        = Kokkos::View<DecisionType*, memory_space>;
  using magnitudeType       = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using boundary_nodes_type = typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type;
  using Level               = MueLu::Level;

  static void runDroppingFunctors_on_A(matrix_type& A,
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
                                       const Factory& factory);
};

}  // namespace MueLu
#endif
