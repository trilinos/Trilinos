// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_VECTORDROPPINGBASE_HPP
#define MUELU_VECTORDROPPINGBASE_HPP

#include "Xpetra_Matrix.hpp"

#include "MueLu_DroppingCommon.hpp"
#include "MueLu_MatrixConstruction.hpp"
#include "MueLu_Factory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class VectorDroppingBase {
 public:
  using matrix_type             = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_matrix_type         = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType               = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type       = typename crs_matrix_type::local_matrix_device_type;
  using local_graph_type        = typename GraphType::local_graph_device_type;
  using rowptr_type             = typename local_graph_type::row_map_type::non_const_type;
  using device_type             = typename Node::device_type;
  using memory_space            = typename device_type::memory_space;
  using execution_space         = typename local_matrix_type::execution_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using Level                   = MueLu::Level;
  using range_type              = Kokkos::RangePolicy<LocalOrdinal, execution_space>;
  using LocalOrdinalVector      = Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  using nnz_count_type          = Kokkos::pair<LocalOrdinal, LocalOrdinal>;
  using block_indices_view_type = typename Kokkos::View<LocalOrdinal*, typename Node::device_type>;

  template <class... Functors>
  static void runDroppingFunctorsImpl(local_matrix_type& lclA, typename matrix_type::local_ordinal_type& blkPartSize, block_indices_view_type& colTranslation, results_view& results, rowptr_type& filtered_rowptr, rowptr_type& graph_rowptr, nnz_count_type& nnz, Functors&... functors) {
    auto numNodes = graph_rowptr.extent(0) - 1;
    auto range    = range_type(0, numNodes);
#if !defined(HAVE_MUELU_DEBUG)
    auto countingFunctor = MatrixConstruction::VectorCountingFunctor(lclA, blkPartSize, colTranslation, results, filtered_rowptr, graph_rowptr, functors...);

#else
    auto debug           = Misc::DebugFunctor(lclA, results);
    auto countingFunctor = MatrixConstruction::VectorCountingFunctor(lclA, blkPartSize, colTranslation, results, filtered_rowptr, graph_rowptr, functors...);
#endif
    Kokkos::parallel_scan("MueLu::CoalesceDrop::CountEntries", range, countingFunctor, nnz);
  }

  template <class... Functors>
  static void runDroppingFunctors(matrix_type& A, matrix_type& mergedA, typename matrix_type::local_ordinal_type& blkPartSize, block_indices_view_type& rowTranslation, block_indices_view_type& colTranslation, results_view& results, rowptr_type& filtered_rowptr, rowptr_type& graph_rowptr, nnz_count_type& nnz, const bool useBlocking, Level& level, const Factory& factory, Functors&... functors) {
    auto lclA = A.getLocalMatrixDevice();
    if (useBlocking) {
      auto BlockNumber       = level.template Get<Teuchos::RCP<LocalOrdinalVector>>("BlockNumber", factory.GetFactory("BlockNumber").get());
      auto block_diagonalize = Misc::BlockDiagonalizeVectorFunctor(A, *BlockNumber, mergedA.getCrsGraph()->getImporter(), results, rowTranslation, colTranslation);

      runDroppingFunctorsImpl(lclA, blkPartSize, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, block_diagonalize, functors...);
    } else
      runDroppingFunctorsImpl(lclA, blkPartSize, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, functors...);
  }
};

}  // namespace MueLu
#endif
