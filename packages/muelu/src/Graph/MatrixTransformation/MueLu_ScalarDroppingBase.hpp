// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SCALARDROPPINGBASE_HPP
#define MUELU_SCALARDROPPINGBASE_HPP

#include "Xpetra_Matrix.hpp"

#include "MueLu_DroppingCommon.hpp"
#include "MueLu_MatrixConstruction.hpp"
#include "MueLu_Factory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ScalarDroppingBase {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_matrix_type    = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType          = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename crs_matrix_type::local_matrix_device_type;
  using local_graph_type   = typename GraphType::local_graph_device_type;
  using rowptr_type        = typename local_graph_type::row_map_type::non_const_type;
  using device_type        = typename Node::device_type;
  using memory_space       = typename device_type::memory_space;
  using execution_space    = typename local_matrix_type::execution_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
  using Level              = MueLu::Level;
  using range_type         = Kokkos::RangePolicy<LocalOrdinal, execution_space>;
  using LocalOrdinalVector = Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;

  template <class... Functors>
  static void runDroppingFunctorsImpl(local_matrix_type& lclA, results_view& results, rowptr_type& filtered_rowptr, LocalOrdinal& nnz_filtered, Functors&... functors) {
    auto range = range_type(0, lclA.numRows());
#if !defined(HAVE_MUELU_DEBUG)
    auto countingFunctor = MatrixConstruction::PointwiseCountingFunctor(lclA, results, filtered_rowptr, functors...);

#else
    auto debug           = Misc::DebugFunctor(lclA, results);
    auto countingFunctor = MatrixConstruction::PointwiseCountingFunctor(lclA, results, filtered_rowptr, functors..., debug);
#endif
    Kokkos::parallel_scan("MueLu::CoalesceDrop::CountEntries", range, countingFunctor, nnz_filtered);
  }

  template <class... Functors>
  static void runDroppingFunctors(matrix_type& A, results_view& results, rowptr_type& filtered_rowptr, LocalOrdinal& nnz_filtered, const bool useBlocking, Level& level, const Factory& factory, Functors&... functors) {
    auto lclA = A.getLocalMatrixDevice();
    if (useBlocking) {
      auto BlockNumber       = level.template Get<Teuchos::RCP<LocalOrdinalVector>>("BlockNumber", factory.GetFactory("BlockNumber").get());
      auto block_diagonalize = Misc::BlockDiagonalizeFunctor(A, *BlockNumber, results);

      runDroppingFunctorsImpl(lclA, results, filtered_rowptr, nnz_filtered, block_diagonalize, functors...);
    } else
      runDroppingFunctorsImpl(lclA, results, filtered_rowptr, nnz_filtered, functors...);
  }
};

}  // namespace MueLu
#endif
