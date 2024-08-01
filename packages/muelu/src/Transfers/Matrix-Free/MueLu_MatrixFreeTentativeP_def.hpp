// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATRIXFREETENTATIVEP_DEF_HPP
#define MUELU_MATRIXFREETENTATIVEP_DEF_HPP

#include "MueLu_MatrixFreeTentativeP_decl.hpp"

#include "MueLu_Aggregates.hpp"

#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "Teuchos_ScalarTraits.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Operator.hpp"

namespace MueLu {

// compute Y = alpha*R*X + beta*Y
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixFreeTentativeP<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector &X,
                                                                            MultiVector &Y,
                                                                            Teuchos::ETransp mode,
                                                                            Scalar alpha,
                                                                            Scalar beta) const {
  using impl_scalar_type     = typename Kokkos::ArithTraits<Scalar>::val_type;
  impl_scalar_type implAlpha = alpha;

  // Step 1: Y*=beta*Y, setup necessary structures
  Y.scale(beta);

  // TODO: probably smarter to sqrt the whole aggSizes once, but may be slower if it's done in a separate kernel launch?
  typename Aggregates::aggregates_sizes_type::const_type aggSizes = aggregates_->ComputeAggregateSizes();

  auto kokkos_view_X = X.getDeviceLocalView(Xpetra::Access::ReadOnly);
  auto kokkos_view_Y = Y.getDeviceLocalView(Xpetra::Access::ReadWrite);
  LO numCols         = kokkos_view_X.extent(1);

  if (mode == Teuchos::TRANS) {  // if we're in restrictor mode
    auto vertex2AggId     = aggregates_->GetVertex2AggId();
    auto vertex2AggIdView = vertex2AggId->getDeviceLocalView(Xpetra::Access::ReadOnly);
    LO numNodes           = kokkos_view_X.extent(0);

    // Step 2: Compute Y=Y+alpha*R*X
    // recall R*X is an average of X over aggregates
    Kokkos::parallel_for(
        "MueLu:MatrixFreeTentativeR_kokkos:apply", md_range_type({0, 0}, {numCols, numNodes}),
        KOKKOS_LAMBDA(const int colIdx, const int NodeIdx) {
          LO aggIdx = vertex2AggIdView(NodeIdx, 0);
          if (aggIdx != -1) {  // depending on maps, vertex2agg
            Kokkos::atomic_add(&kokkos_view_Y(aggIdx, colIdx), implAlpha * kokkos_view_X(NodeIdx, colIdx) / Kokkos::sqrt(aggSizes(aggIdx)));
          }
        });
  } else {  // if we're in prolongator mode
    const auto vertex2Agg = aggregates_->GetVertex2AggId();
    auto vertex2AggView   = vertex2Agg->getDeviceLocalView(Xpetra::Access::ReadOnly);
    LO numNodes           = kokkos_view_Y.extent(0);

    // Step 2: Compute Y=Y+alpha*P*X
    // recall P*X is essentially injection of X, but sum if a node belongs to multiple aggregates
    Kokkos::parallel_for(
        "MueLu:MatrixFreeTentativeP:apply", md_range_type({0, 0}, {numCols, numNodes}),
        KOKKOS_LAMBDA(const int colIdx, const int fineIdx) {
          LO aggIdx = vertex2AggView(fineIdx, 0);
          kokkos_view_Y(fineIdx, colIdx) += implAlpha * kokkos_view_X(aggIdx, colIdx) / Kokkos::sqrt(aggSizes(aggIdx));
        });
  }
}

// I don't care
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixFreeTentativeP<Scalar, LocalOrdinal, GlobalOrdinal, Node>::residual(const MultiVector &X, const MultiVector &B, MultiVector &R) const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MatrixFreeTentativeP residual would make no sense as the operator is not square!");
}

}  // namespace MueLu

#define MUELU_MATRIXFREETENTATIVEP_SHORT
#endif  // MUELU_MATRIXFREETENTATIVEP_DEF_HPP
