// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_VECTORDROPPINGDISTANCELAPLACIAN_DEF_HPP
#define MUELU_VECTORDROPPINGDISTANCELAPLACIAN_DEF_HPP

#include "MueLu_VectorDroppingDistanceLaplacian_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure SoC>
void VectorDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, SoC>::runDroppingFunctors_on_dlap(matrix_type& A,
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
                                                                                                                  const std::string& distanceLaplacianMetric,
                                                                                                                  Teuchos::Array<double>& dlap_weights,
                                                                                                                  LocalOrdinal interleaved_blocksize,
                                                                                                                  Level& level,
                                                                                                                  const Factory& factory) {
  using doubleMultiVector = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
  auto coords             = level.template Get<Teuchos::RCP<doubleMultiVector>>("Coordinates", factory.GetFactory("Coordinates").get());
  if (distanceLaplacianMetric == "unweighted") {
    auto dist2 = DistanceLaplacian::UnweightedDistanceFunctor(mergedA, coords);
    runDroppingFunctors_on_dlap_inner(A, mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
  } else if (distanceLaplacianMetric == "weighted") {
    auto k_dlap_weights_host = Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(&dlap_weights[0], dlap_weights.size());
    auto k_dlap_weights      = Kokkos::View<double*>("dlap_weights", k_dlap_weights_host.extent(0));
    Kokkos::deep_copy(k_dlap_weights, k_dlap_weights_host);
    auto dist2 = DistanceLaplacian::WeightedDistanceFunctor(mergedA, coords, k_dlap_weights);
    runDroppingFunctors_on_dlap_inner(A, mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
  } else if (distanceLaplacianMetric == "block weighted") {
    auto k_dlap_weights_host = Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(&dlap_weights[0], dlap_weights.size());
    auto k_dlap_weights      = Kokkos::View<double*>("dlap_weights", k_dlap_weights_host.extent(0));
    Kokkos::deep_copy(k_dlap_weights, k_dlap_weights_host);
    auto dist2 = DistanceLaplacian::BlockWeightedDistanceFunctor(mergedA, coords, k_dlap_weights, interleaved_blocksize);
    runDroppingFunctors_on_dlap_inner(A, mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
  } else if (distanceLaplacianMetric == "material") {
    auto material = level.template Get<Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>("Material", factory.GetFactory("Material").get());

    if (factory.IsPrint(Runtime0)) {
      auto spatialDim = coords->getNumVectors();
      if (material->getNumVectors() == 1) {
        factory.GetOStream(Runtime0) << "material scalar mean = " << material->getVector(0)->meanValue() << std::endl;
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(spatialDim * spatialDim != material->getNumVectors(), Exceptions::RuntimeError, "Need \"Material\" to have spatialDim^2 vectors.");
        {
          Teuchos::Array<Scalar> means(material->getNumVectors());
          material->meanValue(means());
          std::stringstream ss;
          ss << "material tensor mean =" << std::endl;
          size_t k = 0;
          for (size_t i = 0; i < spatialDim; ++i) {
            ss << "   ";
            for (size_t j = 0; j < spatialDim; ++j) {
              ss << means[k] << " ";
              ++k;
            }
            ss << std::endl;
          }
          factory.GetOStream(Runtime0) << ss.str();
        }
      }
    }

    if (material->getNumVectors() == 1) {
      auto dist2 = DistanceLaplacian::ScalarMaterialDistanceFunctor(mergedA, coords, material);
      runDroppingFunctors_on_dlap_inner(A, mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
    } else {
      auto dist2 = DistanceLaplacian::TensorMaterialDistanceFunctor(mergedA, coords, material);
      runDroppingFunctors_on_dlap_inner(A, mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
    }
  }
}
}  // namespace MueLu

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  MUELU_ETI_SLGN_SoC(MueLu::VectorDroppingDistanceLaplacian, SC, LO, GO, NO)

#endif
