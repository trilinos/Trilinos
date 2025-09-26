// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SCALARDROPPINGDISTANCELAPLACIAN_DEF_HPP
#define MUELU_SCALARDROPPINGDISTANCELAPLACIAN_DEF_HPP

#include "MueLu_ScalarDroppingDistanceLaplacian_decl.hpp"

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure SoC>
void ScalarDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, SoC>::runDroppingFunctors_on_dlap(matrix_type& A,
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
                                                                                                                  const Factory& factory) {
  using doubleMultiVector = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
  auto coords             = level.template Get<Teuchos::RCP<doubleMultiVector>>("Coordinates", factory.GetFactory("Coordinates").get());
  if (distanceLaplacianMetric == "unweighted") {
    auto dist2 = DistanceLaplacian::UnweightedDistanceFunctor(A, coords);
    runDroppingFunctors_on_dlap_inner(A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
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
      auto dist2 = DistanceLaplacian::ScalarMaterialDistanceFunctor(A, coords, material);
      runDroppingFunctors_on_dlap_inner(A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
    } else {
      auto dist2 = DistanceLaplacian::TensorMaterialDistanceFunctor(A, coords, material);
      runDroppingFunctors_on_dlap_inner(A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, dist2, level, factory);
    }
  }
}
}  // namespace MueLu

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  MUELU_ETI_SLGN_SoC(MueLu::ScalarDroppingDistanceLaplacian, SC, LO, GO, NO)

#endif
