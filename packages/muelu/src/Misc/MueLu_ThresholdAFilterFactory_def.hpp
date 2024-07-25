// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_THRESHOLDAFILTERFACTORY_DEF_HPP
#define MUELU_THRESHOLDAFILTERFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_ThresholdAFilterFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ThresholdAFilterFactory(const std::string& ename, const magnitudeType threshold, const bool keepDiagonal, const GlobalOrdinal expectedNNZperRow)
  : varName_(ename)
  , threshold_(threshold)
  , keepDiagonal_(keepDiagonal)
  , expectedNNZperRow_(expectedNNZperRow) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, varName_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "A filter (thresholding)", currentLevel);

  RCP<Matrix> Ain = Get<RCP<Matrix> >(currentLevel, varName_);
  RCP<CrsMatrixWrap> Aout =
      MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetThresholdedMatrix(Ain, threshold_, keepDiagonal_, expectedNNZperRow_);

  GetOStream(Statistics0) << "Nonzeros in " << varName_ << "(input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering " << varName_ << " (parameter: " << threshold_ << "): " << Aout->getGlobalNumEntries() << std::endl;
  currentLevel.Set(varName_, Teuchos::rcp_dynamic_cast<Matrix>(Aout), this);
}

}  // namespace MueLu

#endif  // MUELU_THRESHOLDAFILTERFACTORY_DEF_HPP
