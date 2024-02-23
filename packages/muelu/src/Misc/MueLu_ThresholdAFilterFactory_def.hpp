// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ThresholdAFilterFactory(const std::string& ename, const Scalar threshold, const bool keepDiagonal, const GlobalOrdinal expectedNNZperRow)
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
