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

#ifndef MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_
#define MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_BlockedCoarseMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_Factory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("Aggregates", Teuchos::null, "Generating factory for aggregates.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Generating factory for null space.");
  validParamList->set<RCP<const FactoryBase>>("CoarseMap", Teuchos::null, "Generating factory of previous coarse map. (must be set by user!).");

  // do we need this?
  validParamList->set<std::string>("Striding info", "{}", "Striding information");
  validParamList->set<LocalOrdinal>("Strided block id", -1, "Strided block id");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  this->Input(currentLevel, "Aggregates");
  this->Input(currentLevel, "Nullspace");

  // Get CoarseMap from previously defined block
  RCP<const FactoryBase> prevCoarseMapFact = this->GetFactory("CoarseMap");
  TEUCHOS_TEST_FOR_EXCEPTION(prevCoarseMapFact == Teuchos::null, Exceptions::RuntimeError, "MueLu::BlockedCoarseMapFactory::getDomainMapOffset: user did not specify CoarseMap of previous block. Do not forget to set the CoarseMap factory.");
  currentLevel.DeclareInput("CoarseMap", prevCoarseMapFact.get(), this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  GlobalOrdinal domainGIDOffset = GetDomainGIDOffset(currentLevel);
  CoarseMapFactory::BuildCoarseMap(currentLevel, domainGIDOffset);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetDomainGIDOffset(
    Level &currentLevel) const {
  RCP<const FactoryBase> prevCoarseMapFact = this->GetFactory("CoarseMap");
  RCP<const Map> subPDomainMap             = currentLevel.Get<RCP<const Map>>("CoarseMap", prevCoarseMapFact.get());
  GlobalOrdinal maxGlobalIndex             = subPDomainMap->getMaxAllGlobalIndex();

  return maxGlobalIndex + Teuchos::ScalarTraits<GlobalOrdinal>::one();
}

}  // namespace MueLu

#endif /* MUELU_BLOCKEDCOARSEMAPFACTORY_DEF_HPP_ */
