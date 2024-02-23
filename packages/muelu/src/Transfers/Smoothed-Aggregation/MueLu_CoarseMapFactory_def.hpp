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
/*
 * MueLu_CoarseMapFactory_def.hpp
 *
 *  Created on: Oct 12, 2012
 *      Author: wiesner
 */

#ifndef MUELU_COARSEMAPFACTORY_DEF_HPP_
#define MUELU_COARSEMAPFACTORY_DEF_HPP_

#include <Teuchos_Array.hpp>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_CoarseMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Generating factory for aggregates.");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory for null space.");

  validParamList->set<std::string>("Striding info", "{}", "Striding information");
  validParamList->set<LocalOrdinal>("Strided block id", -1, "Strided block id");

  // Domain GID offsets: list of offset values for domain gids (coarse gids) of tentative prolongator  (default = 0).
  // For each multigrid level we can define a different offset value, ie for the prolongator between
  // level 0 and level 1 we use the offset entry domainGidOffsets_[0], for the prolongator between
  // level 1 and level 2 we use domainGidOffsets_[1]...
  // If the vector domainGidOffsets_ is too short and does not contain a sufficient number of gid offset
  // values for all levels, a gid offset of 0 is used for all the remaining levels!
  // The GIDs for the domain dofs of Ptent start with domainGidOffset, are contiguous and distributed
  // equally over the procs (unless some reordering is done).
  validParamList->set<std::string>("Domain GID offsets", "{0}", "vector with offsets for GIDs for each level. If no offset GID value is given for the level we use 0 as default.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "Aggregates");
  Input(currentLevel, "Nullspace");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setStridingData(std::vector<size_t> stridingInfo) {
  // store striding map in internal variable
  stridingInfo_ = stridingInfo;

  // try to remove string "Striding info" from parameter list to make sure,
  // that stridingInfo_ is not replaced in the Build call.
  std::string strStridingInfo;
  strStridingInfo.clear();
  SetParameter("Striding info", ParameterEntry(strStridingInfo));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  GlobalOrdinal domainGIDOffset = GetDomainGIDOffset(currentLevel);
  BuildCoarseMap(currentLevel, domainGIDOffset);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCoarseMap(
    Level& currentLevel, const GlobalOrdinal domainGIDOffset) const {
  RCP<Aggregates> aggregates = Get<RCP<Aggregates> >(currentLevel, "Aggregates");
  GlobalOrdinal numAggs      = aggregates->GetNumAggregates();
  RCP<const Map> aggMap      = aggregates->GetMap();

  RCP<MultiVector> nullspace = Get<RCP<MultiVector> >(currentLevel, "Nullspace");

  const size_t NSDim                  = nullspace->getNumVectors();
  RCP<const Teuchos::Comm<int> > comm = aggMap->getComm();
  const ParameterList& pL             = GetParameterList();

  LocalOrdinal stridedBlockId = pL.get<LocalOrdinal>("Strided block id");

  // read in stridingInfo from parameter list and fill the internal member variable
  // read the data only if the parameter "Striding info" exists and is non-empty
  if (pL.isParameter("Striding info")) {
    std::string strStridingInfo = pL.get<std::string>("Striding info");
    if (strStridingInfo.empty() == false) {
      Teuchos::Array<size_t> arrayVal = Teuchos::fromStringToArray<size_t>(strStridingInfo);
      stridingInfo_                   = Teuchos::createVector(arrayVal);
    }
  }

  CheckForConsistentStridingInformation(stridedBlockId, NSDim);

  GetOStream(Statistics2) << "domainGIDOffset: " << domainGIDOffset << " block size: " << getFixedBlockSize() << " stridedBlockId: " << stridedBlockId << std::endl;

  // number of coarse level dofs (fixed by number of aggregates and blocksize data)
  GlobalOrdinal nCoarseDofs = numAggs * getFixedBlockSize();
  GlobalOrdinal indexBase   = aggMap->getIndexBase();

  RCP<const Map> coarseMap = StridedMapFactory::Build(aggMap->lib(),
                                                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                      nCoarseDofs,
                                                      indexBase,
                                                      stridingInfo_,
                                                      comm,
                                                      stridedBlockId,
                                                      domainGIDOffset);

  Set(currentLevel, "CoarseMap", coarseMap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetDomainGIDOffset(
    Level& currentLevel) const {
  GlobalOrdinal domainGidOffset = Teuchos::ScalarTraits<GlobalOrdinal>::zero();

  std::vector<GlobalOrdinal> domainGidOffsets;
  domainGidOffsets.clear();
  const ParameterList& pL = GetParameterList();
  if (pL.isParameter("Domain GID offsets")) {
    std::string strDomainGIDs = pL.get<std::string>("Domain GID offsets");
    if (strDomainGIDs.empty() == false) {
      Teuchos::Array<GlobalOrdinal> arrayVal = Teuchos::fromStringToArray<GlobalOrdinal>(strDomainGIDs);
      domainGidOffsets                       = Teuchos::createVector(arrayVal);
      if (currentLevel.GetLevelID() < Teuchos::as<int>(domainGidOffsets.size())) {
        domainGidOffset = domainGidOffsets[currentLevel.GetLevelID()];
      }
    }
  }

  return domainGidOffset;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckForConsistentStridingInformation(
    const LocalOrdinal stridedBlockId, const size_t nullspaceDimension) const {
  // check for consistency of striding information with NSDim and nCoarseDofs
  if (stridedBlockId == -1) {
    // this means we have no real strided map but only a block map with constant blockSize "nullspaceDimension"
    TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() > 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");
    stridingInfo_.clear();
    stridingInfo_.push_back(nullspaceDimension);
    TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo_.size() != 1, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): stridingInfo_.size() but must be one");

  } else {
    // stridedBlockId > -1, set by user
    TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId > Teuchos::as<LO>(stridingInfo_.size() - 1), Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): it is stridingInfo_.size() <= stridedBlockId_. error.");
    size_t stridedBlockSize = stridingInfo_[stridedBlockId];
    TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockSize != nullspaceDimension, Exceptions::RuntimeError, "MueLu::CoarseMapFactory::Build(): dimension of strided block != nullspaceDimension. error.");
  }
}

}  // namespace MueLu

#endif /* MUELU_COARSEMAPFACTORY_DEF_HPP_ */
