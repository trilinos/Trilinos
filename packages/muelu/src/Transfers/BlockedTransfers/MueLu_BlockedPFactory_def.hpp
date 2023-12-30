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

#ifndef MUELU_BLOCKEDPFACTORY_DEF_HPP_
#define MUELU_BLOCKEDPFACTORY_DEF_HPP_

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>

#include "MueLu_BlockedPFactory_decl.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlockedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A (block matrix)");
  validParamList->set<bool>("backwards", false, "Forward/backward order");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "A");

  const ParameterList& pL = GetParameterList();
  const bool backwards    = pL.get<bool>("backwards");

  const int numFactManagers = FactManager_.size();
  for (int k = 0; k < numFactManagers; k++) {
    int i                                            = (backwards ? numFactManagers - 1 - k : k);
    const RCP<const FactoryManagerBase>& factManager = FactManager_[i];

    SetFactoryManager fineSFM(rcpFromRef(fineLevel), factManager);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), factManager);

    if (!restrictionMode_)
      coarseLevel.DeclareInput("P", factManager->GetFactory("P").get(), this);
    else
      coarseLevel.DeclareInput("R", factManager->GetFactory("R").get(), this);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& /* fineLevel */, Level& /* coarseLevel */) const {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  RCP<Matrix> Ain = Get<RCP<Matrix> >(fineLevel, "A");

  RCP<BlockedCrsMatrix> A = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);
  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), Exceptions::BadCast, "Input matrix A is not a BlockedCrsMatrix.");

  const int numFactManagers = FactManager_.size();

  // Plausibility check
  TEUCHOS_TEST_FOR_EXCEPTION(A->Rows() != as<size_t>(numFactManagers), Exceptions::RuntimeError,
                             "Number of block rows [" << A->Rows() << "] does not match the number of SubFactorManagers [" << numFactManagers << "]");
  TEUCHOS_TEST_FOR_EXCEPTION(A->Cols() != as<size_t>(numFactManagers), Exceptions::RuntimeError,
                             "Number of block cols [" << A->Cols() << "] does not match the number of SubFactorManagers [" << numFactManagers << "]");

  // Build blocked prolongator
  std::vector<RCP<Matrix> > subBlockP(numFactManagers);
  std::vector<RCP<const Map> > subBlockPRangeMaps(numFactManagers);
  std::vector<RCP<const Map> > subBlockPDomainMaps(numFactManagers);

  std::vector<GO> fullRangeMapVector;
  std::vector<GO> fullDomainMapVector;

  RCP<const MapExtractor> rangeAMapExtractor  = A->getRangeMapExtractor();
  RCP<const MapExtractor> domainAMapExtractor = A->getDomainMapExtractor();

  const ParameterList& pL = GetParameterList();
  const bool backwards    = pL.get<bool>("backwards");

  // Build and store the subblocks and the corresponding range and domain
  // maps.  Since we put together the full range and domain map from the
  // submaps, we do not have to use the maps from blocked A
  for (int k = 0; k < numFactManagers; k++) {
    int i = (backwards ? numFactManagers - 1 - k : k);
    if (restrictionMode_)
      GetOStream(Runtime1) << "Generating R for block " << k << "/" << numFactManagers << std::endl;
    else
      GetOStream(Runtime1) << "Generating P for block " << k << "/" << numFactManagers << std::endl;

    const RCP<const FactoryManagerBase>& factManager = FactManager_[i];

    SetFactoryManager fineSFM(rcpFromRef(fineLevel), factManager);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), factManager);

    if (!restrictionMode_)
      subBlockP[i] = coarseLevel.Get<RCP<Matrix> >("P", factManager->GetFactory("P").get());
    else
      subBlockP[i] = coarseLevel.Get<RCP<Matrix> >("R", factManager->GetFactory("R").get());

    // Check if prolongator/restrictor operators have strided maps
    TEUCHOS_TEST_FOR_EXCEPTION(subBlockP[i]->IsView("stridedMaps") == false, Exceptions::BadCast,
                               "subBlock P operator [" << i << "] has no strided map information.");

    // Append strided row map (= range map) to list of range maps.
    // TAW the row map GIDs extracted here represent the actual row map GIDs.
    // No support for nested operators! (at least if Thyra style gids are used)
    RCP<const StridedMap> strPartialMap = Teuchos::rcp_dynamic_cast<const StridedMap>(subBlockP[i]->getRowMap("stridedMaps"));
    std::vector<size_t> stridedRgData   = strPartialMap->getStridingData();
    subBlockPRangeMaps[i]               = StridedMapFactory::Build(
                      subBlockP[i]->getRangeMap(),  // actual global IDs (Thyra or Xpetra)
                      stridedRgData,
                      strPartialMap->getStridedBlockId(),
                      strPartialMap->getOffset());
    // subBlockPRangeMaps[i] = subBlockP[i]->getRowMap("stridedMaps");

    // Use plain range map to determine the DOF ids
    ArrayView<const GO> nodeRangeMap = subBlockPRangeMaps[i]->getLocalElementList();
    fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap.begin(), nodeRangeMap.end());

    // Append strided col map (= domain map) to list of range maps.
    RCP<const StridedMap> strPartialMap2 = Teuchos::rcp_dynamic_cast<const StridedMap>(subBlockP[i]->getColMap("stridedMaps"));
    std::vector<size_t> stridedRgData2   = strPartialMap2->getStridingData();
    subBlockPDomainMaps[i]               = StridedMapFactory::Build(
                      subBlockP[i]->getDomainMap(),  // actual global IDs (Thyra or Xpetra)
                      stridedRgData2,
                      strPartialMap2->getStridedBlockId(),
                      strPartialMap2->getOffset());
    // subBlockPDomainMaps[i] = subBlockP[i]->getColMap("stridedMaps");

    // Use plain domain map to determine the DOF ids
    ArrayView<const GO> nodeDomainMap = subBlockPDomainMaps[i]->getLocalElementList();
    fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap.begin(), nodeDomainMap.end());
  }

  // check if GIDs for full maps have to be sorted:
  // For the Thyra mode ordering they do not have to be sorted since the GIDs are
  // numbered as 0...n1,0...,n2 (starting with zero for each subblock). The MapExtractor
  // generates unique GIDs during the construction.
  // For Xpetra style, the GIDs have to be reordered. Such that one obtains a ordered
  // list of GIDs in an increasing ordering. In Xpetra, the GIDs are all unique through
  // out all submaps.
  bool bDoSortRangeGIDs  = rangeAMapExtractor->getThyraMode() ? false : true;
  bool bDoSortDomainGIDs = domainAMapExtractor->getThyraMode() ? false : true;

  if (bDoSortRangeGIDs) {
    sort(fullRangeMapVector.begin(), fullRangeMapVector.end());
    fullRangeMapVector.erase(std::unique(fullRangeMapVector.begin(), fullRangeMapVector.end()), fullRangeMapVector.end());
  }
  if (bDoSortDomainGIDs) {
    sort(fullDomainMapVector.begin(), fullDomainMapVector.end());
    fullDomainMapVector.erase(std::unique(fullDomainMapVector.begin(), fullDomainMapVector.end()), fullDomainMapVector.end());
  }

  // extract map index base from maps of blocked A
  GO rangeIndexBase  = 0;
  GO domainIndexBase = 0;
  if (!restrictionMode_) {
    // Prolongation mode: just use index base of range and domain map of A
    rangeIndexBase  = A->getRangeMap()->getIndexBase();
    domainIndexBase = A->getDomainMap()->getIndexBase();

  } else {
    // Restriction mode: switch range and domain map for blocked restriction operator
    rangeIndexBase  = A->getDomainMap()->getIndexBase();
    domainIndexBase = A->getRangeMap()->getIndexBase();
  }

  // Build full range map.
  // If original range map has striding information, then transfer it to the
  // new range map
  RCP<const StridedMap> stridedRgFullMap = rcp_dynamic_cast<const StridedMap>(rangeAMapExtractor->getFullMap());
  RCP<const Map> fullRangeMap            = Teuchos::null;

  ArrayView<GO> fullRangeMapGIDs(fullRangeMapVector.size() ? &fullRangeMapVector[0] : 0, fullRangeMapVector.size());
  if (stridedRgFullMap != Teuchos::null) {
    std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
    fullRangeMap                    = StridedMapFactory::Build(
                           A->getRangeMap()->lib(),
                           Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                           fullRangeMapGIDs,
                           rangeIndexBase,
                           stridedData,
                           A->getRangeMap()->getComm(),
                           -1, /* the full map vector should always have strided block id -1! */
                           stridedRgFullMap->getOffset());
  } else {
    fullRangeMap = MapFactory::Build(
        A->getRangeMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        fullRangeMapGIDs,
        rangeIndexBase,
        A->getRangeMap()->getComm());
  }

  ArrayView<GO> fullDomainMapGIDs(fullDomainMapVector.size() ? &fullDomainMapVector[0] : 0, fullDomainMapVector.size());
  RCP<const StridedMap> stridedDoFullMap = rcp_dynamic_cast<const StridedMap>(domainAMapExtractor->getFullMap());
  RCP<const Map> fullDomainMap           = null;
  if (stridedDoFullMap != null) {
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoFullMap == Teuchos::null, Exceptions::BadCast,
                               "MueLu::BlockedPFactory::Build: full map in domain map extractor has no striding information!");

    std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
    fullDomainMap                    = StridedMapFactory::Build(
                           A->getDomainMap()->lib(),
                           Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                           fullDomainMapGIDs,
                           domainIndexBase,
                           stridedData2,
                           A->getDomainMap()->getComm(),
                           -1, /* the full map vector should always have strided block id -1! */
                           stridedDoFullMap->getOffset());
  } else {
    fullDomainMap = MapFactory::Build(
        A->getDomainMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        fullDomainMapGIDs,
        domainIndexBase,
        A->getDomainMap()->getComm());
  }

  // Build map extractors
  RCP<const MapExtractor> rangeMapExtractor  = MapExtractorFactory::Build(fullRangeMap, subBlockPRangeMaps, rangeAMapExtractor->getThyraMode());
  RCP<const MapExtractor> domainMapExtractor = MapExtractorFactory::Build(fullDomainMap, subBlockPDomainMaps, domainAMapExtractor->getThyraMode());

  RCP<BlockedCrsMatrix> P = rcp(new BlockedCrsMatrix(rangeMapExtractor, domainMapExtractor, 10));
  for (size_t i = 0; i < subBlockPRangeMaps.size(); i++)
    for (size_t j = 0; j < subBlockPRangeMaps.size(); j++)
      if (i == j) {
        RCP<CrsMatrixWrap> crsOpii = rcp_dynamic_cast<CrsMatrixWrap>(subBlockP[i]);
        TEUCHOS_TEST_FOR_EXCEPTION(crsOpii == Teuchos::null, Xpetra::Exceptions::BadCast,
                                   "Block [" << i << "," << j << "] must be of type CrsMatrixWrap.");
        P->setMatrix(i, i, crsOpii);
      } else {
        P->setMatrix(i, j, Teuchos::null);
      }

  P->fillComplete();

  // Level Set
  if (!restrictionMode_) {
    // Prolongation mode
    coarseLevel.Set("P", rcp_dynamic_cast<Matrix>(P), this);
    // Stick the CoarseMap on the level if somebody wants it (useful for repartitioning)
    coarseLevel.Set("CoarseMap", P->getBlockedDomainMap(), this);
  } else {
    // Restriction mode
    // We do not have to transpose the blocked R operator since the subblocks
    // on the diagonal are already valid R subblocks
    coarseLevel.Set("R", Teuchos::rcp_dynamic_cast<Matrix>(P), this);
  }
}

}  // namespace MueLu

#endif /* MUELU_BLOCKEDPFACTORY_DEF_HPP_ */
