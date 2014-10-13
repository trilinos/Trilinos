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
#ifndef MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_DEF_HPP_

#ifdef HAVE_MUELU_EXPERIMENTAL

#include <Teuchos_Tuple.hpp>

#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_RebalanceBlockInterpolationFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");
  //validParamList->set< RCP<const FactoryBase> >("R",              Teuchos::null, "Factory of the restriction operator that need to be rebalanced (only used if type=Restriction)");
  //validParamList->set< RCP<const FactoryBase> >("Nullspace",      Teuchos::null, "Factory of the nullspace that need to be rebalanced (only used if type=Restriction)");
  /*validParamList->set< RCP<const FactoryBase> >("Coordinates",    Teuchos::null, "Factory of the coordinates that need to be rebalanced (only used if type=Restriction)");
    validParamList->set< RCP<const FactoryBase> >("Importer",       Teuchos::null, "Factory of the importer object used for the rebalancing");
    // The value of "useSubcomm" paramter here must be the same as in RebalanceAcFactory
    validParamList->set< bool >                  ("useSubcomm",              true, "Construct subcommunicators");
    validParamList->set< int >                   ("write start",    -1, "first level at which coordinates should be written to file");
    validParamList->set< int >                   ("write end",      -1, "last level at which coordinates should be written to file");*/

  // TODO validation: "P" parameter valid only for type="Interpolation" and "R" valid only for type="Restriction". Like so:
  // if (paramList.isEntry("type") && paramList.get("type) == "Interpolation) {
  //     validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");

  return validParamList;
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
  FactManager_.push_back(FactManager);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {

  Input(coarseLevel, "P");

  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
    SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);
    coarseLevel.DeclareInput("Importer",(*it)->GetFactory("Importer").get(), this);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::RCP<Matrix> originalTransferOp = Teuchos::null;
  originalTransferOp = Get< RCP<Matrix> >(coarseLevel, "P");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bOriginalTransferOp =
    Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > (originalTransferOp);
  TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp==Teuchos::null, Exceptions::BadCast, "MueLu::RebalanceBlockTransferFactory::Build: input matrix P or R is not of type BlockedCrsMatrix! error.");

  // plausibility check
  TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp->Rows() != 2,Exceptions::RuntimeError, "MueLu::RebalanceBlockTransferFactory::Build: number of block rows of transfer operator is not equal 2. error.");
  TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp->Cols() != 2,Exceptions::RuntimeError, "MueLu::RebalanceBlockTransferFactory::Build: number of block columns of transfer operator is not equal 2. error.");

  // declare variables for maps of blocked rebalanced prolongation operator
  std::vector<GO> fullRangeMapVector;
  std::vector<GO> fullDomainMapVector;
  std::vector<RCP<const Map> > subBlockPRangeMaps;
  std::vector<RCP<const Map> > subBlockPDomainMaps;
  subBlockPRangeMaps.reserve(bOriginalTransferOp->Rows());       // reserve size for block P operators
  subBlockPDomainMaps.reserve(bOriginalTransferOp->Cols());       // reserve size for block P operators

  std::vector<Teuchos::RCP<Matrix> > subBlockRebP;
  subBlockRebP.reserve(bOriginalTransferOp->Rows());

  int curBlockId = 0;
  Teuchos::RCP<const Import> rebalanceImporter = Teuchos::null;
  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    // begin SubFactoryManager environment
    SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());

    // extract diagonal matrix block
    Teuchos::RCP<CrsMatrix> Pmii = bOriginalTransferOp->getMatrix(curBlockId, curBlockId);
    Teuchos::RCP<CrsMatrixWrap> Pwii = Teuchos::rcp(new CrsMatrixWrap(Pmii));
    Teuchos::RCP<Matrix> Pii = Teuchos::rcp_dynamic_cast<Matrix>(Pwii);

    // rebalance P11
    if(rebalanceImporter != Teuchos::null) {
      std::stringstream ss; ss << "Rebalancing prolongator block P(" << curBlockId << "," << curBlockId << ")";
      SubFactoryMonitor m1(*this, ss.str(), coarseLevel);

      // P is the transfer operator from the coarse grid to the fine grid.
      // P must transfer the data from the newly reordered coarse A to the (unchanged) fine A.
      // This means that the domain map (coarse) of P must be changed according to the new partition. The range map (fine) is kept unchanged.
      //
      // The domain map of P must match the range map of R.
      // See also note below about domain/range map of R and its implications for P.
      //
      // To change the domain map of P, P needs to be fillCompleted again with the new domain map.
      // To achieve this, P is copied into a new matrix that is not fill-completed.
      // The doImport() operation is just used here to make a copy of P: the importer is trivial and there is no data movement involved.
      // The reordering actually happens during the fillComplete() with domainMap == rebalanceImporter->getTargetMap().

      RCP<const Import> newImporter;
      {
        SubFactoryMonitor subM(*this, "Rebalancing prolongator  -- fast map replacement", coarseLevel);
        newImporter = ImportFactory::Build(rebalanceImporter->getTargetMap(), Pii->getColMap());
        Pmii->replaceDomainMapAndImporter(rebalanceImporter->getTargetMap(), newImporter);
      }

      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      std::stringstream ss2; ss2 << "P(" << curBlockId << "," << curBlockId << ") rebalanced:";
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*Pii, ss2.str(), params);

      // store rebalanced P block
      subBlockRebP.push_back(Pii);
    } // end rebalance P(1,1)
    else {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      std::stringstream ss; ss << "P(" << curBlockId << "," << curBlockId << ") not rebalanced:";
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*Pii, ss.str(), params);
      // store rebalanced P block
      subBlockRebP.push_back(Pii);
    }

    // fix striding information for rebalanced diagonal block Pii
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgPMapExtractor = bOriginalTransferOp->getRangeMapExtractor(); // original map extractor
    Teuchos::RCP<const StridedMap> orig_stridedRgMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rgPMapExtractor->getMap(Teuchos::as<size_t>(curBlockId)));
    Teuchos::RCP<const Map> stridedRgMap = Teuchos::null;
    if(orig_stridedRgMap != Teuchos::null) {
      std::vector<size_t> stridingData = orig_stridedRgMap->getStridingData();
      Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMapii = Pii->getRangeMap()->getNodeElementList();
      stridedRgMap = StridedMapFactory::Build(
          originalTransferOp->getRangeMap()->lib(),
          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
          nodeRangeMapii,
          Pii->getRangeMap()->getIndexBase(),
          stridingData,
          originalTransferOp->getRangeMap()->getComm(),
          orig_stridedRgMap->getStridedBlockId(),
          orig_stridedRgMap->getOffset());
    } else stridedRgMap = Pii->getRangeMap();
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > doPMapExtractor = bOriginalTransferOp->getDomainMapExtractor(); // original map extractor
    Teuchos::RCP<const StridedMap> orig_stridedDoMap = Teuchos::rcp_dynamic_cast<const StridedMap>(doPMapExtractor->getMap(Teuchos::as<size_t>(curBlockId)));
    Teuchos::RCP<const Map> stridedDoMap = Teuchos::null;
    if(orig_stridedDoMap != Teuchos::null) {
      std::vector<size_t> stridingData = orig_stridedDoMap->getStridingData();
      Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMapii = Pii->getDomainMap()->getNodeElementList();
      stridedDoMap = StridedMapFactory::Build(originalTransferOp->getDomainMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        nodeDomainMapii,
        Pii->getDomainMap()->getIndexBase(),
        stridingData,
        originalTransferOp->getDomainMap()->getComm(),
        orig_stridedDoMap->getStridedBlockId(),
        orig_stridedDoMap->getOffset());
    } else stridedDoMap = Pii->getDomainMap();

    TEUCHOS_TEST_FOR_EXCEPTION(stridedRgMap == Teuchos::null,Exceptions::RuntimeError, "MueLu::RebalanceBlockInterpolationFactory::Build: failed to generate striding information. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoMap == Teuchos::null,Exceptions::RuntimeError, "MueLu::RebalanceBlockInterpolationFactory::Build: failed to generate striding information. error.");

    // replace stridedMaps view in diagonal sub block
    if(Pii->IsView("stridedMaps")) Pii->RemoveView("stridedMaps");
      Pii->CreateView("stridedMaps", stridedRgMap, stridedDoMap);

    // append strided row map (= range map) to list of range maps.
    Teuchos::RCP<const Map> rangeMapii = Pii->getRowMap("stridedMaps"); //Pii->getRangeMap();
    subBlockPRangeMaps.push_back(rangeMapii);
    Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMapii = Pii->getRangeMap()->getNodeElementList();
    fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMapii.begin(), nodeRangeMapii.end());
    sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

    // append strided col map (= domain map) to list of range maps.
    Teuchos::RCP<const Map> domainMapii = Pii->getColMap("stridedMaps"); //Pii->getDomainMap();
    subBlockPDomainMaps.push_back(domainMapii);
    Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMapii = Pii->getDomainMap()->getNodeElementList();
    fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMapii.begin(), nodeDomainMapii.end());
    sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

    curBlockId++; // increase block id index

  } // end SubFactoryManager environment

  // extract map index base from maps of blocked P
  GO rangeIndexBase = originalTransferOp->getRangeMap()->getIndexBase();
  GO domainIndexBase= originalTransferOp->getDomainMap()->getIndexBase();

  RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangePMapExtractor = bOriginalTransferOp->getRangeMapExtractor(); // original map extractor
  Teuchos::ArrayView<GO> fullRangeMapGIDs(fullRangeMapVector.size() ? &fullRangeMapVector[0] : 0,fullRangeMapVector.size());
  Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangePMapExtractor->getFullMap());
  Teuchos::RCP<const Map > fullRangeMap = Teuchos::null;
  if(stridedRgFullMap != Teuchos::null) {
    std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
    fullRangeMap =
        StridedMapFactory::Build(
            originalTransferOp->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            stridedData,
            originalTransferOp->getRangeMap()->getComm(),
            stridedRgFullMap->getStridedBlockId(),
            stridedRgFullMap->getOffset());
  } else {
    fullRangeMap =
        MapFactory::Build(
            originalTransferOp->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            originalTransferOp->getRangeMap()->getComm());
  }

  RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainAMapExtractor = bOriginalTransferOp->getDomainMapExtractor();
  Teuchos::ArrayView<GO> fullDomainMapGIDs(fullDomainMapVector.size() ? &fullDomainMapVector[0] : 0,fullDomainMapVector.size());
  Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainAMapExtractor->getFullMap());
  Teuchos::RCP<const Map > fullDomainMap = Teuchos::null;
  if(stridedDoFullMap != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoFullMap==Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: full map in domain map extractor has no striding information! error.");
    std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
    fullDomainMap =
        StridedMapFactory::Build(
            originalTransferOp->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            stridedData2,
            originalTransferOp->getDomainMap()->getComm(),
            stridedDoFullMap->getStridedBlockId(),
            stridedDoFullMap->getOffset());
  } else {

    fullDomainMap =
        MapFactory::Build(
            originalTransferOp->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            originalTransferOp->getDomainMap()->getComm());
  }

  // build map extractors
  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeMapExtractor  =
      Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullRangeMap,  subBlockPRangeMaps);
  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainMapExtractor =
      Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullDomainMap, subBlockPDomainMaps);

  Teuchos::RCP<BlockedCrsMatrix> bRebP = Teuchos::rcp(new BlockedCrsMatrix(rangeMapExtractor,domainMapExtractor,10));
  for(size_t i = 0; i<subBlockPRangeMaps.size(); i++) {
     Teuchos::RCP<const CrsMatrixWrap> crsOpii = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(subBlockRebP[i]);
     Teuchos::RCP<CrsMatrix> crsMatii = crsOpii->getCrsMatrix();
     bRebP->setMatrix(i,i,crsMatii);
   }
  bRebP->fillComplete();

  Set(coarseLevel, "P", Teuchos::rcp_dynamic_cast<Matrix>(bRebP));


} // Build

} // namespace MueLu

#endif /* HAVE_MUELU_EXPERIMENTAL */
#endif /* MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_DEF_HPP_ */
