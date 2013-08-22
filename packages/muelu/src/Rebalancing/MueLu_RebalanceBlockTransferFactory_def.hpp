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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_REBALANCEBLOCKTRANSFERFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKTRANSFERFACTORY_DEF_HPP_

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

#include "MueLu_RebalanceBlockTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> RebalanceBlockTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    {
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
      RCP<validatorType> typeValidator = rcp (new validatorType(Teuchos::tuple<std::string>("Interpolation", "Restriction"), "type"));
      validParamList->set("type", "Interpolation", "Type of the transfer operator that need to be rebalanced (Interpolation or Restriction)", typeValidator);
    }

    validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");
    validParamList->set< RCP<const FactoryBase> >("R",              Teuchos::null, "Factory of the restriction operator that need to be rebalanced (only used if type=Restriction)");
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
    FactManager_.push_back(FactManager);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    const Teuchos::ParameterList & pL = GetParameterList();

    if (pL.get<std::string>("type") == "Interpolation") {
      Input(coarseLevel, "P");
    } else {
      Input(coarseLevel, "R");
      /*Input(coarseLevel, "Nullspace");
      if(pL.isParameter("Coordinates") && pL.get< RCP<const FactoryBase> >("Coordinates") != Teuchos::null)
      Input(coarseLevel, "Coordinates");*/
    }

    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
      SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

      coarseLevel.DeclareInput("Importer",(*it)->GetFactory("Importer").get(), this);

      if (pL.get<std::string>("type") == "Interpolation") {
        // no null space necessary
      } else {
        //Input(coarseLevel, "Nullspace");
        coarseLevel.DeclareInput("Nullspace",(*it)->GetFactory("Nullspace").get(), this);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);
    const Teuchos::ParameterList & pL = GetParameterList();

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    //std::cout << "\n This is just a test \n" << std::endl;

    Teuchos::RCP<Matrix> originalTransferOp = Teuchos::null;
    if (pL.get<std::string>("type") == "Interpolation") {
      originalTransferOp = Get< RCP<Matrix> >(coarseLevel, "P");
    } else {
      originalTransferOp = Get< RCP<Matrix> >(coarseLevel, "R");
    }

    RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > bOriginalTransferOp = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(originalTransferOp);
    TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp==Teuchos::null, Exceptions::BadCast, "MueLu::RebalanceBlockTransferFactory::Build: input matrix P or R is not of type BlockedCrsMatrix! error.");

    // plausibility check
    TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp->Rows() != 2,Exceptions::RuntimeError, "MueLu::RebalanceBlockTransferFactory::Build: number of block rows of transfer operator is not equal 2. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp->Cols() != 2,Exceptions::RuntimeError, "MueLu::RebalanceBlockTransferFactory::Build: number of block columns of transfer operator is not equal 2. error.");

    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    it = FactManager_.begin(); // we're only interested in the first factory manager

    if (pL.get<std::string>("type") == "Interpolation") {
      Teuchos::RCP<const Import> rebalanceImporter = Teuchos::null;
      { // begin SubFactoryManager environment
        SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
        SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

        rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());
        if(rebalanceImporter == Teuchos::null) {
            Set(coarseLevel, "P", originalTransferOp); // do nothing  // TODO remove this!
        } else {

          // extract matrix block
          Teuchos::RCP<CrsMatrix> Pm00 = bOriginalTransferOp->getMatrix(0, 0);
          Teuchos::RCP<CrsMatrixWrap> Pw00 = Teuchos::rcp(new CrsMatrixWrap(Pm00));
          Teuchos::RCP<Matrix> P00 = Teuchos::rcp_dynamic_cast<Matrix>(Pw00);

          // rebalance P11
          Teuchos::RCP<Matrix> rebP00;
          {
            SubFactoryMonitor m1(*this, "Rebalancing prolongator P(0,0)", coarseLevel);

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
              newImporter = ImportFactory::Build(rebalanceImporter->getTargetMap(), P00->getColMap());
              Pm00->replaceDomainMapAndImporter(rebalanceImporter->getTargetMap(), newImporter);
            }

            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*P00, "P(0,0) rebalanced:", params);

            ///////////////////////// EXPERIMENTAL
            // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
            // That is probably something for an external permutation factory
            //if(originalP->IsView("stridedMaps")) rebalancedP->CreateView("stridedMaps", originalP);
            ///////////////////////// EXPERIMENTAL
          } // end rebalance P(1,1)

          // extract matrix block
          Teuchos::RCP<CrsMatrix> Pm11 = bOriginalTransferOp->getMatrix(1, 1);
          Teuchos::RCP<CrsMatrixWrap> Pw11 = Teuchos::rcp(new CrsMatrixWrap(Pm11));
          Teuchos::RCP<Matrix> P11 = Teuchos::rcp_dynamic_cast<Matrix>(Pw11);

          // TODO do rebalancing of P(1,1)

          // rebuild rebalanced blocked P operator
          std::vector<GO> fullRangeMapVector;
          std::vector<GO> fullDomainMapVector;
          std::vector<RCP<const Map> > subBlockPRangeMaps;
          std::vector<RCP<const Map> > subBlockPDomainMaps;
          subBlockPRangeMaps.reserve(bOriginalTransferOp->Rows());       // reserve size for block P operators
          subBlockPDomainMaps.reserve(bOriginalTransferOp->Cols());       // reserve size for block P operators

          // append strided row map (= range map) to list of range maps.
          Teuchos::RCP<const Map> rangeMap00 = Teuchos::null;
          if(rebP00 != Teuchos::null) rangeMap00 = rebP00->getRangeMap(); //getRowMap("stridedMaps");
          else                        rangeMap00 = P00->getRangeMap();
          subBlockPRangeMaps.push_back(rangeMap00);
          Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap00 = rangeMap00->getNodeElementList(); //subBlockPRangeMaps.back()->getNodeElementList();
          fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap00.begin(), nodeRangeMap00.end());
          sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

          // append strided row map (= range map) to list of range maps.
          Teuchos::RCP<const Map> rangeMap11 = P11->getRangeMap(); //getRowMap("stridedMaps");
          subBlockPRangeMaps.push_back(rangeMap11);
          Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap11 = rangeMap11->getNodeElementList(); //subBlockPRangeMaps.back()->getNodeElementList();
          fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap11.begin(), nodeRangeMap11.end());
          sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

          // append strided col map (= domain map) to list of range maps.
          Teuchos::RCP<const Map> domainMap00 = Teuchos::null;
          if(rebP00 != Teuchos::null) domainMap00 = rebP00->getDomainMap(); //getColMap("stridedMaps");
          else                        domainMap00 = P00->getDomainMap();
          subBlockPDomainMaps.push_back(domainMap00);
          Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap00 = domainMap00->getNodeElementList(); //subBlockPDomainMaps.back()->getNodeElementList();
          fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap00.begin(), nodeDomainMap00.end());
          sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

          // append strided col map (= domain map) to list of range maps.
          Teuchos::RCP<const Map> domainMap11 = P11->getDomainMap(); //getColMap("stridedMaps");
          subBlockPDomainMaps.push_back(domainMap11);
          Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap11 = domainMap11->getNodeElementList(); //subBlockPDomainMaps.back()->getNodeElementList();
          fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap11.begin(), nodeDomainMap11.end());
          sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

          // extract map index base from maps of blocked P
          GO rangeIndexBase = originalTransferOp->getRangeMap()->getIndexBase();
          GO domainIndexBase= originalTransferOp->getDomainMap()->getIndexBase();

          RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangePMapExtractor = bOriginalTransferOp->getRangeMapExtractor(); // original map extractor
          Teuchos::ArrayView<GO> fullRangeMapGIDs(&fullRangeMapVector[0],fullRangeMapVector.size());
          Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangePMapExtractor->getFullMap());
          std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();

          Teuchos::RCP<const StridedMap > fullRangeMap =
              StridedMapFactory::Build(
                  originalTransferOp->getRangeMap()->lib(),
                  Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                  fullRangeMapGIDs,
                  rangeIndexBase,
                  stridedData,
                  originalTransferOp->getRangeMap()->getComm(),
                  stridedRgFullMap->getStridedBlockId(),
                  stridedRgFullMap->getOffset());

          RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainAMapExtractor = bOriginalTransferOp->getDomainMapExtractor();
          Teuchos::ArrayView<GO> fullDomainMapGIDs(&fullDomainMapVector[0],fullDomainMapVector.size());
          Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainAMapExtractor->getFullMap());
          std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();

          Teuchos::RCP<const StridedMap > fullDomainMap =
              StridedMapFactory::Build(
                  originalTransferOp->getDomainMap()->lib(),
                  Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                  fullDomainMapGIDs,
                  domainIndexBase,
                  stridedData2,
                  originalTransferOp->getDomainMap()->getComm(),
                  stridedDoFullMap->getStridedBlockId(),
                  stridedDoFullMap->getOffset());

          // build map extractors
          Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeMapExtractor  =
              Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullRangeMap,  subBlockPRangeMaps);
          Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainMapExtractor =
              Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullDomainMap, subBlockPDomainMaps);

          Teuchos::RCP<BlockedCrsMatrix> bRebP = Teuchos::rcp(new BlockedCrsMatrix(rangeMapExtractor,domainMapExtractor,10));
          /*for(size_t i = 0; i<subBlockPRangeMaps.size(); i++) {
            Teuchos::RCP<CrsMatrixWrapClass> crsOpii = Teuchos::rcp_dynamic_cast<CrsMatrixWrapClass>(subBlockP[i]);
            Teuchos::RCP<CrsMatrixClass> crsMatii = crsOpii->getCrsMatrix();
            bP->setMatrix(i,i,crsMatii);
          }*/
          if(!rebP00.is_null()) {
            Teuchos::RCP<const CrsMatrixWrap> rebPw00 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebP00);
            Teuchos::RCP<CrsMatrix> rebPm00 = rebPw00->getCrsMatrix();
            bRebP->setMatrix(0,0,rebPm00);
          } else {
            Teuchos::RCP<const CrsMatrixWrap> rebPw00 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(P00);
            Teuchos::RCP<CrsMatrix> rebPm00 = rebPw00->getCrsMatrix();
            bRebP->setMatrix(0,0,rebPm00);
          }

          {
            Teuchos::RCP<const CrsMatrixWrap> rebPw11 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(P11);
            Teuchos::RCP<CrsMatrix> rebPm11 = rebPw11->getCrsMatrix();
            bRebP->setMatrix(1,1,rebPm11);
          }

          bRebP->fillComplete();

          Set(coarseLevel, "P", Teuchos::rcp_dynamic_cast<Matrix>(bRebP)); // do nothing  // TODO remove this!

          //Set(coarseLevel, "P", originalTransferOp); // do nothing  // TODO remove this!
        }
      } // end SubFactoryManager environment
    } // end interpolation environment

    // begin restriction environment
    else {
      Teuchos::RCP<const Import> rebalanceImporter = Teuchos::null;
      { // begin SubFactoryManager environment
        SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
        SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

        rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());
        if(rebalanceImporter == Teuchos::null) {
            Set(coarseLevel, "R", originalTransferOp); // do nothing  // TODO remove this!

            // do nothing!
            RCP<MultiVector> nullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", (*it)->GetFactory("Nullspace").get());
            //std::cout << "input: nullspace generated by factory " << (*it)->GetFactory("Nullspace").get() << std::endl;
            coarseLevel.Set<RCP<MultiVector> >("Nullspace", nullspace, (*it)->GetFactory("Nullspace").get());
            //std::cout << "output: nullspace generated by factory " << (*it)->GetFactory("Nullspace").get() << std::endl;

        } else /* do rebalancing */
        {

          // TODO rebalance R
          //Set(coarseLevel, "R", originalTransferOp); // do nothing  // TODO remove this!

          // extract matrix block
          Teuchos::RCP<CrsMatrix> Rm00 = bOriginalTransferOp->getMatrix(0, 0);
          Teuchos::RCP<CrsMatrixWrap> Rw00 = Teuchos::rcp(new CrsMatrixWrap(Rm00));
          Teuchos::RCP<Matrix> R00 = Teuchos::rcp_dynamic_cast<Matrix>(Rw00);

          Teuchos::RCP<Matrix> rebR00;
          {
            SubFactoryMonitor m1(*this, "Rebalancing restriction R(0,0)", coarseLevel);
            {
              SubFactoryMonitor subM(*this, "Rebalancing restriction -- fusedImport", coarseLevel);
              // Note: The 3rd argument says to use originalR's domain map.

              RCP<Map> dummy;
              rebR00 = MatrixFactory::Build(R00,*rebalanceImporter,dummy,rebalanceImporter->getTargetMap());
            }
            //Set(coarseLevel, "R", rebR00);

            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebR00, "R(0,0) rebalanced:", params);
          }

          // extract matrix block
          Teuchos::RCP<CrsMatrix> Rm11 = bOriginalTransferOp->getMatrix(1, 1);
          Teuchos::RCP<CrsMatrixWrap> Rw11 = Teuchos::rcp(new CrsMatrixWrap(Rm11));
          Teuchos::RCP<Matrix> R11 = Teuchos::rcp_dynamic_cast<Matrix>(Rw11);

          //////////////////////////////

          // rebuild rebalanced blocked P operator
          std::vector<GO> fullRangeMapVector;
          std::vector<GO> fullDomainMapVector;
          std::vector<RCP<const Map> > subBlockRRangeMaps;
          std::vector<RCP<const Map> > subBlockRDomainMaps;
          subBlockRRangeMaps.reserve(bOriginalTransferOp->Rows());       // reserve size for block P operators
          subBlockRDomainMaps.reserve(bOriginalTransferOp->Cols());       // reserve size for block P operators

          // append strided row map (= range map) to list of range maps.
          Teuchos::RCP<const Map> rangeMap00 = Teuchos::null;
          if(rebR00 != Teuchos::null) rangeMap00 = rebR00->getRangeMap(); //getRowMap("stridedMaps");
          else                        rangeMap00 = R00->getRangeMap();
          subBlockRRangeMaps.push_back(rangeMap00);
          Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap00 = rangeMap00->getNodeElementList(); //subBlockPRangeMaps.back()->getNodeElementList();
          fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap00.begin(), nodeRangeMap00.end());
          sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

          // append strided row map (= range map) to list of range maps.
          Teuchos::RCP<const Map> rangeMap11 = R11->getRangeMap(); //getRowMap("stridedMaps");
          subBlockRRangeMaps.push_back(rangeMap11);
          Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap11 = rangeMap11->getNodeElementList(); //subBlockPRangeMaps.back()->getNodeElementList();
          fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap11.begin(), nodeRangeMap11.end());
          sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

          // append strided col map (= domain map) to list of range maps.
          Teuchos::RCP<const Map> domainMap00 = Teuchos::null;
          if(rebR00 != Teuchos::null) domainMap00 = rebR00->getDomainMap(); //getColMap("stridedMaps");
          else                        domainMap00 = R00->getDomainMap();
          subBlockRDomainMaps.push_back(domainMap00);
          Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap00 = domainMap00->getNodeElementList(); //subBlockPDomainMaps.back()->getNodeElementList();
          fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap00.begin(), nodeDomainMap00.end());
          sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

          // append strided col map (= domain map) to list of range maps.
          Teuchos::RCP<const Map> domainMap11 = R11->getDomainMap(); //getColMap("stridedMaps");
          subBlockRDomainMaps.push_back(domainMap11);
          Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap11 = domainMap11->getNodeElementList(); //subBlockPDomainMaps.back()->getNodeElementList();
          fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap11.begin(), nodeDomainMap11.end());
          sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

          // extract map index base from maps of blocked P
          GO rangeIndexBase = originalTransferOp->getRangeMap()->getIndexBase();  // ! use transposed of originalTransferOp!
          GO domainIndexBase= originalTransferOp->getDomainMap()->getIndexBase();

          // check this
          RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeRMapExtractor = bOriginalTransferOp->getRangeMapExtractor(); // original map extractor
          Teuchos::ArrayView<GO> fullRangeMapGIDs(&fullRangeMapVector[0],fullRangeMapVector.size());
          Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeRMapExtractor->getFullMap());
          std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();

          Teuchos::RCP<const StridedMap > fullRangeMap =
              StridedMapFactory::Build(
                  originalTransferOp->getRangeMap()->lib(),
                  Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                  fullRangeMapGIDs,
                  rangeIndexBase,
                  stridedData,
                  originalTransferOp->getRangeMap()->getComm(),
                  stridedRgFullMap->getStridedBlockId(),
                  stridedRgFullMap->getOffset());

          RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainAMapExtractor = bOriginalTransferOp->getDomainMapExtractor();
          Teuchos::ArrayView<GO> fullDomainMapGIDs(&fullDomainMapVector[0],fullDomainMapVector.size());
          Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainAMapExtractor->getFullMap());
          std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();

          Teuchos::RCP<const StridedMap > fullDomainMap =
              StridedMapFactory::Build(
                  originalTransferOp->getDomainMap()->lib(),
                  Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                  fullDomainMapGIDs,
                  domainIndexBase,
                  stridedData2,
                  originalTransferOp->getDomainMap()->getComm(),
                  stridedDoFullMap->getStridedBlockId(),
                  stridedDoFullMap->getOffset());

          // build map extractors
          Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeMapExtractor  =
              Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullRangeMap,  subBlockRRangeMaps);
          Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainMapExtractor =
              Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullDomainMap, subBlockRDomainMaps);

          Teuchos::RCP<BlockedCrsMatrix> bRebR = Teuchos::rcp(new BlockedCrsMatrix(rangeMapExtractor,domainMapExtractor,10));
          /*for(size_t i = 0; i<subBlockPRangeMaps.size(); i++) {
            Teuchos::RCP<CrsMatrixWrapClass> crsOpii = Teuchos::rcp_dynamic_cast<CrsMatrixWrapClass>(subBlockP[i]);
            Teuchos::RCP<CrsMatrixClass> crsMatii = crsOpii->getCrsMatrix();
            bP->setMatrix(i,i,crsMatii);
          }*/
          if(!rebR00.is_null()) {
            Teuchos::RCP<const CrsMatrixWrap> rebRw00 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebR00);
            Teuchos::RCP<CrsMatrix> rebRm00 = rebRw00->getCrsMatrix();
            bRebR->setMatrix(0,0,rebRm00);
          } else {
            Teuchos::RCP<const CrsMatrixWrap> rebRw00 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(R00);
            Teuchos::RCP<CrsMatrix> rebRm00 = rebRw00->getCrsMatrix();
            bRebR->setMatrix(0,0,rebRm00);
          }

          {
            Teuchos::RCP<const CrsMatrixWrap> rebRw11 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(R11);
            Teuchos::RCP<CrsMatrix> rebRm11 = rebRw11->getCrsMatrix();
            bRebR->setMatrix(1,1,rebRm11);
          }

          bRebR->fillComplete();

          Set(coarseLevel, "R", Teuchos::rcp_dynamic_cast<Matrix>(bRebR)); // do nothing  // TODO remove this!

          ////////////////////////////

          { // rebalance null space
            SubFactoryMonitor subM(*this, "Rebalancing nullspace", coarseLevel);

            //RCP<MultiVector>         nullspace = Get< RCP<MultiVector> >(coarseLevel, "Nullspace");
            RCP<MultiVector> nullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", (*it)->GetFactory("Nullspace").get());
            RCP<MultiVector> permutedNullspace = MultiVectorFactory::Build(rebalanceImporter->getTargetMap(), nullspace->getNumVectors());
            permutedNullspace->doImport(*nullspace, *rebalanceImporter, Xpetra::INSERT);

            // TODO think about this
            //if (pL.get<bool>("useSubcomm") == true) // TODO either useSubcomm is enabled everywhere or nowhere
            //permutedNullspace->replaceMap(permutedNullspace->getMap()->removeEmptyProcesses());

            coarseLevel.Set<RCP<MultiVector> >("Nullspace", permutedNullspace, (*it)->GetFactory("Nullspace").get());

          } // end rebalance null space

        }

      } // end SubFactoryManager environment

    } // end restriction environment



  } // Build

} // namespace MueLu


#endif /* MUELU_REBALANCEBLOCKTRANSFERFACTORY_DEF_HPP_ */
