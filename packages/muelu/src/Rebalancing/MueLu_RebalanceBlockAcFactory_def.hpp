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
/*
 * MueLu_RebalanceBlockAcFactory_def.hpp
 *
 *  Created on: Aug 15, 2013
 *      Author: tobias
 */

#ifndef MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_RAPFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
#include "MueLu_RepartitionFactory.hpp"
//#include "MueLu_RebalanceTransferFactory.hpp"
//#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_IsorropiaInterface.hpp"
//#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#endif

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RebalanceBlockAcFactory() {  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",         Teuchos::null, "Generating factory of the matrix A for rebalancing");
    //validParamList->set< bool >                  ("useSubcomm",         true, "Construct subcommunicators");
    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
    FactManager_.push_back(FactManager);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(coarseLevel, "A");

    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
      SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

      coarseLevel.DeclareInput("Importer",(*it)->GetFactory("Importer").get(), this);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Computing blocked Ac", coarseLevel);

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    RCP<Matrix> originalAc = Get< RCP<Matrix> >(coarseLevel, "A");

    RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > bA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(originalAc);
    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

    // plausibility check
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != 2,Exceptions::RuntimeError, "MueLu::BlockedPFactory::Build: number of block rows of A is not equal 2. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols() != 2,Exceptions::RuntimeError, "MueLu::BlockedPFactory::Build: number of block columns of A is not equal 2. error.");

    // store map extractors
    Teuchos::RCP<const MapExtractorClass> rangeMapExtractor  = bA->getRangeMapExtractor();
    Teuchos::RCP<const MapExtractorClass> domainMapExtractor = bA->getDomainMapExtractor();


    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;

#if 1

    std::vector<GO> fullRangeMapVector;
    std::vector<GO> fullDomainMapVector;

    std::vector<RCP<const Map> > subBlockARangeMaps;
    std::vector<RCP<const Map> > subBlockADomainMaps;
    subBlockARangeMaps.reserve(bA->Rows());
    subBlockADomainMaps.reserve(bA->Cols());

    std::vector<Teuchos::RCP<Matrix> > subBlockRebA;
    subBlockRebA.reserve(bA->Cols() * bA->Rows());

    for(int i=0; i<bA->Rows(); i++) {
      for(int j=0; j<bA->Cols(); j++) {
        // extract matrix block
        Teuchos::RCP<CrsMatrix> Amij = bA->getMatrix(i, j);
        Teuchos::RCP<CrsMatrixWrap> Awij = Teuchos::rcp(new CrsMatrixWrap(Amij));
        Teuchos::RCP<Matrix> Aij = Teuchos::rcp_dynamic_cast<Matrix>(Awij);
        //subBlockRebA[i*bA->Cols() + j] = Aij;
        subBlockRebA.push_back(Aij);
      }
    }

    int curBlockId = 0;

    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
      SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

      Teuchos::RCP<const Import> rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());

      // rebalance diagonal block

      // extract matrix block
      Teuchos::RCP<Matrix> Aii = subBlockRebA[curBlockId*bA->Cols() + curBlockId];

      Teuchos::RCP<Matrix> rebAii;
      if(rebalanceImporter != Teuchos::null) {
        std::stringstream ss; ss << "Rebalancing matrix block A(" << curBlockId << "," << curBlockId << ")";
        SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
        RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

        const ParameterList & pL = GetParameterList();

        ParameterList XpetraList;
        //if (pL.get<bool>("useSubcomm") == true) {
          //GetOStream(Runtime0,0) << "Replacing maps with a subcommunicator" << std::endl;
          XpetraList.set("Restrict Communicator",false /*true*/ /*XXX*/);
        //}
        // NOTE: If the communicator is restricted away, Build returns Teuchos::null.
        rebAii = MatrixFactory::Build(Aii, *rebalanceImporter, targetMap, targetMap, rcp(&XpetraList,false));

        if (!rebAii.is_null())
          rebAii->SetFixedBlockSize(Aii->GetFixedBlockSize());

        if (!rebAii.is_null()) {
          RCP<ParameterList> params = rcp(new ParameterList());
          params->set("printLoadBalancingInfo", true);
          std::stringstream ss2; ss2 << "A(" << curBlockId << "," << curBlockId << ") rebalanced:";
          GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAii, ss2.str(), params);
        }
      }  // rebalance matrix block A(i,i)
      else {
        rebAii = Aii;
        /*RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        std::stringstream ss2; ss2 << "A(" << curBlockId << "," << curBlockId << ") not rebalanced:";
        GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAii, ss2.str(), params);*/
      }

      subBlockRebA[curBlockId*bA->Cols() + curBlockId] = rebAii;

      // rebalance off-diagonal matrix blocks in same row
      for(int j=0; j<bA->Cols(); j++) {
        if(j==curBlockId) continue;  // jump over block diagonal matrix block

        // extract matrix block
        Teuchos::RCP<Matrix> Aij = subBlockRebA[curBlockId*bA->Cols() + j];

        Teuchos::RCP<Matrix> rebAij;
        if(rebalanceImporter!=Teuchos::null) {
          std::stringstream ss3; ss3 << "Rebalancing matrix block A(" << curBlockId << "," << j << ")";
          SubFactoryMonitor subM(*this, ss3.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          const ParameterList & pL = GetParameterList();

          Teuchos::RCP<Map> dummy; // The 3rd argument says to use the original domain map
          rebAij = MatrixFactory::Build(Aij, *rebalanceImporter, dummy, targetMap);

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebAij.is_null() && Aij->IsView("stridedMaps"))
            rebAij->CreateView("stridedMaps", Aij);

          if (!rebAij.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            std::stringstream ss4; ss4 << "A(" << curBlockId << "," << j << ") rebalanced:";
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss4.str(), params);
          }
        } // rebalance matrix block A(i,j)
        else {
          rebAij = Aij;
          /*RCP<ParameterList> params = rcp(new ParameterList());
          params->set("printLoadBalancingInfo", true);
          std::stringstream ss2; ss2 << "A(" << curBlockId << "," << j << ") not rebalanced:";
          GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss2.str(), params);*/
        }

        subBlockRebA[curBlockId*bA->Cols() + j] = rebAij;
      } // end loop over all columns

      // rebalance off-diagonal matrix blocks in same column
      for(int i=0; i<bA->Rows(); i++) {
        if(i==curBlockId) continue;  // jump over block diagonal matrix block

        // extract matrix block
        Teuchos::RCP<Matrix> Aij = subBlockRebA[i*bA->Cols() + curBlockId];

        Teuchos::RCP<Matrix> rebAij;
        if(rebalanceImporter!=Teuchos::null) {
          std::stringstream ss; ss << "Rebalancing matrix block (" << i << "," << curBlockId << ")";
          SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          rebAij = Aij; // just a copy
          Teuchos::RCP<const CrsMatrixWrap> rebAwij = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebAij);
          Teuchos::RCP<CrsMatrix> rebAmij = rebAwij->getCrsMatrix();
          Teuchos::RCP<const Import> rebAijImport = ImportFactory::Build(targetMap,Aij->getColMap());
          rebAmij->replaceDomainMapAndImporter(targetMap,rebAijImport);

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebAij.is_null() && Aij->IsView("stridedMaps"))
            rebAij->CreateView("stridedMaps", Aij);

          if (!rebAij.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            std::stringstream ss2; ss2 << "A(" << i << "," << curBlockId << ") rebalanced:";
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss2.str(), params);
          }
        } // rebalance matrix block A(1,0)
        else {
          rebAij = Aij;
          /*RCP<ParameterList> params = rcp(new ParameterList());
          params->set("printLoadBalancingInfo", true);
          std::stringstream ss2; ss2 << "A(" << i << "," << curBlockId << ") not rebalanced:";
          GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebAij, ss2.str(), params);*/
        }

        subBlockRebA[i*bA->Cols() + curBlockId] = rebAij;
      } // end loop over all rows


      // build full range and full domain map (strided)
      subBlockARangeMaps.push_back(rebAii->getRangeMap());
      Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap = subBlockARangeMaps.back()->getNodeElementList();
      fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap.begin(), nodeRangeMap.end());
      sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

      subBlockADomainMaps.push_back(rebAii->getDomainMap());
      Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap = subBlockADomainMaps.back()->getNodeElementList();
      fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap.begin(), nodeDomainMap.end());
      sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

      curBlockId++;
    } // end loop over all block rows

    // now, subBlockRebA contains all rebalanced matrix blocks

    // extract map index base from maps of blocked A
    GO rangeIndexBase  = 0;
    GO domainIndexBase = 0;
    rangeIndexBase = bA->getRangeMap()->getIndexBase();
    domainIndexBase= bA->getDomainMap()->getIndexBase();

    Teuchos::ArrayView<GO> fullRangeMapGIDs(&fullRangeMapVector[0],fullRangeMapVector.size());
    Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getFullMap());
    std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
    Teuchos::RCP<const StridedMap > fullRangeMap =
        StridedMapFactory::Build(
            bA->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            stridedData,
            bA->getRangeMap()->getComm(),
            stridedRgFullMap->getStridedBlockId(),
            stridedRgFullMap->getOffset());

    Teuchos::ArrayView<GO> fullDomainMapGIDs(&fullDomainMapVector[0],fullDomainMapVector.size());

    Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getFullMap());
    std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
    Teuchos::RCP<const StridedMap > fullDomainMap =
        StridedMapFactory::Build(
            bA->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            stridedData2,
            bA->getDomainMap()->getComm(),
            stridedDoFullMap->getStridedBlockId(),
            stridedDoFullMap->getOffset());

    // build map extractors
    Teuchos::RCP<const MapExtractorClass> rebRangeMapExtractor  = MapExtractorFactoryClass::Build(fullRangeMap, subBlockARangeMaps);
    Teuchos::RCP<const MapExtractorClass> rebDomainMapExtractor = MapExtractorFactoryClass::Build(fullDomainMap, subBlockADomainMaps);

    Teuchos::RCP<BlockedCrsMatrix> reb_bA = Teuchos::rcp(new BlockedCrsMatrix(rebRangeMapExtractor,rebDomainMapExtractor,10));
    for(int i=0; i<bA->Rows(); i++) {
      for(int j=0; j<bA->Cols(); j++) {
       Teuchos::RCP<const CrsMatrixWrap> crsOpij = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(subBlockRebA[i*bA->Cols() + j]);
       Teuchos::RCP<CrsMatrix> crsMatij = crsOpij->getCrsMatrix();
       reb_bA->setMatrix(i,j,crsMatij);
      }
    }
    reb_bA->fillComplete();
    //reb_bA->describe(*out,Teuchos::VERB_EXTREME);
    coarseLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(reb_bA), this);

#else
    it = FactManager_.begin(); // we're only interested in the first factory manager

    Teuchos::RCP<const Import> rebalanceImporter = Teuchos::null;
    {
      SetFactoryManager fineSFM  (rcpFromRef(fineLevel),   *it);
      SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

      rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());
      if(rebalanceImporter != Teuchos::null) {

        // extract matrix block
        Teuchos::RCP<CrsMatrix> Am00 = bA->getMatrix(0, 0);
        Teuchos::RCP<CrsMatrixWrap> Aw00 = Teuchos::rcp(new CrsMatrixWrap(Am00));
        Teuchos::RCP<Matrix> A00 = Teuchos::rcp_dynamic_cast<Matrix>(Aw00);

        Teuchos::RCP<Matrix> rebA00;
        {
          std::stringstream ss;
          ss << "Rebalancing matrix block (" << 0 << "," << 0 << ")";
          SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          const ParameterList & pL = GetParameterList();

          ParameterList XpetraList;
          //if (pL.get<bool>("useSubcomm") == true) {
            //GetOStream(Runtime0,0) << "Replacing maps with a subcommunicator" << std::endl;
            XpetraList.set("Restrict Communicator",false /*true*/ /*XXX*/);
          //}
          // NOTE: If the communicator is restricted away, Build returns Teuchos::null.
          rebA00 = MatrixFactory::Build(A00, *rebalanceImporter, targetMap, targetMap, rcp(&XpetraList,false));

          if (!rebA00.is_null())
            rebA00->SetFixedBlockSize(A00->GetFixedBlockSize());

          // TODO reset block matrix
          //Set(coarseLevel, "A", rebalancedAc);

          if (!rebA00.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebA00, "A(0,0) rebalanced:", params);
          }
        }  // rebalance matrix block A(0,0)

        // extract matrix block
        Teuchos::RCP<CrsMatrix> Am01 = bA->getMatrix(0, 1);
        Teuchos::RCP<CrsMatrixWrap> Aw01 = Teuchos::rcp(new CrsMatrixWrap(Am01));
        Teuchos::RCP<Matrix> A01 = Teuchos::rcp_dynamic_cast<Matrix>(Aw01);

        Teuchos::RCP<Matrix> rebA01;
        {
          std::stringstream ss;
          ss << "Rebalancing matrix block (" << 0 << "," << 1 << ")";
          SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          const ParameterList & pL = GetParameterList();

          /*ParameterList XpetraList;
          if (pL.get<bool>("useSubcomm") == true) {
            GetOStream(Runtime0,0) << "Replacing maps with a subcommunicator" << std::endl;
            XpetraList.set("Restrict Communicator",true);
          }*/
          // NOTE: If the communicator is restricted away, Build returns Teuchos::null.
          Teuchos::RCP<Map> dummy; // The 3rd argument says to use the original domain map
          rebA01 = MatrixFactory::Build(A01, *rebalanceImporter, dummy, targetMap/*, rcp(&XpetraList,false)*/);

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebA01.is_null() && A01->IsView("stridedMaps"))
            rebA01->CreateView("stridedMaps", A01);

          if (!rebA01.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebA01, "A(0,1) rebalanced:", params);
          }
        } // rebalance matrix block A(0,1)

        // extract matrix block
        Teuchos::RCP<CrsMatrix> Am10 = bA->getMatrix(1, 0);
        Teuchos::RCP<CrsMatrixWrap> Aw10 = Teuchos::rcp(new CrsMatrixWrap(Am10));
        Teuchos::RCP<Matrix> A10 = Teuchos::rcp_dynamic_cast<Matrix>(Aw10);

        Teuchos::RCP<Matrix> rebA10;
        {
          std::stringstream ss;
          ss << "Rebalancing matrix block (" << 1 << "," << 0 << ")";
          SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          rebA10 = A10; // just a copy
          Teuchos::RCP<const CrsMatrixWrap> rebAw10 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebA10);
          Teuchos::RCP<CrsMatrix> rebAm10 = rebAw10->getCrsMatrix();
          Teuchos::RCP<const Import> rebA10Import = ImportFactory::Build(targetMap,A10->getColMap());
          rebAm10->replaceDomainMapAndImporter(targetMap,rebA10Import);

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebA10.is_null() && A10->IsView("stridedMaps"))
            rebA10->CreateView("stridedMaps", A10);

          if (!rebA10.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebA10, "A(1,0) rebalanced:", params);
          }
        } // rebalance matrix block A(1,0)

        // extract matrix block
        Teuchos::RCP<CrsMatrix> Am11 = bA->getMatrix(1, 1);
        Teuchos::RCP<CrsMatrixWrap> Aw11 = Teuchos::rcp(new CrsMatrixWrap(Am11));
        Teuchos::RCP<Matrix> A11 = Teuchos::rcp_dynamic_cast<Matrix>(Aw11);

        Teuchos::RCP<Matrix> rebA11;
        {
          std::stringstream ss;
          ss << "Rebalancing matrix block (" << 1 << "," << 1 << ")";
          SubFactoryMonitor subM(*this, ss.str(), coarseLevel);
          RCP<const Map> targetMap = rebalanceImporter->getTargetMap();

          rebA11 = A11; // just a copy
          /*Teuchos::RCP<const CrsMatrixWrap> rebAw10 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebA10);
          Teuchos::RCP<CrsMatrix> rebAm10 = rebAw10->getCrsMatrix();
          Teuchos::RCP<const Import> rebA10Import = ImportFactory::Build(targetMap,A10->getColMap());
          rebAm10->replaceDomainMapAndImporter(targetMap,rebA10Import);*/

          // copy strided map info from non-rebalanced to rebalanced matrix
          if (!rebA11.is_null() && A11->IsView("stridedMaps"))
            rebA11->CreateView("stridedMaps", A11);

          /*if (!rebA10.is_null()) {
            RCP<ParameterList> params = rcp(new ParameterList());
            params->set("printLoadBalancingInfo", true);
            GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*rebA10, "A(1,0) rebalanced:", params);
          }*/
        } // rebalance matrix block A(1,1)



        std::vector<GO> fullRangeMapVector;
        std::vector<GO> fullDomainMapVector;

        std::vector<RCP<const Map> > subBlockARangeMaps;
        std::vector<RCP<const Map> > subBlockADomainMaps;
        subBlockARangeMaps.reserve(bA->Rows());
        subBlockADomainMaps.reserve(bA->Cols());

        // build full range and full domain map (strided)
        if(!rebA00.is_null()){
          subBlockARangeMaps.push_back(rebA00->getRangeMap());
        } else {
          subBlockARangeMaps.push_back(A00->getRangeMap());
        }
        Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap = subBlockARangeMaps.back()->getNodeElementList();
        fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap.begin(), nodeRangeMap.end());
        sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

        if(!rebA10.is_null()){
          subBlockARangeMaps.push_back(rebA10->getRangeMap());
        } else {
          subBlockARangeMaps.push_back(A10->getRangeMap());
        }
        Teuchos::ArrayView< const GlobalOrdinal > nodeRangeMap2 = subBlockARangeMaps.back()->getNodeElementList();
        fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap2.begin(), nodeRangeMap2.end());
        sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

        if(!rebA00.is_null()){
          subBlockADomainMaps.push_back(rebA00->getDomainMap());
        } else {
          subBlockADomainMaps.push_back(A00->getDomainMap());
        }
        Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap = subBlockADomainMaps.back()->getNodeElementList();
        fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap.begin(), nodeDomainMap.end());
        sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

        if(!rebA01.is_null()){
          subBlockADomainMaps.push_back(rebA01->getDomainMap());
        } else {
          subBlockADomainMaps.push_back(A01->getDomainMap());
        }
        Teuchos::ArrayView< const GlobalOrdinal > nodeDomainMap2 = subBlockADomainMaps.back()->getNodeElementList();
        fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap2.begin(), nodeDomainMap2.end());
        sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

        // build full row (=range) map and full domain map

        // extract map index base from maps of blocked A
        GO rangeIndexBase  = 0;
        GO domainIndexBase = 0;
        rangeIndexBase = bA->getRangeMap()->getIndexBase();
        domainIndexBase= bA->getDomainMap()->getIndexBase();


        Teuchos::ArrayView<GO> fullRangeMapGIDs(&fullRangeMapVector[0],fullRangeMapVector.size());
        Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getFullMap());
        std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
        Teuchos::RCP<const StridedMap > fullRangeMap =
            StridedMapFactory::Build(
                bA->getRangeMap()->lib(),
                Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                fullRangeMapGIDs,
                rangeIndexBase,
                stridedData,
                bA->getRangeMap()->getComm(),
                stridedRgFullMap->getStridedBlockId(),
                stridedRgFullMap->getOffset());



        Teuchos::ArrayView<GO> fullDomainMapGIDs(&fullDomainMapVector[0],fullDomainMapVector.size());

        Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getFullMap());
        std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
        Teuchos::RCP<const StridedMap > fullDomainMap =
            StridedMapFactory::Build(
                bA->getDomainMap()->lib(),
                Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                fullDomainMapGIDs,
                domainIndexBase,
                stridedData2,
                bA->getDomainMap()->getComm(),
                stridedDoFullMap->getStridedBlockId(),
                stridedDoFullMap->getOffset());

        // build map extractors
        Teuchos::RCP<const MapExtractorClass> rebRangeMapExtractor  = MapExtractorFactoryClass::Build(fullRangeMap, subBlockARangeMaps);
        Teuchos::RCP<const MapExtractorClass> rebDomainMapExtractor = MapExtractorFactoryClass::Build(fullDomainMap, subBlockADomainMaps);

        Teuchos::RCP<BlockedCrsMatrix> reb_bA = Teuchos::rcp(new BlockedCrsMatrix(rebRangeMapExtractor,rebDomainMapExtractor,10));
        if(!rebA00.is_null()) {
          Teuchos::RCP<const CrsMatrixWrap> rebAw00 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebA00);
          Teuchos::RCP<CrsMatrix> rebAm00 = rebAw00->getCrsMatrix();
          reb_bA->setMatrix(0,0,rebAm00);
        } else {
          Teuchos::RCP<const CrsMatrixWrap> rebAw00 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A00);
          Teuchos::RCP<CrsMatrix> rebAm00 = rebAw00->getCrsMatrix();
          reb_bA->setMatrix(0,0,rebAm00);
        }
        if(!rebA01.is_null()) {
          Teuchos::RCP<const CrsMatrixWrap> rebAw01 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebA01);
          Teuchos::RCP<CrsMatrix> rebAm01 = rebAw01->getCrsMatrix();
          reb_bA->setMatrix(0,1,rebAm01);
        } else {
          Teuchos::RCP<const CrsMatrixWrap> rebAw01 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A01);
          Teuchos::RCP<CrsMatrix> rebAm01 = rebAw01->getCrsMatrix();
          reb_bA->setMatrix(0,1,rebAm01);
        }
        if(!rebA10.is_null()) {
          Teuchos::RCP<const CrsMatrixWrap> rebAw10 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebA10);
          Teuchos::RCP<CrsMatrix> rebAm10 = rebAw10->getCrsMatrix();
          reb_bA->setMatrix(1,0,rebAm10);
        } else {
          Teuchos::RCP<const CrsMatrixWrap> rebAw10 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A10);
          Teuchos::RCP<CrsMatrix> rebAm10 = rebAw10->getCrsMatrix();
          reb_bA->setMatrix(1,0,rebAm10);
        }
        if(!rebA11.is_null()) {
          Teuchos::RCP<const CrsMatrixWrap> rebAw11 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rebA11);
          Teuchos::RCP<CrsMatrix> rebAm11 = rebAw11->getCrsMatrix();
          reb_bA->setMatrix(1,1,rebAm11);
        } else {
          Teuchos::RCP<const CrsMatrixWrap> rebAw11 = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A11);
          Teuchos::RCP<CrsMatrix> rebAm11 = rebAw11->getCrsMatrix();
          reb_bA->setMatrix(1,1,rebAm11);
        }
        reb_bA->fillComplete();
        //reb_bA->describe(*out,Teuchos::VERB_EXTREME);

        coarseLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(reb_bA), this);
      } // end if rebalanceImporter != Teuchos::null

      // no rebalancing necessary
      else {
        // Ac already built by the load balancing process and no load balancing needed
        GetOStream(Warnings0, 0) << "No rebalancing" << std::endl;
        GetOStream(Warnings0, 0) << "Jamming blocked A into Level " << coarseLevel.GetLevelID() << " w/ generating factory "
                                 << this << std::endl;
        Set(coarseLevel, "A", originalAc);
      } // end "no rebalancing"

      // TODO remove me
      //Set(coarseLevel, "A", originalAc);
    } // end "do rebalancing"

#endif

#if 0
    if (rebalanceFacts_.begin() != rebalanceFacts_.end()) {
      SubFactoryMonitor m2(*this, "Rebalance additional data", coarseLevel);

      // call Build of all user-given transfer factories
      for (std::vector<RCP<const FactoryBase> >::const_iterator it = rebalanceFacts_.begin(); it != rebalanceFacts_.end(); ++it) {
        GetOStream(Runtime0, 0) << "RebalanceAc: call rebalance factory " << (*it).get() << ": " << (*it)->description() << std::endl;
        (*it)->CallBuild(coarseLevel);
      }
    }
#endif
  } //Build()

#if 0
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceBlockedAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddRebalanceFactory(const RCP<const FactoryBase>& factory) {

    /*TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                               "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                               "This is very strange. (Note: you can remove this exception if there's a good reason for)");
    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");*/
    rebalanceFacts_.push_back(factory);
  } //AddRebalanceFactory()
#endif
} //namespace MueLu


#endif /* MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_ */
