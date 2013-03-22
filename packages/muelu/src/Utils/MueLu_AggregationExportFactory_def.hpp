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
 * MueLu_AggregationExportFactory_def.hpp
 *
 *  Created on: Feb 10, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_
#define MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_AggregationExportFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("Aggregates",          Teuchos::null,                        "Generating factory for aggregates");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode",         Teuchos::null,                        "Generating factory for number of dofs per node");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo",  Teuchos::null,                        "Generating factory for amalgamation");

    validParamList->set< std::string >           ("Output filename",     "aggs_level%LEVELID_proc%PROCID.out", "Output filename template (%LEVELID is replaced level id, %PROCID is replaced by processor id)");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "Aggregates");         //< factory which created aggregates
    Input(fineLevel, "DofsPerNode");        //< CoalesceAndDropFactory (needed for DofsPerNode variable)
    Input(fineLevel, "UnAmalgamationInfo"); //< AmalgamationFactory (needed for UnAmalgamationInfo variable)
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "AggregationExportFactory", coarseLevel);

    Teuchos::RCP<Aggregates> aggregates      = Get< Teuchos::RCP<Aggregates> >(fineLevel,"Aggregates");
    LocalOrdinal DofsPerNode                 = Get< LocalOrdinal >            (fineLevel,"DofsPerNode");
    Teuchos::RCP<AmalgamationInfo> amalgInfo = Get< RCP<AmalgamationInfo> >   (fineLevel,"UnAmalgamationInfo");

    GetOStream(Runtime0, 0) << "AggregationExportFactory: DofsPerNode: " << DofsPerNode << std::endl;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = aggregates->GetMap()->getComm();

    Teuchos::RCP<LocalOrdinalVector>vertex2AggId_vector = aggregates->GetVertex2AggId();
    Teuchos::RCP<LocalOrdinalVector>procWinner_vector = aggregates->GetProcWinner();
    Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
    Teuchos::ArrayRCP<LocalOrdinal> procWinner = aggregates->GetProcWinner()->getDataNonConst(0);

    // prepare for calculating global aggregate ids
    std::vector<GlobalOrdinal> numAggsGlobal(comm->getSize(),0);
    std::vector<GlobalOrdinal> numAggsLocal(comm->getSize(),0);
    std::vector<GlobalOrdinal> minGlobalAggId(comm->getSize(),0);
    numAggsLocal[comm->getRank()] = aggregates->GetNumAggregates();
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, comm->getSize(),&numAggsLocal[0], &numAggsGlobal[0]);
    for(int i=1; i<Teuchos::as<int>(numAggsGlobal.size()); ++i) {
      numAggsGlobal[i] += numAggsGlobal[i-1];
      minGlobalAggId[i] = numAggsGlobal[i-1];
    }

    ArrayRCP<LO> aggStart;
    ArrayRCP<GlobalOrdinal> aggToRowMap;
    AmalgamationFactory::UnamalgamateAggregates(*aggregates, *amalgInfo, aggStart, aggToRowMap);

    // write to file
    //std::string outFile = outputFileName_;
    const ParameterList & pL = GetParameterList();
    std::string outFile = pL.get<std::string> ("Output filename");
    std::stringstream streamLevel; streamLevel << fineLevel.GetLevelID();
    outFile = replaceAll(outFile,"%LEVELID", streamLevel.str());
    std::stringstream streamProc; streamProc << comm->getRank();
    outFile = replaceAll(outFile,"%PROCID", streamProc.str());

    GetOStream(Runtime0, 0) << "AggregationExportFactory: outputfilel \"" << outFile << "\"" << std::endl;
    std::ofstream fout(outFile.c_str());

    GO numAggs = aggregates->GetNumAggregates();
    std::vector<GlobalOrdinal> nodeIds;
    for (int i=0; i< numAggs; ++i) {
      fout << "Agg " << minGlobalAggId[comm->getRank()] + i << " Proc " << comm->getRank() << ":";
      for (int k=aggStart[i]; k< aggStart[i+1]; ++k) {
        /*std::cout << "proc: " << comm->getRank() << "\t aggToRowMap[" << i << "][" << k << "]=" <<aggToRowMap[i][k] << "\t node GID: " << aggToRowMap[i][k]/DofsPerNode << "\t GID in colMap=" << aggToRowMap[i][k];
          if(colMap->isNodeGlobalElement(aggToRowMap[i][k])==false)
          std::cout << " NOT ON CUR PROC!";
          std::cout << std::endl;*/

        nodeIds.push_back(aggToRowMap[k]/DofsPerNode);
      }

      // remove duplicate entries from nodeids
      std::sort(nodeIds.begin(),nodeIds.end());
      typename std::vector<GlobalOrdinal>::iterator endLocation = std::unique(nodeIds.begin(),nodeIds.end());
      nodeIds.erase(endLocation,nodeIds.end());

      // print out nodeids
      for(typename std::vector<GlobalOrdinal>::iterator printIt = nodeIds.begin(); printIt!=nodeIds.end(); printIt++) {
        fout << " " << *printIt;
      }
      nodeIds.clear();
      fout << std::endl;
    }

    fout.close();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const {
    while(1)
    {
      const int pos = result.find(replaceWhat);
      if (pos==-1) break;
      result.replace(pos,replaceWhat.size(),replaceWithWhat);
    }
    return result;
  }

} // namespace MueLu

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_ */
