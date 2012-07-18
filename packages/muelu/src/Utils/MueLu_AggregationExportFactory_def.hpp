/*
 * MueLu_AggregationExportFactory_def.hpp
 *
 *  Created on: Feb 10, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_
#define MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>

#include "MueLu_AggregationExportFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AggregationExportFactory(const std::string outputFileName, const FactoryBase* AggFact, const FactoryBase* CoalesceDropFact)
    : outputFileName_(outputFileName), AggFact_(AggFact), CoalesceDropFact_(CoalesceDropFact)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~AggregationExportFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("Aggregates",AggFact_,this);
    fineLevel.DeclareInput("DofsPerNode",CoalesceDropFact_,this);
    fineLevel.DeclareInput("UnAmalgamationInfo", CoalesceDropFact_, this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "AggregationExportFactory", coarseLevel);

    Teuchos::RCP<Aggregates> aggregates = fineLevel.Get< Teuchos::RCP<Aggregates> >("Aggregates",AggFact_);
    LocalOrdinal DofsPerNode =                 fineLevel.Get< LocalOrdinal > ("DofsPerNode", CoalesceDropFact_);
    Teuchos::RCP<AmalgamationInfo> amalgInfo = fineLevel.Get< RCP<AmalgamationInfo> >("UnAmalgamationInfo", CoalesceDropFact_);
    
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

    GO numAggs = aggregates->GetNumAggregates();
    ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
    ComputeAggregateSizes(*aggregates, *amalgInfo, aggSizes);
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > aggToRowMap(aggSizes.size());
    ComputeAggregateToRowMap(*aggregates, *amalgInfo, aggSizes, aggToRowMap);
        
    // write to file
    std::string outFile = outputFileName_;
    std::stringstream streamLevel; streamLevel << fineLevel.GetLevelID();
    outFile = replaceAll(outFile,"%LEVELID", streamLevel.str());
    std::stringstream streamProc; streamProc << comm->getRank();
    outFile = replaceAll(outFile,"%PROCID", streamProc.str());

    std::ofstream fout(outFile.c_str());


      std::vector<GlobalOrdinal> nodeIds;
      for (int i=0; i< aggToRowMap.size(); ++i) {
        fout << "Agg " << minGlobalAggId[comm->getRank()] + i << " Proc " << comm->getRank() << ":";
        for (int k=0; k< aggToRowMap[i].size(); ++k) {
          /*std::cout << "proc: " << comm->getRank() << "\t aggToRowMap[" << i << "][" << k << "]=" <<aggToRowMap[i][k] << "\t node GID: " << aggToRowMap[i][k]/DofsPerNode << "\t GID in colMap=" << aggToRowMap[i][k];
          if(colMap->isNodeGlobalElement(aggToRowMap[i][k])==false)
            std::cout << " NOT ON CUR PROC!";
          std::cout << std::endl;*/

          nodeIds.push_back(aggToRowMap[i][k]/DofsPerNode);
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
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateSizes(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, Teuchos::ArrayRCP<LocalOrdinal> & aggSizes) const {
    // TODO: change routines in TentativePFactory to static and use the TentativePFactory routines
    // we expect the aggSizes array to be initialized as follows
    // aggSizes = Teuchos::ArrayRCP<LO>(nAggregates_);
    // furthermore we suppose the (un)amalgamation info to be set (even for 1 dof per node examples)

    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();

    for (LO i = 0; i< aggregates.GetNumAggregates(); ++i) aggSizes[i] = 0;
    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];
      if (procWinner[lnode] == myPid) {
        GO gnodeid = aggregates.GetMap()->getGlobalElement(lnode);

        std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
        aggSizes[myAgg] += Teuchos::as<LO>(gDofIds.size());
      }
    }    
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateToRowMap(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes, Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > & aggToRowMap) const {
    // TODO: change routines in TentativePFactory to static and use the TentativePFactory routines
    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();

    // initialize array aggToRowMap with empty arrays for each aggregate (with correct aggSize)
    LO t = 0;
    for (typename ArrayRCP<ArrayRCP<GO> >::iterator a2r = aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r) {
      *a2r = ArrayRCP<GO>(aggSizes[t++]);
    }

    // count, how many dofs have been recorded for each aggregate
    ArrayRCP<LO> numDofs(aggregates.GetNumAggregates(),0); // empty array with number of Dofs for each aggregate

    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];
      if (procWinner[lnode] == myPid) {
        GO gnodeid = aggregates.GetMap()->getGlobalElement(lnode);
        std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
        for (LO gDofId=0; gDofId < Teuchos::as<LO>(gDofIds.size()); gDofId++) {
          aggToRowMap[ myAgg ][ numDofs[myAgg] ] = gDofIds[gDofId]; // fill aggToRowMap structure
          ++(numDofs[myAgg]);
        }
      }
    }    
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
