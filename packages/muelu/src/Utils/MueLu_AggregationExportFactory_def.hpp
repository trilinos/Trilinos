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
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AggregationExportFactory(const std::string outputFileName, const FactoryBase* AggFact, const FactoryBase* CoalesceDropFact, const FactoryBase* AFact)
    : outputFileName_(outputFileName), AggFact_(AggFact), CoalesceDropFact_(CoalesceDropFact), AFact_(AFact)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~AggregationExportFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A",AFact_,this);
    fineLevel.DeclareInput("Aggregates",AggFact_,this);
    fineLevel.DeclareInput("DofsPerNode",CoalesceDropFact_,this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
#include "MueLu_UseShortNames.hpp" // needed because some classes are not forward declared in _decl.hpp
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO

    Teuchos::RCP<OOperator> A =           fineLevel.Get< Teuchos::RCP<OOperator> >("A", AFact_);
    Teuchos::RCP<Aggregates> aggregates = fineLevel.Get< Teuchos::RCP<Aggregates> >("Aggregates",AggFact_);
    LocalOrdinal DofsPerNode =            fineLevel.Get< LocalOrdinal > ("DofsPerNode", CoalesceDropFact_);

    Monitor m(*this, "AggregationExportFactory");

    GetOStream(Runtime0, 0) << "AggregationExportFactory: DofsPerNode: " << DofsPerNode << std::endl;

    ExportAggregates(Teuchos::rcpFromRef(fineLevel), aggregates, A, DofsPerNode);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void  AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ExportAggregates(const Teuchos::RCP<Level>& level, const Teuchos::RCP<Aggregates>& aggregates, const Teuchos::RCP<Operator>& Op, LocalOrdinal DofsPerNode) const {
#include "MueLu_UseShortNames.hpp" // needed because some classes are not forward declared in _decl.hpp

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Op->getRowMap()->getComm();

    Teuchos::RCP<LOVector>vertex2AggId_vector = aggregates->GetVertex2AggId();
    Teuchos::RCP<LOVector>procWinner_vector = aggregates->GetProcWinner();
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> procWinner = aggregates->GetProcWinner()->getDataNonConst(0);

    // 1.) prepare for calculating global aggregate ids
    std::vector<GlobalOrdinal> numAggsGlobal(comm->getSize(),0);
    std::vector<GlobalOrdinal> numAggsLocal(comm->getSize(),0);
    std::vector<GlobalOrdinal> minGlobalAggId(comm->getSize(),0);
    numAggsLocal[comm->getRank()] = aggregates->GetNumAggregates();
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, comm->getSize(),&numAggsLocal[0], &numAggsGlobal[0]);
    for(int i=1; i<Teuchos::as<int>(numAggsGlobal.size()); ++i) {
      numAggsGlobal[i] += numAggsGlobal[i-1];
      minGlobalAggId[i] = numAggsGlobal[i-1];
    }

    Teuchos::ArrayRCP<LocalOrdinal> aggSizes = aggregates->ComputeAggregateSizesDofs();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> > aggToRowMap(aggSizes.size());
    aggregates->ComputeAggregateToRowMap(aggToRowMap);

    // 4.) write to file
    std::string outFile = outputFileName_;
    std::stringstream streamLevel; streamLevel << level->GetLevelID();
    outFile = replaceAll(outFile,"%LEVELID", streamLevel.str());
    std::stringstream streamProc; streamProc << comm->getRank();
    outFile = replaceAll(outFile,"%PROCID", streamProc.str());

    std::ofstream fout(outFile.c_str());

      Teuchos::RCP<const Map> colMap = Op->getColMap();
      std::vector<GlobalOrdinal> nodeIds;
      for (int i=0; i< aggToRowMap.size(); ++i) {
        fout << "Agg " << minGlobalAggId[comm->getRank()] + i /*+ 1*/ << " Proc " << comm->getRank() << ":";
        for (int k=0; k< aggToRowMap[i].size(); ++k) {
          nodeIds.push_back(/*1+*/colMap()->getGlobalElement(aggToRowMap[i][k])/DofsPerNode);
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
