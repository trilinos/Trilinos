#ifndef MUELU_AGGREGATES_DEF_HPP
#define MUELU_AGGREGATES_DEF_HPP

#include "MueLu_Aggregates_decl.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Utilities_decl.hpp" // sumAll

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Aggregates(const Graph & graph) {
    nAggregates_  = 0;
    
    vertex2AggId_ = LOVectorFactory::Build(graph.GetImportMap());
    vertex2AggId_->putScalar(MUELU_UNAGGREGATED);
    
    procWinner_ = LOVectorFactory::Build(graph.GetImportMap());
    procWinner_->putScalar(MUELU_UNASSIGNED);
    
    isRoot_ = Teuchos::ArrayRCP<bool>(graph.GetImportMap()->getNodeNumElements());
    for (size_t i=0; i < graph.GetImportMap()->getNodeNumElements(); i++)
      isRoot_[i] = false;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  Teuchos::ArrayRCP<LocalOrdinal>  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateSizes() const {
    if (aggregateSizes_ == Teuchos::null)
      {
        aggregateSizes_ = Teuchos::ArrayRCP<LO>(nAggregates_);
        int myPid = vertex2AggId_->getMap()->getComm()->getRank();
        Teuchos::ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
        Teuchos::ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);
        LO size = procWinner.size();

        for (LO i = 0; i < nAggregates_; ++i) aggregateSizes_[i] = 0;
        for (LO k = 0; k < size; ++k ) {
          if (procWinner[k] == myPid) aggregateSizes_[vertex2AggId[k]]++;
        }
      }

    return aggregateSizes_;
  } //ComputeAggSizes

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateToRowMap(Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> > &aggToRowMap) const {
    int myPid = vertex2AggId_->getMap()->getComm()->getRank();
    ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
    ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);

    ArrayRCP<LO> aggSizes = ComputeAggregateSizes();
    LO t=0;
    for (typename ArrayRCP<ArrayRCP<LO> >::iterator a2r=aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r)
      *a2r = ArrayRCP<LO>(aggSizes[t++]);
    ArrayRCP< LO > numDofs(nAggregates_,0);  //Track how many DOFS have been recorded so far
    //for each each aggregate in aggToRowMap.
    LO size = procWinner.size();
    for (LO k = 0; k < size; ++k ) {
      LO myAgg = vertex2AggId[k];
      if (procWinner[k] == myPid) {
        aggToRowMap[ myAgg ][ numDofs[myAgg] ] = k;
        ++(numDofs[myAgg]);
      }
    }

  } //AggregateToRowMap

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  std::string Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << BaseClass::description();
    out << "{nGlobalAggregates = " << GetNumGlobalAggregates() << "}";
    return out.str();
  }
     
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    MUELU_DESCRIBE;
      
    if (verbLevel & Statistics0) {
      out0 << "Global number of aggregates: " << GetNumGlobalAggregates() << std::endl;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetNumGlobalAggregates() const {
    LO nAggregates = GetNumAggregates();
    GO nGlobalAggregates; sumAll(vertex2AggId_->getMap()->getComm(), (GO)nAggregates, nGlobalAggregates);
    return nGlobalAggregates;
  }
    
} //namespace MueLu

#endif // MUELU_AGGREGATES_DEF_HPP
