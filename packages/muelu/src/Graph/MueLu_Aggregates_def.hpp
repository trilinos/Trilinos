#ifndef MUELU_AGGREGATES_DEF_HPP
#define MUELU_AGGREGATES_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_Aggregates_decl.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationInfo.hpp"
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

    // create an empty container for amalgamation information
    // transfer amalgamation data from graph to AmalgamationInfo container
    // TODO: move this?
    amalgamationData_ = rcp(new AmalgamationInfo());
    amalgamationData_->SetAmalgamationParams(graph.GetMyAmalgamationParams(),graph.GetGlobalAmalgamationParams());

    GenerateImportDofMap(); // amalgamation parameters have to be set before!

  }

  // special copy constructor, does not handle amalgamation information!
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Aggregates(Aggregates & a) {
    nAggregates_  = a.GetNumAggregates();

    vertex2AggId_ = a.GetVertex2AggId();
    procWinner_   = a.GetProcWinner();

    isRoot_ = Teuchos::ArrayRCP<bool>(vertex2AggId_->getLocalLength());
    for (size_t i=0; i < vertex2AggId_->getLocalLength(); i++)
          isRoot_[i] = a.IsRoot(i);

    amalgamationData_ = rcp(new AmalgamationInfo());

    importDofMap_ = Teuchos::null;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateToRowMap(Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> > &aggToRowMap) const {
    // decide whether we need the DOF version (for problems with amalgamated matrix)
    // or just the Node version (for problems with 1 DOF per node)
    if(GetAmalgamationInfo()->GetMyAmalgamationParams() == Teuchos::null) {
      ComputeAggregateToRowMapNodes(aggToRowMap);
    }
    else {
      ComputeAggregateToRowMapDofs(aggToRowMap);
    }
  }

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateToRowMapNodes(Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> > &aggToRowMap) const {
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

  } //AggregateToRowMapNodes

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateToRowMapDofs(Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> > &aggToRowMap) const {

    int myPid = vertex2AggId_->getMap()->getComm()->getRank();
    ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
    ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0); // vector: node 2 aggid

    ArrayRCP<LO> aggSizes = ComputeAggregateSizesDofs(); // size of aggregats (# of nodes in agg)

    // length of aggToRowMap should be the number of DOFs (not the number of nodes)
    LO t=0;
    for (typename ArrayRCP<ArrayRCP<LO> >::iterator a2r=aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r) {
      *a2r = ArrayRCP<LO>(aggSizes[t++]);
    }

    // track how many dofs have been recorded so far
    ArrayRCP< LO > numDofs(nAggregates_,0);

    //for each each aggregate in aggToRowMap.
    LO size = procWinner.size();

    // loop over local node ids
    for (LO lnode = 0; lnode < size; ++lnode ) {
      LO myAgg = vertex2AggId[lnode]; // local agg id for current local node id

      if (procWinner[lnode] == myPid) {
        // for loop over all local row ids for current block id = global node id?
        GlobalOrdinal gblockid = vertex2AggId_->getMap()->getGlobalElement(lnode);

        // unumalgamate graph-based information to dof-based information
        std::vector<LocalOrdinal> blockdofs = (*(GetAmalgamationInfo()->GetMyAmalgamationParams()))[gblockid];
        //std::cout << blockdofs.size() << ": ";
        for (LocalOrdinal blockdof=0; blockdof<Teuchos::as<LocalOrdinal>(blockdofs.size()); blockdof++) {
          //std::cout << blockdofs[blockdof] << " " ;
          aggToRowMap[ myAgg ][ numDofs[myAgg] ] = blockdofs[blockdof];  // add DOF to current aggregate
          ++(numDofs[myAgg]);
        }
        //std::cout << std::endl;
      }
    }

  }

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<LocalOrdinal>  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateSizes() const {
    if (aggregateSizes_ == Teuchos::null)
    {
      // decide whether we need the DOF version (for problems with amalgamated matrix)
      // or just the Node version (for problems with 1 DOF per node)
      //if(globalamalblockid2myrowid_ == Teuchos::null) {
      if (GetAmalgamationInfo()->GetMyAmalgamationParams() == Teuchos::null) {
        ComputeAggregateSizesNodes();
      }
      else {
        ComputeAggregateSizesDofs();
      }
    }

    return aggregateSizes_;
  } //ComputeAggSizesNodes

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<LocalOrdinal>  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateSizesNodes() const {
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
  } //ComputeAggSizesNodes

  ///////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<LocalOrdinal>  Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeAggregateSizesDofs() const {
    if (aggregateSizes_ == Teuchos::null)
    {
        aggregateSizes_ = Teuchos::ArrayRCP<LO>(nAggregates_);
        int myPid = vertex2AggId_->getMap()->getComm()->getRank();
        Teuchos::ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
        Teuchos::ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);
        LO size = procWinner.size();

        for (LO i = 0; i < nAggregates_; ++i) aggregateSizes_[i] = 0;
        for (LO lnode = 0; lnode < size; ++lnode ) {
          LO myAgg = vertex2AggId[lnode];
          if (procWinner[lnode] == myPid) {
            GlobalOrdinal gblockid = vertex2AggId_->getMap()->getGlobalElement(lnode);

            // unumalgamate graph-based information to dof-based information
            std::vector<LocalOrdinal> blockdofs = (*(GetAmalgamationInfo()->GetMyAmalgamationParams()))[gblockid];
            aggregateSizes_[myAgg] += Teuchos::as<LocalOrdinal>(blockdofs.size());
          }
        }
    }
    return aggregateSizes_;
  } //ComputeAggSizesDofs

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

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GenerateImportDofMap() const
  {
    RCP<const Map> nodeMap = vertex2AggId_->getMap(); // use import node map from graph

    // special case: 1 dof per node
    if(GetAmalgamationInfo()->GetMyAmalgamationParams() == Teuchos::null &&
        GetAmalgamationInfo()->GetGlobalAmalgamationParams() == Teuchos::null) {
      GetOStream(Debug, 0) << "MueLu::Aggregates::GenerateImportDofMap: 1 dof per node -> skip reconstruction of import DOF map!" << std::endl;
      // TODO: add debug statement

      // no amalgamation information -> we can assume that we have 1 dof per node
      // just use nodeMap as DOFMap!
      importDofMap_ = nodeMap;
    
      return;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(GetAmalgamationInfo()->GetGlobalAmalgamationParams()==Teuchos::null, Exceptions::RuntimeError, "MueLu::Aggregates::GenerateImportDofMap: insufficient amalgamation information. Error");
    TEUCHOS_TEST_FOR_EXCEPTION(GetAmalgamationInfo()->GetMyAmalgamationParams()==Teuchos::null    , Exceptions::RuntimeError, "MueLu::Aggregates::GenerateImportDofMap: insufficient amalgamation information. Error");

    // build dof map from node map
    RCP<std::vector<GlobalOrdinal> > myDofGIDs = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    for(LocalOrdinal n=0; n<Teuchos::as<LocalOrdinal>(nodeMap->getNodeNumElements()); n++) {
      GlobalOrdinal globalblockid = (GlobalOrdinal) nodeMap->getGlobalElement(n);

      TEUCHOS_TEST_FOR_EXCEPTION(GetAmalgamationInfo()->GetGlobalAmalgamationParams()->count(globalblockid)<=0, Exceptions::RuntimeError, "MueLu::Aggregates::GenerateImportDofMap: empty global block? Error.");
      std::vector<GlobalOrdinal> myrowGIDs = (*(GetAmalgamationInfo()->GetGlobalAmalgamationParams()))[globalblockid];
      TEUCHOS_TEST_FOR_EXCEPTION(myrowGIDs.size()==0, Exceptions::RuntimeError, "MueLu::Aggregates::GenerateImportDofMap: no amalgamation information! Error.");

      typename std::vector<GlobalOrdinal>::iterator gidIt;
      for(gidIt = myrowGIDs.begin(); gidIt!=myrowGIDs.end(); gidIt++) {
        myDofGIDs->push_back(*gidIt); // append local row ids
      }
    }

    // generate row dof map for amalgamated matrix with same distribution over all procs as row node map
    Teuchos::ArrayRCP<GlobalOrdinal> arr_myDofGIDs = Teuchos::arcp( myDofGIDs );

    importDofMap_ = MapFactory::Build(nodeMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), arr_myDofGIDs(), nodeMap->getIndexBase(), nodeMap->getComm());
  }

} //namespace MueLu

#endif // MUELU_AGGREGATES_DEF_HPP
