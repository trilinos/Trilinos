#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include "MueLu_CoalesceDropFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  static const std::string color_esc = "\x1b[";
  static const std::string color_std = "39;49;00m";
  static const std::string color_purple = "35m";

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoalesceDropFactory(RCP<const FactoryBase> AFact, RCP<const FactoryBase> nullspaceFact)
    : AFact_(AFact), nullspaceFact_(nullspaceFact), blksize_(1), fixedBlkSize_(true)
  {
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get(), this);
    currentLevel.DeclareInput("Nullspace", nullspaceFact_.get(), this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetFixedBlockSize(LocalOrdinal blksize) {
    blksize_ = blksize;
    fixedBlkSize_ = true;
    GetOStream(Debug, 0) << color_esc << color_purple << "CoalesceDropFactory::SetFixedBlockSize()" << color_esc << color_std << std::endl;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &predrop) { predrop_ = predrop; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    Monitor m(*this, "CoalesceDropFactory");

    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());
    RCP<MultiVector> nullspace  = currentLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());

    //RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    // pre-dropping
    RCP<Graph> graph;
    //if (predrop_ == Teuchos::null) {
      //graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));
      LocalOrdinal blockdim = 1;
      if(currentLevel.GetLevelID() == 0) {
        blockdim = blksize_;
      } else {
        blockdim = Teuchos::as<LocalOrdinal>(nullspace->getNumVectors());
      }

      if (blockdim > 1) {
        Amalgamate(A, blockdim, graph);
      } else {
        graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));
      }
    /*} else {
      //FIXME predropping does not fit to amalgamation routine
      graph = predrop_->Drop(A);
    }*/

    // coalesce

    currentLevel.Set("Graph", graph, this);

    // post-dropping?

  } // Build

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Amalgamate(const RCP<Operator>& A, const LocalOrdinal blocksize, RCP<Graph>& graph) const {

    GetOStream(Runtime0, 0) << color_esc << color_purple << "CoalesceDropFactory::Amalgamate()" << color_esc << color_std << " constant blocksize=" << blocksize << std::endl;

    // do amalgamation
    int nUnamalgamatedBlockSize = blocksize;

    // map: global block id of amalagamated matrix -> vector of local row ids of unamalgamated matrix (only for global block ids of current proc)
    RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid;
    RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid;
    RCP<std::vector<GlobalOrdinal> > globalamalblockids;
    globalamalblockid2myrowid = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<LocalOrdinal> >);
    globalamalblockid2globalrowid = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);
    globalamalblockids = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    globalamalblockids->empty();

    // extract information from overlapping column map of A
    GlobalOrdinal cnt_amalRows = 0;
    RCP<std::map<GlobalOrdinal,GlobalOrdinal> > globalrowid2globalamalblockid  = Teuchos::rcp(new std::map<GlobalOrdinal, GlobalOrdinal>);
    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getColMap()->getNodeNumElements());i++) {
      GlobalOrdinal gDofId = A->getColMap()->getGlobalElement(i);
      // fixme for variable block size

      GlobalOrdinal globalblockid = (GlobalOrdinal) gDofId / nUnamalgamatedBlockSize;

      // gDofId -> gblockid
      (*globalrowid2globalamalblockid)[gDofId] = globalblockid;

      // gblockid -> gDofId/lDofId
      if(globalamalblockid2myrowid->count(globalblockid) > 0) {
        globalamalblockid2myrowid->find(globalblockid)->second.push_back(i);
        globalamalblockid2globalrowid->find(globalblockid)->second.push_back(gDofId);
      } else {
        (*globalamalblockid2myrowid)[globalblockid] = std::vector<LocalOrdinal>(1,i);
        (*globalamalblockid2globalrowid)[globalblockid] = std::vector<GlobalOrdinal>(1,gDofId);
        if(A->getRowMap()->isNodeGlobalElement(gDofId)) {
          globalamalblockids->push_back(globalblockid);
          cnt_amalRows++; // new local block row in amalgamated matrix graph
        }
      }
    }

    // clean up DofVectors (remove duplicate entries)
    typename std::map<GlobalOrdinal,std::vector<LocalOrdinal> >::iterator lit;
    typename std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >::iterator git;
    for (lit=globalamalblockid2myrowid->begin(); lit!=globalamalblockid2myrowid->end(); lit++) {
      std::vector<LocalOrdinal> lrowids = lit->second;
      sort(lrowids.begin(), lrowids.end());
      typename std::vector<LocalOrdinal>::iterator lendLocation;
      lendLocation = std::unique(lrowids.begin(), lrowids.end());
      lrowids.erase(lendLocation,lrowids.end());
    }
    for (git=globalamalblockid2globalrowid->begin(); git!=globalamalblockid2globalrowid->end(); git++) {
      std::vector<GlobalOrdinal> growids = git->second;
      sort(growids.begin(), growids.end());
      typename std::vector<GlobalOrdinal>::iterator gendLocation;
      gendLocation = std::unique(growids.begin(), growids.end());
      growids.erase(gendLocation,growids.end());
    }


    // inter processor communication: sum up number of block ids
    GlobalOrdinal num_blockids = 0;
    Teuchos::reduceAll<int,GlobalOrdinal>(*(A->getRowMap()->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, Teuchos::ptr(&num_blockids) );
    // TODO: check me: is num_blockids = map->getGlobalNumElements()/nUnamalgamatedBlockSize???
    // for constant block size we can avoid the communication and just use above formula!
    // for variable block size, this information has to be provided

    GetOStream(Debug, 0) << color_esc << color_purple << "CoalesceDropFactory::Amalgamate()" << color_esc << color_std << " # of amalgamated blocks=" << num_blockids << std::endl;

    // generate row map for amalgamated matrix with same distribution over all procs as row map of A

    Teuchos::ArrayRCP<GlobalOrdinal> arr_amalGIDs = Teuchos::arcp( globalamalblockids );
    Teuchos::RCP<Map> amal_map = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_amalGIDs(), A->getRowMap()->getIndexBase(), A->getRowMap()->getComm());
    GetOStream(Debug, 0) << "CoalesceDropFactory: amal_Amap " << amal_map->getNodeNumElements() << "/" << amal_map->getGlobalNumElements() << " elements" << std::endl;

    // create new CrsGraph for amalgamated matrix (TODO: no shortcut for CrsGraphFactory?)
    RCP<CrsGraph> crsGraph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(amal_map, 10, Xpetra::DynamicProfile);

    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements());i++) {
      GlobalOrdinal gDofId = A->getRowMap()->getGlobalElement(i);
      // fixme for variable block size
      GlobalOrdinal globalblockid = (GlobalOrdinal) gDofId / nUnamalgamatedBlockSize;

      size_t nnz = A->getNumEntriesInLocalRow(i);
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      A->getLocalRowView(i, indices, vals);
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: number of nonzeros not equal to number of indices? Error.");

      RCP<std::vector<GlobalOrdinal> > colblocks = Teuchos::rcp(new std::vector<GlobalOrdinal>);  // global column block ids
      LocalOrdinal realnnz = 0;
      for(LocalOrdinal k=0; k<Teuchos::as<LocalOrdinal>(nnz); k++) {
        TEUCHOS_TEST_FOR_EXCEPTION(A->getColMap()->isNodeLocalElement(indices[k])==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: Problem with columns. Error.");
        GlobalOrdinal gcid = A->getColMap()->getGlobalElement(indices[k]); // global column id
        // TODO: decide whether to add or skip a matrix entry in resulting graph
        if(vals[k]!=0.0) {  // avoid zeros
          colblocks->push_back(globalrowid2globalamalblockid->find(gcid)->second); // add column block id to column ids of amalgamated matrix
          realnnz++; // increment number of nnz in matrix row
        }
      }

      Teuchos::ArrayRCP<GlobalOrdinal> arr_colblocks = Teuchos::arcp( colblocks );

      // fill matrix graph
      TEUCHOS_TEST_FOR_EXCEPTION(crsGraph->getRowMap()->isNodeGlobalElement(globalblockid)==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: global row id does not belong to current proc. Error.");
      crsGraph->insertGlobalIndices(globalblockid, arr_colblocks());
    }

    crsGraph->fillComplete(amal_map,amal_map);

    // create MueLu::Graph object
    graph = rcp(new Graph(crsGraph, "amalgamated graph of A"));

    // store information in Graph object for unamalgamation of vectors
    graph->SetAmalgamationParams(globalamalblockid2myrowid, globalamalblockid2globalrowid);
  }

} //namespace MueLu

#endif // MUELU_COALESCEDROPFACTORY_DEF_HPP
