#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include "MueLu_CoalesceDropFactory_decl.hpp"

namespace MueLu {

  static const std::string color_esc = "\x1b[";
  static const std::string color_std = "39;49;00m";
  static const std::string color_purple = "35m";

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoalesceDropFactory(RCP<const FactoryBase> AFact)
    : AFact_(AFact), blksize_(1), fixedBlkSize_(true)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
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
    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

    // pre-dropping
    RCP<Graph> graph;
    //if (predrop_ == Teuchos::null) {
      //graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));
      Amalgamate(A, graph);
    /*} else {
      //FIXME predropping does not fit to amalgamation routine
      graph = predrop_->Drop(A);
    }*/

    // coalesce

    currentLevel.Set("Graph", graph, this);

    // post-dropping?

  } // Build

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Amalgamate(const RCP<Operator>& A, RCP<Graph>& graph) const {

#if 1

    // do amalgamation
    int nUnamalgamatedBlockSize = blksize_;

    std::cout << "CoalesceDropFactory: nUnamalgamatedBlockSize=" << blksize_ << std::endl;

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
    Teuchos::reduceAll<int,GlobalOrdinal>(*(A->getRowMap()->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, &num_blockids );
    // TODO: check me: is num_blockids = map->getGlobalNumElements()/nUnamalgamatedBlockSize???
    // for constant block size we can avoid the communication and just use above formula!
    // for variable block size, this information has to be provided

    std::cout << "CoalesceDropFactory: PROC:" << A->getRowMap()->getComm()->getRank() << " " << cnt_amalRows << "/" << num_blockids << " blocks" << std::endl;

    // generate row map for amalgamated matrix with same distribution over all procs as row map of A

    Teuchos::ArrayRCP<GlobalOrdinal> arr_amalGIDs = Teuchos::arcp( globalamalblockids );
    Teuchos::RCP<Map> amal_map = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_amalGIDs(), A->getRowMap()->getIndexBase(), A->getRowMap()->getComm());
    std::cout << "CoalesceDropFactory: PROC:" << A->getRowMap()->getComm()->getRank() << " amal_Amap " << amal_map->getNodeNumElements() << "/" << amal_map->getGlobalNumElements() << " elements" << std::endl;

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
    std::cout << "CoalesceDropFactory: PROC:" << A->getRowMap()->getComm()->getRank() << " " << crsGraph->getNodeNumRows() << "/" << crsGraph->getGlobalNumRows() << " rows " << std::endl;
    std::cout << "CoalesceDropFactory: PROC:" << A->getRowMap()->getComm()->getRank() << " " << crsGraph->getNodeNumCols() << "/" << crsGraph->getGlobalNumCols() << " cols " << std::endl;
    graph = rcp(new Graph(/*A->getCrsGraph()*/crsGraph, "amalgamated graph of A")); // TODO change me

    // store information in Graph object for unamalgamation of vectors
    graph->SetAmalgamationParams(globalamalblockid2myrowid, globalamalblockid2globalrowid);

    std::cout << "CoalesceDropFactory: amalgamation complete" << std::endl;
#else
    // old amalgamation routine

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate" << std::endl;

    // do amalgamation
    int nUnamalgamatedBlockSize = blksize_; // TODO fix me, change me to BlkSize???
    //int nPDE = 2; // TODO fix me, find something smart to switch BlkSize from 2 (finest level) to 3 (all other levels...)

    // map: global block id of amalagamated matrix -> vector of local row ids of unamalgamated matrix (only for global block ids of current proc)
    RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid;
    // TODO: add map with GlobalIds!!!
    RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid;
    // vector with global block ids for amalagamated matrix on current proc
    RCP<std::vector<GlobalOrdinal> > globalamalblockids; // TODO remove me

    globalamalblockid2myrowid = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<LocalOrdinal> >);
    globalamalblockid2globalrowid = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);
    /// vector with global block ids for amalagamated matrix on current proc
    globalamalblockids = Teuchos::rcp(new std::vector<GlobalOrdinal>);

    const RCP<const Map> map = A->getRowMap();

    // fill maps for communication of row map GIDs and global block ids of amalgamated matrix
    GlobalOrdinal cnt_amalRows = 0;
    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(map->getNodeNumElements()); i++) {
      // note: this is for constant block size //fixme
      GlobalOrdinal globalblockid = (GlobalOrdinal) map->getGlobalElement(i)/nUnamalgamatedBlockSize;

      if(globalamalblockid2myrowid->count(globalblockid) > 0) {
        globalamalblockid2myrowid->find(globalblockid)->second.push_back(i);
        globalamalblockid2globalrowid->find(globalblockid)->second.push_back(map->getGlobalElement(i));
      } else {
        (*globalamalblockid2myrowid)[globalblockid] = std::vector<LocalOrdinal>(1,i);
        (*globalamalblockid2globalrowid)[globalblockid] = std::vector<LocalOrdinal>(1,map->getGlobalElement(i));
      }

      if(cnt_amalRows>0) {
        if((*globalamalblockids)[cnt_amalRows-1]!=globalblockid) {
          // add "globalblockid" as new row in amalgamated matrix
          globalamalblockids->push_back(globalblockid);
          cnt_amalRows++;
        }
      } else {
        // no amalgamated rows available: just add first one
        globalamalblockids->push_back(globalblockid);
        cnt_amalRows++;
      }
    }

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate: 1" << std::endl;

    // inter processor communication: sum up number of block ids
    GlobalOrdinal num_blockids = 0;
    Teuchos::reduceAll<int,GlobalOrdinal>(*(map->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, &num_blockids );
    // TODO: check me: is num_blockids = map->getGlobalNumElements()/nUnamalgamatedBlockSize???
    // for constant block size we can avoid the communication and just use above formula!
    // for variable block size, this information has to be provided

    // vector for block sizes
    std::vector<int> myblockid2blocksize;
    for(GlobalOrdinal i=0; i<cnt_amalRows; i++) {
      myblockid2blocksize.push_back(nUnamalgamatedBlockSize); // fixme for variable block size
    }
    // now we know, how many blocks the amalagamated matrix has

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate: 2" << std::endl;

    // generate row map for amalgamated matrix with same distribution over all procs as row map of A
    Teuchos::ArrayRCP<GlobalOrdinal> arr_amalGIDs = Teuchos::arcp( globalamalblockids );
    Teuchos::RCP<Map> amal_map = MapFactory::Build(map->lib(), num_blockids, arr_amalGIDs(), map->getIndexBase(), map->getComm());

    // TODO print out number of blocks on current proc

    // extract information from overlapping column map of A
    std::map<GlobalOrdinal,GlobalOrdinal> globalcolid2globalamalblockid;
    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getColMap()->getNodeNumElements());i++) {
      GlobalOrdinal gcolid = A->getColMap()->getGlobalElement(i);
      // fixme for variable block size
      GlobalOrdinal globalblockid = (GlobalOrdinal) gcolid / nUnamalgamatedBlockSize;
      globalcolid2globalamalblockid[gcolid] = globalblockid;
    }

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate: 3" << std::endl;

    // create new CrsGraph for amalgamated matrix (TODO: no shortcut for CrsGraphFactory?)
    RCP<CrsGraph> crsGraph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(amal_map, 10, Xpetra::DynamicProfile);

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate: 3.5" << std::endl;

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " number of local rows: " << A->getRowMap()->getNodeNumElements() << std::endl;

    // fill graph for amalgamated matrix
    for(GlobalOrdinal i=0; i<cnt_amalRows; i++) {
      RCP<std::vector<GlobalOrdinal> > colblocks = Teuchos::rcp(new std::vector<GlobalOrdinal>);  // global column block ids

      // loop over blocksize of block[i]
      // loop over all rows of current block i
      for(int bi=0; bi<myblockid2blocksize[i]; bi++) {
        // extract local row information
        LocalOrdinal row = (LocalOrdinal)i*myblockid2blocksize[i] + bi; // current local row index (TODO: change me for variable block size)

        size_t nnz = A->getNumEntriesInLocalRow(row);
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        A->getLocalRowView(row, indices, vals);
        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: number of nonzeros not equal to number of indices? Error.");

        LocalOrdinal realnnz = 0;
        for(LocalOrdinal k=0; k<Teuchos::as<LocalOrdinal>(nnz); k++) {
          TEUCHOS_TEST_FOR_EXCEPTION(A->getColMap()->isNodeLocalElement(indices[k])==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: Problem with columns. Error.");
          GlobalOrdinal gcid = A->getColMap()->getGlobalElement(indices[k]); // global column id
          if(vals[k]!=0.0) {  // avoid zeros
            colblocks->push_back(globalcolid2globalamalblockid.find(gcid)->second); // add column block id to column ids of amalgamated matrix
            realnnz++; // increment number of nnz in matrix row
          }
        }

      }

      // eliminate duplicate values in colblocks
      sort(colblocks->begin(), colblocks->end());
      typename std::vector<GlobalOrdinal>::iterator endLocation;
      endLocation = std::unique(colblocks->begin(), colblocks->end());
      colblocks->erase(endLocation,colblocks->end());
      Teuchos::ArrayRCP<GlobalOrdinal> arr_colblocks = Teuchos::arcp( colblocks );

      // colblocks now contains the global column ids of the amalgamated matrix
      // now we have all information for one block row
      // fill matrix graph
      GlobalOrdinal grid = amal_map->getGlobalElement(i);
      TEUCHOS_TEST_FOR_EXCEPTION(crsGraph->getRowMap()->isNodeGlobalElement(grid)==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: global row id does not belong to current proc. Error.");
      crsGraph->insertGlobalIndices(grid, arr_colblocks());
    }

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate: 4" << std::endl;

    crsGraph->fillComplete(amal_map,amal_map);

    std::cout << "PROC: " << A->getRowMap()->getComm()->getRank() << " enter CoalesceDropFactory::Amalgamate: 5" << std::endl;

    std::cout << "amalgamated graph: rows=" << crsGraph->getRowMap()->getGlobalNumElements() << std::endl;

    graph = rcp(new Graph(/*A->getCrsGraph()*/crsGraph, "amalgamated graph of A")); // TODO change me

    // store information in Graph object for unamalgamation of vectors
    std::cout << "GLOBALAMALBLOCKID2MYROWID" << std::endl;
    std::cout << globalamalblockid2myrowid << std::endl;
    graph->SetAmalgamationParams(globalamalblockid2myrowid, globalamalblockid2globalrowid, globalamalblockids);

#endif
  }

} //namespace MueLu

#endif // MUELU_COALESCEDROPFACTORY_DEF_HPP
