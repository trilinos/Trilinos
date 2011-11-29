#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include "MueLu_CoalesceDropFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_Monitor.hpp"

#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_ImportFactory.hpp"

namespace MueLu {

  static const std::string color_esc = "\x1b[";
  static const std::string color_std = "39;49;00m";
  static const std::string color_purple = "35m";

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoalesceDropFactory(RCP<const FactoryBase> AFact, RCP<const FactoryBase> nullspaceFact)
    : AFact_(AFact), nullspaceFact_(nullspaceFact), blksize_(1), fixedBlkSize_(true), blkSizeInfo_(Teuchos::null)
  {
    predrop_ = Teuchos::null;  // no pre-dropping filter
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get(), this);
    currentLevel.DeclareInput("Nullspace", nullspaceFact_.get(), this);
    if(fixedBlkSize_ == false && currentLevel.GetLevelID() == 0)
      currentLevel.DeclareInput("VariableBlockSizeInfo", MueLu::NoFactory::get(), this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetFixedBlockSize(LocalOrdinal blksize) {
    blksize_ = blksize;
    fixedBlkSize_ = true;
    GetOStream(Debug, 0) << color_esc << color_purple << "CoalesceDropFactory::SetFixedBlockSize()" << color_esc << color_std << std::endl;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetVariableBlockSize() {
    fixedBlkSize_ = false;
    GetOStream(Debug, 0) << color_esc << color_purple << "CoalesceDropFactory::SetVariableBlockSize()" << color_esc << color_std << std::endl;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &predrop) { predrop_ = predrop; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    Monitor m(*this, "CoalesceDropFactory");

    //RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());
    RCP<MultiVector> nullspace  = currentLevel.Get< RCP<MultiVector> >("Nullspace", nullspaceFact_.get());

    LocalOrdinal blockdim = 1; // block dim for fixed size blocks

    if(currentLevel.GetLevelID() == 0) {
      // switch between constant and variable block size on finest level
      if(fixedBlkSize_ == false) {
        // variable block size
        blockdim = -1; // no constant block size
        // read in and transform variable block size information
        RCP<Vector> blkSizeInfo = currentLevel.Get<RCP<Vector> >("VariableBlockSizeInfo", MueLu::NoFactory::get());
        TEUCHOS_TEST_FOR_EXCEPTION(blkSizeInfo->getMap()->isSameAs(*(A->getRowMap()))==false, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Build: map of blkSizeInfo does not match the row map of A. Error.");
        RCP<const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(A->getRowMap(),A->getColMap());
        blkSizeInfo_ = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A->getColMap());
        blkSizeInfo_->doImport(*blkSizeInfo,*importer,Xpetra::INSERT);
        TEUCHOS_TEST_FOR_EXCEPTION(blkSizeInfo_->getMap()->isSameAs(*(A->getColMap()))==false, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Build: map of blkSizeInfo does not match the column map of A. Error.");
      } else {
        // constant block size
        blockdim = blksize_;
        blkSizeInfo_ = Teuchos::null;
      }
    } else {
      // switch to constant blocksize on intermediate and coarse levels
      blockdim = Teuchos::as<LocalOrdinal>(nullspace->getNumVectors());
      fixedBlkSize_ = true; // TODO: what about this??
    }

    // pre-dropping
    RCP<Graph> graph;

    if (blockdim > 1 || fixedBlkSize_ == false) {
      Amalgamate(A, blockdim, graph);
    } else {
      graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));
    }


    currentLevel.Set("Graph", graph, this);

    // post-dropping?

  } // Build

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GlobalId2GlobalAmalBlockId(GlobalOrdinal gid, const RCP<Operator>& A, const RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& globalgid2globalamalblockid_vector, LocalOrdinal blockSize) const {

    if(fixedBlkSize_ == true) {
      //GetOStream(Runtime0, 0) << "fixed block size..." << std::endl;
      GlobalOrdinal globalblockid = (GlobalOrdinal) gid / blockSize;
      return globalblockid;
    } else {
      //GetOStream(Runtime0, 0) << "variable block size..." << std::endl;
      Teuchos::ArrayRCP< Scalar > ovamalblockid_data = globalgid2globalamalblockid_vector->getDataNonConst(0);

      Teuchos::RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > overlappingMap = globalgid2globalamalblockid_vector->getMap();

      // check if map of globalgid2globalamalblockid_vector is same as overlapping column map of A
      TEUCHOS_TEST_FOR_EXCEPTION(overlappingMap->isSameAs(*(A->getColMap()))==false, Exceptions::RuntimeError, "MueLu::CoalesceFactory::GlobalId2GlobalAmalBlockId: map of globalgid2globalamalblockid_vector in GlobalId2GlobalAmalBlockId must be same as column map of A. Error.");

      LocalOrdinal lid = overlappingMap->getLocalElement(gid);
      GlobalOrdinal amalgid = ovamalblockid_data[lid];
      return amalgid;
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupAmalgamationData(const RCP<Operator>& A, const RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& globalgid2globalamalblockid_vector, LocalOrdinal blockSize) const {

    GetOStream(Runtime0, 0) << color_esc << color_purple << "CoalesceDropFactory::SetupAmalgamationData()" << color_esc << color_std << std::endl;

    globalamalblockid2myrowid_ = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<LocalOrdinal> >);
    globalamalblockid2globalrowid_ = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);

    RCP<std::vector<GlobalOrdinal> > globalamalblockids;
    globalamalblockids = Teuchos::rcp(new std::vector<GlobalOrdinal>); // vector of global amal block ids on current processor
    globalamalblockids->empty();

    // extract information from overlapping column map of A
    GlobalOrdinal cnt_amalRows = 0;
    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getColMap()->getNodeNumElements());i++) {
      // get global DOF id
      GlobalOrdinal gDofId = A->getColMap()->getGlobalElement(i);

      // translate gDofId to global amal block id
      GlobalOrdinal globalblockid = GlobalId2GlobalAmalBlockId(gDofId, A, globalgid2globalamalblockid_vector, blockSize);

      // gblockid -> gDofId/lDofId
      if(globalamalblockid2myrowid_->count(globalblockid) > 0) {
        globalamalblockid2myrowid_->find(globalblockid)->second.push_back(i);
        globalamalblockid2globalrowid_->find(globalblockid)->second.push_back(gDofId);
      } else {
        (*globalamalblockid2myrowid_)[globalblockid] = std::vector<LocalOrdinal>(1,i);
        (*globalamalblockid2globalrowid_)[globalblockid] = std::vector<GlobalOrdinal>(1,gDofId);
        if(A->getRowMap()->isNodeGlobalElement(gDofId)) {
          globalamalblockids->push_back(globalblockid);
          cnt_amalRows++; // new local block row in amalgamated matrix graph
        }
      }
    }

    // clean up DofVectors (remove duplicate entries)
    typename std::map<GlobalOrdinal,std::vector<LocalOrdinal> >::iterator lit;
    typename std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >::iterator git;
    for (lit=globalamalblockid2myrowid_->begin(); lit!=globalamalblockid2myrowid_->end(); lit++) {
      std::vector<LocalOrdinal> lrowids = lit->second;
      sort(lrowids.begin(), lrowids.end());
      typename std::vector<LocalOrdinal>::iterator lendLocation;
      lendLocation = std::unique(lrowids.begin(), lrowids.end());
      lrowids.erase(lendLocation,lrowids.end());
    }
    for (git=globalamalblockid2globalrowid_->begin(); git!=globalamalblockid2globalrowid_->end(); git++) {
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

    GetOStream(Debug, 0) << color_esc << color_purple << "CoalesceDropFactory::SetupAmalgamationData()" << color_esc << color_std << " # of amalgamated blocks=" << num_blockids << std::endl;

    // generate row map for amalgamated matrix with same distribution over all procs as row map of A

    Teuchos::ArrayRCP<GlobalOrdinal> arr_amalGIDs = Teuchos::arcp( globalamalblockids );
    Teuchos::RCP<Map> amal_map = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_amalGIDs(), A->getRowMap()->getIndexBase(), A->getRowMap()->getComm());
    GetOStream(Debug, 0) << "CoalesceDropFactory: amal_map " << amal_map->getNodeNumElements() << "/" << amal_map->getGlobalNumElements() << " elements" << std::endl;

    return amal_map;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Amalgamate(const RCP<Operator>& A, const LocalOrdinal blockSize, RCP<Graph>& graph) const {

    GetOStream(Runtime0, 0) << color_esc << color_purple << "CoalesceDropFactory::Amalgamate()" << color_esc << color_std << " constant blocksize=" << blockSize << std::endl;

    // do amalgamation


    // setup amalgamation information (will be stored in Graph in the end of the routine)
    RCP<Map> amal_map = SetupAmalgamationData(A, blkSizeInfo_, blockSize);


    // create new CrsGraph for amalgamated matrix (TODO: no shortcut for CrsGraphFactory?)
    RCP<CrsGraph> crsGraph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(amal_map, 10, Xpetra::DynamicProfile);

    for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements());i++) {
      // get global DOF id
      GlobalOrdinal gDofId = A->getRowMap()->getGlobalElement(i);

      // translate to global block id
      GlobalOrdinal globalblockid = GlobalId2GlobalAmalBlockId(gDofId, A, blkSizeInfo_, blockSize);

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
        //if(vals[k]!=0.0) {  // avoid zeros
        if((predrop_ == Teuchos::null && vals[k]!=0.0) ||
           (predrop_ != Teuchos::null && predrop_->Drop(i,gDofId, k,indices[k],gcid,indices,vals) == false)) {
          //colblocks->push_back(globalrowid2globalamalblockid->find(gcid)->second); // add column block id to column ids of amalgamated matrix
          GlobalOrdinal globalcolblockid = GlobalId2GlobalAmalBlockId(gcid, A, blkSizeInfo_, blockSize);
          colblocks->push_back(globalcolblockid);
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
    graph->SetAmalgamationParams(globalamalblockid2myrowid_, globalamalblockid2globalrowid_);
  }

} //namespace MueLu

#endif // MUELU_COALESCEDROPFACTORY_DEF_HPP
