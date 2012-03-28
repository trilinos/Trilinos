#ifndef MUELU_GRAPH_DEF_HPP
#define MUELU_GRAPH_DEF_HPP

#include <Xpetra_MapFactory.hpp>

#include "MueLu_Graph_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  Teuchos::ArrayView<const LocalOrdinal> Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getNeighborVertices(LocalOrdinal v) const { 
    Teuchos::ArrayView<const LocalOrdinal> neighborVertices;
    graph_->getLocalRowView(v, neighborVertices); 
    return neighborVertices;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetAmalgamationParams(RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid,RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid) const {
    globalamalblockid2myrowid_ = globalamalblockid2myrowid;
    globalamalblockid2globalrowid_ = globalamalblockid2globalrowid;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMyAmalgamationParams() const {
    return globalamalblockid2myrowid_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetGlobalAmalgamationParams() const {
    return globalamalblockid2globalrowid_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetImportDofMap() const
  {
    //TEUCHOS_TEST_FOR_EXCEPTION(globalamalblockid2globalrowid_==Teuchos::null, Exceptions::RuntimeError, "MueLu::Graph::GetImportDofMap: no amalgamation information! Error.");

    RCP<const Map> nodeMap = GetImportMap();  // import node map

    // special case: 1 dof per node
    if(globalamalblockid2myrowid_ == Teuchos::null &&
       globalamalblockid2globalrowid_ == Teuchos::null) {
      GetOStream(Debug, 0) << "MueLu::Graph::GetImportDofMap: 1 dof per node -> skip reconstruction of import DOF map!" << std::endl;
      // no amalgamation information -> we can assume that we have 1 dof per node
      // just return nodeMap as DOFMap!
      return nodeMap;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(globalamalblockid2globalrowid_==Teuchos::null, Exceptions::RuntimeError, "MueLu::Graph::GetImportDofMap: insufficient amalgamation information. Error");
    TEUCHOS_TEST_FOR_EXCEPTION(globalamalblockid2myrowid_==Teuchos::null    , Exceptions::RuntimeError, "MueLu::Graph::GetImportDofMap: insufficient amalgamation information. Error");

    // build dof map from node map
    RCP<std::vector<GlobalOrdinal> > myDofGIDs = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    for(LocalOrdinal n=0; n<Teuchos::as<LocalOrdinal>(nodeMap->getNodeNumElements()); n++) {
      GlobalOrdinal globalblockid = (GlobalOrdinal) nodeMap->getGlobalElement(n);

      TEUCHOS_TEST_FOR_EXCEPTION(globalamalblockid2globalrowid_->count(globalblockid)<=0, Exceptions::RuntimeError, "MueLu::Graph::GetImportDofMap: empty global block? Error.");
      std::vector<GlobalOrdinal> myrowGIDs = (*globalamalblockid2globalrowid_)[globalblockid];
      TEUCHOS_TEST_FOR_EXCEPTION(myrowGIDs.size()==0, Exceptions::RuntimeError, "MueLu::Graph::GetImportDofMap: no amalgamation information! Error.");

      typename std::vector<GlobalOrdinal>::iterator gidIt;
      for(gidIt = myrowGIDs.begin(); gidIt!=myrowGIDs.end(); gidIt++) {
        myDofGIDs->push_back(*gidIt); // append local row ids
      }
    }

    // generate row dof map for amalgamated matrix with same distribution over all procs as row node map
    Teuchos::ArrayRCP<GlobalOrdinal> arr_myDofGIDs = Teuchos::arcp( myDofGIDs );
    Teuchos::RCP<Map> ImportDofMap = MapFactory::Build(nodeMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), arr_myDofGIDs(), nodeMap->getIndexBase(), nodeMap->getComm());
    return ImportDofMap;
  }


#ifdef MUELU_UNUSED
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  size_t Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetNodeNumGhost() const { 
    /*
      Ray's comments about nGhost:
      Graph->NGhost == graph_->RowMatrixColMap()->NumMyElements() - graph_->OperatorDomainMap()->NumMyElements()
      is basically right. But we've had some issues about how epetra handles empty columns.
      Probably worth discussing this with Jonathan and Chris to see if this is ALWAYS right. 
    */
    size_t nGhost = graph_->getColMap()->getNodeNumElements() - graph_->getDomainMap()->getNodeNumElements();
    if (nGhost < 0) nGhost = 0; // FIXME: size_t is unsigned.
      
    return nGhost;
  }
#endif

  /// Return a simple one-line description of this object.
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  std::string Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    return "MueLu.description()";
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  //using MueLu::Describable::describe; // overloading, not hiding
  //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      //out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1) {
      //out0 << "Linear Algebra: " << toString(lib_) << std::endl;
      //out0 << "PrecType: " << type_ << std::endl;
      //out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      //out0 << "Overlap: " << overlap_ << std::endl;
    }

    if (verbLevel & Debug) {
      graph_->describe(out0, Teuchos::VERB_EXTREME);
    }
  }

}

#endif // MUELU_GRAPH_DEF_HPP
