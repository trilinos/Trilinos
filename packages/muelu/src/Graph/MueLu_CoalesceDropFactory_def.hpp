#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include "Xpetra_Operator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"

namespace MueLu {

  static const std::string color_esc = "\x1b[";
  static const std::string color_std = "39;49;00m";
  static const std::string color_purple = "35m";

  /*!
    @class CoalesceDropFactory
    @brief Factory for creating a graph base on a given matrix.

    Factory for creating graphs from matrices with entries selectively dropped.
  
    - TODO This factory is very incomplete.
    - TODO The Build method simply builds the matrix graph with no dropping.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CoalesceDropFactory : public SingleLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    CoalesceDropFactory(RCP<const FactoryBase> AFact = Teuchos::null)
      : AFact_(AFact), fixedBlkSize_(true)
    { }

    //! Destructor
    virtual ~CoalesceDropFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      currentLevel.DeclareInput("A", AFact_.get());
    }

    /// set fixed block size
    void SetFixedBlockSize(GO blksize) {
    	blksize_ = blksize;
    	fixedBlkSize_ = true;
    	GetOStream(Debug, 0) << color_esc << color_purple << "CoalesceDropFactory::SetFixedBlockSize()" << color_esc << color_std << std::endl;
    }

    /// set predrop function
    void SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &predrop) { predrop_ = predrop; }

    // todo: method that takes a block map...

    //@}

    void Build(Level &currentLevel) const {
      RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

      // pre-dropping
      RCP<Graph> graph;
      if (predrop_ == Teuchos::null) {
        graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));
      } else {
        graph = predrop_->Drop(A);
      }

      // coalesce

      currentLevel.Set("Graph", graph, this);

      // post-dropping?

    } // Build

    // amalgamation routine
    // input: RCP<Operator> A: input matrix (unamalgamated)
    // input/output: RCP<Graph>: input: empty graph, output: graph of amalgamated matrix for aggregation routine
    void Amalgamate(const RCP<Operator>& A, RCP<Graph>& graph) const {

      int nUnamalgamatedBlockSize = 3; // TODO fix me, change me to BlkSize???
      int nPDE = 2; // TODO fix me, find something smart to switch BlkSize from 2 (finest level) to 3 (all other levels...)

      // map: global block id of amalagamated matrix -> vector of local row ids of unamalgamated matrix (only for global block ids of current proc)
      RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid;
      // vector with global block ids for amalagamated matrix on current proc
      RCP<std::vector<GlobalOrdinal> > globalamalblockids; // TODO remove me

      globalamalblockid2myrowid = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<LocalOrdinal> >);
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
        } else {
          (*globalamalblockid2myrowid)[globalblockid] = std::vector<LocalOrdinal>(1,i);
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

      // inter processor communication: sum up number of block ids
      GlobalOrdinal num_blockids = 0;
      Teuchos::reduceAll<int,GlobalOrdinal>(*(map->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, &num_blockids );
      // TODO: check me: is num_blockids = map->getGlobalNumElements()/nUnamalgamatedBlockSize???
      // for constant block size we can avoid the communication and just use above formula!
      // for variable block size, this information has to be provided

      // vector for block sizes
      std::vector<int> myblockid2blocksize;
      for(GlobalOrdinal i=0; i<cnt_amalRows; i++) {
        myblockid2blocksize.push_back(nPDE); // fixme for variable block size
      }
      // now we know, how many blocks the amalagamated matrix has

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

      // create new CrsGraph for amalgamated matrix (TODO: no shortcut for CrsGraphFactory?)
      RCP<CrsGraph> crsGraph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(amal_map, 10, Xpetra::DynamicProfile);

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

      crsGraph->fillComplete(amal_map,amal_map);

      graph = rcp(new Graph(A->getCrsGraph()/*crsGraph*/, "amalgamated graph of A")); // TODO change me

      // store information in Graph object for unamalgamation of vectors
      graph->SetAmalgamationParams(globalamalblockid2myrowid, globalamalblockids);
      // TODO return graph...
      //graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));

    }

  private:
    //! A Factory
    RCP<const FactoryBase> AFact_;

    /// blocksize for fixed blocksize setup
    GO blksize_;

    /// are we doing fixed or variable blocks
    bool fixedBlkSize_;

    /// pre-drop function
    RCP<PreDropFunctionBaseClass> predrop_;

  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT
#endif // MUELU_COALESCEDROPFACTORY_DEF_HPP
