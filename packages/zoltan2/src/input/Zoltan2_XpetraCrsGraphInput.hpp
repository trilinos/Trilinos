// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsGraphInput.hpp
    \brief Defines XpetraCrsGraphInput class.
*/

#ifndef _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_

#include <Zoltan2_GraphInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Util.hpp>

#include <Xpetra_CrsGraph.hpp>

namespace Zoltan2 {

/*!  \brief Provides access for Zoltan2 to Xpetra::CrsGraph data.

    \todo test for memory alloc failure when we resize a vector
    \todo we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.

    The template parameter is the user's input object - an Epetra
    graph or a templated Tpetra graph 
    or a templated Xpetra::CrsGraph.
*/

template <typename User>
class XpetraCrsGraphInput : public GraphInput<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
  typedef GraphInput<User>       base_adapter_t;
  typedef User user_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraCrsGraphInput() { }

  /*! \brief Constructor
   */
  XpetraCrsGraphInput(const RCP<const User> &ingraph):
    ingraph_(ingraph), graph_(), comm_() ,
    offs_(), eids_()
#if 0
    ,rowMap_(), colMap_(), edgeOffsets_(),
    vtxWeightDim_(0), edgeWeightDim_(0), coordinateDim_(0),
    edgeWgt_(), vertexWgt_(), xyz_()
#endif
  {
    graph_ = XpetraTraits<User>::convertToXpetra(ingraph);
    comm_ = graph_->getComm();
    size_t nvtx = graph_->getNodeNumRows();
    size_t nedges = graph_->getNodeNumEntries();
    Environment env;

    size_t n = nvtx + 1;
    lno_t *offs = new lno_t [n];
    env.localMemoryAssertion(__FILE__, __LINE__, n, offs);

    gid_t *eids = NULL;
    if (nedges){
      eids = new gid_t [nedges];
      env.localMemoryAssertion(__FILE__, __LINE__, nedges, eids);
    }

    offs[0] = 0;
    for (lno_t v=0; v < nvtx; v++){
      ArrayView<const lno_t> nbors;
      graph_->getLocalRowView(v, nbors);
      offs[v+1] = offs[v] + nbors.size();
      for (lno_t e=offs[v], i=0; e < offs[v+1]; e++)
        eids[e] = graph_->getColMap()->getGlobalElement(nbors[i++]);
    }

    offs_ = arcp(offs, 0, n, true);
    eids_ = arcp(eids, 0, nedges, true);
#if 0
    makeOffsets();
#endif
  }

#if 0
  /* \brief Provide optional vertex coordinates.
   *  \param lid  The vertex local id.
   *  \param xyz The coordinates(s) associated with the corresponding vertex
   *    local id.  They should be ordered by vertex by coordinate axis.
   */
  void setVertexCoordinates(std::vector<lid_t> &lid, std::vector<scalar_t> &xyz)
  {
    size_t veclen = xyz.size();
    if (veclen == 0) return;
    
    size_t numIds = lid.size();
    int dim = veclen / numIds;
    if (numIds * dim != veclen)
      throw std::runtime_error("invalid number of coordinates");

    if (coordinateDim_){
      if (dim != coordinateDim_)
        throw std::runtime_error("inconsistent number of coordinates");
    }
    else{
      if (dim > 3)
        throw std::runtime_error("coordinates exceed 3");
      coordinateDim_ = dim;
      xyz_.clear();
      xyz_.resize(veclen,0);  // TODO need an "unset" value
    }

    // TODO - they're always consecutive, right?
    lid_t min = rowMap_->getMinLocalIndex();
    lid_t max = rowMap_->getMaxLocalIndex();

    for (size_t i = 0; i < numIds; i++){
      if ( (lid[i] < min) || (lid[i] > max))
        throw std::runtime_error("invalid vertex local id");
      lid_t to_pos = coordinateDim_ * (lid[i] - min);
      lid_t from_pos = coordinateDim_ * i;
      for (int j=0; j < coordinateDim_; j++){
        xyz_[to_pos++] = xyz[from_pos++];
      }
    }
  }

  /* \brief Provide optional vertex weights.
   *  \param lid  The vertex local id.
   *  \param wgts The weight(s) associated with the corresponding vertex
   *    local id.  Weights should be ordered by vertex by weight coordinate.
   */
  void setVertexWeights(std::vector<lid_t> &lid, std::vector<scalar_t> &wgts)
  {
    size_t veclen = wgts.size();
    if (veclen == 0) return;
    
    size_t numIds = lid.size();
    int dim = veclen / numIds;
    if (numIds * dim != veclen)
      throw std::runtime_error("invalid number of weights");

    if (vtxWeightDim_){
      if (dim != vtxWeightDim_)
        throw std::runtime_error("inconsistent number of weights");
    }
    else{
      vtxWeightDim_ = dim;
      vertexWgt_.clear();
      vertexWgt_.resize(veclen,0);
    }

    // TODO - they're always consecutive, right?
    lid_t min = rowMap_->getMinLocalIndex();
    lid_t max = rowMap_->getMaxLocalIndex();

    for (size_t i = 0; i < numIds; i++){
      if ( (lid[i] < min) || (lid[i] > max))
        throw std::runtime_error("invalid vertex local id");
      lid_t to_pos = vtxWeightDim_ * (lid[i] - min);
      lid_t from_pos = vtxWeightDim_ * i;
      for (int j=0; j < vtxWeightDim_; j++){
        vertexWgt_[to_pos++] = wgts[from_pos++];
      }
    }
  }

  /* \brief Provide optional edge weights.
   *  \param vertexLid  The vertex local id.
   *  \param numNbors   The number of edge weights provided.
   *  \param nborGid    The global vertex id of the neighbor.
   *  \param wgts The weight(s) associated with the corresponding edge.
   *    Weights should be ordered by edge by weight coordinate.
   */
  void setEdgeWeights(std::vector<lid_t> &vertexLid, 
    std::vector<lid_t> &numNbors,
    std::vector<gid_t> &nborGid, std::vector<scalar_t> &wgts )
  {
    lno_t nvtx = vertexLid.size();

    if ((nvtx==0) || (nborGid.size()==0) || (wgts.size()==0))
      return;

    if (edgeWeightDim_ == 0){
      edgeWeightDim_ = wgts.size() / nborGid.size();
      if (edgeWeightDim_ * nborGid.size() != wgts.size())
        throw std::runtime_error("Invalid number of edge weights");
      edgeWgt_.resize(edgeWeightDim_ * getLocalNumEdges(), scalar_t(1));
    }
    else if ((nborGid.size() * edgeWeightDim_) != wgts.size()){
      throw std::runtime_error("Invalid number of edge weights");
    }

    int nextNbor=0, nextWgt=0;

    for (lno_t v=0; v < nvtx; v++){
      int nnbors = numNbors[v];

      if (nnbors < 1)
        continue;

      lid_t lid = vertexLid[v];
      gid_t gid = rowMap_->getGlobalElement(lid);
      std::vector<gid_t> edges;
      std::vector<scalar_t> ewgts;
      getVertexEdgeCopy(gid, lid, edges, ewgts); 

      if (nnbors > edges.size())
        throw std::runtime_error("invalid number of neighbors");

      std::vector<gid_t> nbors(nnbors);
      std::vector<gid_t> idx(nnbors);
      for (int i=0; i < nnbors; i++){
        nbors[i] = nborGid[nextNbor++];
        idx[i] = i;
      }

      if (edges != nbors){
        // TODO make it more efficient to match up edge IDs with their index
        for (int i=0; i < nnbors; i++){
          typename std::vector<gid_t>::iterator loc = std::find(edges.begin(), edges.end(),nbors[i]);
          if (loc == edges.end())
            throw std::runtime_error("Invalid edge global id");
          idx[i] = loc - edges.begin();
        }
      }

      for (int i=0; i < nnbors; i++){
        int toOffset = (edgeOffsets_[lid] + idx[i]) * edgeWeightDim_;
        int fromOffset = nextWgt + (i * edgeWeightDim_);
        for (int j=0; j < edgeWeightDim_; j++)
          edgeWgt_[toOffset+j] = wgts[fromOffset+j];
      }
      nextWgt += nnbors * edgeWeightDim_;
    }
  }
#endif

  /*! \brief Access to xpetra graph 
   */ 
   
  RCP<const xgraph_t> getXpetraGraph() const
  {
    return graph_;
  }

  /*! \brief Access to user's graph 
   */ 
  RCP<const User> getUserGraph() const
  {
    return ingraph_;
  }

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  std::string inputAdapterName()const { return std::string("XpetraCrsGraph");}

  size_t getLocalNumberOfObjects() const { return getLocalNumVertices();}

  int getNumberOfWeightsPerObject() const { return 0;}

  ////////////////////////////////////////////////////
  // The GraphInput interface.
  ////////////////////////////////////////////////////

  /*! \brief Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { 
    return graph_->getNodeNumRows(); 
  }

  /*! \brief Returns the number vertices in the entire graph.
   */
  global_size_t getGlobalNumVertices() const { 
    return graph_->getGlobalNumRows(); 
  }

  /*! \brief Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const { 
    return graph_->getNodeNumEntries();
  }

  /*! \brief Returns the number edges on this entire graph.
   *    what about directional edges, count twice?
   */
  global_size_t getGlobalNumEdges() const { 
    return graph_->getGlobalNumEntries();
  }

#if 0
  /*! Returns the number weights supplied for each vertex.
   */
  int getVertexWeightDim() const { 
    return vtxWeightDim_;
  }

  /*! Returns the number weights supplied for each edge.
   */
  int getEdgeWeightDim() const { 
    return edgeWeightDim_;
  }

  /*! Returns the number of coordinates per vertex
   */
  int getCoordinateDim() const { 
    return coordinateDim_;
  }
#endif

  /*! \brief Return a read only view of the data.
   */
  size_t getVertexListView(const gid_t *&ids,
    const lno_t *&offsets, const gid_t *& edgeId) const
  {
    size_t nvtx = getLocalNumVertices();
    ids = edgeId = NULL;
    offsets = NULL;

    if (nvtx){
      ids = graph_->getRowMap()->getNodeElementList().getRawPtr();
      offsets = offs_.getRawPtr();
      edgeId = eids_.getRawPtr();
    }
    
    return nvtx;
  }

  /*! \brief Repartition a graph that has the same structure as
   *   the graph that instantiated this input adapter.
   */
  template<typename User2>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User2> &solution)
  {
    // Get an import list

    Zoltan2::Environment env;
    size_t len = solution.getNumberOfIds();
    const gid_t *gids = solution.getGlobalIdList();
    const size_t *parts = solution.getPartList();
    ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false);
    ArrayRCP<size_t> partList = arcp(const_cast<size_t *>(parts), 0, len, false);

    ArrayRCP<lno_t> dummyIn;
    ArrayRCP<gid_t> importList;
    ArrayRCP<lno_t> dummyOut;
    size_t numNewVtx;
    const RCP<const Comm<int> > comm = graph_->getComm();

    try{
      numNewVtx = convertPartListToImportList<gid_t, lno_t, lno_t>(
        *comm, partList, gidList, dummyIn, importList, dummyOut);
    }
    Z2_FORWARD_EXCEPTIONS;

    gno_t lsum = numNewVtx;
    gno_t gsum = 0;
    reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum, &gsum);

    RCP<const User> inPtr = rcp(&in, false);

    RCP<const User> outPtr = XpetraTraits<User>::doMigration(
     inPtr, lsum, importList.getRawPtr());

    out = const_cast<User *>(outPtr.get());
    outPtr.release();
    return numNewVtx;
  }

private:

  RCP<const User > ingraph_;
  RCP<const xgraph_t > graph_;
  RCP<const Comm<int> > comm_;

  // FOR NOW  TODO - how to manage these buffers?
  ArrayRCP<const lno_t> offs_;
  ArrayRCP<const gid_t> eids_;

#if 0
  RCP<const Xpetra::Map<lid_t, gid_t, node_t> > rowMap_;
  RCP<const Xpetra::Map<lid_t, gid_t, node_t> > colMap_;
  std::vector<int> edgeOffsets_; 

  int vtxWeightDim_;
  int edgeWeightDim_;
  int coordinateDim_;
  std::vector<scalar_t> edgeWgt_;
  std::vector<scalar_t> vertexWgt_;
  std::vector<scalar_t> xyz_;

  void makeOffsets()
  {
    rowMap_ = graph_->getRowMap();
    colMap_ = graph_->getColMap();
    int numV = rowMap_->getNodeNumElements();
    edgeOffsets_.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      edgeOffsets_[i+1] = edgeOffsets_[i] + graph_->getNumEntriesInLocalRow(i);
    }
  }
#endif
};
  
}  //namespace Zoltan2
  
#endif
