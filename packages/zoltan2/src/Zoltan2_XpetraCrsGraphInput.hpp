// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsGraphInput.hpp

    \brief An input adapter for a Xpetra::CrsGraph.
*/

#ifndef _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_

#include <Xpetra_EpetraCrsGraph.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Zoltan2_GraphInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>

namespace Zoltan2 {

/////////////////////////////////////////////////////////////////////////////
/*! Zoltan2::XpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Xpetra::CrsGraph data.

    TODO -test for memory alloc failure when we resize a vector
    TODO: we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.

    The template parameter is the user's input object - an Epetra
    graph or a templated Tpetra graph 
    or a templated Xpetra::CrsGraph.
*/

template <typename User>
class XpetraCrsGraphInput : public GraphInput<User> {

public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
  typedef Xpetra::TpetraCrsGraph<lno_t, gno_t, node_t> xtgraph_t;
  typedef Xpetra::EpetraCrsGraph xegraph_t;

  /*! Name of input adapter type   TODO make this a trait
   */
  std::string inputAdapterName()const {return std::string("XpetraCrsGraph");}

  /*! Destructor
   */
  ~XpetraCrsGraphInput() { }

  /*! Constructor
   */
  XpetraCrsGraphInput(const RCP<const User> &ingraph):
    ingraph_(ingraph),
    graph_(),
    rowMap_(), colMap_(), edgeOffsets_(),
    vtxWeightDim_(0), edgeWeightDim_(0), coordinateDim_(0),
    edgeWgt_(), vertexWgt_(), xyz_()
  {
    cout << __func__ << " getting Traits from "
         << InputTraits<User>::name() << endl;

    graph_ = XpetraTraits<User>::convertToXpetra(ingraph);
    makeOffsets();
  }

  /* Provide optional vertex coordinates.
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

  /* Provide optional vertex weights.
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

  /* Provide optional edge weights.
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
        int toOffset = (edgeOffsets_[lid-base_] + idx[i]) * edgeWeightDim_;
        int fromOffset = nextWgt + (i * edgeWeightDim_);
        for (int j=0; j < edgeWeightDim_; j++)
          edgeWgt_[toOffset+j] = wgts[fromOffset+j];
      }
      nextWgt += nnbors * edgeWeightDim_;
    }
  }

  ////////////////////////////////////////////////////
  // The GraphInput interface.
  ////////////////////////////////////////////////////

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { 
    return rowMap_->getNodeNumElements(); 
  }

  /*! Returns the number vertices in the entire graph.
   */
  global_size_t getGlobalNumVertices() const { 
    return rowMap_->getGlobalNumElements(); 
  }

  /*! Return whether input adapter wants to use local IDs.
   */

  bool haveLocalIds() const {return true;}

  /*! Return whether local ids are consecutive and if so the base.
   */

  bool haveConsecutiveLocalIds (size_t &base) const
  {
    base = static_cast<size_t>(base_);
    return true;
  }

  /*! Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const { 
    return graph_->getNodeNumEntries();
  }

  /*! Returns the number edges on this entire graph.
   *    what about directional edges, count twice?
   */
  global_size_t getGlobalNumEdges() const { 
    return graph_->getGlobalNumEntries();
  }

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

  /*! Get the list of vertex IDs and their weights.
   */
  void getVertexListCopy(std::vector<gid_t> &ids, 
    std::vector<lid_t> &localIDs, std::vector<scalar_t> &xyz,
    std::vector<scalar_t> &wgt) const
  {
    // Global IDs are in local ID order, so we omit localIDs
    // TODO: For Tpetra and Epetra maps, are the GIDs always
    //    in local ID order?  From looking at the source it
    //    seems to be true.  But not positive.
    localIDs.clear();

    size_t numVtx = this->getLocalNumVertices();
    int nweights = vtxWeightDim_ * numVtx;
    wgt.resize(nweights); 

    if (nweights){
      scalar_t *wTo = &wgt[0];
      const scalar_t *wFrom = &vertexWgt_[0];
      memcpy(wTo, wFrom, sizeof(scalar_t) * nweights);
    }

    ids.resize(numVtx);
    for (unsigned i=0; i < rowMap_->getNodeNumElements(); i++){
      ids[i] =  rowMap_->getGlobalElement(i);
    }

    int ncoords = coordinateDim_ * numVtx;
    xyz.resize(ncoords);

    if (ncoords){
      scalar_t *cTo = &xyz[0];
      const scalar_t *cFrom = &xyz_[0];
      memcpy(cTo, cFrom, sizeof(scalar_t) * ncoords);
    }
  }

  /*! Return a read only view of the data.
   */
  lid_t getVertexListView(const gid_t *&ids, const lid_t *& localIDs,
      const scalar_t *& xyz, const scalar_t *&wgts)
  {
    // TODO we need to verify that gids are actually stored
    //   in lid order
    int nvtx = this->getLocalNumVertices();
    ids = rowMap_->getNodeElementList().getRawPtr();
    localIDs = NULL;   // because haveConsecutiveLocalIds() == true
    xyz = &xyz_[0];
    wgts = &vertexWgt_[0];
    return nvtx;
  }

  /*! Return a copy of the edge IDs and edge weights for a vertex.
   */
  void getVertexEdgeCopy(gid_t vtxId, lid_t localId, 
    std::vector<gid_t> &edgeId, std::vector<scalar_t> &wgts) const
  {
    size_t nvtx = this->getLocalNumVertices();

    if (localId < base_ || localId >= base_+nvtx)
      throw std::runtime_error("invalid local vertex ID");

    edgeId.clear();
    wgts.clear();

    ArrayView<const lno_t> nbors;
    graph_->getLocalRowView(localId, nbors);
    size_t nedges = nbors.size();

    if (nedges > 0){
      edgeId.resize(nedges);
      for (unsigned i=0; i < nedges; i++){
        edgeId[i] = colMap_->getGlobalElement(nbors[i]);
      }

      if (edgeWeightDim_ > 0){
        int offset = edgeOffsets_[localId-base_] * edgeWeightDim_;
        const scalar_t *fromWgt = &edgeWgt_[offset];
        wgts.resize(edgeWeightDim_ * nedges);
        scalar_t *toWgt = &wgts[0];
        memcpy(toWgt, fromWgt, sizeof(scalar_t) * edgeWeightDim_ * nedges);
      }
      else{
        wgts.clear();
      }
    }
  }

  /*! Access to xpetra graph 
   */ 
   
  RCP<const xgraph_t> getXpetraGraph() const
  {
    return graph_;
  }

  RCP<const User> getUserGraph() const
  {
    return ingraph_;
  }


  /*! Return pointers to the edge IDs and edge weights for a vertex.
   *      The edges are available as local IDs only at this point
   *      so this is not defined.  TODO explain better.
   */
  //int getVertexEdgeView(gid_t vtxId, lid_t localId, 
  //  const gid_t *&edgeId, const scalar_t *&wgts) const{}

private:


  RCP<const User > ingraph_;
  RCP<const xgraph_t > graph_;
  RCP<const Xpetra::Map<lid_t, gid_t, node_t> > rowMap_;
  RCP<const Xpetra::Map<lid_t, gid_t, node_t> > colMap_;
  std::vector<int> edgeOffsets_; 
  lid_t base_;

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
    base_ = rowMap_->getMinLocalIndex();
    int numV = rowMap_->getNodeNumElements();
    edgeOffsets_.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      edgeOffsets_[i+1] = edgeOffsets_[i] + graph_->getNumEntriesInLocalRow(i);
    }
  }
};
  
}  //namespace Zoltan2
  
#endif
