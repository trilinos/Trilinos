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

namespace Zoltan2 {

/*! Zoltan2::XpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Xpetra::CrsGraph data.

    The template parameter is the weight type.  Xpetra local and global IDs 
    are ints.
*/
    // TODO -test for memory alloc failure when we resize a vector

CONSISTENT_CLASS_TEMPLATE_LINE
class XpetraCrsGraphInput : public GraphInput<CONSISTENT_TEMPLATE_PARAMS> {
private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// local typedefs should go here
#endif

  RCP<Xpetra::CrsGraph<LID, GID, Node> > _graph;
  RCP<const Xpetra::Map<LID, GID, Node> > _rowMap;
  RCP<const Xpetra::Map<LID, GID, Node> > _colMap;
  std::vector<int> _edgeOffsets; 

  int _vtxWeightDim;
  int _edgeWeightDim;
  int _coordinateDim;
  int _base;
  std::vector<Scalar> _edgeWgt;
  std::vector<Scalar> _vertexWgt;
  std::vector<Scalar> _xyz;

  void makeOffsets()
  {
    _rowMap = _graph->getRowMap();
    _colMap = _graph->getColMap();
    _base = _rowMap->getMinLocalIndex();
    int numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
  }

public:

  // TODO: should this be part of InputAdapter interface?
  std::string inputAdapterName()const {return std::string("XpetraCrsGraph");}

  ~XpetraCrsGraphInput() { }

  /*! Default constructor - can't build a valid object this way
   *    TODO - remove?
   */
  XpetraCrsGraphInput(): _graph(), _rowMap(), _colMap(), _edgeOffsets(),
    _vtxWeightDim(0), _edgeWeightDim(0), _coordinateDim(0),
    _edgeWgt(), _vertexWgt(), _xyz() {}

  /*! Constructor with an Xpetra::CrsGraph
   */
  XpetraCrsGraphInput(
    RCP<Xpetra::CrsGraph<LID, GID, Node> > graph):
    _graph(graph), _rowMap(), _colMap(), _edgeOffsets(),
    _vtxWeightDim(0), _edgeWeightDim(0), _coordinateDim(0),
    _edgeWgt(), _vertexWgt(), _xyz()
  {
    makeOffsets();
  }

  /*! Constructor with a Tpetra::CrsGraph
   */
  XpetraCrsGraphInput(
    RCP<Tpetra::CrsGraph<LID, GID, Node> > graph):
    _graph(), _rowMap(), _colMap(), _edgeOffsets(),
    _vtxWeightDim(0), _edgeWeightDim(0), _coordinateDim(0),
    _edgeWgt(), _vertexWgt(), _xyz()
  {
     Xpetra::TpetraCrsGraph<LID, GID, Node> *xgraph =
       new Xpetra::TpetraCrsGraph<LID, GID, Node>(graph);

    _graph = Teuchos::rcp_implicit_cast<Xpetra::CrsGraph<LID, GID, Node> >(
      Teuchos::rcp(xgraph));
    makeOffsets();
  }

  /*! Constructor with an Epetra_CrsGraph
   */
  XpetraCrsGraphInput(RCP<Epetra_CrsGraph> graph):
    _graph(), _rowMap(), _colMap(), _edgeOffsets(),
    _vtxWeightDim(0), _edgeWeightDim(0), _coordinateDim(0),
    _edgeWgt(), _vertexWgt(), _xyz()
  {
     Xpetra::EpetraCrsGraph *xgraph = new Xpetra::EpetraCrsGraph(graph);

    _graph = Teuchos::rcp_implicit_cast<Xpetra::CrsGraph<LID, GID, Node> >(
      Teuchos::rcp(xgraph));
    makeOffsets();
  }

  /* Provide optional vertex coordinates.
   *  \param lid  The vertex local id.
   *  \param xyz The coordinates(s) associated with the corresponding vertex
   *    local id.  They should be ordered by vertex by coordinate axis.
   */
  void setVertexCoordinates(std::vector<LID> &lid, std::vector<Scalar> &xyz)
  {
    size_t veclen = xyz.size();
    if (veclen == 0) return;
    
    size_t numIds = lid.size();
    int dim = veclen / numIds;
    if (numIds * dim != veclen)
      throw std::runtime_error("invalid number of coordinates");

    if (_coordinateDim){
      if (dim != _coordinateDim)
        throw std::runtime_error("inconsistent number of coordinates");
    }
    else{
      if (dim > 3)
        throw std::runtime_error("coordinates exceed 3");
      _coordinateDim = dim;
      _xyz.clear();
      _xyz.resize(veclen,0);  // TODO need an "unset" value
    }

    // TODO - they're always consecutive, right?
    LID min = _rowMap->getMinLocalIndex();
    LID max = _rowMap->getMaxLocalIndex();

    for (size_t i = 0; i < numIds; i++){
      if ( (lid[i] < min) || (lid[i] > max))
        throw std::runtime_error("invalid vertex local id");
      LID to_pos = _coordinateDim * (lid[i] - min);
      LID from_pos = _coordinateDim * i;
      for (int j=0; j < _coordinateDim; j++){
        _xyz[to_pos++] = xyz[from_pos++];
      }
    }
  }

  /* Provide optional vertex weights.
   *  \param lid  The vertex local id.
   *  \param wgts The weight(s) associated with the corresponding vertex
   *    local id.  Weights should be ordered by vertex by weight coordinate.
   */
  void setVertexWeights(std::vector<LID> &lid, std::vector<Scalar> &wgts)
  {
    size_t veclen = wgts.size();
    if (veclen == 0) return;
    
    size_t numIds = lid.size();
    int dim = veclen / numIds;
    if (numIds * dim != veclen)
      throw std::runtime_error("invalid number of weights");

    if (_vtxWeightDim){
      if (dim != _vtxWeightDim)
        throw std::runtime_error("inconsistent number of weights");
    }
    else{
      _vtxWeightDim = dim;
      _vertexWgt.clear();
      _vertexWgt.resize(veclen,0);
    }

    // TODO - they're always consecutive, right?
    LID min = _rowMap->getMinLocalIndex();
    LID max = _rowMap->getMaxLocalIndex();

    for (size_t i = 0; i < numIds; i++){
      if ( (lid[i] < min) || (lid[i] > max))
        throw std::runtime_error("invalid vertex local id");
      LID to_pos = _vtxWeightDim * (lid[i] - min);
      LID from_pos = _vtxWeightDim * i;
      for (int j=0; j < _vtxWeightDim; j++){
        _vertexWgt[to_pos++] = wgts[from_pos++];
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
  void setEdgeWeights(std::vector<LID> &vertexLid, 
    std::vector<LID> &numNbors,
    std::vector<GID> &nborGid, std::vector<Scalar> &wgts )
  {
    LNO nvtx = vertexLid.size();

    if ((nvtx==0) || (nborGid.size()==0) || (wgts.size()==0))
      return;

    if (_edgeWeightDim == 0){
      _edgeWeightDim = wgts.size() / nborGid.size();
      if (_edgeWeightDim * nborGid.size() != wgts.size())
        throw std::runtime_error("Invalid number of edge weights");
      _edgeWgt.resize(_edgeWeightDim * getLocalNumEdges(), Scalar(1));
    }
    else if ((nborGid.size() * _edgeWeightDim) != wgts.size()){
      throw std::runtime_error("Invalid number of edge weights");
    }

    int nextNbor=0, nextWgt=0;

    for (LNO v=0; v < nvtx; v++){
      int nnbors = numNbors[v];

      if (nnbors < 1)
        continue;

      LID lid = vertexLid[v];
      GID gid = _rowMap->getGlobalElement(lid);
      std::vector<GID> edges;
      std::vector<Scalar> ewgts;
      getVertexEdgeCopy(gid, lid, edges, ewgts); 

      if (nnbors > edges.size())
        throw std::runtime_error("invalid number of neighbors");

      std::vector<GID> nbors(nnbors);
      std::vector<GID> idx(nnbors);
      for (int i=0; i < nnbors; i++){
        nbors[i] = nborGid[nextNbor++];
        idx[i] = i;
      }

      if (edges != nbors){
        // TODO make it more efficient to match up edge IDs with their index
        for (int i=0; i < nnbors; i++){
          typename std::vector<GID>::iterator loc = std::find(edges.begin(), edges.end(),nbors[i]);
          if (loc == edges.end())
            throw std::runtime_error("Invalid edge global id");
          idx[i] = loc - edges.begin();
        }
      }

      for (int i=0; i < nnbors; i++){
        int toOffset = (_edgeOffsets[lid-_base] + idx[i]) * _edgeWeightDim;
        int fromOffset = nextWgt + (i * _edgeWeightDim);
        for (int j=0; j < _edgeWeightDim; j++)
          _edgeWgt[toOffset+j] = wgts[fromOffset+j];
      }
      nextWgt += nnbors * _edgeWeightDim;
    }
  }

  ////////////////////////////////////////////////////
  // The GraphInput interface.
  ////////////////////////////////////////////////////

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { 
    return _rowMap->getNodeNumElements(); 
  }

  /*! Returns the number vertices in the entire graph.
   */
  global_size_t getGlobalNumVertices() const { 
    return _rowMap->getGlobalNumElements(); 
  }

  /*! Return whether input adapter wants to use local IDs.
   */

  bool haveLocalIds() const {return true;}

  /*! Return whether local ids are consecutive and if so the base.
   */

  bool haveConsecutiveLocalIds (size_t &base) const
  {
    base = static_cast<size_t>(_base);
    return true;
  }

  /*! Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const { 
    return _graph->getNodeNumEntries();
  }

  /*! Returns the number edges on this entire graph.
   *    what about directional edges, count twice?
   */
  global_size_t getGlobalNumEdges() const { 
    return _graph->getGlobalNumEntries();
  }

  /*! Returns the number weights supplied for each vertex.
   */
  int getVertexWeightDim() const { 
    return _vtxWeightDim;
  }

  /*! Returns the number weights supplied for each edge.
   */
  int getEdgeWeightDim() const { 
    return _edgeWeightDim;
  }

  /*! Returns the number of coordinates per vertex
   */
  int getCoordinateDim() const { 
    return _coordinateDim;
  }

  /*! Get the list of vertex IDs and their weights.
   */
  void getVertexListCopy(std::vector<GID> &ids, 
    std::vector<LID> &localIDs, std::vector<Scalar> &xyz,
    std::vector<Scalar> &wgt) const
  {
    // Global IDs are in local ID order, so we omit localIDs
    // TODO: For Tpetra and Epetra maps, are the GIDs always
    //    in local ID order?  From looking at the source it
    //    seems to be true.  But not positive.
    localIDs.clear();

    size_t numVtx = this->getLocalNumVertices();
    int nweights = _vtxWeightDim * numVtx;
    wgt.resize(nweights); 

    if (nweights){
      Scalar *wTo = &wgt[0];
      const Scalar *wFrom = &_vertexWgt[0];
      memcpy(wTo, wFrom, sizeof(Scalar) * nweights);
    }

    ids.resize(numVtx);
    for (unsigned i=0; i < _rowMap->getNodeNumElements(); i++){
      ids[i] =  _rowMap->getGlobalElement(i);
    }

    int ncoords = _coordinateDim * numVtx;
    xyz.resize(ncoords);

    if (ncoords){
      Scalar *cTo = &xyz[0];
      const Scalar *cFrom = &_xyz[0];
      memcpy(cTo, cFrom, sizeof(Scalar) * ncoords);
    }
  }

  /*! Return a read only view of the data.
   */
  LID getVertexListView(const GID *&ids, const LID *& localIDs,
      const Scalar *& xyz, const Scalar *&wgts)
  {
    // TODO we need to verify that gids are actually stored
    //   in lid order
    int nvtx = this->getLocalNumVertices();
    ids = _rowMap->getNodeElementList().getRawPtr();
    localIDs = NULL;   // because haveConsecutiveLocalIds() == true
    xyz = &_xyz[0];
    wgts = &_vertexWgt[0];
    return nvtx;
  }

  /*! Return a copy of the edge IDs and edge weights for a vertex.
   */
  void getVertexEdgeCopy(GID vtxId, LID localId, 
    std::vector<GID> &edgeId, std::vector<Scalar> &wgts) const
  {
    size_t nvtx = this->getLocalNumVertices();

    if (localId < _base || localId >= _base+nvtx)
      throw std::runtime_error("invalid local vertex ID");

    edgeId.clear();
    wgts.clear();

    ArrayView<const LNO> nbors;
    _graph->getLocalRowView(localId, nbors);
    size_t nedges = nbors.size();

    if (nedges > 0){
      edgeId.resize(nedges);
      for (unsigned i=0; i < nedges; i++){
        edgeId[i] = _colMap->getGlobalElement(nbors[i]);
      }

      if (_edgeWeightDim > 0){
        int offset = _edgeOffsets[localId-_base] * _edgeWeightDim;
        const Scalar *fromWgt = &_edgeWgt[offset];
        wgts.resize(_edgeWeightDim * nedges);
        Scalar *toWgt = &wgts[0];
        memcpy(toWgt, fromWgt, sizeof(Scalar) * _edgeWeightDim * nedges);
      }
      else{
        wgts.clear();
      }
    }
  }

  /*! Return pointers to the edge IDs and edge weights for a vertex.
   *      The edges are available as local IDs only at this point
   *      so this is not defined.  TODO explain better.
   */
  //int getVertexEdgeView(GID vtxId, LID localId, 
  //  const GID *&edgeId, const Scalar *&wgts) const{}
};
  
}  //namespace Zoltan2
  
#endif
