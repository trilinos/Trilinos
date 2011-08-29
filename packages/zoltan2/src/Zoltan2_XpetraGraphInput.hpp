// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraGraphInput.hpp

    \brief The methods used in common for Tpetra and Epetra graphs.
*/


#ifndef _ZOLTAN2_XPETRAGRAPHINPUT_HPP_
#define _ZOLTAN2_XPETRAGRAPHINPUT_HPP_

#include <Zoltan2_Environment.hpp>
#include <Zoltan2_GraphInput.hpp>

namespace Zoltan2 {

/*! Zoltan2::XpetraGraphInput
    \brief XpetraGraphInput defines methods used by both the Tpetra and Epetra graph adapters.

   To simplify working with both Tpetra and Epetra graphs, we use one interface for both.

   TODO: For now we assume the graph (represented as a matrix) is either lower or upper
    triangular.  We should check to see if it is not, and adjust our answers accordingly.
   (i.e. number of edges is 1/2 number of nonzeros).

*/

template<typename Scalar, typename LNO, typename GNO>
  class XpetraGraphInput : public GraphInput<Scalar, LNO, GNO>{
private:

  Xpetra::CrsGraph<Scalar, LNO, GNO> &_graph;
  Teuchos::RCP<const Xpetra::Map<LNO, GNO> > _rowMap;
  Teuchos::RCP<const Teuchos::Comm<int> > _comm;
  Teuchos::ArrayRCP<Scalar> _vtxWgt;
  Teuchos::ArrayRCP<Scalar> _edgeWgt;
  int _vtxWeightDim;
  int _edgeWeightDim;
  std::vector<LNO> _edgeOffsets;   // TODO do we need these?

public:

  /*! Set the graph, the vertex weights and the edge weights.
   *
   *  Tpetra and Epetra Graphs have a concept of local index, which
   *  is contiguous and begins at zero on each process.  The vertexWeights
   *  should be supplied in vertex (row) local index order.  (For multiple
   *  weights the order is by vertex by weight component.)
   *
   *  When the Xpetra::CrsGraph is queried for the entries in a row (the edges
   *  of a vertex), they are listed in a given order.  The edgeWeights
   *  supplied here should be listed in that order.
   */

  void setXpetraCrsGraph(Xpetra::CrsGraph<Scalar, LNO, GNO> &g,
                         Teuchos::ArrayRCP<Scalar> vertexWeights,
                         Teuchos::ArrayRCP<Scalar> edgeWeights){
    _graph=g; 
    _comm=g.getComm();
    _rowMap=g.getRowMap();
    _vtxWeightDim = 0;
    _edgeWeightDim = 0;

    Z2_GLOBAL_INPUT_ASSERTION(_comm, _env, "broke assumption about input", 
        g.isLowerTriangular || g.isUpperTriangular, Z2_BASIC_ASSERTION);

    LNO numVwgts = vertexWeights.size();
    LNO numEwgts = edgeWeights.size();
    LNO numV = this->getLocalNumVertices();
    LNO numE = this->getLocalNumEdges();

    if (numVwgts > 0){
      _vtxWeightDim = numVwgts / numV;
      Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "bad number of weights", 
         _vtxWeightDim * numV == numVwgts, Z2_BASIC_ASSERTION);
    }

    if (numEwgts > 0){
      _edgeWeightDim = numEwgts / numE;
      Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "bad number of weights", 
         _edgeWeightDim * numE == numEwgts, Z2_BASIC_ASSERTION);
    }

    _edgeOffsets.clear();
    _edgeOffsets.resize(numV+1, 0);

    for (LNO i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph.getNumEntriesInLocalRow(i);
    }
  }

  /*! Returns the global number of vertices in the graph.
   */
  GNO getGlobalNumVertices() const { return _rowMap->getGlobalNumElements();}

  /*! Returns the global number edges in the graph.
   */
  GNO getGlobalNumEdges() const { return _graph.getGlobalNumEntries();}

  /*! Returns the number vertices on this process.
   */
  LNO getLocalNumVertices() const { return _rowMap->getNodeNumElements(); }

  /*! Returns the number edges on this process.
   */
  LNO getLocalNumEdges() const { return _graph.getNodeNumEntries();}

  /*! Returns the number weights supplied for each vertex.
   */
  LNO getVertexWeightDim() const { return _vtxWeightDim;}

  /*! Returns the number weights supplied for each edge.
   */
  LNO getEdgeWeightDim() const { return _edgeWeightDim;}

  /*! Returns the minimum global vertex ID on this process.
   */
  GNO getMinGlobalVertex() const {return _rowMap->getMinGlobalIndex(); }

  /*! Returns the maximum global vertex ID on this process.
   */
  GNO getMaxGlobalVertex() const  {return _rowMap->getMaxGlobalIndex(); }

  /*! Returns the minimum local global vertex ID everywhere.
   */
  GNO getMinAllGlobalVertex() const {return _rowMap->getMinAllGlobalIndex(); }

  /*! Returns the maximum local global vertex ID everywhere.
   */
  GNO getMaxAllGlobalVertex() const {return _rowMap->getMaxAllGlobalIndex(); }

  /*! Translate vertex ID from global to local.
      \param vtxId The global ID of a vertex on this process.
      \result the local ID for the specified vertex.
   */
  LNO getLocalVertex(GNO vtxId) const {return _rowMap->getLocalElement(vtxId);}

  /*! Translate vertex ID from local to global.
      \param vtxId The local ID of a vertex on this process.
      \result the global ID for the specified vertex.
   */
  GNO getGlobalVertex(LNO vtxId) const {return _rowMap->getGlobalElement(vtxId);}

  /*! Find the owner of and local ID of a list vertices.
      \param vtxId   A list of global vertex IDs
      \param nodeId  A view of memory allocated by the caller
      \param nodeVtxLid A view of memory allocated by the caller
      \result nodeId  A list of the process rank for the owner corresponding to each vertex in vtxId
      \result nodeVtxLid  A list of the local ID at the owning process for each vertex vtxId

       This method must be called by all processes in the communicator.
   */
  void getRemoteVertexList(const Teuchos::ArrayView<const GNO> &vtxId,
    const Teuchos::ArrayView<const int> &nodeId,
    const Teuchos::ArrayView<const LNO> &nodeVtxLid) const {

      Xpetra::LookupStatus s = getRemoteIndexList(vtxId, nodeId, nodeVtxLid);
      Z2_GLOBAL_INPUT_ASSERTION(_comm, _env, "getRemoteVertexList invalid IDs", 
         s == Xpetra::AllIDsPresent, Z2_BASIC_ASSERTION);
    }

  /*! Find the process rank for the owner of each vertex in a list.
      \param vtxId   A list of global vertex IDs.
      \param nodeId  A view of memory allocated by the caller.
      \result nodeId  A list of the process rank for the owner corresponding to each vertex in vtxId.

       This method must be called by all processes in the communicator.
   */
  void getRemoteVertexList(const Teuchos::ArrayView<const GNO> &vtxId,
    const Teuchos::ArrayView<const int> &nodeId) const{

      Xpetra::LookupStatus s = getRemoteIndexList(vtxId, nodeId);
      Z2_GLOBAL_INPUT_ASSERTION(_comm, _env, "getRemoteVertexList invalid IDs", 
         s == Xpetra::AllIDsPresent, Z2_BASIC_ASSERTION);
    }

  /*! Get the list of vertex global IDs and their weights.
      \param ids a view of an array allocated by the caller.
      \param wgts a view of an array allocated by the caller.
      \result ids a list of the global IDs of each vertex on this process.
      \result wgts a list of the weight or weights associated with each vertex in the ids list.  Weights are
                     listed by vertex by weight component.

      TODO - we should support a view instead of copy when possible
   */
  void getGlobalVertexList(Teuchos::ArrayView<const GN0> &ids, 
    Teuchos::ArrayView<const Scalar> &wgt) const{

    LNO numVtx = this->getLocalNumVertices();

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid input array size", 
      (wgt.size() >= numVtx*_vtxWeightDim) && (ids.size() >= numVtx), 
      Z2_BASIC_ASSERTION);

    LNO nweights = _vtxWeightDim * numVtx;
    if (nweights){
      Scalar *wTo = wgt.getRawPtr();
      Scalar *wFrom = _vtxWgt.getRawPtr();
      memcpy(wTo, wFrom, sizeof(Scalar) * nweights);
    }

    getGlobalVertexList(ids);
  }

  /*! Get the list of vertex global IDs.
      \param ids a view of an array allocated by the caller.
      \result ids a list of the global IDs of each vertex on this process.
   */
  void getGlobalVertexList(Teuchos::ArrayView<const GN0> &ids) const{

    LNO numVtx = this->getLocalNumVertices();

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid input array size", 
      ids.size() >= numVtx, Z2_BASIC_ASSERTION);

    GNO *idTo= ids.getRawPtr();
    for (LNO i=0; i < _rowMap->getNodeNumElements(); i++){
      *idTo++ = _rowMap->GetGlobalElement(i);
    }
  }

  /*! Does the supplied local vertex ID exist on this process.
      \result true if the ID is a vertex local ID on this process.
      \result false if the ID is not a vertex local ID on this process.
   */
  bool isNodeLocalVertex(LNO vtxId) const { return _rowMap->isNodeLocalElement(vtxId); }

  /*! Does the supplied global vertex ID exist on this process.
      \result true if the ID is a vertex global ID on this process.
      \result false if the ID is not a vertex global ID on this process.
   */
  bool isNodeGlobalVertex(GNO vtxId) const  { return _rowMap->isNodeGlobalElement(vtxId); }

  /*! Are the vertex global IDs contiguous across process?
      \result true if global IDs process contiguously from one process to the next
      \result false if global IDs are not contiguous
   */
  bool isContiguous() const { return _rowMap->isContiguous(); }

  /*! Are the graph distributed across more than this one process.
      \result true if the graph is distributed.
      \result false if the graph is not distributed.
   */
  bool isDistributed() const {return _rowMap->isDistributed(); }

  /*! Returns the sum of the number of neighbors of all vertices
   */
  GNO getGlobalNumberOfNbors() const { return _graph.getGlobalNumEntries(); }

  /*! Returns the sum of the number of neighbors of vertices on this process
   */
  LNO getLocalNumberOfNbors() const { return _graph.getNodeNumEntries(); }

  /*! Returns the number of neighbors for a given vertex.
      \param vtxId The global ID of a vertex owned by this process.
   */
  LNO getNumNborsOfGlobalVertex(GNO vtxId) const { return _graph.getNumEntriesInGlobalRow(vtxId); }

  /*! Returns the number of neighbors for a given vertex.
      \param vtxId The local ID of a vertex.
   */
  LNO getNumNborsOfLocalVertex(LNO vtxId) const { return _graph.getNumEntriesInLocalRow(vtxId); }

  /*! Returns the global number of self-edges.  TODO: change the name?
   */
  GNO getGlobalNumDiags() const { return _graph.getGlobalNumDiags(); }

  /*! Returns the local number of self-edges.
   */
  LNO getLocalNumDiags() const { return _graph.getNodeNumDiags(); }

  /*! Returns the global maximum vertex degree.
   */
  LNO getGlobalMaxNumVertexNbors() const { return _graph.getGlobalMaxNumRowEntries(); }

  /*! Returns the local maximum vertex degree.
   */
  LNO getLocalMaxNumVertexNbors() const { return _graph.getNodeMaxNumRowEntries(); }

  /*! Edge IDs may or may not be supplied in increasing order in queries.
      \return true if edge IDs are supplied in increasing order.
      \return false if edge IDs may not be supplied in increasing order.

   TODO - move this to Tpetra::CrsGraph adapter.  Not available from Xpetra
  bool isSorted() const { return _graph.
   */

  /*! Vertex weights are optional in a user supplied graph.
      \return true if vertex weights are available
      \return false if vertex weights are not available
   */
  bool hasVertexWeights() const {return _vtxWeightDim > 0; }

  /*! Edge weights are optional in a user supplied graph.
      \return true if edge weights are available
      \return false if edge weights are not available
   */
  bool hasEdgeWeights() const {return _edgeWeightDim > 0; }

  /*! Obtain a copy of the edge IDs of the input vertex global ID
      \param vtxId  global ID for a vertex on this process
      \param edgeId user allocated array for edge IDs
      \result edgeId global edge IDs are copied to user's array
      \result numEdges  the number of edgeIds written to the array
   */
  void getGlobalVertexCopy(GNO vtxId, const Teuchos::ArrayView<GNO> &edgeId, 
    LNO &numEdges) const{

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid global vertex ID", 
     this->isNodeGlobalVertex(vtxId), Z2_BASIC_ASSERTION);

    numEdges = this->getNumNborsOfGlobalVertex(vtxId)

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid global vertex ID", 
     edgeId.size() >= numEdges, ASIC_ASSERTION);

    GNO *toId = edgeId.getRawPtr();
    Teuchos::ArrayView<GNO> nbors;
    graph.getGlobalRowView(vtxId, nbors);    // TODO trap errors
    GNO *fromId = nbors.getRawPtr();

    memcpy(toId, fromId, sizeof(GNO) * numEdges);
  }

  /*! Obtain a copy of the edge IDs of the input vertex local ID
      \param vtxId  a vertex local ID 
      \param edgeId user allocated array for edge IDs
      \result edgeId global edge IDs are copied to user's array
      \result numEdges  the number of edgeIds written to the array
   */
  void getLocalVertexCopy(LNO vtxId, const Teuchos::ArrayView<GNO> &edgeId, 
    LNO &numEdges) const{

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid local vertex ID", 
     this->isNodeLocalVertex(vtxId), Z2_BASIC_ASSERTION);

    GNO globalVtxId = this->getGlobalVertex(vtxId);

    this->getGlobalVertexCopy(globalVtxId, edgeId, numEdges);
  }

  /*! Obtain a read-only view of the edge IDs of the input vertex
      \param vtxId  global ID for a vertex on this process
      \param edgeId user supplied ArrayView object
      \result edgeId the ArrayView is set to a read-only view of the edge Ids
   */
  void getGlobalVertexView(GNO vtxId, Teuchos::ArrayView<const GNO> &edgeId) const{

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid global vertex ID", 
     this->isNodeGlobalVertex(vtxId), Z2_BASIC_ASSERTION);

    _graph.getGlobalRowView(vtxId, edgeId);
  }

  /*! Obtain a read-only view of the edge IDs of the input vertex
      \param vtxId  a local vertex ID
      \param edgeId user supplied ArrayView object
      \result edgeId the ArrayView is set to a read-only view of the edge GLOBAL IDs

      Note: This behavior differs from Xpetra::CrsGraph which returns edge local IDs.
   */
  void getLocalVertexView(LNO vtxId, Teuchos::ArrayView<const GNO> &edgeId) const{

    Z2_LOCAL_INPUT_ASSERTION(_comm, _env, "invalid local vertex ID", 
     this->isNodeLocalVertex(vtxId), Z2_BASIC_ASSERTION);

    GNO globalVtxId = this->getGlobalVertex(vtxId);

    this->getGlobalVertexView(globalVtxId, edgeId);
  }


};
  
  
}  //namespace Zoltan2
  
#endif
