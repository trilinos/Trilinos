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

#include <Zoltan2_InputAdapter.hpp>
#include <Xpetra_EpetraCrsGraph.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Teuchos_RCP.hpp>


namespace Zoltan2 {

/*! Zoltan2::XpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Xpetra::CrsGraph data.

    The template parameter is the weight type.  Xpetra local and global IDs 
    are ints.
*/

CONSISTENT_CLASS_TEMPLATE_LINE
class XpetraCrsGraphInput : public InputAdapter {
private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#endif

  bool _valid;
  Teuchos::RCP<Xpetra::CrsGraph> _graph;
  Teuchos::RCP<const Xpetra::Map<int, int> > _rowMap;
  std::vector<int> _edgeOffsets; 

  int _vtxWeightDim;
  int _edgeWeightDim;
  int _coordinateDim;
  std::vector<Scalar> _edgeWgt;
  std::vector<Scalar> _vertexWgt;
  std::vector<Scalar> _xyz;

public:

  // TODO: should this be part of InputAdapter interface?
  std::string inputAdapterName()const {return std::string("XpetraCrsGraph");}

  //~XpetraCrsGraphInput() { }

  /*! Default constructor - can't build a valid object this way
   */
  XpetraCrsGraphInput():_valid(false), _graph(), _rowMap(), _edgeOffsets() {}

  /*! Constructor with an Xpetra::CrsGraph
   */
  XpetraCrsGraphInput(
    Teuchos::RCP<Xpetra::CrsGraph<CONSISTENT_TEMPLATE_PARAMS> > graph):
    _valid(false), _graph(graph) 
  {
    _rowMap = _graph->getRowMap();
    int numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    _valid=true;
  }

  /*! Constructor with a Tpetra::CrsGraph
   */
  XpetraCrsGraphInput(
    Teuchos::RCP<Tpetra::CrsGraph<CONSISTENT_TEMPLATE_PARAMS> > graph):
    _valid(false)
  {
     Xpetra::TpetraCrsGraph<CONSISTENT_TEMPLATE_PARAMS> *xgraph =
       new Xpetra::TpetraCrsGraph<CONSISTENT_TEMPLATE_PARAMS>(graph);

    _graph = Teuchos::rcp_implicit_cast<Xpetra::CrsGraph<CONSISTENT_TEMPLATE_PARAMS> >(Teuchos::rcp(xgraph));

    _rowMap = _graph->getRowMap();
    int numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    _valid=true;
  }

  /*! Constructor with an Epetra_CrsGraph
   */
  XpetraCrsGraphInput(Teuchos::RCP<Epetra_CrsGraph> graph):
    _valid(false)
  {
     Xpetra::EpetraCrsGraph *xgraph = new Xpetra::EpetraCrsGraph(graph);

    _graph = Teuchos::rcp_implicit_cast<Xpetra::CrsGraph>(Teuchos::rcp(xgraph));

    _rowMap = _graph->getRowMap();
    int numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    _valid=true;
  }

  /* Provide optional vertex coordinates.
   *  \param lid  The vertex local id.
   *  \param xyz The coordinates(s) associated with the corresponding vertex
   *    local id.  They should be ordered by vertex by coordinate axis.
   */
  void setVertexCoordinates(std::vector<AppLID> &lid, std::vector<Scalar> &xyz)
  {
    if (!_valid)
      throw std::runtime_error("improperly constructed adapter");

    size_t veclen = xyz.size();
    if (veclen == 0) return;
    
    size_t numIds = lid.size();
    int dim = veclen / numIds;
    if (numIds * dim != veclen)
      throw std::runtime_error("invalid number of coordinates");

    if (_coordinateDim){
      if (dim != _coordinateDIm))
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
    LNO min = _rowMap->getMinLocalIndex();
    LNO max = _rowMap->getMaxLocalIndex();

    for (size_t i = 0; i < numIds; i++){
      if ( (lid[i] < min) || (lid[i] > max))
        throw std::runtime_error("invalid vertex local id");
      LNO to_pos = _coordinateDim * (lid[i] - min);
      LNO from_pos = _coordinateDim * i;
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
  void setVertexWeights(std::vector<AppLID> &lid, std::vector<Scalar> &wgts)
  {
    if (!_valid)
      throw std::runtime_error("improperly constructed adapter");

    size_t veclen = wgts.size();
    if (veclen == 0) return;
    
    size_t numIds = lid.size();
    int dim = veclen / numIds;
    if (numIds * dim != veclen)
      throw std::runtime_error("invalid number of weights");

    if (_vtxWeightDim){
      if (dim != _vtxWeightDIm))
        throw std::runtime_error("inconsistent number of weights");
    }
    else{
      _vtxWeightDim = dim;
      _vertexWgt.clear();
      _vertexWgt.resize(veclen,0);
    }

    // TODO - they're always consecutive, right?
    LNO min = _rowMap->getMinLocalIndex();
    LNO max = _rowMap->getMaxLocalIndex();

    for (size_t i = 0; i < numIds; i++){
      if ( (lid[i] < min) || (lid[i] > max))
        throw std::runtime_error("invalid vertex local id");
      LNO to_pos = _vtxWeightDim * (lid[i] - min);
      LNO from_pos = _vtxWeightDim * i;
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
  void setEdgeWeights(std::vector<AppLID> &vertexLid, 
    std::vector<LNO> &numNbors,
    std::vector<AppGID> &nborGid, std::vertex<Scalar> &wgts )
  {
    if (!_valid)
      throw std::runtime_error("improperly constructed adapter");

    // TODO write this

  }

  // TODO: should this be part of InputAdapter interface?
  bool validInput() { return _valid;}

  // The GraphInput interface.

  /*! Returns the number vertices on this process.
   */
  int getLocalNumVertices() const { 
    if (!_valid)
      throw std::runtime_error("invalid input object");
    return _rowMap->getNodeNumElements(); 
  }

  /*! Returns the number edges on this process.
   */
  int getLocalNumEdges() const { 
    if (!_valid)
      throw std::runtime_error("invalid input object");
    return _graph->getNodeNumEntries();
  }

  /*! Returns the number weights supplied for each vertex.
   */
  int getVertexWeightDim() const { 
    if (!_valid)
      throw std::runtime_error("invalid input object");
    
    return _vtxWeightDim;
  }

  /*! Returns the number weights supplied for each edge.
   */
  int getEdgeWeightDim() const { 
    if (!_valid)
      throw std::runtime_error("invalid input object");
  
    return _edgeWeightDim;
  }

  /*! Returns the number of coordinates per vertex
   */
  int getCoordinateDim() const { 
    if (!_valid)
      throw std::runtime_error("invalid input object");
 
    return 0;
  }

  /*! Get the list of vertex IDs and their weights.
   */
  void getVertexListCopy(std::vector<int> &ids, 
    std::vector<int> &localIDs, std::vector<Scalar> &xyz,
    std::vector<Scalar> &wgt) const
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    xyz.clear();

    int numVtx = this->getLocalNumVertices();
    int nweights = _vtxWeightDim * numVtx;
    // TODO -test for memory alloc failure when we resize
    wgt.resize(nweights); 

    if (nweights){
      Scalar *wTo = &wgt[0];
      Scalar *wFrom = _vtxWgt.getRawPtr();
      memcpy(wTo, wFrom, sizeof(Scalar) * nweights);
    }

    ids.resize(numVtx);
    localIDs.resize(numVtx);
    for (unsigned i=0; i < _rowMap->getNodeNumElements(); i++){
      ids[i] =  _rowMap->getGlobalElement(i);
      localIDs[i] = i;
    }
  }

  /*! Return a read only view of the data.
   */
  int getVertexListView(int const *&ids,
    std::vector<int> &localIDs, Scalar const *& xyz, Scalar const *&wgts)
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");
    
    int nvtx = this->getLocalNumVertices();
    ids = _rowMap->getNodeElementList().getRawPtr();
    xyz = NULL;
    wgts = _vtxWgt.getRawPtr();
    localIDs.resize(nvtx);
    for (int i=0; i < nvtx; i++)
      localIDs[i] = i;
    return nvtx;
  }

  /*! Return a copy of the edge IDs and edge weights for a vertex.
   */
  void getVertexEdgeCopy(int vtxId, int localId, 
    std::vector<int> &edgeId, std::vector<Scalar> &wgts) const
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");
    
    int nvtx = this->getLocalNumVertices();

    if (localId < 0 || localId >= nvtx)
      throw std::runtime_error("invalid local vertex ID");

    edgeId.clear();
    wgts.clear();

    Teuchos::ArrayView<const int> nbors;
    _graph->getLocalRowView(localId, nbors);
    size_t nedges = nbors.size();

    if (nedges > 0){
      const int *fromId = nbors.getRawPtr();
      edgeId.resize(nedges);
      int *toId = &edgeId[0];
      memcpy(toId, fromId, sizeof(int) * nedges);

      if (_edgeWeightDim > 0){
        int offset = _edgeOffsets[localId] * _edgeWeightDim;
        Scalar *fromWgt = &_edgeWgt[offset];
        wgts.resize(_edgeWeightDim);
        Scalar *toWgt = &wgts[0];
        memcpy(toWgt, fromWgt, sizeof(Scalar) * _edgeWeightDim);
      }
      else{
        wgts.clear();
      }
    }
  }

  /*! Return pointers to the edge IDs and edge weights for a vertex.
   */
  int getVertexEdgeView(int vtxId, int localId, 
    int const *&edgeId, Scalar const *&wgts) const
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");
    
    int nvtx = this->getLocalNumVertices();

    if (localId < 0 || localId >= nvtx)
      throw std::runtime_error("invalid local vertex ID");

    edgeId = NULL;
    wgts = NULL;

    Teuchos::ArrayView<const int> nbors;
    _graph->getLocalRowView(localId, nbors);
    size_t nedges = nbors.size();
    edgeId = &nbors[0];

    if ((nedges > 0) && (_edgeWeightDim > 0)){
      int offset = _edgeOffsets[localId] * _edgeWeightDim;
      wgts = &_edgeWgt[offset];
    }
    return nedges;
  }
};
  
}  //namespace Zoltan2
  
#endif
