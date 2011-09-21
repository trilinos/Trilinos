// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsGraphInput.hpp

    \brief An input adapter for a Epetra::CrsGraph.
*/

#ifndef _ZOLTAN2_TPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_TPETRACRSGRAPHINPUT_HPP_

#include <Xpetra_TpetraCrsGraph.hpp>
#include <Zoltan2_GraphInput.hpp>
#include <Teuchos_RCP.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Tpetra::CrsGraph data plus weights.

    TpetraCrsGraphInput provides access for Zoltan2 to a
    Tpetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 are in the XpetraGraphInput parent class.  The methods used
    by the Zoltan2 caller are in this class.

    The template parameter is the weight type.

    TODO: include NODE
*/

template <typename Scalar, typename LNO, typename GNO>
class TpetraCrsGraphInput : public GraphInput<Scalar, LNO, GNO>{
private:

  bool _valid;
  Teuchos::RCP<Xpetra::TpetraCrsGraph<LNO, GNO> > _graph;
  Teuchos::RCP<const Xpetra::Map<LNO, GNO> > _rowMap;
  Teuchos::ArrayRCP<Scalar> _vtxWgt;
  Teuchos::ArrayRCP<Scalar> _edgeWgt;
  int _vtxWeightDim;
  int _edgeWeightDim;
  std::vector<LNO> _edgeOffsets;

public:
  std::string inputAdapterName()const {return std::string("TpetraCrsGraph");}

  //~TpetraCrsGraphInput() { }

  /*! Default constructor 
   */
  TpetraCrsGraphInput():_valid(false), _graph(), _rowMap(),
    _vtxWgt(), _edgeWgt(),
    _vtxWeightDim(0), _edgeWeightDim(0), _edgeOffsets() {}

  /*! Constructor with graph only
   */
  TpetraCrsGraphInput(Teuchos::RCP<Tpetra::CrsGraph<LNO, GNO> > graph):
   _valid(false), _vtxWgt(), _edgeWgt(), _vtxWeightDim(0), _edgeWeightDim(0)
  {
    _graph = Teuchos::rcp(new Xpetra::TpetraCrsGraph<LNO, GNO>(graph));
    _rowMap = _graph->getRowMap();
    LNO numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (LNO i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    _valid=true;
  }

  /*! Constructor with weights
   */
  TpetraCrsGraphInput(Teuchos::RCP<Tpetra::CrsGraph<LNO, GNO> > graph,
                      Teuchos::ArrayRCP<Scalar> vertexWeights,
                      Teuchos::ArrayRCP<Scalar> edgeWeights)_valid(false)
  {
    _graph = Teuchos::rcp(new Xpetra::TpetraCrsGraph<LNO, GNO>(graph));
    _rowMap = _graph->getRowMap();
    LNO numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (LNO i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    LNO numVwgts = vertexWeights.size();
    LNO numEwgts = edgeWeights.size();

    if (numVwgts > 0){
      _vtxWeightDim = numVwgts / numV;
      if (_vtxWeightDim * numV != numVwgts)
        throw std::runtime_error("bad number of vertex weights");
    }

    if (numEwgts > 0){
      LNO numE = _graph->getNodeNumEntries();
      _edgeWeightDim = numEwgts / numE;
      if (_edgeWeightDim * numE != numEwgts)
        throw std::runtime_error("bad number of edge weights");
    }
    _valid=true;
  }
  bool validInput() { return _valid;}

  /*! Returns the number vertices on this process.
   */
  LNO getLocalNumVertices() const {
    if (!_valid)
      throw std::runtime_error("invalid input object");
    return _rowMap->getNodeNumElements();
  }

  /*! Returns the number edges on this process.
   */
  LNO getLocalNumEdges() const {
    if (!_valid)
      throw std::runtime_error("invalid input object");
    return _graph->getNodeNumEntries();
  }

  /*! Returns the number weights supplied for each vertex.
   */
  LNO getVertexWeightDim() const {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    return _vtxWeightDim;
  }

  /*! Returns the number weights supplied for each edge.
   */
  LNO getEdgeWeightDim() const {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    return _edgeWeightDim;
  }

  /*! Returns the number of coordinates per vertex
   */
  LNO getCoordinateDim() const {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    return 0;
  }

  /*! Get the list of vertex IDs and their weights.
   */
  void getVertexListCopy(std::vector<GNO> &ids,
    std::vector<LNO> &localIDs, std::vector<Scalar> &xyz,
    std::vector<Scalar> &wgt) const
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    xyz.clear();

    LNO numVtx = this->getLocalNumVertices();
    LNO nweights = _vtxWeightDim * numVtx;
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
  LNO getVertexListView(GNO const *&ids,
    std::vector<LNO> &localIDs, Scalar const *& xyz, Scalar const *&wgts)
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    LNO nvtx = this->getLocalNumVertices();
    ids = _rowMap->getNodeElementList().getRawPtr();
    xyz = NULL;
    wgts = _vtxWgt.getRawPtr();
    localIDs.resize(nvtx);
    for (LNO i=0; i < nvtx; i++)
      localIDs[i] = i;
    return nvtx;
  }

  /*! Return a copy of the edge IDs and edge weights for a vertex.
   */
  void getVertexEdgeCopy(GNO vtxId, LNO localId,
    std::vector<GNO> &edgeId, std::vector<Scalar> &wgts) const
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");
   
    LNO nvtx = this->getLocalNumVertices();
  
    if (localId < 0 || localId >= nvtx)
      throw std::runtime_error("invalid local vertex ID");
    
    edgeId.clear();
    wgts.clear();
    
    Teuchos::ArrayView<const GNO> nbors;
    _graph->getLocalRowView(localId, nbors);
    size_t nedges = nbors.size();

    if (nedges > 0){
      const GNO *fromId = nbors.getRawPtr();
      edgeId.resize(nedges); 
      GNO *toId = &edgeId[0];
      memcpy(toId, fromId, sizeof(GNO) * nedges);

      if (_edgeWeightDim > 0){
        LNO offset = _edgeOffsets[localId] * _edgeWeightDim;
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
  int getVertexEdgeView(GNO vtxId, LNO localId,
    GNO const *&edgeId, Scalar const *&wgts) const
  {
    if (!_valid)
      throw std::runtime_error("invalid input object");

    LNO nvtx = this->getLocalNumVertices();
    
    if (localId < 0 || localId >= nvtx)
      throw std::runtime_error("invalid local vertex ID");
    
    edgeId = NULL;
    wgts = NULL;
    
    Teuchos::ArrayView<const GNO> nbors;
    _graph->getLocalRowView(localId, nbors);
    size_t nedges = nbors.size();
    edgeId = &nbors[0];

    if ((nedges > 0) && (_edgeWeightDim > 0)){
      LNO offset = _edgeOffsets[localId] * _edgeWeightDim;
      wgts = &_edgeWgt[offset];
    }
    return nedges;
  }
};
  
  
}  //namespace Zoltan2
  
#endif
