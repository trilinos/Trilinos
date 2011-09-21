// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsGraphInput.hpp

    \brief An input adapter for a Epetra_CrsGraph.
*/

#ifndef _ZOLTAN2_EPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_EPETRACRSGRAPHINPUT_HPP_

#include <Zoltan2_GraphInput.hpp>
#include <Epetra_CrsGraph.h>
#include <Xpetra_EpetraCrsGraph.hpp>
#include <Teuchos_RCP.hpp>

namespace Zoltan2 {

/*! Zoltan2::EpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Epetra_CrsGraph data and weights.

    EpetraCrsGraphInput provides access for Zoltan2 to an Epetra_CrsGraph and
    vertex and edges weights.  Most of the methods used by Zoltan2 are in the 
    XpetraGraphInput parent class.  The methods used by the Zoltan2 caller 
    to initialize the objects are in this class.

    The template parameter is the weight type.  Epetra local and global IDs 
    are ints.
*/

template <typename Scalar>
class EpetraCrsGraphInput : public GraphInput<Scalar, int, int>{
private:

  bool _valid;
  Teuchos::RCP<Xpetra::EpetraCrsGraph> _graph;
  Teuchos::RCP<const Xpetra::Map<int, int> > _rowMap;
  Teuchos::ArrayRCP<Scalar> _vtxWgt;
  Teuchos::ArrayRCP<Scalar> _edgeWgt;
  int _vtxWeightDim;
  int _edgeWeightDim;
  std::vector<int> _edgeOffsets; 

public:

  std::string inputAdapterName()const {return std::string("EpetraCrsGraph");}

  //~EpetraCrsGraphInput() { }

  /*! Default constructor
   */
  EpetraCrsGraphInput():_valid(false), _graph(), _rowMap(),
    _vtxWgt(), _edgeWgt(),
    _vtxWeightDim(0), _edgeWeightDim(0), _edgeOffsets() {}

  /*! Constructor with graph only
      \param graph  the graph represented by this input adapter
   */
  EpetraCrsGraphInput(Teuchos::RCP<Epetra_CrsGraph> graph):
    _valid(false), _vtxWgt(), _edgeWgt(), _vtxWeightDim(0), _edgeWeightDim(0)
  {
    _graph = Teuchos::rcp(new Xpetra::EpetraCrsGraph(graph));
    _rowMap = _graph->getRowMap();
    int numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    _valid=true;
  }

  /*! Constructor with weights
      \param graph  the graph represented by this input adapter
      \param vertexWeights are given in vertex local number order
      \param edgeWeights when queried, the Epetra_CrsGraph gives edges in a 
               certain order.  The edgeWeights follow this order, omitting
               any self edges.
   */
  EpetraCrsGraphInput(Teuchos::RCP<Epetra_CrsGraph> graph,
                      Teuchos::ArrayRCP<Scalar> vertexWeights,
                      Teuchos::ArrayRCP<Scalar> edgeWeights):_valid(false)
  { 
    _graph = Teuchos::rcp(new Xpetra::EpetraCrsGraph(graph));
    _rowMap = _graph->getRowMap();
    int numV = _rowMap->getNodeNumElements();
    _edgeOffsets.resize(numV+1, 0);
    for (int i=0; i < numV; i++){
      _edgeOffsets[i+1] = _edgeOffsets[i] + _graph->getNumEntriesInLocalRow(i);
    }
    int numVwgts = vertexWeights.size();
    int numEwgts = edgeWeights.size();

    if (numVwgts > 0){
      _vtxWeightDim = numVwgts / numV;
      if (_vtxWeightDim * numV != numVwgts)
        throw std::runtime_error("bad number of vertex weights");
    }

    if (numEwgts > 0){
      int numE = _graph->getNodeNumEntries();
      _edgeWeightDim = numEwgts / numE;
      if (_edgeWeightDim * numE != numEwgts)
        throw std::runtime_error("bad number of edge weights");
    }
    _valid=true;
  }

  bool validInput() { return _valid;}

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
