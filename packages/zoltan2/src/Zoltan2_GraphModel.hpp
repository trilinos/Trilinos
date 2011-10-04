// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphModel.hpp

    \brief The abstract interface for a graph model.
*/


#ifndef _ZOLTAN2_GRAPHMODEL_HPP_
#define _ZOLTAN2_GRAPHMODEL_HPP_

#include <vector>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Hashtable.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_IdentifierMap.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphModel
    \brief GraphModel defines the interface required for
            graph models.  

    GraphModels are templated on an input adapter.

    Your concrete implementation can require that all processes must
    call the constructor.  The remaining methods must be able
    to be called asynchronously.
*/
template <
  template <typename, typename, typename, typename, typename, typename> 
  class AdapterType,
  typename Scalar, typename LNO, typename GNO, typename LID, typename GID, 
  typename Node
>
    struct GraphModel {
private:

public:

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return 0; }

  /*! Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const { return 0; }

  /*! Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const { return 0; }

  /*! Returns the global number edges.
   */
  global_size_t getGlobalNumEdges() const { return 0; }

  /*! Returns the dimension (0 or greater) of vertex weights.
   */
  int getVertexWeightDim() const { return 0; }

  /*! Returns the dimension (0 or greater) of edge weights.
   */
  int getEdgeWeightDim() const { return 0; }

  /*! Returns the dimension (0 to 3) of vertex coordinates.
   */
  int getCoordinateDim() const { return 0; }

  /*! Sets pointers to this process' vertex Ids and their weights.
      \param Ids will on return point to the list of the global Ids for 
        each vertex on this process.
      \param xyz will on return point to a list coordinates for
         each vertex in the Ids list.  Coordinates are listed by 
         vertex by component.  
      \param wgts will on return point to a list of the weight or weights 
         associated with each vertex in the Ids list.  Weights are listed by 
         vertex by weight component.  
       \return The number of ids in the Ids list.
   */

  size_t getVertexList( ArrayView<const GNO> &Ids,
    ArrayView<const Scalar> &xyz, ArrayView<const Scalar> &wgts) const { return 0; }

  /*! Obtain a read-only view of the edge Ids of the input vertex.
      \param Id  is the global Id for a vertex on this process.
      \param edgeId on return will point to the list of edge neighbors.
      \param procId on return holds the list of each process owning the 
        corresponding neighbor vertex.
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.
   */
  size_t getVertexGlobalEdge( GNO Id, 
    ArrayView<const GNO> &edgeId, ArrayView<const int> &procId, 
    ArrayView<const Scalar> *&wgts) const { return 0; }

  /*! Obtain a read-only view of the edge Ids of the input vertex.
      \param localRef  is the local id associated with vertex.  Local ids
        are consecutive, begin at 0 and follow the order of vertices returned
        by getVertexList.
      \param edgeId on return will point to the list of edge neighbor global 
         Ids.
      \param procId on return holds the list of each process owning the 
        corresponding neighbor vertex.
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.
   */
  size_t getVertexLocalEdge( LNO localRef, 
    ArrayView<const GNO> &edgeId, ArrayView<const int> &procId, 
    ArrayView<const Scalar> &wgts) const { return 0; }
};

////////////////////////////////////////////////////////////////
// Graph model derived from matrix input.
////////////////////////////////////////////////////////////////

/*! Zoltan2::GraphModel<MatrixInput>
    \brief A (partial) specialization of GraphModel 
           for a Zoltan2::MatrixInput object.
          
         TODO test for memory allocation failure when allocating arrays
     TODO - default template arguments are not allowed in partial specializations
*/

template <typename Adapter, Z2FN_TEMPLATE>
struct GraphModel<MatrixInput, Z2PARAM_TEMPLATE>
{

private:

  RCP<const MatrixInput<Scalar, LNO, GNO, LID, GID, Node> > _input;
  RCP<const Teuchos::Comm<int> > _comm;
  RCP<const Environment > _env;

  ArrayRCP<GID> _gids;
  ArrayRCP<GNO> _gnos;
  ArrayRCP<LNO> _lids;
  ArrayRCP<LNO> _nedges;
  ArrayRCP<GID> _edgeGids;
  ArrayRCP<GNO> _edgeGnos;
  ArrayRCP<int> _procIds;
  Array<LNO> _offsets;

  RCP<Teuchos::Hashtable<GNO, LNO> > _gnoToLno;

  // Transpose is only required if vertices are columns.
  RCP<Tpetra::CrsMatrix< Scalar, LNO, GNO, Node> > _inputTranspose; 

  RCP<IdentifierMap< LID, GID, LNO, GNO> > _idMap; 

  global_size_t _numGlobalEdges;
  bool _gnosAreGids;

  void makeOffsets()
  {
    if (_offsets.size() > 0)
      return;
    _offsets = Array<LNO>(_gid.size() + 1);
    _offsets[0] = 0;

    for (size_t i=1; i <= _gid.size(); i++)
      _offsets[i] = _offsets[i-1] + _nedges[i-1];
  }

  void makeGnoToLno()
  {
    if (!_gnoToLno.is_null())
      return;

    _gnoToLno = rcp(new Teuchos::Hashtable<GNO, LNO>(_gid.size()+1));

    if (_gid.size() == 0)
      return;

    if (_gnosAreGids)
      for (size_t i=0; i < _gid.size(); i++)
        _gnoToLno.put(_gid[i], i);
    else
      for (size_t i=0; i < _gid.size(); i++)
        _gnoToLno.put(_gno[i], i);
  }

public:

  /*! Constructor
   *  All processes in the communicator must call the constructor.
   */
  GraphModel(
    RCP<const MatrixInput<Scalar, LNO, GNO, LID, GNO, Node> > inputAdapter,
    RCP<const Comm<int> > comm, RCP<const Environment <int> > env) : 
      _input(inputAdapter), _comm(comm), _env(env),
      _gids(), _gnos(), _lids(), _nedges(), _edgeGids(), _edgeGnos(),
      _procIds(), _offsets(), _gnoToLno(),
      _inputTranspose(), _idMap(), _numGlobalEdges(0), _gnosAreGids(false)
  {
    // TODO: interpretation of matrix should be given in constructor,
    //   not in later set methods. For now we assume vertices are matrix rows,
    //   rows are uniquely owned, and the matrix is symmetric.

    // Get graph from matrix input

    try{
      typedef std::vector<GID> gid_vector_t;
      typedef std::vector<LID> lid_vector_t;
      typedef std::vector<LNO> lno_vector_t;

      RCP<gid_vector_t> gids = rcp(new std::vector<gid_vector_t>);
      RCP<lid_vector_t> lids = rcp(new std::vector<lid_vector_t>);
      RCP<lno_vector_t> numEdges = rcp(new std::vector<lno_vector_t>);
      RCP<gid_vector_t> edgeGids = rcp(new std::vector<gid_vector_t>);

      _input->getRowListCopy(gids.get(), lids.get(), numEdges.get(), 
        edgeGids.get());
    }
    catch (std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(_env, e);

    _gids = arcp(gids);
    _lids = arcp(lids);
    _nEdges = arcp(numEdges);
    _edgeIds = arcp(edgeGids);

    // Translate user's global Ids if necessary

    _idMap = rcp(new IdentifierMap<LID, GID, LNO, GNO>
        (_comm, _env, _gids, _lids));

    _gnosAreGids = _idMap->gnosAreGids();

    makeOffsets();    // offsets into edge Ids

    global_size_t numEdges = _offsets[_gids.size()];
    Teuchos::reduceAll<size_t>(*comm, Teuchos::REDUCE_SUM, 1,
      &numEdges, &_numGlobalEdges);

    _procIds = arcp<int>(numEdges);

    if (!_gnosAreGids){
      _gnos = arcp<GNO>(_gids.size());
      _idMap.gidTranslate(_gids, _gnos, TRANSLATE_GID_TO_GNO);

      _edgeGnos = arcp<GNO>(numEdges);
    }
    else{
      _gnos = arcp<GNO>(0);
      _edgeGnos = arcp<GNO>(0);
    }

    // TODO is there a short cut for creating a view that
    //    is actually the whole array.

    _idMap.gidGlobalTranslate(_edgeIds.view(0,numEdges),
       _edgeGnos.view(0, _edgeGnos.size()),
       _procIds.view(0, numEdges));

    // process owning each edge Id (assuming vertices are rows,
    //   and rows are uniquely owned, and matrix is symmetric)
    //   TODO which assumptions must be changed?
  }

  RCP<const MatrixInput> getMatrixInput() { return _input;}

  // // // // // // // // // // // // // // // // // // // // // /
  // Configuration methods for a graph derived from matrix input.
  // // // // // // // // // // // // // // // // // // // // // /

#if 0
  void setVerticesAreRows() 
  {
    _verticesAreRows = true;
    _verticesAreCols = false;
  }

  void setVerticesAreCols() 
  {
    throw std::runtime("setVerticesAreCols not yet implemented");
    _verticesAreRows = false;
    _verticesAreCols = true;
    // TODO build the transpose
  }

  void setVertexWeightsAreUnit()
  {
    _vertexWeightsAreUnit = true;
    _vertexWeightsAreNumEdges = false;
  }

  void setVertexWeightsAreNumEdges()
  {
    _vertexWeightsAreUnit = false;
    _vertexWeightsAreNumEdges = true;
  }

  bool getVerticesAreRows() { return _verticesAreRows; } 

  bool getVerticesAreCols() { return _verticesAreCols; } 

  bool getVertexWeightsAreUnit() { return _vertexWeightsAreUnit; }

  bool getVertexWeightsAreNumEdges() { return _vertexWeightsAreNumEdges; }
#endif

  // // // // // // // // // // // // // // // // // // // // // /
  // The GraphModel interface.
  // // // // // // // // // // // // // // // // // // // // // /

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const
  {
    return _input->getLocalNumRows();
  }

  /*! Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const
  {
    return _input->getGlobalNumRows();
  }

  /*! Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const 
  {
    return static_cast<size_t>(_offsets[_gids.size()]);
  }

  /*! Returns the global number edges.
   */
  global_size_t getGlobalNumEdges() const
  {
    return _numGlobalEdges;
  }

  /*! Returns the dimension (0 or greater) of vertex weights.
   */
  int getVertexWeightDim() const
  {
    return 0;   // TODO
  }

  /*! Returns the dimension (0 or greater) of edge weights.
   */
  int getEdgeWeightDim() const
  {
    return 0;   // TODO
  }

  /*! Returns the dimension (0 to 3) of vertex coordinates.
   */
  int getCoordinateDim() const
  {
    return 0;   // TODO
  }

  /*! Sets pointers to this process' vertex Ids and their weights.
      \param Ids will on return point to the list of the global Ids for 
        each vertex on this process.
      \param xyz will on return point to a list coordinates for
         each vertex in the Ids list.  Coordinates are listed by 
         vertex by component.  
      \param wgts will on return point to a list of the weight or weights 
         associated with each vertex in the Ids list.  Weights are listed by 
         vertex by weight component.  
       \return The number of ids in the Ids list.
      
      TODO - do we want to pass around ArrayRCP instead? 
   */

  size_t getVertexList( ArrayView<const GNO> &Ids,
    ArrayView<const Scalar> &xyz, ArrayView<const Scalar> &wgts) const
  {
    if (_gnosAreGids){
      Ids = _gids.view(0, _gids.size());
    }
    else{
      Ids = _gnos.view(0, _gnos.size());
    }
    // TODO coordinates and weights
    return _gids.size();
  }

  /*! Obtain a read-only view of the edge Ids of the input vertex.
      \param Id  is the global Id for a vertex on this process.
      \param edgeId on return will point to the list of edge neighbors.
      \param procId on return holds the list of each process owning the 
        corresponding neighbor vertex.
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.

      It is more efficient to query with getVertexLocalEdge.
   */
  size_t getVertexGlobalEdge( GNO Id, ArrayView<const GNO> &edgeId, 
    ArrayView<const int> &procId, ArrayView<const Scalar> *&wgts) const
  {
    LNO lno(0);

    try {
      lno = _gnoToLno(Id);
    }
    catch (std::exception &e){
      Z2_THROW_OUTSIDE_ERROR(_env, e);
    }

    if (_gnosAreGids){
      edgeId = _edgeIds.view(_offsets[lno], _offsets[lno+1] - _offsets[lno]);
    }
    else{
      edgeId = _edgeGnos.view(_offsets[lno], _offsets[lno+1] - _offsets[lno]);
    }
    // TODO edge weights
  }

  /*! Obtain a read-only view of the edge Ids of the input vertex.
      \param lno is the local id associated with local vertex. Local ids
         are consecutive beginning with 0, and ordered corresponding to the
         global Ids returned in getVertexList.
      \param edgeId on return will point to the list of edge neighbor global 
         Ids.
      \param procId on return holds the list of each process owning the 
        corresponding neighbor vertex.
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.
   */
  size_t getVertexLocalEdge( LNO lno, ArrayView<const GNO> &edgeId, 
    ArrayView<const int> &procId, ArrayView<const Scalar> *&wgts) const
  {
    Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, "invalid local id",
      localRef >= 0 && localRef < gids.size(), Z2_BASIC_ASSERTION);

    size_t firstId =  _offsets[lno];
    size_t nEdges = _offsets[lno+1] - _offsets[lno];

    if (_gnosAreGids){
      edgeId = _edgeIds.view(firstIdx, nEdges);
    }
    else{
      edgeId = _edgeGnos.view(firstIdx, nEdges);
    }

    procId = _procIds.view(firstIdx, nEdges);
    // TODO edge weights

    return nEdges;
  }
};

#endif
