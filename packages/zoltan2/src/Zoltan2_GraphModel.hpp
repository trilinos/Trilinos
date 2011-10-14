// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphModel.hpp

    \brief The interface and implementations of a graph model.
*/


#ifndef _ZOLTAN2_GRAPHMODEL_HPP_
#define _ZOLTAN2_GRAPHMODEL_HPP_

#include <vector>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Hashtable.hpp>
#include <Zoltan2_Model.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_IdentifierMap.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphModel
    \brief GraphModel defines the interface required for graph models.  

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    method should be local methods.
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
private:

public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef typename Adapter::gid_t     gid_t;
  typedef typename Adapter::lid_t     lid_t;
  typedef typename Adapter::node_t    node_t;
  typedef typename Adapter::user_t    user_t;

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

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const {
      return 0; }

  /*! Sets pointers to this process' edge (neighbor) global Ids.
      \param edgeIds This is the list of global neighbor Ids corresponding
        to the vertices listed in getVertexList.
      \param procIds lists the process owning each neighbor in the edgeIds
         list.
      \param offsets offsets[i] is the offset into edgeIds to the start
        of neighbors for ith vertex.
      \param wgts will on return point to a list of the weight or weights
         associated with each edge in the edgeIds list.  Weights are listed by
         edge by weight component.
       \return The number of ids in the edgeIds list.
   */

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const { return 0; }

  /*! Obtain a view of the edge Ids of the input vertex.
      \param Id  is the global Id for a vertex on this process.
      \param edgeIds on return will point to the list of edge neighbors.
      \param procIds on return holds the list of each process owning the
        corresponding neighbor vertex. 
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.
  
      This method is defined for convenience when obtaining the
      neighbors of a vertex.  It is not efficient to call this method
      many times in a loop, due to the construction and destruction of
      ArrayViews.  Call getEdgeList instead.
   */
  size_t getVertexGlobalEdge( gno_t Id, 
    ArrayView<const gno_t> &edgeIds, ArrayView<const int> &procIds,
    ArrayView<const scalar_t> *&wgts) const { return 0; }
   
  /*! Obtain a view of the edge Ids of the input vertex.
      \param localRef  is the local id associated with vertex.  Local ids
        are consecutive, begin at 0 and follow the order of vertices returned
        by getVertexList.
      \param edgeIds on return will point to the list of edge neighbor global
         Ids.
      \param procIds on return holds the list of each process owning the
        corresponding neighbor vertex.
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.

      This method is defined for convenience when obtaining the 
      neighbors of a vertex.  It is not efficient to call this method 
      many times in a loop, due to the construction and destruction of
      ArrayViews.  Call getEdgeList instead.
   */
  size_t getVertexLocalEdge( lno_t localRef,
    ArrayView<const gno_t> &edgeIds, ArrayView<const int> &procIds,
    ArrayView<const scalar_t> &wgts) const { return 0; }
};

////////////////////////////////////////////////////////////////
// Graph model derived from XpetraCrsMatrixInput.
//    We know that Xpetra input does not need an IdentifierMap
////////////////////////////////////////////////////////////////

/*! Zoltan2::GraphModel<XpetraCrsMatrixInput>
    \brief A (partial) specialization of GraphModel
           for a Zoltan2::XpetraCrsMatrixInput object.
*/

template <typename User>
class GraphModel<XpetraCrsMatrixInput<User> >
{
public:

  typedef typename XpetraCrsMatrixInput<User>::scalar_t  scalar_t;
  typedef typename XpetraCrsMatrixInput<User>::gno_t     gno_t;
  typedef typename XpetraCrsMatrixInput<User>::lno_t     lno_t;
  typedef typename XpetraCrsMatrixInput<User>::gid_t     gid_t;
  typedef typename XpetraCrsMatrixInput<User>::lid_t     lid_t;
  typedef typename XpetraCrsMatrixInput<User>::node_t    node_t;

  /*! Constructor
   *  All processes in the communicator must call the constructor.
   */
  GraphModel(
    const RCP<const XpetraCrsMatrixInput<User> > &inputAdapter,
    const RCP<const Comm<int> > &comm, const RCP<const Environment> &env) :
      input_(inputAdapter), comm_(comm), env_(env),
      gnos_(), edgeGnos_(), procIds_(), offsets_(),
      numLocalEdges_(), numGlobalEdges_(0)
  {
    gno_t const *vtxIds=NULL, *nborIds=NULL;
    lno_t const  *offsets=NULL, *lids=NULL; 
    lno_t numVtx;
    try{
      numVtx = input_->getRowListView(vtxIds, lids, offsets, nborIds);
    }
    catch (std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(env_, e);

    gnos_ = arcp(vtxIds, 0, numVtx, false);   // non-owning ArrayRCPs
    offsets_ = arcp(offsets, 0, numVtx+1, false);

    numLocalEdges_ = offsets_[numVtx];

    edgeGnos_ = arcp(nborIds, 0, numLocalEdges_, false);

    Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1,
      &numLocalEdges_, &numGlobalEdges_);

    procIds_ = arcp<int>(numLocalEdges_);

    RCP<const Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> > xmatrix =
      input_->getMatrix();

    xmatrix->getRowMap()->getRemoteIndexList(edgeGnos_(), procIds_());
  }

  // // // // // // // // // // // // // // // // // // // // // /
  // The GraphModel interface.
  // // // // // // // // // // // // // // // // // // // // // /

  size_t getLocalNumVertices() const
  {
    return input_->getLocalNumRows();
  }

  global_size_t getGlobalNumVertices() const
  {
    return input_->getGlobalNumRows();
  }

  size_t getLocalNumEdges() const
  {
    return numLocalEdges_;
  }
   
  global_size_t getGlobalNumEdges() const
  { 
    return numGlobalEdges_;
  } 

  int getVertexWeightDim() const
  {
    return 0;   // TODO
  } 

  int getEdgeWeightDim() const
  {
    return 0;   // TODO
  } 
    
  int getCoordinateDim() const
  {
    return 0;   // TODO
  } 

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const
  {
    Ids = gnos_();    // () operator - an ArrayView
    return gnos_.size();
  }

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds, 
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const
  {
    edgeIds = edgeGnos_();
    procIds = procIds_();
    offsets = offsets_();

    return numLocalEdges_;
  }

  size_t getVertexGlobalEdge( gno_t Id, ArrayView<const gno_t> &edgeId,
    ArrayView<const int> &procId, ArrayView<const scalar_t> *&wgts) const
  {
    lno_t lno(0);
    // TODO map lno to gno
    throw std::runtime_error("not implemented");
    return getVertexLocalEdge(lno, edgeId, procId, wgts);
  }

  size_t getVertexLocalEdge( lno_t lno, ArrayView<const gno_t> &edgeId,
    ArrayView<const int> &procId, ArrayView<const scalar_t> *&wgts) const
  { 
    Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, "invalid local id",
      lno >= 0 && lno < gnos_.size(), BASIC_ASSERTION);

    lno_t thisVtx =  offsets_[lno];
    lno_t nextVtx = (lno < gnos_.size()-1) ? offsets_[lno+1] : numLocalEdges_;
    size_t nEdges = nextVtx - thisVtx;

    edgeId = edgeGnos_.view(thisVtx, nEdges);
    procId = procIds_.view(thisVtx, nEdges);
    return nEdges;
  }

private:

  RCP<const XpetraCrsMatrixInput<User> > input_;
  RCP<const Teuchos::Comm<int> > comm_;
  RCP<const Environment > env_;

  ArrayRCP<const gno_t> gnos_;
  ArrayRCP<const gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  // Transpose is required only if vertices are columns.
  // KDDKDD ??  We won't form an actual transpose, will we?
  // KDDKDD RCP<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> > inputTranspose;

  global_size_t numLocalEdges_;
  global_size_t numGlobalEdges_;

};


#if 0
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

  RCP<const MatrixInput<Scalar, LNO, GNO, LID, GID, Node> > input_;
  RCP<const Teuchos::Comm<int> > comm_;
  RCP<const Environment > env_;

  ArrayRCP<GID> _gids;
  ArrayRCP<GNO> gnos_;
  ArrayRCP<LNO> _lids;
  ArrayRCP<LNO> _nedges;
  ArrayRCP<GID> _edgeGids;
  ArrayRCP<GNO> edgeGnos_;
  ArrayRCP<int> procIds_;
  Array<LNO> offsets_;

  RCP<Teuchos::Hashtable<GNO, LNO> > _gnoToLno;

  // Transpose is only required if vertices are columns.
  RCP<Tpetra::CrsMatrix< Scalar, LNO, GNO, Node> > input_Transpose; 

  RCP<IdentifierMap< LID, GID, LNO, GNO> > _idMap; 

  global_size_t numGlobalEdges_;
  bool gnos_AreGids;

  void makeOffsets()
  {
    if (offsets_.size() > 0)
      return;
    offsets_ = Array<LNO>(_gid.size() + 1);
    offsets_[0] = 0;

    for (size_t i=1; i <= _gid.size(); i++)
      offsets_[i] = offsets_[i-1] + _nedges[i-1];
  }

  void makeGnoToLno()
  {
    if (!_gnoToLno.is_null())
      return;

    _gnoToLno = rcp(new Teuchos::Hashtable<GNO, LNO>(_gid.size()+1));

    if (_gid.size() == 0)
      return;

    if (gnos_AreGids)
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
    const RCP<const MatrixInput<Scalar, LNO, GNO, LID, GNO, Node> > &inputAdapter,
    const RCP<const Comm<int> > &comm, const RCP<const Environment <int> > &env) : 
      input_(inputAdapter), comm_(comm), env_(env),
      _gids(), gnos_(), _lids(), _nedges(), _edgeGids(), edgeGnos_(),
      procIds_(), offsets_(), _gnoToLno(),
      input_Transpose(), _idMap(), numGlobalEdges_(0), gnos_AreGids(false)
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

      input_->getRowListCopy(gids.get(), lids.get(), numEdges.get(), 
        edgeGids.get());
    }
    catch (std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(env_, e);

    _gids = arcp(gids);
    _lids = arcp(lids);
    _nEdges = arcp(numEdges);
    _edgeIds = arcp(edgeGids);

    // Translate user's global Ids if necessary

    _idMap = rcp(new IdentifierMap<LID, GID, LNO, GNO>
        (comm_, env_, _gids, _lids));

    gnos_AreGids = _idMap->gnosAreGids();

    makeOffsets();    // offsets into edge Ids

    global_size_t numEdges = offsets_[_gids.size()];
    Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1,
      &numEdges, &numGlobalEdges_);

    procIds_ = arcp<int>(numEdges);

    if (!gnos_AreGids){
      gnos_ = arcp<GNO>(_gids.size());
      _idMap.gidTranslate(_gids, gnos_, TRANSLATE_GID_TO_GNO);

      edgeGnos_ = arcp<GNO>(numEdges);
    }
    else{
      gnos_ = arcp<GNO>(0);
      edgeGnos_ = arcp<GNO>(0);
    }

    // TODO is there a short cut for creating a view that
    //    is actually the whole array.

    _idMap.gidGlobalTranslate(_edgeIds.view(0,numEdges),
       edgeGnos_.view(0, edgeGnos_.size()),
       procIds_.view(0, numEdges));

    // process owning each edge Id (assuming vertices are rows,
    //   and rows are uniquely owned, and matrix is symmetric)
    //   TODO which assumptions must be changed?
  }

  RCP<const MatrixInput> getMatrixInput() { return input_;}

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
    return input_->getLocalNumRows();
  }

  /*! Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const
  {
    return input_->getGlobalNumRows();
  }

  /*! Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const 
  {
    return static_cast<size_t>(offsets_[_gids.size()]);
  }

  /*! Returns the global number edges.
   */
  global_size_t getGlobalNumEdges() const
  {
    return numGlobalEdges_;
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
    if (gnos_AreGids){
      Ids = _gids.view(0, _gids.size());
    }
    else{
      Ids = gnos_.view(0, gnos_.size());
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
      Z2_THROW_OUTSIDE_ERROR(env_, e);
    }

    if (gnos_AreGids){
      edgeId = _edgeIds.view(offsets_[lno], offsets_[lno+1] - offsets_[lno]);
    }
    else{
      edgeId = edgeGnos_.view(offsets_[lno], offsets_[lno+1] - offsets_[lno]);
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
    Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, "invalid local id",
      localRef >= 0 && localRef < gids.size(), BASIC_ASSERTION);

    size_t firstId =  offsets_[lno];
    size_t nEdges = offsets_[lno+1] - offsets_[lno];

    if (gnos_AreGids){
      edgeId = _edgeIds.view(firstIdx, nEdges);
    }
    else{
      edgeId = edgeGnos_.view(firstIdx, nEdges);
    }

    procId = procIds_.view(firstIdx, nEdges);
    // TODO edge weights

    return nEdges;
  }
};
#endif

}   // namespace Zoltan2

#endif
