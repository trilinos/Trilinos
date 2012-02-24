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

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_IdentifierInput.hpp>

#include <vector>
#include <Teuchos_Hashtable.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphModel
    \brief GraphModel defines the interface required for graph models.  

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an Input Adapter.  Input adapters are
    templated on the basic user input type.
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  
  GraphModel(){
    throw std::logic_error("in non-specialized GraphModel");
  }

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return 0; }

  /*! Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const { return 0; }

  /*! Returns the number of global edges on this process.
   *  Includes remote edges.
   */
  size_t getLocalNumGlobalEdges() const { return 0; }

  /*! Returns the number of local edges on this process.
   *  Does not include remote edges.
   */
  size_t getLocalNumLocalEdges() const { return 0; }

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
      KDDKDD This function should not return coordinates, as coordinates are
      KDDKDD not needed in most graph algorithms.  There should be a separate
      KDDKDD function for coordinates.  Or there should be an option to 
      KDDKDD specify whether or not coordinates are needed.
   */
  // GNOs are returned.
  // LNOs are implied by the order that the vertices have been returned.

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const {
      return 0; }

  /*! getEdgeList:
      Sets pointers to this process' edge (neighbor) global Ids, including
      off-process edges.
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
  // Implied Vertex LNOs from getVertexList are used as indices to offsets 
  // array.
  // Vertex GNOs are returned as neighbors in edgeIds.

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const { return 0; }

  /*! getLocalEdgeList:
      Sets pointers to this process' local-only edge (neighbor) LNOs, using
      the same implied vertex LNOs returned in getVertexList.
      \param edgeIds lists the only neighbors of the vertices in getVertexList
        which are on this process.  The Id returned is not the neighbor's 
        global Id, but rather the index of the neighbor in the list 
        returned by getVertexList.
      \param offsets offsets[i] is the offset into edgeIds to the start
        of neighbors for ith vertex returned in getVertexList.
      \param wgts will on return point to a list of the weight or weights
         associated with each edge in the edgeIds list.  Weights are listed by
         edge by weight component.
       \return The number of ids in the edgeIds list.

       This method is not const, because the local edge list is not created
       unless this method is called.
   */

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) { return 0; }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumVertices();
  }

  global_size_t getGlobalNumObjects() const
  {
    return getGlobalNumVertices();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const { return; }

  int getNumWeights() const { return 0; }
};

////////////////////////////////////////////////////////////////
// Graph model derived from MatrixInput.
//
//   TODO: support a flag that says the vertices are columns or
//           non-zeros rather than rows
////////////////////////////////////////////////////////////////

/*!  \brief A specialization of GraphModel for Zoltan2::MatrixInput.
*/

template <typename User>
class GraphModel<MatrixInput<User> > : public Model<MatrixInput<User> >
{
public:

  typedef typename MatrixInput<User>::scalar_t  scalar_t;
  typedef typename MatrixInput<User>::gno_t     gno_t;
  typedef typename MatrixInput<User>::lno_t     lno_t;
  typedef typename MatrixInput<User>::gid_t     gid_t;
  typedef typename MatrixInput<User>::node_t    node_t;
  typedef IdentifierMap<User>     idmap_t;

  /*! Constructor
   *  All processes in the communicator must call the constructor.
   *
   *  \param  inputAdapter  an encapsulation of the user data
   *  \param  env           environment (library configuration settings)
   *  \param  comm       communicator for the problem
   *  \param  modelFlags  a bit map of GraphModelFlags
   */

  GraphModel(const MatrixInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    unsigned int modelFlags):
     input_(ia), env_(env), comm_(comm),
     gids_(), gnos_(), edgeGids_(), edgeGnos_(), procIds_(), 
     offsets_(), gnosConst_(), edgeGnosConst_(), procIdsConst_(), 
     numLocalEdges_(0), numGlobalEdges_(0), numLocalVtx_(0), 
     gidsAreGnos_(false), nearEdgeLnos_(), nearEdgeOffsets_(), 
     numNearLocalEdges_(0)
  {
    initializeData(modelFlags);
  }

  //!  Destructor
  ~GraphModel() { }

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

  // TODO - move these out of definition
  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const
  {
    size_t n = getLocalNumVertices();
    Ids = ArrayView<const gno_t>(Teuchos::null);

    if (n){
      if (gidsAreGnos_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
    // KDDKDD  Is it dangerous that xyz is ignored here?  Perhaps coordinates
    // KDDKDD  should be separate.
    // LRIESEN  Coordinates and weights are not yet implemented.
  }

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const
  {
    edgeIds = ArrayView<const gno_t>(Teuchos::null);
    procIds = procIdsConst_(0, numLocalEdges_);
    offsets = offsets_(0, numLocalVtx_+1);
    wgts = ArrayView<const scalar_t>(Teuchos::null);

    if (numLocalEdges_){
      if (gidsAreGnos_)
        edgeIds = edgeGids_(0, numLocalEdges_);
      else
        edgeIds = edgeGnosConst_(0, numLocalEdges_);
    }

    return numLocalEdges_;
  }

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts)
  {
    int nvtx = numLocalVtx_;
    
    edgeIds = ArrayView<const lno_t>(Teuchos::null);
    wgts = ArrayView<const scalar_t>(Teuchos::null);

    if (nearEdgeOffsets_.size() > 0){
      if (numNearLocalEdges_ > 0)  // Teuchos bug: view(0,0) crashes
        edgeIds = nearEdgeLnos_.view(0, numNearLocalEdges_);

      offsets = nearEdgeOffsets_.view(0, nvtx + 1);
    }
    else{
      lno_t *offs = new lno_t [nvtx + 1];
      Z2_LOCAL_MEMORY_ASSERTION(*env_, nvtx+1, offs);
      numNearLocalEdges_ = 0;

      offs[0] = 0;
      if (nvtx > 0){
        for (lno_t i=0; i < nvtx; i++){
          offs[i+1] = 0;
          for (lno_t j=offsets_[i]; j < offsets_[i+1]; j++){
            if (procIds_[j] == env_->myRank_){
              offs[i+1]++;
              numNearLocalEdges_++;
            }
          }
        }

        if (numNearLocalEdges_ > 0){
          gno_t *gnos = new gno_t [numNearLocalEdges_];
          Z2_LOCAL_MEMORY_ASSERTION(*env_, numNearLocalEdges_, gnos);
          ArrayView<gno_t> gnoList(gnos, numNearLocalEdges_);

          for (lno_t i=2; i < nvtx; i++){
            offs[i] += offs[i-1];
          }

          const gno_t *graphGnos = NULL;
          if (gidsAreGnos_)
            graphGnos = reinterpret_cast<const gno_t *>(edgeGids_.getRawPtr());
          else
            graphGnos = edgeGnosConst_.getRawPtr();

          for (lno_t i=0; i < nvtx; i++){
            for (lno_t j=offsets_[i]; j < offsets_[i+1]; j++){
              if (procIds_[j] == env_->myRank_){
                gnoList[offs[i]] = graphGnos[j];
                offs[i]++;
              }
            }
          }

          lno_t *lnos = new lno_t [numNearLocalEdges_];
          Z2_LOCAL_MEMORY_ASSERTION(*env_, numNearLocalEdges_, lnos);
          ArrayRCP<lno_t> lnoList(lnos, 0, numNearLocalEdges_, true);
          RCP<const idmap_t > idMap = this->getIdentifierMap();

          idMap->lnoTranslate(lnoList.view(0, numNearLocalEdges_), gnoList, 
             TRANSLATE_LIB_TO_APP);

          delete [] gnos;
          nearEdgeLnos_ = arcp_const_cast<const lno_t>(lnoList);
          edgeIds = nearEdgeLnos_.view(0, numNearLocalEdges_);
          
          for (lno_t i = nvtx; i > 0; i--){
            offs[i] = offs[i-1];
          }
          offs[0] = 0;
        }
      }

      nearEdgeOffsets_ = arcp(offs, 0, nvtx+1);
      offsets = nearEdgeOffsets_.view(0, nvtx + 1);
    }
    return numNearLocalEdges_;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumVertices();
  }

  global_size_t getGlobalNumObjects() const
  {
    return getGlobalNumVertices();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<const scalar_t> xyz, wgts;
    getVertexList(gnos, xyz, wgts);
  }

  int getNumWeights() const { return 0; }

private:

  void initializeData(unsigned int);

  const MatrixInput<User> *input_;
  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  ArrayRCP<const gid_t> gids_;
  ArrayRCP<gno_t> gnos_;

  // Note: in case of graph subsetting, size of these arrays
  // may be larger than numLocalEdges_.  So do not use .size().

  ArrayRCP<const gid_t> edgeGids_;
  ArrayRCP<gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  ArrayRCP<const gno_t> gnosConst_;
  ArrayRCP<const gno_t> edgeGnosConst_;
  ArrayRCP<const int> procIdsConst_;

  size_t numLocalEdges_;
  global_size_t numGlobalEdges_;
  size_t numLocalVtx_;

  bool gidsAreGnos_;

  // For local graphs (graph restricted to local process).  We only
  // create these arrays if required by the algorithm.

  ArrayRCP<const lno_t> nearEdgeLnos_;
  ArrayRCP<const lno_t> nearEdgeOffsets_;
  size_t numNearLocalEdges_;
};

template <typename User>
  void GraphModel<MatrixInput<User> >::initializeData(unsigned int modelFlags)
{
  // Model creation flags

  bool symTranspose = modelFlags & SYMMETRIZE_INPUT_TRANSPOSE;
  bool symBipartite = modelFlags & SYMMETRIZE_INPUT_BIPARTITE;
  bool vertexRows = modelFlags & VERTICES_ARE_MATRIX_ROWS;
  bool vertexCols = modelFlags & VERTICES_ARE_MATRIX_COLUMNS;
  bool vertexNz = modelFlags & VERTICES_ARE_MATRIX_NONZEROS;
  bool consecutiveIdsRequired = modelFlags & IDS_MUST_BE_GLOBALLY_CONSECUTIVE;
  bool removeSelfEdges = modelFlags & SELF_EDGES_MUST_BE_REMOVED;
  bool subsetGraph = modelFlags & GRAPH_IS_A_SUBSET_GRAPH;

  if (symTranspose || symBipartite || vertexCols || vertexNz){
    throw std::runtime_error("graph build option not yet implemented");
  }

  // Get the matrix from the input adapter

  gid_t const *vtxIds=NULL, *nborIds=NULL;
  lno_t const  *offsets=NULL;
  try{
    numLocalVtx_ = input_->getRowListView(vtxIds, offsets, nborIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  gids_ = arcp(vtxIds, 0, numLocalVtx_, false);

  numLocalEdges_ = offsets[numLocalVtx_];

  // Remove self edges if necessary.

  ArrayRCP<lno_t> tmpOffsets;
  ArrayRCP<gid_t> tmpEdges;
  lno_t nSelfEdges = 0;

  size_t numOffsets = numLocalVtx_ + 1;

  if (removeSelfEdges && input_->diagonalEntriesMayBePresent()) {

    lno_t *offArray = new lno_t [numOffsets];
    Z2_LOCAL_MEMORY_ASSERTION(*env_, numOffsets, offArray);
    gid_t *edArray = new gid_t [numLocalEdges_];
    Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalEdges_, !numLocalEdges_||edArray);

    for (lno_t i=0; i < numLocalVtx_; i++){

      offArray[i] = offsets[i] - nSelfEdges;

      for (lno_t j = offsets[i]; j < offsets[i+1]; j++) {
        if (gids_[i] == nborIds[j]) { // self edge; remove it
          nSelfEdges++;
        }
        else {  // Not a self-edge; keep it.
          edArray[j-nSelfEdges] = nborIds[j];
        }
      }
    }
    numLocalEdges_ -= nSelfEdges;
    offArray[numLocalVtx_] = numLocalEdges_;

    if (nSelfEdges > 0){
      tmpOffsets = arcp(offArray, 0, numLocalVtx_+1);
      tmpEdges = arcp(edArray, 0, numLocalEdges_);
    }
    else{
      delete [] offArray;
      if (numLocalEdges_) delete [] edArray;
    }
  }

  if (nSelfEdges == 0){
    offsets_ = arcp(const_cast<lno_t *>(offsets), 0, numOffsets, false);
    edgeGids_ = arcp(const_cast<gno_t *>(nborIds), 0, numLocalEdges_, false);
  }
  else{
    offsets_ = tmpOffsets;
    edgeGids_ =  tmpEdges;
  }

  // Create an IdentifierMap, which will map the user's global IDs to
  //   Zoltan2 internal global numbers if neccesary.
  //   The map can also give us owners of our vertex neighbors.

  RCP<const idmap_t> idMap;

  try{
    idMap = rcp(new idmap_t(env_, comm_, gids_, consecutiveIdsRequired));
  }
  Z2_FORWARD_EXCEPTIONS;

  gidsAreGnos_ = idMap->gnosAreGids();

  if (numLocalVtx_ && !gidsAreGnos_){
    gno_t *tmp = new gno_t [numLocalVtx_];
    Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalVtx_, tmp)
    gnos_ = arcp(tmp, 0, numLocalVtx_);

    try{
      // Because gidTranslate can translate gids to gnos or
      // gnos to gids depending on a flag, neither the gids nor
      // the gnos are declared to be const.
      ArrayRCP<gid_t> tmpGids = arcp_const_cast<gid_t>(gids_);

      idMap->gidTranslate(tmpGids(0,numLocalVtx_), gnos_(0,numLocalVtx_), 
        TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;

    if (numLocalEdges_){
      tmp = new gno_t [numLocalEdges_];
      Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalEdges_, tmp)
      edgeGnos_ = arcp(tmp, 0, numLocalEdges_);
    }
  }

  if (numLocalEdges_){
    int *p = new int [numLocalEdges_];
    Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalEdges_, p)
    procIds_ = arcp(p, 0, numLocalEdges_);
  }

  ArrayView<const gid_t> gidArray(Teuchos::null);
  ArrayView<gno_t> gnoArray(Teuchos::null);
  ArrayView<int> procArray(Teuchos::null);

  if (numLocalEdges_){
    gidArray = edgeGids_.view(0, numLocalEdges_);
    procArray = procIds_.view(0, numLocalEdges_);
    if (!gidsAreGnos_)
      gnoArray = edgeGnos_.view(0, numLocalEdges_);
  }
    
  try{
    // All processes must make this call.
    idMap->gidGlobalTranslate(gidArray, gnoArray, procArray);
  }
  Z2_FORWARD_EXCEPTIONS;

  this->setIdentifierMap(idMap);   // Zoltan2::Model method

  // Check for edges that are not in our list of global Ids.  If
  // we are subsetting a graph that is not an error, otherwise
  // it is an error.

  int *nborProc = procIds_.getRawPtr();
  size_t numRemoteEdges = 0;
  for (size_t i=0; i < numLocalEdges_; i++){
    if (nborProc[i] < 0)
      numRemoteEdges++;
  }

  if (numRemoteEdges > 0){

    if (!subsetGraph){
      Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid edge ids", 1, BASIC_ASSERTION)
    }
    else{ // Remove edges that are not in the sub graph

      size_t numNewEdges = numLocalEdges_ - numRemoteEdges;

      const lno_t *offFrom = offsets_.getRawPtr();
      lno_t *offTo = const_cast<lno_t *>(offFrom);

      if (offFrom == offsets){  // can't overwrite user's data
        offTo = new lno_t [numLocalVtx_ + 1];
        Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalVtx_+1, offTo);
      }
  
      const gid_t *egidFrom = edgeGids_.getRawPtr();
      gid_t *egidTo = const_cast<gid_t *>(egidFrom);

      if (egidFrom == nborIds){ // can't overwrite user's data
        egidTo = new gid_t [numNewEdges];
        Z2_LOCAL_MEMORY_ASSERTION(*env_, numNewEdges, egidTo);
      }

      gno_t *egno = NULL;
      if (!gidsAreGnos_)
        egno = edgeGnos_.getRawPtr();

      offTo[0] = 0;

      for (lno_t i=0, idx=0; i < numLocalVtx_; i++){
        for (lno_t j=offFrom[i]; j < offFrom[i+1]; j++){
          if (nborProc[j] >= 0){
            egidTo[idx] = egidFrom[j];
            nborProc[idx] = nborProc[j];
            if (egno)
              egno[idx] = egno[j];
            idx++;
          }
        }
        offTo[i+1] = idx;
      }

      if (offTo != offFrom)
        offsets_ = arcp(offTo, 0, numLocalVtx_+1, true);

      if (egidTo != egidFrom)
        edgeGids_ = arcp(egidTo, 0, numNewEdges, true);

      numLocalEdges_ = numNewEdges;
    }
  }

  reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
    &numLocalEdges_, &numGlobalEdges_);

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);
  edgeGnosConst_ = arcp_const_cast<const gno_t>(edgeGnos_);
  procIdsConst_ = arcp_const_cast<const int>(procIds_);
}

////////////////////////////////////////////////////////////////
// Graph model derived from IdentifierInput.
//
//  We do not build a graph model from identifiers.  We include
//  this definition so that other code will compile.
////////////////////////////////////////////////////////////////

/*!  \brief An empty specialization of GraphModel for IdentifierInput

    We do not build graphs from identifier lists, but this definition
    must exist in order for other code to compile.
*/


template <typename User>
class GraphModel<IdentifierInput<User> > : public Model<IdentifierInput<User> >
{
public:

  typedef typename IdentifierInput<User>::scalar_t  scalar_t;
  typedef typename IdentifierInput<User>::gno_t     gno_t;
  typedef typename IdentifierInput<User>::lno_t     lno_t;

  GraphModel(const IdentifierInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    bool graphVerticesAreMatrixRows=true,
    bool consecutiveIdsRequired=false, bool removeSelfEdges=false,
    bool subsetGraph=false) 

  {
    throw std::runtime_error("may not build a graph with identifiers");
  }

  // GraphModel interface

  size_t getLocalNumVertices() const { return 0;}
  global_size_t getGlobalNumVertices() const { return 0;}
  size_t getLocalNumEdges() const { return 0;}
  global_size_t getGlobalNumEdges() const {return 0;}
  int getVertexWeightDim() const { return 0; }
  int getEdgeWeightDim() const { return 0; }
  int getCoordinateDim() const { return 0; }
  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, 
    ArrayView<const scalar_t> &wgts) const { return 0; }
  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const { return 0; }
  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) { return 0; }

  // Model interface

  size_t getLocalNumObjects() const { return 0; }
  global_size_t getGlobalNumObjects() const { return 0; }
  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {}
  int getNumWeights() const { return 0; }

};

}   // namespace Zoltan2

#endif

