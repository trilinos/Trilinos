// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_GraphModel.hpp
    \brief Defines the GraphModel interface.
*/

#ifndef _ZOLTAN2_GRAPHMODEL_HPP_
#define _ZOLTAN2_GRAPHMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_IdentifierInput.hpp>
#include <Zoltan2_CoordinateInput.hpp>
#include <Zoltan2_VectorInput.hpp>

#include <vector>
#include <Teuchos_Hashtable.hpp>

namespace Zoltan2 {

/*!  \brief GraphModel defines the interface required for graph models.  

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.

    Explicit instantiations of the GraphModel exist for:
      \li MatrixInput

    \todo instantiations for GraphInput, MeshInput
    \todo use StridedData objects for coordinates and weights, call
             base model class setWeightArrayLengths
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
#endif

  /*! \brief Constructor
   *
   *  \param  inputAdapter  a pointer to the user's data
   *  \param  env           object containing the parameters
   *  \param  comm          communicator for the problem
   *  \param  modelFlags    a bit map of Zoltan2::GraphModelFlags
   *
   *  All processes in the communicator must call the constructor.
   *  \todo document the model flags that might be set
   */

  GraphModel(const Adapter *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags)
  {
    throw std::logic_error("in non-specialized GraphModel");
  }

  /*! \brief Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return 0; }

  /*! \brief Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const { return 0; }

  /*! \brief Returns the number of global edges on this process.
   *  Includes remote edges.
   */
  size_t getLocalNumGlobalEdges() const { return 0; }

  /*! \brief Returns the number of local edges on this process.
   *  Does not include remote edges.
   */
  size_t getLocalNumLocalEdges() const { return 0; }

  /*! \brief Returns the global number edges.
   */
  global_size_t getGlobalNumEdges() const { return 0; }

  /*! \brief Returns the dimension (0 or greater) of vertex weights.
   */
  int getVertexWeightDim() const { return 0; }

  /*! \brief Returns the dimension (0 or greater) of edge weights.
   */
  int getEdgeWeightDim() const { return 0; }

  /*! \brief Returns the dimension (0 to 3) of vertex coordinates.
   */
  int getCoordinateDim() const { return 0; }

  /*! \brief Sets pointers to this process' vertex Ids and their weights.

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

  /*! \brief Sets pointers to this process' edge (neighbor) global Ids, including
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

  /*! \brief Sets pointers to this process' local-only edge (neighbor) LNOs, using
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
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

////////////////////////////////////////////////////////////////
// Graph model derived from MatrixInput.
////////////////////////////////////////////////////////////////

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
   *  \param  modelFlags  a bit map of Zoltan2::GraphModelFlags
   */

  GraphModel(const MatrixInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags):
     input_(ia), env_(env), comm_(comm),
     gids_(), gnos_(), numGlobalVertices_(), 
     edgeGids_(), edgeGnos_(), procIds_(), 
     offsets_(), gnosConst_(), edgeGnosConst_(), procIdsConst_(), 
     numLocalEdges_(0), numGlobalEdges_(0), numLocalVtx_(0), 
     gidsAreGnos_(false), nearEdgeLnos_(), nearEdgeOffsets_(), 
     numNearLocalEdges_(0)
  {
    initializeData(modelFlags);
    env_->memory("After construction of graph model");
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
    return numGlobalVertices_;
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
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const;

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const;

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts);

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

  void initializeData(modelFlag_t &);

  const MatrixInput<User> *input_;
  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  ArrayRCP<const gid_t> gids_;
  ArrayRCP<gno_t> gnos_;
  global_size_t numGlobalVertices_;

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
  void GraphModel<MatrixInput<User> >::initializeData(modelFlag_t &modelFlags)
{
  // Model creation flags

  bool symTranspose = modelFlags.test(SYMMETRIZE_INPUT_TRANSPOSE);
  bool symBipartite = modelFlags.test(SYMMETRIZE_INPUT_BIPARTITE);
  bool vertexCols = modelFlags.test(VERTICES_ARE_MATRIX_COLUMNS);
  bool vertexNz = modelFlags.test(VERTICES_ARE_MATRIX_NONZEROS);
  bool consecutiveIdsRequired = 
    modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  bool removeSelfEdges = modelFlags.test(SELF_EDGES_MUST_BE_REMOVED);
  bool subsetGraph = modelFlags.test(GRAPH_IS_A_SUBSET_GRAPH);

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
    env_->localMemoryAssertion(__FILE__, __LINE__, numOffsets, offArray);
    gid_t *edArray = new gid_t [numLocalEdges_];
    env_->localMemoryAssertion(__FILE__, __LINE__, numLocalEdges_,
      !numLocalEdges_||edArray);

    for (size_t i=0; i < numLocalVtx_; i++){

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

  numGlobalVertices_ = idMap->getGlobalNumberOfIds();
  gidsAreGnos_ = idMap->gnosAreGids();

  if (numLocalVtx_ && !gidsAreGnos_){
    gno_t *tmp = new gno_t [numLocalVtx_];
    env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVtx_, tmp);
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
      env_->localMemoryAssertion(__FILE__, __LINE__, numLocalEdges_, tmp);
      edgeGnos_ = arcp(tmp, 0, numLocalEdges_);
    }
  }

  if (numLocalEdges_){
    int *p = new int [numLocalEdges_];
    env_->localMemoryAssertion(__FILE__, __LINE__, numLocalEdges_, p);
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

  Array<lno_t> weightArrayLengths;
  this->setWeightArrayLengths(weightArrayLengths, *comm_);
  this->setIdentifierMap(idMap);   // Zoltan2::Model method

  // Check for edges that are not in our list of global Ids.  If
  // we are subsetting a graph that is not an error, otherwise
  // it is an error.

  int *nborProc = NULL;
  if (numLocalEdges_)
    nborProc = procIds_.getRawPtr();

  size_t numRemoteEdges = 0;
  for (size_t i=0; i < numLocalEdges_; i++){
    if (nborProc[i] < 0)
      numRemoteEdges++;
  }

  if (numRemoteEdges > 0){

    if (!subsetGraph){
      env_->localInputAssertion(__FILE__, __LINE__, "invalid edge ids", 1, 
        BASIC_ASSERTION);
    }
    else{ // Remove edges that are not in the sub graph

      size_t numNewEdges = numLocalEdges_ - numRemoteEdges;

      const lno_t *offFrom = offsets_.getRawPtr();
      lno_t *offTo = const_cast<lno_t *>(offFrom);

      if (offFrom == offsets){  // can't overwrite user's data
        offTo = new lno_t [numLocalVtx_ + 1];
        env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVtx_+1, offTo);
      }
  
      const gid_t *egidFrom = edgeGids_.getRawPtr();
      gid_t *egidTo = const_cast<gid_t *>(egidFrom);

      if (egidFrom == nborIds){ // can't overwrite user's data
        egidTo = new gid_t [numNewEdges];
        env_->localMemoryAssertion(__FILE__, __LINE__, numNewEdges, egidTo);
      }

      gno_t *egno = NULL;
      if (!gidsAreGnos_)
        egno = edgeGnos_.getRawPtr();

      offTo[0] = 0;

      for (size_t i=0, idx=0; i < numLocalVtx_; i++){
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

template <typename User>
  size_t GraphModel<MatrixInput<User> >::getVertexList( 
    ArrayView<const gno_t> &Ids,
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

template <typename User>
  size_t GraphModel<MatrixInput<User> >::getEdgeList( 
    ArrayView<const gno_t> &edgeIds,
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

template <typename User>
  size_t GraphModel<MatrixInput<User> >::getLocalEdgeList( 
    ArrayView<const lno_t> &edgeIds,
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
    env_->localMemoryAssertion(__FILE__, __LINE__, nvtx+1, offs);
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
        env_->localMemoryAssertion(__FILE__, __LINE__, numNearLocalEdges_, 
          gnos);
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
        env_->localMemoryAssertion(__FILE__, __LINE__, numNearLocalEdges_, 
          lnos);
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

////////////////////////////////////////////////////////////////
// Graph model derived from CoordinateInput.
//
//  We do not build a graph model from coordinates.  We include
//  this definition so that other code will compile.
////////////////////////////////////////////////////////////////

template <typename User>
class GraphModel<CoordinateInput<User> > : public Model<CoordinateInput<User> >
{
public:

  typedef typename CoordinateInput<User>::scalar_t  scalar_t;
  typedef typename CoordinateInput<User>::gno_t     gno_t;
  typedef typename CoordinateInput<User>::lno_t     lno_t;

  GraphModel(const CoordinateInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &flags)
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

////////////////////////////////////////////////////////////////
// Graph model derived from VectorInput.
//
//  We do not build a graph model from a vector.  We include
//  this definition so that other code will compile.
////////////////////////////////////////////////////////////////

template <typename User>
class GraphModel<VectorInput<User> > : public Model<VectorInput<User> >
{
public:

  typedef typename VectorInput<User>::scalar_t  scalar_t;
  typedef typename VectorInput<User>::gno_t     gno_t;
  typedef typename VectorInput<User>::lno_t     lno_t;

  GraphModel(const VectorInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &flags)
  {
    throw std::runtime_error("can not build a graph from a vector");
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

////////////////////////////////////////////////////////////////
// Graph model derived from IdentifierInput.
//
//  We do not build a graph model from identifiers.  We include
//  this definition so that other code will compile.
////////////////////////////////////////////////////////////////

template <typename User>
class GraphModel<IdentifierInput<User> > : public Model<IdentifierInput<User> >
{
public:

  typedef typename IdentifierInput<User>::scalar_t  scalar_t;
  typedef typename IdentifierInput<User>::gno_t     gno_t;
  typedef typename IdentifierInput<User>::lno_t     lno_t;

  GraphModel(const IdentifierInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &flags)
  {
    throw std::runtime_error("can not build a graph with identifiers");
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

#endif    // DOXYGEN_SHOULD_SKIP_THIS

}   // namespace Zoltan2

#endif

