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
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_MeshAdapter.hpp>

#include <vector>
#include <Teuchos_Hashtable.hpp>

namespace Zoltan2 {


//////////////////////////////////////////////////////////////////////////
/*! \brief Helper function to remove undesired edges from a graph.
 *
 *  \param env the environment
 *  \param myRank is my rank in the problem communicator
 *  \param removeSelfEdges true if self-edges (edges such that both
 *     vertices are the same) should be removed. If true then
 *     gids must be set.
 *  \param removeOffProcessEdges true if edges belonging to processes
 *        other than my process should be removed (in which case
 *       \c procIds must be set).
 *  \param removeOffGroupEdges true if edges belonging to processes
 *    outside of our communicator should be removed (in which case
 *       \c procIds must be set).
 *  \param  gids  vertex global Id list
 *  \param gidNbors list of vertex neighbor global ids (edges)
 *  \param procIds is the list of processes owning the vertices in
 *             the \c gidNbors list.
 *  \param edgeWeights weights for edges in \c gidNbors list
 *  \param offsets offset into above lists for each vertex in \c gids.
 *  \param newGidNbors  on return a list of the desired neighbors
 *  \param newWeights if \c wdim is the number of weights per edge
 *       then on return this points to \c wdim pointers to arrays
 *       of weights for the desired edges.  If it is NULL on return,
 *       then one of these must be true:
 *
 *         - \c wdim is zero
 *         - none of the edges are desired
 *         - all of the edges are desired, so you don't need new lists
 *
 *  \param newOffsets  on return a list of offsets into the above list
 *        for the start of the neighbors for each vertex
 *  \return the number of edges left after removal of undesired edges
 *
 *  The template parameter is an Adapter type.
 */

template <typename User>
size_t removeUndesiredEdges(
  const RCP<const Environment> &env,
  int myRank,
  bool removeSelfEdges,
  bool removeOffProcessEdges,
  bool removeOffGroupEdges,
  ArrayView<const typename InputTraits<User>::zgid_t> &gids,
  ArrayView<const typename InputTraits<User>::zgid_t> &gidNbors,
  ArrayView<const int> &procIds,
  ArrayView<StridedData<typename InputTraits<User>::lno_t,
                        typename InputTraits<User>::scalar_t> > &edgeWeights,
  ArrayView<const typename InputTraits<User>::lno_t> &offsets,
  ArrayRCP<const typename InputTraits<User>::zgid_t> &newGidNbors, // out
  typename InputTraits<User>::scalar_t **&newWeights,             // out
  ArrayRCP<const typename InputTraits<User>::lno_t> &newOffsets)  // out
{
  typedef typename InputTraits<User>::zgid_t zgid_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t lno_t;
  size_t numKeep = 0;

  size_t numVtx = offsets.size() - 1;
  size_t numNbors = gidNbors.size();

  env->localInputAssertion(__FILE__, __LINE__, "need more input",
    (!removeSelfEdges ||
      gids.size() >=
       static_cast<typename ArrayView<const zgid_t>::size_type>(numVtx))
      &&
    (!removeOffProcessEdges ||
      procIds.size() >=
       static_cast<typename ArrayView<const int>::size_type>(numNbors)) &&
    (!removeOffGroupEdges ||
      procIds.size() >=
       static_cast<typename ArrayView<const int>::size_type>(numNbors)),
    BASIC_ASSERTION);

  // initialize edge weight array

  newWeights = NULL;
  int eDim = edgeWeights.size();

  // count desired edges

  lno_t *offs = new lno_t [numVtx + 1];
  env->localMemoryAssertion(__FILE__, __LINE__, numVtx+1, offs);
  for (size_t i = 0; i < numVtx+1; i++) offs[i] = 0;
  ArrayRCP<const lno_t> offArray = arcp(offs, 0, numVtx+1, true);

  const lno_t *allOffs = offsets.getRawPtr();
  const zgid_t *allIds = gidNbors.getRawPtr();

  const zgid_t *vtx = NULL;
  const int *proc = NULL;

  if (gids.size() > 0)
    vtx = gids.getRawPtr();

  if (procIds.size() > 0)
    proc = procIds.getRawPtr();

  offs[0] = 0;
  for (size_t i=0; i < numVtx; i++){
    offs[i+1] = 0;
    zgid_t vid = vtx ? vtx[i] : zgid_t(0);
    for (lno_t j=allOffs[i]; j < allOffs[i+1]; j++){
      int owner = proc ? proc[j] : 0;
      bool keep = (!removeSelfEdges || vid != allIds[j]) &&
               (!removeOffProcessEdges || owner == myRank) &&
               (!removeOffGroupEdges || owner >= 0);

      if (keep)
        offs[i+1]++;
    }
  }

  // from counters to offsets

  for (size_t i=1; i < numVtx; i++)
    offs[i+1] += offs[i];

  numKeep = offs[numVtx];

  // do we need a new neighbor list?

  if (numNbors == numKeep){
    newGidNbors = Teuchos::arcpFromArrayView(gidNbors);
    newOffsets = Teuchos::arcpFromArrayView(offsets);
    return numNbors;
  }
  else if (numKeep == 0){
    newGidNbors = ArrayRCP<const zgid_t>(Teuchos::null);
    newOffsets = offArray;
    return 0;
  }

  // Build the subset neighbor lists (id, weight, and offset).

  zgid_t *newGids = new zgid_t [numKeep];
  env->localMemoryAssertion(__FILE__, __LINE__, numKeep, newGids);

  newGidNbors = arcp(newGids, 0, numKeep, true);
  newOffsets = offArray;

  if (eDim > 0){
    newWeights = new scalar_t * [eDim];
    env->localMemoryAssertion(__FILE__, __LINE__, eDim, newWeights);

    if (numKeep) {
      for (int w=0; w < eDim; w++){
        newWeights[w] = new scalar_t [numKeep];
        env->localMemoryAssertion(__FILE__, __LINE__, numKeep, newWeights[w]);
      }
    }
    else {
      for (int w=0; w < eDim; w++)
        newWeights[w] = NULL;
    }
  }

  size_t next = 0;
  for (size_t i=0; i < numVtx && next < numKeep; i++){
    zgid_t vid = vtx ? vtx[i] : zgid_t(0);
    for (lno_t j=allOffs[i]; j < allOffs[i+1]; j++){
      int owner = proc ? proc[j] : 0;
      bool keep = (!removeSelfEdges || vid != allIds[j]) &&
               (!removeOffProcessEdges || owner == myRank) &&
               (!removeOffGroupEdges || owner >= 0);

      if (keep){
        newGids[next] = allIds[j];
        for (int w=0; w < eDim; w++){
          newWeights[w][next] = edgeWeights[w][j];
        }
        next++;
        if (next == numKeep)
          break;

      }  // if (keep)
    }
  }

  return numKeep;
}

//////////////////////////////////////////////////////////////////////////
/*! \brief Helper function to create new edges lists containing
     only edges connected to a neighbor on this process.
 */

template <typename User>
size_t computeLocalEdgeList(
  const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
  size_t numLocalEdges,           // local edges
  size_t numLocalGraphEdges,      // edges in "local" graph
  RCP<const IdentifierMap<User> > &idMap,
  ArrayRCP<const typename InputTraits<User>::zgid_t> &allEdgeIds, // in
  ArrayRCP<const typename InputTraits<User>::gno_t> &allEdgeGnos, // in
  ArrayRCP<int> &allProcs,                                 // in
  ArrayRCP<const typename InputTraits<User>::lno_t> &allOffs,    // in
  ArrayRCP<StridedData<typename InputTraits<User>::lno_t,
                       typename InputTraits<User>::scalar_t> > &allWeights,// in
  ArrayRCP<const typename InputTraits<User>::lno_t> &edgeLocalIds, //
  ArrayRCP<const typename InputTraits<User>::lno_t> &offsets,      // out
  ArrayRCP<StridedData<typename InputTraits<User>::lno_t,
    typename InputTraits<User>::scalar_t> > &eWeights)             // out
{
  typedef typename InputTraits<User>::zgid_t zgid_t;
  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  int rank = comm->getRank();

  bool gnosAreGids = idMap->gnosAreGids();

  edgeLocalIds = ArrayRCP<const lno_t>(Teuchos::null);
  eWeights = ArrayRCP<input_t>(Teuchos::null);
  offsets = ArrayRCP<const lno_t>(Teuchos::null);

  if (numLocalGraphEdges == 0) {
    // Set the offsets array and return
    size_t allOffsSize = allOffs.size();
    lno_t *offs = new lno_t [allOffsSize];
    env->localMemoryAssertion(__FILE__, __LINE__, allOffsSize, offs);
    for (size_t i = 0; i < allOffsSize; i++) offs[i] = 0;
    offsets = arcp(offs, 0, allOffsSize, true);
    return 0;
  }

  if (numLocalGraphEdges == numLocalEdges){

    // Entire graph is local.

    lno_t *lnos = new lno_t [numLocalGraphEdges];
    env->localMemoryAssertion(__FILE__, __LINE__, numLocalGraphEdges, lnos);
    if (comm->getSize() == 1) {
      // With one rank, Can use gnos as local index.
      if (gnosAreGids)
        for (size_t i=0; i < numLocalEdges; i++) lnos[i] = allEdgeIds[i];
      else
        for (size_t i=0; i < numLocalEdges; i++) lnos[i] = allEdgeGnos[i];
    }
    else {
      ArrayRCP<gno_t> gnoArray;

      if (gnosAreGids){
        ArrayRCP<const gno_t> gnosConst =
                 arcp_reinterpret_cast<const gno_t>(allEdgeIds);
        gnoArray = arcp_const_cast<gno_t>(gnosConst);
      }
      else {
        gnoArray = arcp_const_cast<gno_t>(allEdgeGnos);
      }

      // Need to translate to gnos to local indexing
      ArrayView<lno_t> lnoView(lnos, numLocalGraphEdges);
      try {
        idMap->lnoTranslate(lnoView,
                            gnoArray.view(0,numLocalGraphEdges),
                            TRANSLATE_LIB_TO_APP);
      }
      Z2_FORWARD_EXCEPTIONS;
    }
    edgeLocalIds = arcp(lnos, 0, numLocalGraphEdges, true);
    offsets = allOffs;
    eWeights = allWeights;

  }
  else{

    // Create subset list of local graph edges, offsets and weights.

    int nWeightsPerEdge = allWeights.size();

    ArrayRCP<const zgid_t> newEgids;
    scalar_t **newWeights = NULL;

    ArrayView<const zgid_t> dummyVtx;
    ArrayView<const zgid_t> nborView= allEdgeIds.view(0, numLocalEdges);
    ArrayView<const int> nborOwner = allProcs.view(0, numLocalEdges);
    ArrayView<input_t> eWgts = allWeights.view(0, nWeightsPerEdge);
    ArrayView<const lno_t> offView = allOffs.view(0, allOffs.size());

    try{
      numLocalEdges = removeUndesiredEdges<User>(env, rank, false, true, false,
                                                 dummyVtx, nborView, nborOwner,
                                                 eWgts, offView, newEgids,
                                                 newWeights, offsets);
    }
    Z2_FORWARD_EXCEPTIONS;

    env->localBugAssertion(__FILE__, __LINE__, "local graph miscalculation",
      numLocalEdges == numLocalGraphEdges, BASIC_ASSERTION);

    // offsets array was set by removeUndesiredEdges.  Create weight array.

    if (nWeightsPerEdge > 0){
      input_t *wgts = new input_t [nWeightsPerEdge];
      for (int w=0; w < nWeightsPerEdge; w++){
        ArrayRCP<const scalar_t> wgtArray(newWeights[w], 0, numLocalGraphEdges,true);
        wgts[w] = input_t(wgtArray, 1);
      }
      eWeights = arcp(wgts, 0, nWeightsPerEdge);
      delete [] newWeights;
    }

    // Create local ID array.  First translate gid to gno.
    ArrayRCP<gno_t> gnoArray;

    if (gnosAreGids){
      ArrayRCP<const gno_t> gnosConst =
        arcp_reinterpret_cast<const gno_t>(newEgids);
      gnoArray = arcp_const_cast<gno_t>(gnosConst);
    }
    else{

      ArrayRCP<zgid_t> gidArray = arcp_const_cast<zgid_t>(newEgids);
      gno_t *gnoList= new gno_t [numLocalGraphEdges];
      env->localMemoryAssertion(__FILE__, __LINE__, numLocalGraphEdges,
        gnoList);
      gnoArray = arcp(gnoList, 0, numLocalGraphEdges, true);

      try {
        idMap->gidTranslate(
          gidArray.view(0,numLocalGraphEdges),
          gnoArray.view(0,numLocalGraphEdges),
          TRANSLATE_APP_TO_LIB);
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    // translate gno to lno

    lno_t *lnoList = new lno_t [numLocalGraphEdges];
    env->localMemoryAssertion(__FILE__, __LINE__, numLocalGraphEdges,
      lnoList);
    ArrayView<lno_t> lnoView(lnoList, numLocalGraphEdges);

    try {
      idMap->lnoTranslate(
        lnoView,
        gnoArray.view(0,numLocalGraphEdges),
        TRANSLATE_LIB_TO_APP);
    }
    Z2_FORWARD_EXCEPTIONS;
    edgeLocalIds = arcp<const lno_t>(lnoList, 0, numLocalGraphEdges, true);
  }

  return numLocalGraphEdges;
}

//////////////////////////////////////////////////////////////////////////
/*!  \brief GraphModel defines the interface required for graph models.

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::zgid_t       zgid_t;
  typedef typename Adapter::node_t      node_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef IdentifierMap<user_t>         idmap_t;
  typedef StridedData<lno_t, scalar_t>  input_t;
#endif

  //!  Destructor
  ~GraphModel() { }

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

  GraphModel(const MatrixAdapter<user_t,userCoord_t> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags);

  GraphModel(const GraphAdapter<user_t,userCoord_t> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags);

  GraphModel(const MeshAdapter<user_t> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelflags);

  GraphModel(const VectorAdapter<userCoord_t> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags)
  {
    throw std::runtime_error("cannot build GraphModel from VectorAdapter");
  }

  GraphModel(const IdentifierAdapter<user_t> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags)
  {
    throw std::runtime_error("cannot build GraphModel from IdentifierAdapter");
  }

  /*! \brief Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return numLocalVertices_; }

  /*! \brief Returns the global number vertices.
   */
  size_t getGlobalNumVertices() const { return numGlobalVertices_; }

  /*! \brief Returns the number of global edges on this process.
   *  Includes remote edges.
   */
  size_t getLocalNumGlobalEdges() const { return numLocalEdges_; }

  /*! \brief Returns the number of local edges on this process.
   *  Does not include edges to off-process vertices.
   */
  size_t getLocalNumLocalEdges() const { return numLocalGraphEdges_; }

  /*! \brief Returns the global number edges.
   */
  size_t getGlobalNumEdges() const { return numGlobalEdges_; }

  /*! \brief Returns the number (0 or greater) of weights per vertex
   */
  int getNumWeightsPerVertex() const { return numWeightsPerVertex_; }

  /*! \brief Returns the number (0 or greater) of weights per edge.
   */
  int getNumWeightsPerEdge() const { return nWeightsPerEdge_; }

  /*! \brief Returns the dimension (0 to 3) of vertex coordinates.
   */
  int getCoordinateDim() const { return vCoordDim_; }

  /*! \brief Sets pointers to this process' vertex Ids and their weights.

      \param Ids will on return point to the list of the global Ids for
        each vertex on this process.
      \param xyz If vertex coordinate data is available, \c xyz
         will on return point to a StridedData object of coordinates.
      \param wgts If vertex weights is available, \c wgts
         will on return point to a StridedData object of weights.
   */

  size_t getVertexList(
    ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const
  {
    Ids = ArrayView<const gno_t>();
    size_t nv = gids_.size();
    if (nv){
      if (gnosAreGids_)
        Ids = ArrayView<const gno_t>(
                        reinterpret_cast<const gno_t*>(gids_.getRawPtr()), nv);
      else
        Ids = gnosConst_(0, nv);
    }

    xyz = vCoords_.view(0, vCoordDim_);
    wgts = vWeights_.view(0, numWeightsPerVertex_);
    return nv;
  }

  /*! \brief Sets pointers to this process' edge (neighbor) global Ids,
      including off-process edges.

      \param edgeIds This is the list of global neighbor Ids corresponding
        to the vertices listed in getVertexList.
      \param procIds lists the process owning each neighbor in the edgeIds
         list.
      \param offsets offsets[i] is the offset into edgeIds to the start
        of neighbors for ith vertex.
      \param wgts If edge weights is available, \c wgts
         will on return point to a StridedData object of weights.
       \return The number of ids in the edgeIds list.
   */
  // Implied Vertex LNOs from getVertexList are used as indices to offsets
  // array.
  // Vertex GNOs are returned as neighbors in edgeIds.

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) const
  {
    edgeIds = ArrayView<const gno_t>();
    if (numLocalEdges_) {
      if (gnosAreGids_)
        edgeIds = ArrayView<const gno_t>(
                        reinterpret_cast<const gno_t*>(edgeGids_.getRawPtr()),
                                                       numLocalEdges_);
      else
        edgeIds = edgeGnosConst_(0, numLocalEdges_);
    }

    procIds = procIdsConst_.view(0, numLocalEdges_);
    offsets = offsets_.view(0, numLocalVertices_+1);
    wgts = eWeights_.view(0, nWeightsPerEdge_);
    return numLocalEdges_;
  }

  /*! \brief Sets pointers to this process' local-only edge (neighbor) LNOs,
      using the same implied vertex LNOs returned in getVertexList.

      Local only means the neighbor vertex is owned by this process.

      \param edgeIds lists the only neighbors of the vertices in getVertexList
        which are on this process.  The Id returned is not the neighbor's
        global Id, but rather the index of the neighbor in the list
        returned by getVertexList.
      \param offsets offsets[i] is the offset into edgeIds to the start
        of neighbors for ith vertex returned in getVertexList.
      \param wgts If edge weights is available, \c wgts
         will on return point to a StridedData object of weights.
       \return The number of ids in the edgeIds list.

       This method is not const, because a local edge list is not created
       unless this method is called.

       Note that if there are no local edges, the
         \c edgeIds, \c offsets and \c wgts are returned
         as empty arrays.
   */

  size_t getLocalEdgeList(
    ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts)
  {
    env_->timerStart(MACRO_TIMERS, "GraphModel::getLocalEdgeList");

    if (localGraphEdgeOffsets_.size() == 0) {
      // Local graph not created yet
      RCP<const IdentifierMap<user_t> > idmap = this->getIdentifierMap();
      computeLocalEdgeList(env_, comm_,
        numLocalEdges_, numLocalGraphEdges_,
        idmap, edgeGids_, edgeGnosConst_, procIds_, offsets_, eWeights_,
        localGraphEdgeLnos_, localGraphEdgeOffsets_, localGraphEdgeWeights_);
    }
    edgeIds = localGraphEdgeLnos_();
    offsets = localGraphEdgeOffsets_();
    wgts = localGraphEdgeWeights_();

    env_->timerStop(MACRO_TIMERS, "GraphModel::getLocalEdgeList");

    return numLocalGraphEdges_;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const { return numLocalVertices_; }

  size_t getGlobalNumObjects() const { return numGlobalVertices_; }

  void get2ndAdjsViewFromAdjs(const Adapter *ia,
			      Zoltan2::MeshEntityType sourcetarget,
			      Zoltan2::MeshEntityType through,
			      const lno_t *&offsets,
			      const zgid_t *&adjacencyIds);
private:
  void shared_constructor(const Adapter *ia, modelFlag_t &modelFlags);

  template <typename AdapterWithCoords>
  void shared_GetVertexCoords(const AdapterWithCoords *ia);

  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  ArrayRCP<const zgid_t> gids_;        // vertices of input graph
  ArrayRCP<gno_t> gnos_;

  int numWeightsPerVertex_;
  ArrayRCP<input_t> vWeights_;

  int vCoordDim_;
  ArrayRCP<input_t> vCoords_;

  // Note: in case of graph subsetting, size of these arrays
  // may be larger than numLocalEdges_.  So do not use .size().

  ArrayRCP<const zgid_t> edgeGids_;
  ArrayRCP<gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  int nWeightsPerEdge_;
  ArrayRCP<input_t> eWeights_;

  ArrayRCP<const gno_t> gnosConst_;
  ArrayRCP<const gno_t> edgeGnosConst_;
  ArrayRCP<const int> procIdsConst_;

  bool gnosAreGids_;

  // For local graphs (graph restricted to local process).  We
  // create these arrays only if required by the algorithm.

  ArrayRCP<const lno_t> localGraphEdgeLnos_;
  ArrayRCP<const lno_t> localGraphEdgeOffsets_;
  ArrayRCP<input_t> localGraphEdgeWeights_;

  // For convenience

  size_t numLocalVertices_;
  size_t numGlobalVertices_;
  size_t numLocalEdges_;
  size_t numGlobalEdges_;
  size_t numLocalGraphEdges_;

  // For debugging
  void print();
};


////////////////////////////////////////////////////////////////
template <typename Adapter>
GraphModel<Adapter>::GraphModel(
  const MatrixAdapter<user_t,userCoord_t> *ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags):
       env_(env),
       comm_(comm),
       gids_(),
       gnos_(),
       numWeightsPerVertex_(0),
       vWeights_(),
       vCoordDim_(0),
       vCoords_(),
       edgeGids_(),
       edgeGnos_(),
       procIds_(),
       offsets_(),
       nWeightsPerEdge_(0),
       eWeights_(),
       gnosConst_(),
       edgeGnosConst_(),
       procIdsConst_(),
       gnosAreGids_(false),
       localGraphEdgeLnos_(),
       localGraphEdgeOffsets_(),
       localGraphEdgeWeights_(),
       numLocalVertices_(0),
       numGlobalVertices_(0),
       numLocalEdges_(0),
       numGlobalEdges_(0),
       numLocalGraphEdges_(0)
{
  // Model creation flags
  bool symTranspose = modelFlags.test(SYMMETRIZE_INPUT_TRANSPOSE);
  bool symBipartite = modelFlags.test(SYMMETRIZE_INPUT_BIPARTITE);
  bool vertexCols = modelFlags.test(VERTICES_ARE_MATRIX_COLUMNS);
  bool vertexNz = modelFlags.test(VERTICES_ARE_MATRIX_NONZEROS);

  if (symTranspose || symBipartite || vertexCols || vertexNz){
    throw std::runtime_error("graph build option not yet implemented");
  }

  // Get the matrix from the input adapter
  zgid_t const *vtxIds=NULL, *nborIds=NULL;
  lno_t const *offsets=NULL;
  try{
    numLocalVertices_ = ia->getLocalNumIDs();
    ia->getIDsView(vtxIds);
  }
  Z2_FORWARD_EXCEPTIONS;
  try{
    if (ia->CRSViewAvailable()) {
      ia->getCRSView(offsets, nborIds);
    }
    else {
      // TODO:  Add support for CCS matrix layout
      throw std::runtime_error("Only MatrixAdapter::getCRSView is supported "
                               "in graph model");
    }
  }
  Z2_FORWARD_EXCEPTIONS;

  numLocalEdges_ = offsets[numLocalVertices_];

  gids_ = arcp<const zgid_t>(vtxIds, 0, numLocalVertices_, false);
  edgeGids_ = arcp<const zgid_t>(nborIds, 0, numLocalEdges_, false);
  offsets_ = arcp<const lno_t>(offsets, 0, numLocalVertices_ + 1, false);

  nWeightsPerEdge_ = 0;   // no edge weights from a matrix yet.
                     // TODO:  use matrix values as edge weights

  shared_constructor(ia, modelFlags);

  // Get vertex coordinates, if available
  if (ia->coordinatesAvailable()) {
    typedef VectorAdapter<userCoord_t> adapterWithCoords_t;
    shared_GetVertexCoords<adapterWithCoords_t>(ia->getCoordinateInput());
  }
  //print();
}


////////////////////////////////////////////////////////////////
template <typename Adapter>
GraphModel<Adapter>::GraphModel(
  const GraphAdapter<user_t,userCoord_t> *ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags):
       env_(env),
       comm_(comm),
       gids_(),
       gnos_(),
       numWeightsPerVertex_(0),
       vWeights_(),
       vCoordDim_(0),
       vCoords_(),
       edgeGids_(),
       edgeGnos_(),
       procIds_(),
       offsets_(),
       nWeightsPerEdge_(0),
       eWeights_(),
       gnosConst_(),
       edgeGnosConst_(),
       procIdsConst_(),
       gnosAreGids_(false),
       localGraphEdgeLnos_(),
       localGraphEdgeOffsets_(),
       localGraphEdgeWeights_(),
       numLocalVertices_(0),
       numGlobalVertices_(0),
       numLocalEdges_(0),
       numGlobalEdges_(0),
       numLocalGraphEdges_(0)
{

  // This GraphModel is built with vertices == GRAPH_VERTEX from GraphAdapter.
  // It is not ready to use vertices == GRAPH_EDGE from GraphAdapter.
  env_->localInputAssertion(__FILE__, __LINE__,
    "GraphModel from GraphAdapter is implemented only for "
    "Graph Vertices as primary object, not for Graph Edges",
    ia->getPrimaryEntityType() == Zoltan2::GRAPH_VERTEX, BASIC_ASSERTION);

  // Get the graph from the input adapter

  zgid_t const *vtxIds=NULL, *nborIds=NULL;
  lno_t const *offsets=NULL;
  try{
    numLocalVertices_ = ia->getLocalNumVertices();
    ia->getVertexIDsView(vtxIds);
    ia->getEdgesView(offsets, nborIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  numLocalEdges_ = offsets[numLocalVertices_];

  gids_ = arcp<const zgid_t>(vtxIds, 0, numLocalVertices_, false);
  edgeGids_ = arcp<const zgid_t>(nborIds, 0, numLocalEdges_, false);
  offsets_ = arcp<const lno_t>(offsets, 0, numLocalVertices_ + 1, false);

  nWeightsPerEdge_ = ia->getNumWeightsPerEdge();

  if (nWeightsPerEdge_ > 0){
    input_t *wgts = new input_t [nWeightsPerEdge_];
    eWeights_ = arcp(wgts, 0, nWeightsPerEdge_, true);
  }

  for (int w=0; w < nWeightsPerEdge_; w++){
    const scalar_t *ewgts=NULL;
    int stride=0;

    ia->getEdgeWeightsView(ewgts, stride, w);

    ArrayRCP<const scalar_t> wgtArray(ewgts, 0, numLocalEdges_, false);
    eWeights_[w] = input_t(wgtArray, stride);
  }

  shared_constructor(ia, modelFlags);

  // Get vertex coordinates, if available
  if (ia->coordinatesAvailable()) {
    typedef VectorAdapter<userCoord_t> adapterWithCoords_t;
    shared_GetVertexCoords<adapterWithCoords_t>(ia->getCoordinateInput());
  }
  //print();
}

////////////////////////////////////////////////////////////////
template <typename Adapter>
GraphModel<Adapter>::GraphModel(
  const MeshAdapter<user_t> *ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags):
       env_(env),
       comm_(comm),
       gids_(),
       gnos_(),
       numWeightsPerVertex_(0),
       vWeights_(),
       vCoordDim_(0),
       vCoords_(),
       edgeGids_(),
       edgeGnos_(),
       procIds_(),
       offsets_(),
       nWeightsPerEdge_(0),
       eWeights_(),
       gnosConst_(),
       edgeGnosConst_(),
       procIdsConst_(),
       gnosAreGids_(false),
       localGraphEdgeLnos_(),
       localGraphEdgeOffsets_(),
       localGraphEdgeWeights_(),
       numLocalVertices_(0),
       numGlobalVertices_(0),
       numLocalEdges_(0),
       numGlobalEdges_(0),
       numLocalGraphEdges_(0)
{
  env_->timerStart(MACRO_TIMERS, "GraphModel constructed from MeshAdapter");

  // This GraphModel is built with vertices == ia->getPrimaryEntityType()
  // from MeshAdapter.

  // Get the graph from the input adapter

  Zoltan2::MeshEntityType primaryEType = ia->getPrimaryEntityType();
  Zoltan2::MeshEntityType secondAdjEType = ia->getSecondAdjacencyEntityType();

  // Get the IDs of the primary entity type; these are graph vertices

  zgid_t const *vtxIds=NULL;
  try {
    numLocalVertices_ = ia->getLocalNumOf(primaryEType);
    ia->getIDsViewOf(primaryEType, vtxIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  gids_ = arcp<const zgid_t>(vtxIds, 0, numLocalVertices_, false);

  // Get the second adjacencies to construct edges of the dual graph.
  // TODO:  Enable building the graph from 1st adjacencies

  zgid_t const *nborIds=NULL;
  lno_t const *offsets=NULL;

  if (!ia->avail2ndAdjs(primaryEType, secondAdjEType)) {

    try {
      get2ndAdjsViewFromAdjs(ia,primaryEType,secondAdjEType,offsets,nborIds);
    }
    Z2_FORWARD_EXCEPTIONS;
    /*throw std::logic_error("MeshAdapter must provide 2nd adjacencies for "
      "graph construction");*/

  }
  else {  // avail2ndAdjs

    // Get the edges
    try {
      ia->get2ndAdjsView(primaryEType, secondAdjEType, offsets, nborIds);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  numLocalEdges_ = offsets[numLocalVertices_];

  edgeGids_ = arcp<const zgid_t>(nborIds, 0, numLocalEdges_, false);
  offsets_ = arcp<const lno_t>(offsets, 0, numLocalVertices_ + 1, false);

  // Get edge weights
  nWeightsPerEdge_ = ia->getNumWeightsPer2ndAdj(primaryEType, secondAdjEType);

  if (nWeightsPerEdge_ > 0){
    input_t *wgts = new input_t [nWeightsPerEdge_];
    eWeights_ = arcp(wgts, 0, nWeightsPerEdge_, true);
  }

  for (int w=0; w < nWeightsPerEdge_; w++){
    const scalar_t *ewgts=NULL;
    int stride=0;

    ia->get2ndAdjWeightsView(primaryEType, secondAdjEType,
			     ewgts, stride, w);

    ArrayRCP<const scalar_t> wgtArray(ewgts, 0, numLocalEdges_, false);
    eWeights_[w] = input_t(wgtArray, stride);
  }

  shared_constructor(ia, modelFlags);

  typedef MeshAdapter<user_t> adapterWithCoords_t;
  shared_GetVertexCoords<adapterWithCoords_t>(ia);

  env_->timerStop(MACRO_TIMERS, "GraphModel constructed from MeshAdapter");
  print();
}

template <typename Adapter>
void GraphModel<Adapter>::get2ndAdjsViewFromAdjs(
  const Adapter *ia,
  Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through,
  const lno_t *&offsets, const zgid_t *&adjacencyIds)
{
  typedef int nonzero_t;  // adjacency matrix doesn't need scalar_t
  typedef Tpetra::CrsMatrix<nonzero_t,lno_t,gno_t,node_t>   sparse_matrix_type;
  typedef Tpetra::Map<lno_t, gno_t, node_t>                 map_type;
  //typedef Tpetra::global_size_t GST;
  //const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  /* Find the adjacency for a nodal based decomposition */
  size_t nadj = 0;
  if (ia->availAdjs(sourcetarget, through)) {
    using Tpetra::DefaultPlatform;
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Get the default communicator and Kokkos Node instance
    // TODO:  Default communicator is not correct here; need to get
    // TODO:  communicator from the problem
    RCP<const Comm<int> > comm =
      DefaultPlatform::getDefaultPlatform ().getComm ();

    // Get node-element connectivity

    offsets=NULL;
    adjacencyIds=NULL;
    ia->getAdjsView(sourcetarget, through, offsets, adjacencyIds);

    zgid_t const *Ids=NULL;
    ia->getIDsViewOf(sourcetarget, Ids);

    zgid_t const *throughIds=NULL;
    ia->getIDsViewOf(through, throughIds);

    size_t LocalNumIDs = ia->getLocalNumOf(sourcetarget);

    /***********************************************************************/
    /************************* BUILD MAPS FOR ADJS *************************/
    /***********************************************************************/

    Array<gno_t> sourcetargetGIDs;
    RCP<const map_type> sourcetargetMapG;
    RCP<const map_type> throughMapG;

    // count owned nodes
    size_t LocalNumOfThrough = ia->getLocalNumOf(through);

    // Build a list of the global sourcetarget ids...
    sourcetargetGIDs.resize (LocalNumIDs);
    gno_t min[2];
    min[0] = as<gno_t> (Ids[0]);
    for (size_t i = 0; i < LocalNumIDs; ++i) {
      sourcetargetGIDs[i] = as<gno_t> (Ids[i]);

      if (sourcetargetGIDs[i] < min[0]) {
	min[0] = sourcetargetGIDs[i];
      }
    }

    // min(throughIds[i])
    min[1] = as<gno_t> (throughIds[0]);
    for (size_t i = 0; i < LocalNumOfThrough; ++i) {
      gno_t tmp = as<gno_t> (throughIds[i]);

      if (tmp < min[1]) {
	min[1] = tmp;
      }
    }

    gno_t gmin[2];
    Teuchos::reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MIN, 2, min, gmin);

    //Generate Map for sourcetarget.
    sourcetargetMapG = rcp(new map_type(ia->getGlobalNumOf(sourcetarget),
					sourcetargetGIDs(), gmin[0], comm));

    //Generate Map for through.
// TODO
// TODO Could check for max through id as well, and if all through ids are
// TODO in gmin to gmax, then default constructors works below.
// TODO Otherwise, may need a constructor that is not one-to-one containing
// TODO all through entities on processor, followed by call to createOneToOne
// TODO

    throughMapG = rcp (new map_type(ia->getGlobalNumOf(through),gmin[1],comm));

    /***********************************************************************/
    /************************* BUILD GRAPH FOR ADJS ************************/
    /***********************************************************************/

    RCP<sparse_matrix_type> adjsMatrix;

    // Construct Tpetra::CrsGraph objects.
    adjsMatrix = rcp (new sparse_matrix_type (sourcetargetMapG, 0));

    nonzero_t justOne = 1;
    ArrayView<nonzero_t> justOneAV = Teuchos::arrayView (&justOne, 1);

    for (size_t localElement=0; localElement<LocalNumIDs; ++localElement){

      //globalRow for Tpetra Graph
      gno_t globalRowT = as<gno_t> (Ids[localElement]);

// KDD can we insert all adjacencies at once instead of one at a time
// (since they are contiguous in adjacencyIds)?
// KDD maybe not until we get rid of zgid_t, as we need the conversion to gno_t.
      for (lno_t j=offsets[localElement]; j<offsets[localElement+1]; ++j){
	gno_t globalCol = as<gno_t> (adjacencyIds[j]);
	//create ArrayView globalCol object for Tpetra
	ArrayView<gno_t> globalColAV = Teuchos::arrayView (&globalCol,1);

	//Update Tpetra adjs Graph
	adjsMatrix->insertGlobalValues(globalRowT,globalColAV,justOneAV);
      }// *** through loop ***
    }// *** source loop ***

    //Fill-complete adjs Graph
    adjsMatrix->fillComplete (throughMapG, adjsMatrix->getRowMap());

    // Form 2ndAdjs
    RCP<sparse_matrix_type> secondAdjs =
      rcp (new sparse_matrix_type(adjsMatrix->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*adjsMatrix,false,*adjsMatrix,
                                     true,*secondAdjs);
    Array<gno_t> Indices;
    Array<nonzero_t> Values;

    /* Allocate memory necessary for the adjacency */
    lno_t *start = new lno_t [LocalNumIDs+1];
    std::vector<gno_t> adj;

    for (size_t localElement=0; localElement<LocalNumIDs; ++localElement){
      start[localElement] = nadj;
      const gno_t globalRow = Ids[localElement];
      size_t NumEntries = secondAdjs->getNumEntriesInGlobalRow (globalRow);
      Indices.resize (NumEntries);
      Values.resize (NumEntries);
      secondAdjs->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

      for (size_t j = 0; j < NumEntries; ++j) {
	if(globalRow != Indices[j]) {
	  adj.push_back(Indices[j]);
	  nadj++;;
	}
      }
    }

    Ids = NULL;
    start[LocalNumIDs] = nadj;

    zgid_t *adj_ = new zgid_t [nadj];

    for (size_t i=0; i < nadj; i++) {
      adj_[i] = adj[i];
    }

    offsets = start;
    adjacencyIds = adj_;
  }

  //return nadj;
}

//////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void GraphModel<Adapter>::shared_constructor(
  const Adapter *ia,
  modelFlag_t &modelFlags)
{
  // Model creation flags
  bool consecutiveIdsRequired =
    modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  bool removeSelfEdges = modelFlags.test(SELF_EDGES_MUST_BE_REMOVED);
  bool subsetGraph = modelFlags.test(GRAPH_IS_A_SUBSET_GRAPH);

  // A subset graph is special in that it may contain neighbor
  // vertices that are not owned by processes in this communicator.
  // We remove these.

  ArrayRCP<const int> nborProcs;

  if (subsetGraph){
    RCP<const idmap_t> idMap;
    try{
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
    }
    Z2_FORWARD_EXCEPTIONS;

    ArrayRCP<int> procArray;

    if (numLocalEdges_ > 0){
      int *pids = new int [numLocalEdges_];
      env_->localMemoryAssertion(__FILE__, __LINE__, numLocalEdges_, pids);
      procArray = arcp(pids, 0, numLocalEdges_, true);
    }

    ArrayView<gno_t> dummyGno;

    try{
      // All processes must make this call.
      // procOwner will be -1 if edge Id is not in our communicator.

      idMap->gidGlobalTranslate( edgeGids_.view(0, numLocalEdges_),
        dummyGno, procArray.view(0, numLocalEdges_));
    }
    Z2_FORWARD_EXCEPTIONS;

    int outOfSubset = 0;
    for (size_t i=0; i < numLocalEdges_; i++){
      if (procArray[i] < 0){
        outOfSubset++;
        break;
      }
    }

    if (outOfSubset == 0){
      procArray.clear();
      subsetGraph = false;
    }
    else{
      nborProcs = arcp_const_cast<const int>(procArray);
    }
  }

  // Now remove undesired edges.

  if (subsetGraph || removeSelfEdges){

    ArrayRCP<const zgid_t> newEdges;
    ArrayRCP<const lno_t> newOffs;
    scalar_t **newWeights = NULL;
    size_t numNewEdges = 0;

    ArrayView<const zgid_t> vtxView = gids_.view(0, numLocalVertices_);
    ArrayView<const zgid_t> nborView= edgeGids_.view(0, numLocalEdges_);
    ArrayView<const int> nborOwner = nborProcs.view(0, nborProcs.size());
    ArrayView<input_t> eWgts = eWeights_.view(0, nWeightsPerEdge_);
    ArrayView<const lno_t> offView = offsets_.view(0, numLocalVertices_ + 1);

    try{
      numNewEdges = removeUndesiredEdges<user_t>(env_, comm_->getRank(),
        removeSelfEdges,
        false,
        subsetGraph,
        vtxView,
        nborView,
        nborOwner,
        eWgts,
        offView,
        newEdges,
        newWeights,
        newOffs);
    }
    Z2_FORWARD_EXCEPTIONS;

    nborProcs.clear();

    if (numNewEdges < numLocalEdges_){
      edgeGids_ = newEdges;
      offsets_ = newOffs;
      numLocalEdges_ = numNewEdges;

      for (int w=0; w < nWeightsPerEdge_; w++){
        ArrayRCP<const scalar_t> wgtArray(newWeights[w], 0, numNewEdges, true);
        eWeights_[w] = input_t(wgtArray, 1);
      }
    }
    delete [] newWeights;
  }

  // Create an IdentifierMap, which maps the user's global IDs to
  //   Zoltan2 internal global numbers if necessary.
  //   The map can also give us owners of our vertex neighbors.

  RCP<const idmap_t> idMap;

  try{
    idMap = rcp(new idmap_t(env_, comm_, gids_, consecutiveIdsRequired));
  }
  Z2_FORWARD_EXCEPTIONS;

  // Model base class needs to have IdentifierMap.

  this->setIdentifierMap(idMap);

  numGlobalVertices_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  // Compute internal global numbers if we can not use the
  // user's global Ids.  Also find the process owning each
  // neighboring vertex.

  ArrayView<const zgid_t> gidArray(Teuchos::null);  // edge gid
  ArrayView<gno_t> gnoArray(Teuchos::null);        // edge gno
  ArrayView<int> procArray(Teuchos::null);         // edge owner

  if (numLocalVertices_){

    if (!gnosAreGids_){   // need vertex global numbers, edge global numbers
      gno_t *tmp = new gno_t [numLocalVertices_];
      env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVertices_, tmp);
      gnos_ = arcp(tmp, 0, numLocalVertices_);

      try{
        ArrayRCP<zgid_t> tmpGids = arcp_const_cast<zgid_t>(gids_);

        idMap->gidTranslate(tmpGids(0,numLocalVertices_),
          gnos_(0,numLocalVertices_), TRANSLATE_APP_TO_LIB);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (numLocalEdges_){
        tmp = new gno_t [numLocalEdges_];
        env_->localMemoryAssertion(__FILE__, __LINE__, numLocalEdges_, tmp);
        edgeGnos_ = arcp(tmp, 0, numLocalEdges_);
        gnoArray = edgeGnos_.view(0, numLocalEdges_);
      }
    }

    if (numLocalEdges_){
      gidArray = edgeGids_.view(0, numLocalEdges_);

      int *p = new int [numLocalEdges_];
      env_->localMemoryAssertion(__FILE__, __LINE__, numLocalEdges_, p);
      procIds_ = arcp(p, 0, numLocalEdges_);
      procArray = procIds_.view(0, numLocalEdges_);
    }
  }

  try{
    // All processes must make this call.
    idMap->gidGlobalTranslate(gidArray, gnoArray, procArray);
  }
  Z2_FORWARD_EXCEPTIONS;

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);
  edgeGnosConst_ = arcp_const_cast<const gno_t>(edgeGnos_);
  procIdsConst_ = arcp_const_cast<const int>(procIds_);

  // Number of edges such that neighbor is on the local process.
  // We only create the list of local graph edges if the user
  // calls getLocalEdgeList().

  numLocalGraphEdges_ = 0;
  int *pids = procArray.getRawPtr();
  int me = comm_->getRank();
  for (size_t i=0; i < numLocalEdges_; i++)
    if (pids[i] == me) numLocalGraphEdges_++;

  // Vertex weights

  numWeightsPerVertex_ = ia->getNumWeightsPerID();

  if (numWeightsPerVertex_ > 0){
    input_t *weightInfo = new input_t [numWeightsPerVertex_];
    env_->localMemoryAssertion(__FILE__, __LINE__, numWeightsPerVertex_,
                               weightInfo);

    for (int idx=0; idx < numWeightsPerVertex_; idx++){
      bool useNumNZ = ia->useDegreeAsWeight(idx);
      if (useNumNZ){
        scalar_t *wgts = new scalar_t [numLocalVertices_];
        env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVertices_, wgts);
        ArrayRCP<const scalar_t> wgtArray =
          arcp(wgts, 0, numLocalVertices_, true);
        for (size_t i=0; i < numLocalVertices_; i++){
          wgts[i] = offsets_[i+1] - offsets_[i];
        }
        weightInfo[idx] = input_t(wgtArray, 1);
      }
      else{
        const scalar_t *weights=NULL;
        int stride=0;
        ia->getWeightsView(weights, stride, idx);
        ArrayRCP<const scalar_t> wgtArray = arcp(weights, 0,
                                                 stride*numLocalVertices_,
                                                 false);
        weightInfo[idx] = input_t(wgtArray, stride);
      }
    }

    vWeights_ = arcp<input_t>(weightInfo, 0, numWeightsPerVertex_, true);
  }


  reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
    &numLocalEdges_, &numGlobalEdges_);

  env_->memory("After construction of graph model");
}

//////////////////////////////////////////////////////////////////////////

template <typename Adapter>
template <typename AdapterWithCoords>
void GraphModel<Adapter>::shared_GetVertexCoords(const AdapterWithCoords *ia)
{
  // get Vertex coordinates from input adapter

  vCoordDim_ = ia->getDimension();

  if (vCoordDim_ > 0){
    input_t *coordInfo = new input_t [vCoordDim_];
    env_->localMemoryAssertion(__FILE__, __LINE__, vCoordDim_, coordInfo);

    for (int dim=0; dim < vCoordDim_; dim++){
      const scalar_t *coords=NULL;
      int stride=0;
      ia->getCoordinatesView(coords, stride, dim);
      ArrayRCP<const scalar_t> coordArray = arcp(coords, 0,
                                                 stride*numLocalVertices_,
                                                 false);
      coordInfo[dim] = input_t(coordArray, stride);
    }

    vCoords_ = arcp<input_t>(coordInfo, 0, vCoordDim_, true);
  }
}

//////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
void GraphModel<Adapter>::print()
{
  if (env_->getDebugLevel() < VERBOSE_DETAILED_STATUS)
    return;

  std::ostream *os = env_->getDebugOStream();
  
  int me = comm_->getRank();
  std::string fn(" ");

  *os << me << fn
      << " Nvtx  " << gids_.size()
      << " Nedge " << edgeGids_.size()
      << " NLocalEdge " << numLocalGraphEdges_
      << " NVWgt " << numWeightsPerVertex_
      << " NEWgt " << nWeightsPerEdge_
      << " CDim  " << vCoordDim_
      << " GnosAreGids " << gnosAreGids_ << std::endl;

  for (lno_t i = 0; i < gids_.size(); i++) {
    *os << me << fn << i << " GID " << gids_[i] << ": ";
    for (lno_t j = offsets_[i]; j < offsets_[i+1]; j++)
      *os << edgeGids_[j] << " " << "(" << procIds_[j] << ") ";
    *os << std::endl;
  }

  if (gnos_.size())
    for (lno_t i = 0; i < gnos_.size(); i++) {
      *os << me << fn << i << " GNO " << gnos_[i] << ": ";
      for (lno_t j = offsets_[i]; j < offsets_[i+1]; j++)
        *os << edgeGnos_[j] << " ";//<< "(" << procIds_[j] << ") ";
      *os << std::endl;
    }
  else
    *os << me << fn << " GNOS NOT AVAILABLE " << std::endl;

  if (comm_->getSize() > 1) {
    // Print local graph, with no off-process edges.
    ArrayView<const lno_t> localEdgeIds;
    ArrayView<const lno_t> localOffsets;
    ArrayView<input_t> localWgts;
    this->getLocalEdgeList(localEdgeIds, localOffsets, localWgts);

    for (lno_t i = 0; i < gids_.size(); i++) {
      *os << me << fn << i << " LGNO " << gids_[i] << ": ";
      for (lno_t j = localOffsets[i]; j < localOffsets[i+1]; j++) 
        *os << localEdgeIds[j] << " ";
      *os << std::endl;
    }
  }
  else
    *os << me << fn 
       << " LOCAL GRAPH IS SAME AS GLOBAL GRAPH ON ONE RANK " << std::endl;

  if (vCoordDim_) {
    for (lno_t i = 0; i < gids_.size(); i++) {
      *os << me << fn << i << " COORDS " << gids_[i] << ": ";
      for (int j = 0; j < vCoordDim_; j++)
         *os << vCoords_[j][i] << " ";
      *os << std::endl;
    }
  }
  else
    *os << me << fn << "NO COORDINATES AVAIL " << std::endl;
}

}   // namespace Zoltan2


#endif

