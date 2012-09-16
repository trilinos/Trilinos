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
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_GraphInput.hpp>
#include <Zoltan2_IdentifierInput.hpp>
#include <Zoltan2_CoordinateInput.hpp>
#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_StridedData.hpp>

#include <vector>
#include <Teuchos_Hashtable.hpp>

namespace Zoltan2 {


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
 *  \param newWeights if \c wdim is the edge weight dimension,
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
 *  The template parameter is an InputAdapter type.
 */

template <typename User> size_t removeUndesiredEdges(
  const RCP<const Environment> &env,
  int myRank,
  bool removeSelfEdges,
  bool removeOffProcessEdges,
  bool removeOffGroupEdges,
  ArrayView<const typename InputTraits<User>::gid_t> &gids,
  ArrayView<const typename InputTraits<User>::gid_t> &gidNbors,
  ArrayView<const int> &procIds,
  ArrayView<StridedData<typename InputTraits<User>::lno_t,
    typename InputTraits<User>::scalar_t> > &edgeWeights,
  ArrayView<const typename InputTraits<User>::lno_t> &offsets,
  ArrayRCP<const typename InputTraits<User>::gid_t> &newGidNbors, // out
  typename InputTraits<User>::scalar_t **&newWeights,       // out
  ArrayRCP<const typename InputTraits<User>::lno_t> &newOffsets)  // out
{
  typedef typename InputTraits<User>::gid_t gid_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  size_t numKeep = 0;

  size_t numVtx = offsets.size() - 1;
  size_t numNbors = gidNbors.size();

  env->localInputAssertion(__FILE__, __LINE__, "need more input",
    (!removeSelfEdges ||
      gids.size() >=
       static_cast<typename ArrayView<const gid_t>::size_type>(numVtx))
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
  std::vector<bool> uniformWeight;
  if (eDim > 0){
    for (int i=0; i < eDim; i++)
      uniformWeight[i] = (edgeWeights[i].size() == 0);
  }

  // count desired edges

  lno_t *offs = new lno_t [numVtx + 1];
  env->localMemoryAssertion(__FILE__, __LINE__, numVtx+1, offs);
  ArrayRCP<const lno_t> offArray = arcp(offs, 0, numVtx+1, true);

  const lno_t *allOffs = offsets.getRawPtr();
  const gid_t *allIds = gidNbors.getRawPtr();

  const gid_t *vtx = NULL;
  const int *proc = NULL;

  if (gids.size() > 0)
    vtx = gids.getRawPtr();

  if (procIds.size() > 0)
    proc = procIds.getRawPtr();

  offs[0] = 0;
  for (size_t i=0; i < numVtx; i++){
    offs[i+1] = 0;
    int vid = vtx ? vtx[i] : 0;
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
    newGidNbors = ArrayRCP<const gid_t>(Teuchos::null);
    newOffsets = ArrayRCP<const lno_t>(Teuchos::null);
    return 0;
  }

  // Build the subset neighbor lists (id, weight, and offset).

  gid_t *newGids = new gid_t [numKeep];
  env->localMemoryAssertion(__FILE__, __LINE__, numKeep, newGids);

  newGidNbors = arcp(newGids, 0, numKeep, true);
  newOffsets = offArray;

  if (eDim > 0){
    newWeights = new scalar_t * [eDim];
    env->localMemoryAssertion(__FILE__, __LINE__, eDim, newWeights);

    for (int w=0; w < eDim; w++){
      if (uniformWeight[w])
        newWeights[w] = NULL;  // implies uniform
      else{
        newWeights[w] = new scalar_t [numKeep];
        env->localMemoryAssertion(__FILE__, __LINE__, numKeep, newWeights[w]);
      }
    }
  }

  size_t next = 0;
  for (size_t i=0; i < numVtx && next < numKeep; i++){
    int vid = vtx ? vtx[i] : 0;
    for (lno_t j=allOffs[i]; j < allOffs[i+1]; j++){
      int owner = proc ? proc[j] : 0;
      bool keep = (!removeSelfEdges || vid != allIds[j]) &&
               (!removeOffProcessEdges || owner == myRank) &&
               (!removeOffGroupEdges || owner >= 0);

      if (keep){
        newGids[next] = allIds[j];
        if (eDim > 0){
          for (int w=0; w < eDim; w++){
            if (!uniformWeight[w])
              newWeights[w][next] = edgeWeights[w][j];
          }
        }
        next++;
        if (next == numKeep)
          break;

      }  // if (keep)
    }
  }

  return numKeep;
}

/*! \brief Helper function to create new edges lists containing
     only edges connected to a neighbor on this process.
 */

template <typename User> size_t computeLocalEdgeList(
  const RCP<const Environment> &env, int rank,
  size_t numLocalEdges,           // local edges
  size_t numLocalGraphEdges,      // edges in "local" graph
  RCP<const IdentifierMap<User> > &idMap,
  ArrayRCP<const typename InputTraits<User>::gid_t> &allEdgeIds, // in
  ArrayRCP<int> &allProcs,                                 // in
  ArrayRCP<const typename InputTraits<User>::lno_t> &allOffs,    // in
  ArrayRCP<StridedData<typename InputTraits<User>::lno_t,
    typename InputTraits<User>::scalar_t> > &allWeights,         // in
  ArrayRCP<const typename InputTraits<User>::lno_t> &edgeLocalIds, //
  ArrayRCP<const typename InputTraits<User>::lno_t> &offsets,      // out
  ArrayRCP<StridedData<typename InputTraits<User>::lno_t,
    typename InputTraits<User>::scalar_t> > &eWeights)             // out
{
  typedef typename InputTraits<User>::gid_t gid_t;
  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  bool gnosAreGids = idMap->gnosAreGids();

  edgeLocalIds = ArrayRCP<const lno_t>(Teuchos::null);
  offsets = ArrayRCP<const lno_t>(Teuchos::null);
  eWeights = ArrayRCP<input_t>(Teuchos::null);

  if (numLocalGraphEdges == 0)
    return 0;

  if (numLocalGraphEdges == numLocalEdges){

    // Entire graph is local.

    lno_t *lnos = new lno_t [numLocalEdges];
    env->localMemoryAssertion(__FILE__, __LINE__,numLocalEdges, lnos);
    for (size_t i=0; i < numLocalEdges; i++)
      lnos[i] = i;
    edgeLocalIds = arcp(lnos, 0, numLocalEdges, true);
    offsets = allOffs;
    eWeights = allWeights;
  }
  else{

    // Create subset list of local graph edges, offsets and weights.

    int eWeightDim = allWeights.size();

    ArrayRCP<const gid_t> newEgids;
    scalar_t **newWeights = NULL;

    ArrayView<const gid_t> dummyVtx;
    ArrayView<const gid_t> nborView= allEdgeIds.view(0, numLocalEdges);
    ArrayView<const int> nborOwner = allProcs.view(0, numLocalEdges);
    ArrayView<input_t> eWgts = allWeights.view(0, eWeightDim);
    ArrayView<const lno_t> offView = allOffs.view(0, allOffs.size());

    try{
      numLocalEdges = removeUndesiredEdges<User>(
        env, rank,
        false, true, false,
        dummyVtx,
        nborView,
        nborOwner,
        eWgts,
        offView,
        newEgids,
        newWeights,
        offsets);
    }
    Z2_FORWARD_EXCEPTIONS;

    env->localBugAssertion(__FILE__, __LINE__, "local graph miscalculation",
      numLocalEdges == numLocalGraphEdges, BASIC_ASSERTION);

    // offsets array was set by removeUndesiredEdges.  Create weight array.

    if (eWeightDim > 0){
      input_t *wgts = new input_t [eWeightDim];
      for (int w=0; w < eWeightDim; w++){
        if (newWeights[w]){
          ArrayRCP<const scalar_t> wgtArray(
            newWeights[w], 0, numLocalGraphEdges, true);
          wgts[w] = input_t(wgtArray, 1);
        }
      }
      eWeights = arcp(wgts, 0, eWeightDim);
    }

    // Create local ID array.  First translate gid to gno.

    ArrayRCP<gno_t> gnoArray;

    if (gnosAreGids){
      ArrayRCP<const gno_t> gnosConst =
        arcp_reinterpret_cast<const gno_t>(newEgids);
      gnoArray = arcp_const_cast<gno_t>(gnosConst);
    }
    else{

      ArrayRCP<gid_t> gidArray = arcp_const_cast<gid_t>(newEgids);
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
  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
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
  size_t getGlobalNumVertices() const { return 0; }

  /*! \brief Returns the number of global edges on this process.
   *  Includes remote edges.
   */
  size_t getLocalNumGlobalEdges() const { return 0; }

  /*! \brief Returns the number of local edges on this process.
   *  Does not include edges to off-process vertices.
   */
  size_t getLocalNumLocalEdges() const { return 0; }

  /*! \brief Returns the global number edges.
   */
  size_t getGlobalNumEdges() const { return 0; }

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
      \param xyz If vertex coordinate data is available, \c xyz
         will on return point to a StridedData object of coordinates.
      \param wgts If vertex weights is available, \c wgts
         will on return point to a StridedData object of weights.
   */

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const {return 0;}

  /*! \brief Sets pointers to this process' edge (neighbor) global Ids, including
      off-process edges.

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
    ArrayView<input_t> &wgts) const {return 0;}

  /*! \brief Sets pointers to this process' local-only edge (neighbor) LNOs, using
      the same implied vertex LNOs returned in getVertexList.

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

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) {return 0;}

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const { return 0;}

  size_t getGlobalNumObjects() const { return 0;}

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

  typedef typename InputTraits<User>::scalar_t  scalar_t;
  typedef typename InputTraits<User>::gno_t     gno_t;
  typedef typename InputTraits<User>::lno_t     lno_t;
  typedef typename InputTraits<User>::gid_t     gid_t;
  typedef StridedData<lno_t, scalar_t> input_t;
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
    modelFlag_t &modelFlags);

  //!  Destructor
  ~GraphModel() { }

  // // // // // // // // // // // // // // // // // // // // // /
  // The GraphModel interface.
  // // // // // // // // // // // // // // // // // // // // // /

  size_t getLocalNumVertices() const { return numLocalVertices_; }

  size_t getGlobalNumVertices() const { return numGlobalVertices_; }

  size_t getLocalNumGlobalEdges() const { return numLocalEdges_; }

  size_t getLocalNumLocalEdges() const { return numLocalGraphEdges_; }

  size_t getGlobalNumEdges() const { return numGlobalEdges_; }

  int getVertexWeightDim() const { return vWeightDim_; }

  int getEdgeWeightDim() const { return eWeightDim_; }

  int getCoordinateDim() const { return vCoordDim_; }

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz, ArrayView<input_t> &wgts) const;

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) const;

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts){

    if (localGraphEdgeLnos_.size() <
        static_cast<typename ArrayRCP<const lno_t>::size_type>(numLocalGraphEdges_)){

      RCP<const IdentifierMap<User> > idmap = this->getIdentifierMap();

      computeLocalEdgeList<User>(env_, comm_->getRank(),
        numLocalEdges_, numLocalGraphEdges_,
        idmap, edgeGids_, procIds_, offsets_, eWeights_,
        localGraphEdgeLnos_, localGraphEdgeOffsets_, localGraphEdgeWeights_);
    }

    edgeIds = localGraphEdgeLnos_();
    offsets = localGraphEdgeOffsets_();
    wgts = localGraphEdgeWeights_();

    return numLocalGraphEdges_;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const { return numLocalVertices_; }

  size_t getGlobalNumObjects() const { return numGlobalVertices_; }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const
  {
    ArrayView<input_t> xyz, wgts;
    getVertexList(gnos, xyz, wgts);
  }

private:

  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  ArrayRCP<const gid_t> gids_;        // rows
  ArrayRCP<gno_t> gnos_;

  int vWeightDim_;
  ArrayRCP<input_t> vWeights_;

  int vCoordDim_;
  ArrayRCP<input_t> vCoords_;

  // Note: in case of graph subsetting, size of these arrays
  // may be larger than numLocalEdges_.  So do not use .size().

  ArrayRCP<const gid_t> edgeGids_;
  ArrayRCP<gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  int eWeightDim_;
  ArrayRCP<input_t> eWeights_;

  ArrayRCP<const gno_t> gnosConst_;
  ArrayRCP<const gno_t> edgeGnosConst_;
  ArrayRCP<const int> procIdsConst_;

  bool gnosAreGids_;

  // For local graphs (graph restricted to local process).  We only
  // create these arrays if required by the algorithm.

  ArrayRCP<const lno_t> localGraphEdgeLnos_;
  ArrayRCP<const lno_t> localGraphEdgeOffsets_;
  ArrayRCP<input_t> localGraphEdgeWeights_;

  // For convenience

  size_t numLocalVertices_;
  size_t numGlobalVertices_;
  size_t numLocalEdges_;
  size_t numGlobalEdges_;
  size_t numLocalGraphEdges_;
};

template <typename User>
  GraphModel<MatrixInput<User> >::GraphModel(const MatrixInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags):
     env_(env), comm_(comm),
     gids_(), gnos_(),
     vWeightDim_(0), vWeights_(),
     vCoordDim_(0), vCoords_(),
     edgeGids_(), edgeGnos_(), procIds_(), offsets_(),
     eWeightDim_(0), eWeights_(),
     gnosConst_(), edgeGnosConst_(), procIdsConst_(),
     gnosAreGids_(false),
     localGraphEdgeLnos_(), localGraphEdgeOffsets_(), localGraphEdgeWeights_(),
     numLocalVertices_(0), numGlobalVertices_(0), numLocalEdges_(0),
     numGlobalEdges_(0), numLocalGraphEdges_(0)
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
    numLocalVertices_ = ia->getRowListView(vtxIds, offsets, nborIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  numLocalEdges_ = offsets[numLocalVertices_];

  gids_ = arcp<const gid_t>(vtxIds, 0, numLocalVertices_, false);
  edgeGids_ = arcp<const gid_t>(nborIds, 0, numLocalEdges_, false);
  offsets_ = arcp<const lno_t>(offsets, 0, numLocalVertices_ + 1, false);

  eWeightDim_ = 0;   // no edge weights from a matrix

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

    ArrayRCP<input_t> noEdgeWeights;
    ArrayRCP<const gid_t> newEdges;
    ArrayRCP<const lno_t> newOffs;
    scalar_t **newWeights = NULL;
    size_t numNewEdges = 0;

    // Compiler complained of an error if gids_.view(0, n), etc
    // was listed directly as a parameter in removeUndesiredEdges.
    // So we have to create the ArraView before before the call.

    ArrayView<const gid_t> vtxView= gids_.view(0, numLocalVertices_);
    ArrayView<const gid_t> nborView= edgeGids_.view(0, numLocalEdges_);
    ArrayView<const int> nborOwner = nborProcs.view(0, nborProcs.size());
    ArrayView<input_t> eWgts = noEdgeWeights.view(0, 0);
    ArrayView<const lno_t> offView = offsets_.view(0, numLocalVertices_ + 1);

    try{
      numNewEdges = Zoltan2::removeUndesiredEdges<User>(env_, comm_->getRank(),
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
    }
  }

  // Create an IdentifierMap, which will map the user's global IDs to
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

  ArrayView<const gid_t> gidArray(Teuchos::null);  // edge gid
  ArrayView<gno_t> gnoArray(Teuchos::null);        // edge gno
  ArrayView<int> procArray(Teuchos::null);         // edge owner

  if (numLocalVertices_){

    if (!gnosAreGids_){   // need vertex global numbers, edge global numbers
      gno_t *tmp = new gno_t [numLocalVertices_];
      env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVertices_, tmp);
      gnos_ = arcp(tmp, 0, numLocalVertices_);

      try{
        ArrayRCP<gid_t> tmpGids = arcp_const_cast<gid_t>(gids_);

        idMap->gidTranslate(tmpGids(0,numLocalVertices_),
          gnos_(0,numLocalVertices_), TRANSLATE_APP_TO_LIB);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (numLocalEdges_){
        gno_t *tmp = new gno_t [numLocalEdges_];
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
  for (size_t i=0; i < numLocalEdges_; i++)
    if (pids[i] == comm_->getRank())
      numLocalGraphEdges_++;

  // Vertex weights

  vWeightDim_ = ia->getRowWeightDimension();

  if (vWeightDim_ > 0){
    input_t *weightInfo = new input_t [vWeightDim_];
    env_->localMemoryAssertion(__FILE__, __LINE__, vWeightDim_, weightInfo);

    for (int dim=0; dim < vWeightDim_; dim++){
      bool useNumNZ = ia->getRowWeightIsNumberOfNonZeros(dim);
      if (useNumNZ){
        scalar_t *wgts = new scalar_t [numLocalVertices_];
        env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVertices_, wgts);
        ArrayRCP<const scalar_t> wgtArray =
          arcp(wgts, 0, numLocalVertices_, true);
        for (size_t i=0; i < numLocalVertices_; i++){
          wgts[i] = offsets_[i+1] - offsets_[i];
        }
        weightInfo[dim] = input_t(wgtArray, 1);
      }
      else{
        const scalar_t *weights=NULL;
        int stride=0;
        size_t len = ia->getRowWeights(dim, weights, stride);
        // If weights is NULL, user wants to use uniform weights
        if (weights != NULL){
          ArrayRCP<const scalar_t> wgtArray = arcp(weights, 0, len, false);
          weightInfo[dim] = input_t(wgtArray, stride);
        }
      }
    }

    vWeights_ = arcp<input_t>(weightInfo, 0, vWeightDim_, true);
  }

  // Model base class needs to know if any weights are uniform.

  Array<lno_t> weightArrayLengths(vWeightDim_);
  for (int dim=0; dim < vWeightDim_; dim++){
    weightArrayLengths[dim] = vWeights_[dim].size();
  }
  this->setWeightArrayLengths(weightArrayLengths, *comm_);

  // Vertex coordinates

  vCoordDim_ = ia->getCoordinateDimension();

  if (vCoordDim_ > 0){
    input_t *coordInfo = new input_t [vCoordDim_];
    env_->localMemoryAssertion(__FILE__, __LINE__, vCoordDim_, coordInfo);

    for (int dim=0; dim < vCoordDim_; dim++){
      const scalar_t *coords=NULL;
      int stride=0;
      size_t len = ia->getRowCoordinates(dim, coords, stride);
      ArrayRCP<const scalar_t> coordArray = arcp(coords, 0, len, false);
      coordInfo[dim] = input_t(coordArray, stride);
    }

    vCoords_ = arcp<input_t>(coordInfo, 0, vCoordDim_, true);
  }

  reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
    &numLocalEdges_, &numGlobalEdges_);

  env_->memory("After construction of graph model");
}

template <typename User>
  size_t GraphModel<MatrixInput<User> >::getVertexList(
    ArrayView<const gno_t> &Ids, ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const
  {
    size_t nv = gids_.size();

    if (gnosAreGids_)
      Ids = gids_.view(0, nv);
    else
      Ids = gnosConst_.view(0, nv);

    xyz = vCoords_.view(0, vWeightDim_);
    wgts = vWeights_.view(0, vCoordDim_);

    return nv;
  }

template <typename User>
  size_t GraphModel<MatrixInput<User> >::getEdgeList(
    ArrayView<const gno_t> &edgeIds, ArrayView<const int> &procIds,
    ArrayView<const lno_t> &offsets, ArrayView<input_t> &wgts) const
{
  if (gnosAreGids_)
    edgeIds = edgeGids_.view(0, numLocalEdges_);
  else
    edgeIds = edgeGnosConst_.view(0, numLocalEdges_);

  procIds = procIdsConst_.view(0, numLocalEdges_);
  offsets = offsets_.view(0, numLocalVertices_+1);
  wgts = eWeights_.view(0, eWeightDim_);

  return numLocalEdges_;
}

////////////////////////////////////////////////////////////////
// Graph model derived from GraphInput.
////////////////////////////////////////////////////////////////

template <typename User>
class GraphModel<GraphInput<User> > : public Model<GraphInput<User> >
{
public:

  typedef typename GraphInput<User>::scalar_t  scalar_t;
  typedef typename GraphInput<User>::gno_t     gno_t;
  typedef typename GraphInput<User>::lno_t     lno_t;
  typedef typename GraphInput<User>::gid_t     gid_t;
  typedef typename GraphInput<User>::node_t    node_t;
  typedef IdentifierMap<User>     idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! Constructor
   *  All processes in the communicator must call the constructor.
   *
   *  \param  inputAdapter  an encapsulation of the user data
   *  \param  env           environment (library configuration settings)
   *  \param  comm       communicator for the problem
   *  \param  modelFlags  a bit map of Zoltan2::GraphModelFlags
   */

  GraphModel(const GraphInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags);

  //!  Destructor
  ~GraphModel() { }

  // // // // // // // // // // // // // // // // // // // // // /
  // The GraphModel interface.
  // // // // // // // // // // // // // // // // // // // // // /

  size_t getLocalNumVertices() const { return numLocalVertices_; }

  size_t getGlobalNumVertices() const { return numGlobalVertices_; }

  size_t getLocalNumGlobalEdges() const { return numLocalEdges_; }

  size_t getLocalNumLocalEdges() const { return numLocalGraphEdges_; }

  size_t getGlobalNumEdges() const { return numGlobalEdges_; }

  int getVertexWeightDim() const { return vWeightDim_; }

  int getEdgeWeightDim() const { return eWeightDim_; }

  int getCoordinateDim() const { return vCoordDim_; }

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz, ArrayView<input_t> &wgts) const;

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) const;

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts){

    if (localGraphEdgeLnos_.size() < numLocalGraphEdges_){

      RCP<const IdentifierMap<User> > idmap = this->getIdentifierMap();

      computeLocalEdgeList(env_, comm_->getRank(),
        numLocalEdges_, numLocalGraphEdges_,
        idmap, edgeGids_, procIds_, offsets_, eWeights_,
        localGraphEdgeLnos_, localGraphEdgeOffsets_, localGraphEdgeWeights_);
    }

    edgeIds = localGraphEdgeLnos_;
    offsets = localGraphEdgeOffsets_;
    wgts = localGraphEdgeWeights_;

    return numLocalGraphEdges_;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const { return numLocalVertices_; }

  size_t getGlobalNumObjects() const { return numGlobalVertices_; }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const
  {
    ArrayView<input_t> xyz, wgts;
    getVertexList(gnos, xyz, wgts);
  }

private:

  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  ArrayRCP<const gid_t> gids_;        // vertices of input graph
  ArrayRCP<gno_t> gnos_;

  int vWeightDim_;
  ArrayRCP<input_t> vWeights_;

  int vCoordDim_;
  ArrayRCP<input_t> vCoords_;

  // Note: in case of graph subsetting, size of these arrays
  // may be larger than numLocalEdges_.  So do not use .size().

  ArrayRCP<const gid_t> edgeGids_;
  ArrayRCP<gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  int eWeightDim_;
  ArrayRCP<input_t> eWeights_;

  ArrayRCP<const gno_t> gnosConst_;
  ArrayRCP<const gno_t> edgeGnosConst_;
  ArrayRCP<const int> procIdsConst_;

  bool gnosAreGids_;

  // For local graphs (graph restricted to local process).  We only
  // create these arrays if required by the algorithm.

  ArrayRCP<const lno_t> localGraphEdgeLnos_;
  ArrayRCP<const lno_t> localGraphEdgeOffsets_;
  ArrayRCP<input_t> localGraphEdgeWeights_;

  // For convenience

  size_t numLocalVertices_;
  size_t numGlobalVertices_;
  size_t numLocalEdges_;
  size_t numGlobalEdges_;
  size_t numLocalGraphEdges_;

};

template <typename User>
  GraphModel<GraphInput<User> >::GraphModel(const GraphInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags):
     env_(env), comm_(comm),
     gids_(), gnos_(),
     vWeightDim_(0), vWeights_(),
     vCoordDim_(0), vCoords_(),
     edgeGids_(), edgeGnos_(), procIds_(), offsets_(),
     eWeightDim_(0), eWeights_(),
     gnosConst_(), edgeGnosConst_(), procIdsConst_(),
     gnosAreGids_(false),
     localGraphEdgeLnos_(), localGraphEdgeOffsets_(), localGraphEdgeWeights_(),
     numLocalVertices_(0), numGlobalVertices_(0), numLocalEdges_(0),
     numGlobalEdges_(0), numLocalGraphEdges_(0)
{
  // Model creation flags

  bool consecutiveIdsRequired =
    modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  bool removeSelfEdges = modelFlags.test(SELF_EDGES_MUST_BE_REMOVED);
  bool subsetGraph = modelFlags.test(GRAPH_IS_A_SUBSET_GRAPH);

  // Get the graph from the input adapter

  gid_t const *vtxIds=NULL, *nborIds=NULL;
  lno_t const  *offsets=NULL;
  try{
    numLocalVertices_ = ia->getVertexListView(vtxIds, offsets, nborIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  numLocalEdges_ = offsets[numLocalVertices_];

  gids_ = arcp<const gid_t>(vtxIds, 0, numLocalVertices_, false);
  edgeGids_ = arcp<const gid_t>(nborIds, 0, numLocalEdges_, false);
  offsets_ = arcp<const lno_t>(offsets, 0, numLocalVertices_ + 1, false);

  eWeightDim_ = ia->getEdgeWeightDimension();

  if (eWeightDim_ > 0){
    input_t *wgts = new input_t [eWeightDim_];
    eWeights_ = arcp(wgts, 0, eWeightDim_, true);
  }

  for (int w=0; w < eWeightDim_; w++){
    const scalar_t *ewgts=NULL;
    int stride=0;

    ia->getEdgeWeights(w, ewgts, stride);

    ArrayRCP<const scalar_t> wgtArray(ewgts, 0, numLocalEdges_, false);
    eWeights_[w] = input_t(wgtArray, stride);
  }

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

    ArrayRCP<const gid_t> newEdges;
    ArrayRCP<const lno_t> newOffs;
    scalar_t **newWeights = NULL;
    size_t numNewEdges = 0;

    ArrayView<const gid_t> vtxView = gids_.view(0, numLocalVertices_);
    ArrayView<const gid_t> nborView= edgeGids_.view(0, numLocalEdges_);
    ArrayView<const int> nborOwner = nborProcs.view(0, nborProcs.size());
    ArrayView<input_t> eWgts = eWeights_.view(0, eWeightDim_);
    ArrayView<const lno_t> offView = offsets_.view(0, numLocalVertices_ + 1);

    try{
      numNewEdges = removeUndesiredEdges<User>(env_, comm_->getRank(),
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

      for (int w=0; w < eWeightDim_; w++){
        if (newWeights[w] != NULL){   // non-uniform weights
          ArrayRCP<const scalar_t> wgtArray(newWeights[w],
            0, numNewEdges, true);
          eWeights_[w] = input_t(wgtArray, 1);
        }
      }
    }
  }

  // Create an IdentifierMap, which will map the user's global IDs to
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

  ArrayView<const gid_t> gidArray(Teuchos::null);  // edge gid
  ArrayView<gno_t> gnoArray(Teuchos::null);        // edge gno
  ArrayView<int> procArray(Teuchos::null);         // edge owner

  if (numLocalVertices_){

    if (!gnosAreGids_){   // need vertex global numbers, edge global numbers
      gno_t *tmp = new gno_t [numLocalVertices_];
      env_->localMemoryAssertion(__FILE__, __LINE__, numLocalVertices_, tmp);
      gnos_ = arcp(tmp, 0, numLocalVertices_);

      try{
        ArrayRCP<gid_t> tmpGids = arcp_const_cast<gid_t>(gids_);

        idMap->gidTranslate(tmpGids(0,numLocalVertices_),
          gnos_(0,numLocalVertices_), TRANSLATE_APP_TO_LIB);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (numLocalEdges_){
        gno_t *tmp = new gno_t [numLocalEdges_];
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
  for (lno_t i=0; i < numLocalEdges_; i++)
    if (pids[i] == comm_->getRank())
      numLocalGraphEdges_++;

  // Vertex weights

  vWeightDim_ = ia->getVertexWeightDimension();

  if (vWeightDim_ > 0){
    input_t *weightInfo = new input_t [vWeightDim_];
    env_->localMemoryAssertion(__FILE__, __LINE__, vWeightDim_, weightInfo);

    for (int dim=0; dim < vWeightDim_; dim++){
      const scalar_t *weights=NULL;
      int stride=0;
      size_t len = ia->getVertexWeights(dim, weights, stride);
      // If weights is NULL, user wants to use uniform weights
      if (weights != NULL){
        ArrayRCP<const scalar_t> wgtArray = arcp(weights, 0, len, false);
        weightInfo[dim] = input_t(wgtArray, stride);
      }
    }

    vWeights_ = arcp<input_t>(weightInfo, 0, vWeightDim_, true);
  }

  // Model base class needs to know if any weights are uniform.

  Array<lno_t> weightArrayLengths(vWeightDim_);
  for (int dim=0; dim < vWeightDim_; dim++){
    weightArrayLengths[dim] = vWeights_[dim].size();
  }
  this->setWeightArrayLengths(weightArrayLengths, *comm_);

  // Vertex coordinates

  vCoordDim_ = ia->getCoordinateDimension();

  if (vCoordDim_ > 0){
    input_t *coordInfo = new input_t [vCoordDim_];
    env_->localMemoryAssertion(__FILE__, __LINE__, vCoordDim_, coordInfo);

    for (int dim=0; dim < vCoordDim_; dim++){
      const scalar_t *coords=NULL;
      int stride=0;
      size_t len = ia->getVertexCoordinates(dim, coords, stride);
      ArrayRCP<const scalar_t> coordArray = arcp(coords, 0, len, false);
      coordInfo[dim] = input_t(coordArray, stride);
    }

    vCoords_ = arcp<input_t>(coordInfo, 0, vCoordDim_, true);
  }

  reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
    &numLocalEdges_, &numGlobalEdges_);

  env_->memory("After construction of graph model");
}

template <typename User>
  size_t GraphModel<GraphInput<User> >::getVertexList(
    ArrayView<const gno_t> &Ids, ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const
  {
    size_t nv = gids_.size();

    if (gnosAreGids_)
      Ids = gids_.view(0, nv);
    else
      Ids = gnosConst_.view(0, nv);

    xyz = vCoords_.view(0, vWeightDim_);
    wgts = vWeights_.view(0, vCoordDim_);

    return nv;
  }

template <typename User>
  size_t GraphModel<GraphInput<User> >::getEdgeList(
    ArrayView<const gno_t> &edgeIds, ArrayView<const int> &procIds,
    ArrayView<const lno_t> &offsets, ArrayView<input_t> &wgts) const
{
  if (gnosAreGids_)
    edgeIds = edgeGids_.view(0, numLocalEdges_);
  else
    edgeIds = edgeGnosConst_.view(0, numLocalEdges_);

  procIds = procIdsConst_.view(0, numLocalEdges_);
  offsets = offsets_.view(0, numLocalVertices_+1);
  wgts = eWeights_.view(0, eWeightDim_);

  return numLocalEdges_;
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
  typedef StridedData<lno_t, scalar_t> input_t;

  GraphModel(const CoordinateInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags)
  {
    throw std::runtime_error("may not build a graph with identifiers");
  }

  // GraphModel interface

  size_t getLocalNumVertices() const { return 0;}
  size_t getGlobalNumVertices() const { return 0;}
  size_t getLocalNumGlobalEdges() const { return 0;}
  size_t getLocalNumLocalEdges() const { return 0;}
  size_t getGlobalNumEdges() const {return 0;}
  int getVertexWeightDim() const { return 0; }
  int getEdgeWeightDim() const { return 0; }
  int getCoordinateDim() const { return 0; }
  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const { return 0; }
  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) const { return 0; }
  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) { return 0; }

  // Model interface

  size_t getLocalNumObjects() const { return 0; }
  size_t getGlobalNumObjects() const { return 0; }
  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {}

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
  typedef StridedData<lno_t, scalar_t> input_t;

  GraphModel(const VectorInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags)
  {
    throw std::runtime_error("can not build a graph from a vector");
  }

  // GraphModel interface

  size_t getLocalNumVertices() const { return 0;}
  size_t getGlobalNumVertices() const { return 0;}
  size_t getLocalNumGlobalEdges() const { return 0;}
  size_t getLocalNumLocalEdges() const { return 0;}
  size_t getGlobalNumEdges() const {return 0;}
  int getVertexWeightDim() const { return 0; }
  int getEdgeWeightDim() const { return 0; }
  int getCoordinateDim() const { return 0; }
  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const { return 0; }
  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) const { return 0; }
  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) { return 0; }

  // Model interface

  size_t getLocalNumObjects() const { return 0; }
  size_t getGlobalNumObjects() const { return 0; }
  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {}

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
  typedef StridedData<lno_t, scalar_t> input_t;

  GraphModel(const IdentifierInput<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags)
  {
    throw std::runtime_error("can not build a graph with identifiers");
  }

  // GraphModel interface

  size_t getLocalNumVertices() const { return 0;}
  size_t getGlobalNumVertices() const { return 0;}
  size_t getLocalNumGlobalEdges() const { return 0;}
  size_t getLocalNumLocalEdges() const { return 0;}
  size_t getGlobalNumEdges() const {return 0;}
  int getVertexWeightDim() const { return 0; }
  int getEdgeWeightDim() const { return 0; }
  int getCoordinateDim() const { return 0; }
  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const { return 0; }
  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) const { return 0; }
  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<input_t> &wgts) { return 0; }

  // Model interface

  size_t getLocalNumObjects() const { return 0; }
  size_t getGlobalNumObjects() const { return 0; }
  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {}

};

#endif    // DOXYGEN_SHOULD_SKIP_THIS

}   // namespace Zoltan2

#endif

