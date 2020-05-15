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
#include <Zoltan2_ModelHelpers.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <unordered_map>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////
/*!  \brief GraphModel defines the interface required for graph models.

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.

    GraphModels may represent a local (on-process) graph or 
    a global (all-communicator) graph.
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::node_t      node_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef StridedData<lno_t, scalar_t>  input_t;
  typedef typename Adapter::offset_t    offset_t;
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

  GraphModel(const RCP<const MatrixAdapter<user_t,userCoord_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags);

  GraphModel(const RCP<const GraphAdapter<user_t,userCoord_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags);

  GraphModel(const RCP<const MeshAdapter<user_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelflags);

  GraphModel(const RCP<const VectorAdapter<userCoord_t> > &/* ia */,
    const RCP<const Environment> &/* env */, const RCP<const Comm<int> > &/* comm */,
    modelFlag_t &/* flags */)
  {
    throw std::runtime_error("cannot build GraphModel from VectorAdapter");
  }

  GraphModel(const RCP<const IdentifierAdapter<user_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags)
  {
    throw std::runtime_error("cannot build GraphModel from IdentifierAdapter");
  }

  /*! \brief Return the communicator used by the model
   */
  const RCP<const Comm<int> > getComm() { return comm_; }

  /*! \brief Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return nLocalVertices_; }

  /*! \brief Returns the global number vertices.
   */
  size_t getGlobalNumVertices() const { return nGlobalVertices_; }

  /*! \brief Returns the number of edges on this process.
   *  In global or subset graphs, includes off-process edges.
   */
  size_t getLocalNumEdges() const { return nLocalEdges_; }

  /*! \brief Returns the global number edges.
   *  For local graphs, the number of global edges is the number of local edges.
   */
  size_t getGlobalNumEdges() const { return nGlobalEdges_; }

  /*! \brief Returns the number (0 or greater) of weights per vertex
   */
  int getNumWeightsPerVertex() const { return nWeightsPerVertex_; }

  /*! \brief Returns the number (0 or greater) of weights per edge.
   */
  int getNumWeightsPerEdge() const { return nWeightsPerEdge_; }

  /*! \brief Returns the dimension (0 to 3) of vertex coordinates.
   */
  int getCoordinateDim() const { return vCoordDim_; }

  /*! \brief Sets pointers to this process' vertex Ids and their weights.

      \param Ids will on return point to the list of the global Ids for
        each vertex on this process.
      \param wgts If vertex weights is available, \c wgts
         will on return point to a StridedData object of weights.
  */

  size_t getVertexList(
    ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &wgts) const
  {
    Ids = vGids_.view(0, nLocalVertices_);
    wgts = vWeights_.view(0, nWeightsPerVertex_);
    return nLocalVertices_;
  }

  /*! \brief Sets pointers to this process' vertex coordinates, if available.
      Order of coordinate info matches that of Ids in getVertexList().

      \param xyz If vertex coordinate data is available, \c xyz
         will on return point to a StridedData object of coordinates.
  */
  size_t getVertexCoords(ArrayView<input_t> &xyz) const
  {
    xyz = vCoords_.view(0, vCoordDim_);
    return nLocalVertices_;
  }

  /*! \brief Sets pointers to this process' edge (neighbor) global Ids,
      including off-process edges.

      \param edgeIds This is the list of global neighbor Ids corresponding
        to the vertices listed in getVertexList.
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
    ArrayView<const offset_t> &offsets,
    ArrayView<input_t> &wgts) const
  {
    edgeIds = eGids_.view(0, nLocalEdges_);
    offsets = eOffsets_.view(0, nLocalVertices_+1);
    wgts = eWeights_.view(0, nWeightsPerEdge_);
    return nLocalEdges_;
  }

  /*! \brief Return the vtxDist array 
   *  Array of size comm->getSize() + 1
   *  Array[n+1] - Array[n] is number of vertices on rank n
   */
  inline void getVertexDist(ArrayView<size_t> &vtxdist) const
  {
    vtxdist = vtxDist_();
    if (vtxDist_.size() == 0) {
      throw std::runtime_error("getVertexDist is available only "
                               "when consecutiveIdsRequired");
    }
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const { return nLocalVertices_; }

  size_t getGlobalNumObjects() const { return nGlobalVertices_; }

private:

  void shared_constructor(const RCP<const Adapter>&ia, modelFlag_t &modelFlags);

  template <typename AdapterWithCoords>
  void shared_GetVertexCoords(const AdapterWithCoords *ia);

  void print(); // For debugging

  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  bool localGraph_;    // Flag indicating whether this graph is 
                       // LOCAL with respect to the process;
                       // if !localGraph_, graph is GLOBAL with respect to
                       // the communicator.
                       

  size_t nLocalVertices_;                // # local vertices in built graph
  size_t nGlobalVertices_;               // # global vertices in built graph
  ArrayRCP<gno_t> vGids_;                  // vertices of graph built in model;
                                           // may be same as adapter's input
                                           // or may be renumbered 0 to (N-1).

  int nWeightsPerVertex_;
  ArrayRCP<input_t> vWeights_;

  int vCoordDim_;
  ArrayRCP<input_t> vCoords_;

  // Note: in some cases, size of these arrays
  // may be larger than nLocalEdges_.  So do not use .size().
  // Use nLocalEdges_, nGlobalEdges_

  size_t nLocalEdges_;                  // # local edges in built graph
  size_t nGlobalEdges_;                 // # global edges in built graph
  ArrayRCP<gno_t> eGids_;                 // edges of graph built in model
  ArrayRCP<offset_t> eOffsets_;           // edge offsets build in model
                                          // May be same as adapter's input 
                                          // or may differ
                                          // due to renumbering, self-edge
                                          // removal, or local graph.

  int nWeightsPerEdge_;
  ArrayRCP<input_t> eWeights_;            // edge weights in built graph
                                          // May be same as adapter's input
                                          // or may differ due to self-edge
                                          // removal, or local graph.

  ArrayRCP<size_t> vtxDist_;              // If consecutiveIdsRequired,
                                          // vtxDist (as needed by ParMETIS
                                          // and Scotch) is also created.
                                          // Otherwise, it is Teuchos::null.
};


////////////////////////////////////////////////////////////////
// GraphModel from MatrixAdapter
template <typename Adapter>
GraphModel<Adapter>::GraphModel(
  const RCP<const MatrixAdapter<user_t,userCoord_t> > &ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags):
        env_(env),
        comm_(comm),
        localGraph_(false),
        nLocalVertices_(0),
        nGlobalVertices_(0),
        vGids_(),
        nWeightsPerVertex_(0),
        vWeights_(),
        vCoordDim_(0),
        vCoords_(),
        nLocalEdges_(0),
        nGlobalEdges_(0),
        eGids_(),
        eOffsets_(),
        nWeightsPerEdge_(0),
        eWeights_(),
        vtxDist_()
{
  // Model creation flags
  localGraph_ = modelFlags.test(BUILD_LOCAL_GRAPH);

  bool symTranspose = modelFlags.test(SYMMETRIZE_INPUT_TRANSPOSE);
  bool symBipartite = modelFlags.test(SYMMETRIZE_INPUT_BIPARTITE);
  bool vertexCols = modelFlags.test(VERTICES_ARE_MATRIX_COLUMNS);
  bool vertexNz = modelFlags.test(VERTICES_ARE_MATRIX_NONZEROS);

  if (symTranspose || symBipartite || vertexCols || vertexNz){
    throw std::runtime_error("graph build option not yet implemented");
  }

  // Get the matrix from the input adapter
  gno_t const *vtxIds=NULL, *nborIds=NULL;
  offset_t const *offsets=NULL;
  try{
    nLocalVertices_ = ia->getLocalNumIDs();
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

  // Save the pointers from the input adapter
  nLocalEdges_ = offsets[nLocalVertices_];
  vGids_ = arcp_const_cast<gno_t>(
                arcp<const gno_t>(vtxIds, 0, nLocalVertices_, false));
  eGids_ = arcp_const_cast<gno_t>(
                arcp<const gno_t>(nborIds, 0, nLocalEdges_, false));
  eOffsets_ = arcp_const_cast<offset_t>(
                   arcp<const offset_t>(offsets, 0, nLocalVertices_+1, false));

  // Edge weights
  nWeightsPerEdge_ = 0;   // no edge weights from a matrix yet.
                          // TODO:  use matrix values as edge weights

  // Do constructor common to all adapters
  shared_constructor(ia, modelFlags);

  // Get vertex coordinates, if available
  if (ia->coordinatesAvailable()) {
    typedef VectorAdapter<userCoord_t> adapterWithCoords_t;
    shared_GetVertexCoords<adapterWithCoords_t>(ia->getCoordinateInput());
  }
  // print();
}


////////////////////////////////////////////////////////////////
// GraphModel from GraphAdapter
template <typename Adapter>
GraphModel<Adapter>::GraphModel(
  const RCP<const GraphAdapter<user_t,userCoord_t> > &ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags):
        env_(env),
        comm_(comm),
        localGraph_(false),
        nLocalVertices_(0),
        nGlobalVertices_(0),
        vGids_(),
        nWeightsPerVertex_(0),
        vWeights_(),
        vCoordDim_(0),
        vCoords_(),
        nLocalEdges_(0),
        nGlobalEdges_(0),
        eGids_(),
        eOffsets_(),
        nWeightsPerEdge_(0),
        eWeights_(),
        vtxDist_()
{
  // Model creation flags
  localGraph_ = modelFlags.test(BUILD_LOCAL_GRAPH);

  // This GraphModel is built with vertices == GRAPH_VERTEX from GraphAdapter.
  // It is not ready to use vertices == GRAPH_EDGE from GraphAdapter.
  env_->localInputAssertion(__FILE__, __LINE__,
    "GraphModel from GraphAdapter is implemented only for "
    "Graph Vertices as primary object, not for Graph Edges",
    ia->getPrimaryEntityType() == Zoltan2::GRAPH_VERTEX, BASIC_ASSERTION);

  // Get the graph from the input adapter

  gno_t const *vtxIds=NULL, *nborIds=NULL;
  offset_t const *offsets=NULL;
  try{
    nLocalVertices_ = ia->getLocalNumVertices();
    ia->getVertexIDsView(vtxIds);
    ia->getEdgesView(offsets, nborIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Save the pointers from the input adapter
  nLocalEdges_ = offsets[nLocalVertices_];
  vGids_ = arcp_const_cast<gno_t>(
                arcp<const gno_t>(vtxIds, 0, nLocalVertices_, false));
  eGids_ = arcp_const_cast<gno_t>(
                arcp<const gno_t>(nborIds, 0, nLocalEdges_, false));
  eOffsets_ = arcp_const_cast<offset_t>(
                   arcp<const offset_t>(offsets, 0, nLocalVertices_+1, false));

  // Edge weights
  nWeightsPerEdge_ = ia->getNumWeightsPerEdge();
  if (nWeightsPerEdge_ > 0){
    input_t *wgts = new input_t [nWeightsPerEdge_];
    eWeights_ = arcp(wgts, 0, nWeightsPerEdge_, true);
  }

  for (int w=0; w < nWeightsPerEdge_; w++){
    const scalar_t *ewgts=NULL;
    int stride=0;

    ia->getEdgeWeightsView(ewgts, stride, w);

    ArrayRCP<const scalar_t> wgtArray(ewgts, 0, nLocalEdges_, false);
    eWeights_[w] = input_t(wgtArray, stride);
  }

  // Do constructor common to all adapters
  shared_constructor(ia, modelFlags);

  // Get vertex coordinates, if available
  if (ia->coordinatesAvailable()) {
    typedef VectorAdapter<userCoord_t> adapterWithCoords_t;
    shared_GetVertexCoords<adapterWithCoords_t>(ia->getCoordinateInput());
  }
  // print();
}

////////////////////////////////////////////////////////////////
// GraphModel from MeshAdapter
template <typename Adapter>
GraphModel<Adapter>::GraphModel(
  const RCP<const MeshAdapter<user_t> > &ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags):
        env_(env),
        comm_(comm),
        localGraph_(false),
        nLocalVertices_(0),
        nGlobalVertices_(0),
        vGids_(),
        nWeightsPerVertex_(0),
        vWeights_(),
        vCoordDim_(0),
        vCoords_(),
        nLocalEdges_(0),
        nGlobalEdges_(0),
        eGids_(),
        eOffsets_(),
        nWeightsPerEdge_(0),
        eWeights_(),
        vtxDist_()
{
  env_->timerStart(MACRO_TIMERS, "GraphModel constructed from MeshAdapter");

  // Model creation flags
  localGraph_ = modelFlags.test(BUILD_LOCAL_GRAPH);

  // This GraphModel is built with vertices == ia->getPrimaryEntityType()
  // from MeshAdapter.

  // Get the graph from the input adapter

  Zoltan2::MeshEntityType primaryEType = ia->getPrimaryEntityType();
  Zoltan2::MeshEntityType secondAdjEType = ia->getSecondAdjacencyEntityType();

  // Get the IDs of the primary entity type; these are graph vertices

  gno_t const *vtxIds=NULL;
  try {
    nLocalVertices_ = ia->getLocalNumOf(primaryEType);
    ia->getIDsViewOf(primaryEType, vtxIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  vGids_ = arcp_const_cast<gno_t>(
                arcp<const gno_t>(vtxIds, 0, nLocalVertices_, false));

  // Get the second adjacencies to construct edges of the dual graph.

  if (!ia->avail2ndAdjs(primaryEType, secondAdjEType)) {
    // KDDKDD TODO Want to do this differently for local and global graphs?
    // KDDKDD TODO Currently getting global 2nd Adjs and filtering them for
    // KDDKDD TODO local graphs.  That approach is consistent with other
    // KDDKDD TODO adapters, but is more expensive -- why build them just to
    // KDDKDD TODO throw them away?  Instead, perhaps should build 
    // KDDKDD TODO only local adjacencies.
    // KDDKDD TODO Does it suffice to pass a serial comm for local graph?
    try {
      get2ndAdjsViewFromAdjs(ia, comm_, primaryEType, secondAdjEType, eOffsets_,
                             eGids_);
      nLocalEdges_ = eOffsets_[nLocalVertices_];
    }
    Z2_FORWARD_EXCEPTIONS;
  }
  else {  // avail2ndAdjs
    // Get the edges
    try {
      gno_t const *nborIds=NULL;
      offset_t const *offsets=NULL;

      ia->get2ndAdjsView(primaryEType, secondAdjEType, offsets, nborIds);
      // Save the pointers from the input adapter; we do not control the
      // offsets and nborIds memory
      nLocalEdges_ = offsets[nLocalVertices_];
      eGids_ = arcp_const_cast<gno_t>(
               arcp<const gno_t>(nborIds, 0, nLocalEdges_, false));
      eOffsets_ = arcp_const_cast<offset_t>(
                  arcp<const offset_t>(offsets, 0, nLocalVertices_+1, false));
    }
    Z2_FORWARD_EXCEPTIONS;
  }


  // Edge weights
  // Cannot specify edge weights if Zoltan2 computes the second adjacencies;
  // there's no way to know the correct order for the adjacencies and weights.
  // InputAdapter must provide 2nd adjs in order for edge weights to be used.
  if (ia->avail2ndAdjs(primaryEType, secondAdjEType)) {
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

      ArrayRCP<const scalar_t> wgtArray(ewgts, 0, 
                                        nLocalEdges_*stride, false);
      eWeights_[w] = input_t(wgtArray, stride);
    }
  }

  // Do constructor common to all adapters
  shared_constructor(ia, modelFlags);

  // Get vertex coordinates
  typedef MeshAdapter<user_t> adapterWithCoords_t;
  shared_GetVertexCoords<adapterWithCoords_t>(&(*ia));

  env_->timerStop(MACRO_TIMERS, "GraphModel constructed from MeshAdapter");
  // print();
}


//////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void GraphModel<Adapter>::shared_constructor(
  const RCP<const Adapter> &ia,
  modelFlag_t &modelFlags)
{
  // Model creation flags
  bool consecutiveIdsRequired = modelFlags.test(GENERATE_CONSECUTIVE_IDS);
  bool removeSelfEdges = modelFlags.test(REMOVE_SELF_EDGES);
  bool subsetGraph = modelFlags.test(BUILD_SUBSET_GRAPH);

  // May modify the graph provided from input adapter; save pointers to 
  // the input adapter's data.
  size_t adapterNLocalEdges = nLocalEdges_;
  ArrayRCP<gno_t> adapterVGids = vGids_;     // vertices of graph from adapter
  ArrayRCP<gno_t> adapterEGids = eGids_;     // edges of graph from adapter
  ArrayRCP<offset_t> adapterEOffsets = eOffsets_;    // edge offsets from adapter
  ArrayRCP<input_t> adapterEWeights = eWeights_;  // edge weights from adapter

  if (localGraph_) {
    // Local graph is requested.
    // Renumber vertices 0 to nLocalVertices_-1
    // Filter out off-process edges
    // Renumber edge neighbors to be in range [0,nLocalVertices_-1]

    // Allocate new space for local graph info
    // Note that eGids_ and eWeights_[w] may be larger than needed; 
    // we would have to pre-count the number of local edges to avoid overalloc
    vGids_ = arcp(new gno_t[nLocalVertices_],
                  0, nLocalVertices_, true);
    eGids_ = arcp(new gno_t[adapterNLocalEdges],
                  0, adapterNLocalEdges, true);  
    eOffsets_ = arcp(new offset_t[nLocalVertices_+1],
                     0, nLocalVertices_+1, true);

    scalar_t **tmpEWeights = NULL;
    if (nWeightsPerEdge_ > 0){
      eWeights_ = arcp(new input_t[nWeightsPerEdge_], 0,
                       nWeightsPerEdge_, true);
      // Need to use temporary array because StridedData has const data
      // so we can't write to it.
      tmpEWeights = new scalar_t*[nWeightsPerEdge_];
      for (int w = 0; w < nWeightsPerEdge_; w++)
        tmpEWeights[w] = new scalar_t[adapterNLocalEdges];
    }

    // Build map between global and local vertex numbers
    std::unordered_map<gno_t, lno_t> globalToLocal(nLocalVertices_);
    for (size_t i = 0; i < nLocalVertices_; i++)
      globalToLocal[adapterVGids[i]] = i;

    // Loop over edges; keep only those that are local (i.e., on-rank)
    eOffsets_[0] = 0;
    lno_t ecnt = 0;
    for (size_t i = 0; i < nLocalVertices_; i++) {
      vGids_[i] = gno_t(i);
      for (offset_t j = adapterEOffsets[i]; j < adapterEOffsets[i+1]; j++) {

        if (removeSelfEdges && (adapterEGids[j] == adapterVGids[i]))
          continue;  // Skipping self edge

        // Determine whether neighbor vertex is local
        typename std::unordered_map<gno_t, lno_t>::iterator localidx;
        if ((localidx = globalToLocal.find(adapterEGids[j])) != 
                        globalToLocal.end()) {
          // neighbor vertex is local
          // Keep the edge and its weights
          eGids_[ecnt] = localidx->second;
          for (int w = 0; w < nWeightsPerEdge_; w++)
            tmpEWeights[w][ecnt] = adapterEWeights[w][j];

          ecnt++;
        }  
      }
      eOffsets_[i+1] = ecnt;
    }
    nLocalEdges_ = eOffsets_[nLocalVertices_];
    if (nWeightsPerEdge_) {
      for (int w = 0; w < nWeightsPerEdge_; w++) {
        ArrayRCP<const scalar_t> wgtArray(tmpEWeights[w],
                                          0, adapterNLocalEdges, true);
        eWeights_[w] = input_t(wgtArray, 0);
      }
      delete [] tmpEWeights;
    }
  }  // localGraph_

  else  if (consecutiveIdsRequired || removeSelfEdges || subsetGraph) {
    // Build a Global graph
    // If we are here, we expect SOMETHING in the graph to change from input:
    // removing self edges, or converting to consecutive IDs, or subsetting
    // the graph.


    // Determine vertex GIDs for the global GraphModel
    if (consecutiveIdsRequired) {
      // Allocate new memory for vertices for consecutiveIds
      vGids_ = arcp(new gno_t[nLocalVertices_], 0, nLocalVertices_, true);

      // Build vtxDist_ array with starting vGid on each rank
      int np = comm_->getSize();
      vtxDist_ = arcp(new size_t[np+1], 0, np+1, true);
      vtxDist_[0] = 0;
      Teuchos::gatherAll(*comm_, 1, &nLocalVertices_, np, &vtxDist_[1]);
      for (int i = 0; i < np; i++)
        vtxDist_[i+1] += vtxDist_[i];
    }

    // Allocate new memory for edges and offsets, as needed
    // Note that eGids_ may or may not be larger than needed; 
    // would have to pre-count number of edges kept otherwise
    eGids_ = arcp(new gno_t[adapterNLocalEdges],
                  0, adapterNLocalEdges, true);

    scalar_t **tmpEWeights = NULL;
    if (subsetGraph || removeSelfEdges) {
      // May change number of edges and, thus, the offsets 
      eOffsets_ = arcp(new offset_t[nLocalVertices_+1],
                       0, nLocalVertices_+1, true);
      eOffsets_[0] = 0;

      // Need to copy weights if remove edges
      if (nWeightsPerEdge_ > 0){
        eWeights_ = arcp(new input_t[nWeightsPerEdge_], 0,
                         nWeightsPerEdge_, true);
        // Need to use temporary array because StridedData has const data
        // so we can't write to it.
        tmpEWeights = new scalar_t*[nWeightsPerEdge_];
        for (int w = 0; w < nWeightsPerEdge_; w++)
          tmpEWeights[w] = new scalar_t[adapterNLocalEdges];
      }
    }

    // If needed, determine the owning ranks and its local index off-proc
    Teuchos::ArrayRCP<int> edgeRemoteRanks;
    Teuchos::ArrayRCP<lno_t> edgeRemoteLids;
    std::unordered_map<gno_t, size_t> edgeRemoteUniqueMap;

    if (subsetGraph || consecutiveIdsRequired) {

      // Find global minGID for map construction
      gno_t myMinGID = std::numeric_limits<gno_t>::max();
      size_t nVtx = adapterVGids.size();
      for (size_t i = 0; i < nVtx; i++)
        if (adapterVGids[i] < myMinGID) myMinGID = adapterVGids[i];

      gno_t minGID;
      reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_MIN, 1,
                             &myMinGID, &minGID);

      gno_t dummy = Teuchos::OrdinalTraits<gno_t>::invalid();
      Tpetra::Map<lno_t,gno_t> vtxMap(dummy, adapterVGids(), minGID, comm_);

      // Need to filter requested edges to make a unique list,
      // as Tpetra::Map does not return correct info for duplicated entries
      // (See bug 6412)  
      // The local filter may be more efficient anyway -- fewer communicated
      // values in the Tpetra directory
      Teuchos::ArrayRCP<gno_t> edgeRemoteUniqueGids = 
               arcp(new gno_t[adapterNLocalEdges], 0, adapterNLocalEdges, true);

      size_t nEdgeUnique = 0;
      for (size_t i = 0; i < adapterNLocalEdges; i++) {
        if (edgeRemoteUniqueMap.find(adapterEGids[i]) == 
            edgeRemoteUniqueMap.end()) {
          edgeRemoteUniqueGids[nEdgeUnique] = adapterEGids[i];
          edgeRemoteUniqueMap[adapterEGids[i]] = nEdgeUnique;
          nEdgeUnique++;
        }
      }

      edgeRemoteRanks = arcp(new int[nEdgeUnique], 0, nEdgeUnique, true);
      edgeRemoteLids = arcp(new lno_t[nEdgeUnique], 0, nEdgeUnique, true);
      vtxMap.getRemoteIndexList(edgeRemoteUniqueGids(0, nEdgeUnique),
                                edgeRemoteRanks(), edgeRemoteLids());
    }

    // Renumber and/or filter the edges and vertices
    lno_t ecnt = 0;
    int me = comm_->getRank();
    for (size_t i = 0; i < nLocalVertices_; i++) {

      if (consecutiveIdsRequired)
        vGids_[i] = vtxDist_[me] + i;

      for (offset_t j = adapterEOffsets[i]; j < adapterEOffsets[i+1]; j++) {

        if (removeSelfEdges && (adapterVGids[i] == adapterEGids[j]))
          continue;  // Skipping self edge

        size_t remoteIdx = edgeRemoteUniqueMap[adapterEGids[j]];

        if (subsetGraph && (edgeRemoteRanks[remoteIdx] == -1)) 
          continue;  // Skipping edge with neighbor vertex that was not found 
                     // in communicator

        if (consecutiveIdsRequired)
          // Renumber edge using local number on remote rank plus first 
          // vtx number for remote rank
          eGids_[ecnt] = edgeRemoteLids[remoteIdx]
                       + vtxDist_[edgeRemoteRanks[remoteIdx]];  
        else
          eGids_[ecnt] = adapterEGids[j];

        if (subsetGraph || removeSelfEdges) {
          // Need to copy weights only if number of edges might change
          for (int w = 0; w < nWeightsPerEdge_; w++)
            tmpEWeights[w][ecnt] = adapterEWeights[w][j];
        }

        ecnt++;
      }
      if (subsetGraph || removeSelfEdges)
        eOffsets_[i+1] = ecnt;
    }
    nLocalEdges_ = ecnt;
    if (nWeightsPerEdge_ && (subsetGraph || removeSelfEdges)) {
      for (int w = 0; w < nWeightsPerEdge_; w++) {
        ArrayRCP<const scalar_t> wgtArray(tmpEWeights[w],
                                          0, nLocalEdges_, true);
        eWeights_[w] = input_t(wgtArray, 1);
      }
      delete [] tmpEWeights;
    }
  }

  // Vertex weights
  nWeightsPerVertex_ = ia->getNumWeightsPerID();

  if (nWeightsPerVertex_ > 0){
    input_t *weightInfo = new input_t [nWeightsPerVertex_];
    env_->localMemoryAssertion(__FILE__, __LINE__, nWeightsPerVertex_,
                               weightInfo);

    for (int idx=0; idx < nWeightsPerVertex_; idx++){
      bool useNumNZ = ia->useDegreeAsWeight(idx);
      if (useNumNZ){
        scalar_t *wgts = new scalar_t [nLocalVertices_];
        env_->localMemoryAssertion(__FILE__, __LINE__, nLocalVertices_, wgts);
        ArrayRCP<const scalar_t> wgtArray = arcp(wgts,
                                                 0, nLocalVertices_, true);
        for (size_t i=0; i < nLocalVertices_; i++)
          wgts[i] = eOffsets_[i+1] - eOffsets_[i];
        weightInfo[idx] = input_t(wgtArray, 1);
      }
      else{
        const scalar_t *weights=NULL;
        int stride=0;
        ia->getWeightsView(weights, stride, idx);
        ArrayRCP<const scalar_t> wgtArray = arcp(weights, 0,
                                                 stride*nLocalVertices_,
                                                 false);
        weightInfo[idx] = input_t(wgtArray, stride);
      }
    }

    vWeights_ = arcp<input_t>(weightInfo, 0, nWeightsPerVertex_, true);
  }

  // Compute global values
  if (localGraph_) {
    nGlobalVertices_ = nLocalVertices_;
    nGlobalEdges_ = nLocalEdges_;
  }
  else {
    reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
                           &nLocalVertices_, &nGlobalVertices_);
    reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
                           &nLocalEdges_, &nGlobalEdges_);
  }

  env_->memory("After construction of graph model");
}

//////////////////////////////////////////////////////////////////////////

template <typename Adapter>
template <typename AdapterWithCoords>
void GraphModel<Adapter>::shared_GetVertexCoords(const AdapterWithCoords *ia)
{
  // get pointers to vertex coordinates from input adapter

  vCoordDim_ = ia->getDimension();

  if (vCoordDim_ > 0){
    input_t *coordInfo = new input_t [vCoordDim_];
    env_->localMemoryAssertion(__FILE__, __LINE__, vCoordDim_, coordInfo);

    for (int dim=0; dim < vCoordDim_; dim++){
      const scalar_t *coords=NULL;
      int stride=0;
      ia->getCoordinatesView(coords, stride, dim);
      ArrayRCP<const scalar_t> coordArray = arcp(coords, 0,
                                                 stride*nLocalVertices_,
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
//  if (env_->getDebugLevel() < VERBOSE_DETAILED_STATUS)
//    return;

  std::ostream *os = env_->getDebugOStream();
  
  int me = comm_->getRank();

  *os << me
      << " " << (localGraph_ ? "LOCAL GRAPH  " : "GLOBAL GRAPH  ")
      << " Nvtx  " << nLocalVertices_
      << " Nedge " << nLocalEdges_
      << " NVWgt " << nWeightsPerVertex_
      << " NEWgt " << nWeightsPerEdge_
      << " CDim  " << vCoordDim_
      << std::endl;

  for (size_t i = 0; i < nLocalVertices_; i++) {
    *os << me << " " << i << " GID " << vGids_[i] << ": ";
    for (offset_t j = eOffsets_[i]; j < eOffsets_[i+1]; j++)
      *os << eGids_[j] << " " ;
    *os << std::endl;
  }

  if (nWeightsPerVertex_) {
    for (size_t i = 0; i < nLocalVertices_; i++) {
      *os << me << " " << i << " VWGTS " << vGids_[i] << ": ";
      for (int j = 0; j < nWeightsPerVertex_; j++)
        *os << vWeights_[j][i] << " ";
      *os << std::endl;
    }
  }

  if (nWeightsPerEdge_) {
    for (size_t i = 0; i < nLocalVertices_; i++) {
      *os << me << " " << i << " EWGTS " << vGids_[i] << ": ";
      for (offset_t j = eOffsets_[i]; j < eOffsets_[i+1]; j++) {
        *os << eGids_[j] << " (";
        for (int w = 0; w < nWeightsPerEdge_; w++)
          *os << eWeights_[w][j] << " ";
        *os << ") ";
      }
      *os << std::endl;
    }
  }

  if (vCoordDim_) {
    for (size_t i = 0; i < nLocalVertices_; i++) {
      *os << me << " " << i << " COORDS " << vGids_[i] << ": ";
      for (int j = 0; j < vCoordDim_; j++)
         *os << vCoords_[j][i] << " ";
      *os << std::endl;
    }
  }
  else
    *os << me << " NO COORDINATES AVAIL " << std::endl;
}

}   // namespace Zoltan2


#endif

