// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! CommGraphModel creates a graph representing the communication topology of
the MPI ranks for a given XpetraCrsGraphAdapter object. If there are n MPI ranks
in the given communicator, then this model contains n vertices so that each vertex
represents an MPI rank. If rank i sends a message to rank j (during the mat-vec on
the matrix corresponding to the given adapter), then there is a directed edge
from vertex i to vertex j in the graph. The size of the edge is the number of
nonzeros that cause that message. The weight of vertex i is the number of nonzeros
currently residing at rank i.

Since the above mentioned graph is too small, we migrate it into a subset of ranks,
which we call activeRanks. nActiveRanks_ denotes the number of active ranks and
is computed as n/threshold_.
*/

#ifndef _ZOLTAN2_COMMGRAPHMODEL_HPP_
#define _ZOLTAN2_COMMGRAPHMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_ModelHelpers.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <unordered_map>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////
/*!  \brief CommGraphModel defines the interface required for communication graph.

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.

    For now, this model only works with GraphAdapter (XpetraCrsGraphAdapter).

*/
template <typename Adapter>
class CommGraphModel : public Model<Adapter>
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
  ~CommGraphModel() { }

  /*! \brief Constructor
   *
   *  \param  inputAdapter  a pointer to the user's data
   *  \param  env           object containing the parameters
   *  \param  comm          communicator for the problem
   *  \param  modelflags    a bit map of Zoltan2::GraphModelFlags
   *
   *  All processes in the communicator must call the constructor.
   */

  CommGraphModel(const RCP<const MatrixAdapter<user_t,userCoord_t> > &/* ia */,
		 const RCP<const Environment> &/* env */, const RCP<const Comm<int> > &/* comm */,
     const modelFlag_t &modelflags = modelFlag_t())
  {
    throw std::runtime_error("CommGraphModel is not implemented for MatrixAdapter yet.");
  }

  CommGraphModel(const RCP<const GraphAdapter<user_t,userCoord_t> > &ia,
		 const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
     const modelFlag_t &modelflags = modelFlag_t());

  CommGraphModel(const RCP<const MeshAdapter<user_t> > &/* ia */,
		 const RCP<const Environment> &/* env */, const RCP<const Comm<int> > &/* comm */,
     const modelFlag_t &modelflags = modelFlag_t())
  {
    throw std::runtime_error("CommGraphModel is not implemented for MeshAdapter yet.");
  }

  CommGraphModel(const RCP<const VectorAdapter<userCoord_t> > &/* ia */,
		 const RCP<const Environment> &/* env */, const RCP<const Comm<int> > &/* comm */,
     const modelFlag_t &modelflags = modelFlag_t())
  {
    throw std::runtime_error("cannot build CommGraphModel from VectorAdapter");
  }

  CommGraphModel(const RCP<const IdentifierAdapter<user_t> > &/* ia */,
		 const RCP<const Environment> &/* env */, const RCP<const Comm<int> > &/* comm */,
     const modelFlag_t &modelflags = modelFlag_t())
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

  ////////////////////////////////////////////////////
  // Migration-related functions.
  ////////////////////////////////////////////////////

  int getNumActiveRanks() const { return nActiveRanks_; }

  int getDestinationRank() const { return destRank_; }

  int getStartRank() const { return startRank_; }

  int getEndRank() const { return endRank_; }

private:

  void print(); // For debugging
  void migrateGraph(); // For debugging

  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  int threshold_;        // threshold on #vertices each rank stores post-migration
  int nActiveRanks_ ;    // # ranks for the small graph to be partitioned on
  int destRank_, startRank_, endRank_;


  size_t nLocalVertices_;                // # local vertices in built graph
  size_t nGlobalVertices_;               // # global vertices in built graph
  ArrayRCP<gno_t> vGids_;                  // vertices of graph built in model;
                                           // may be same as adapter's input
                                           // or may be renumbered 0 to (N-1).

  int nWeightsPerVertex_;
  ArrayRCP<input_t> vWeights_;

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
// GraphModel from GraphAdapter
template <typename Adapter>
CommGraphModel<Adapter>::CommGraphModel(
  const RCP<const GraphAdapter<user_t,userCoord_t> > &bia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  const modelFlag_t &/* modelflags */):
        env_(env),
        comm_(comm),
        nLocalVertices_(0),
        nGlobalVertices_(0),
        vGids_(),
        nWeightsPerVertex_(0),
        vWeights_(),
        nLocalEdges_(0),
        nGlobalEdges_(0),
        eGids_(),
        eOffsets_(),
        nWeightsPerEdge_(0),
        eWeights_(),
        vtxDist_()
{
  int commSize = comm_->getSize();

  // Get XpetraCrsGraphAdapter from GraphAdapter
  RCP<XpetraCrsGraphAdapter<user_t, userCoord_t>> ia;
  try{
    RCP<GraphAdapter<user_t, userCoord_t>> tmp =
       rcp_const_cast<GraphAdapter<user_t, userCoord_t>>(bia);
    ia = rcp_dynamic_cast<XpetraCrsGraphAdapter<user_t, userCoord_t>>(tmp);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Get the graph from the input adapter
  auto inGraph = ia->getXpetraGraph();

  // Get the importer of the graph
  auto imp = inGraph->getImporter();

  // Identify nbor PIDs and number of entries sent per PID
  std::map<int,double> exportpidmap;
  auto exportpids = imp->getExportPIDs();
  size_t nexportpids = imp->getNumExportIDs();
  for (size_t i = 0; i < nexportpids; i++) {
    int k = exportpids[i];
    if (exportpidmap.find(k) != exportpidmap.end())
      exportpidmap[k] = exportpidmap[k] + 1.;
    else
      exportpidmap[k] = 1.;
  }

  // Set sizes
  // There is only one vertex in each rank
  nLocalVertices_ = 1;
  nLocalEdges_ = exportpidmap.size();

  // Allocate space
  vGids_ = arcp(new gno_t[nLocalVertices_],
		0, nLocalVertices_, true);
  eGids_ = arcp(new gno_t[nLocalEdges_],
		0, nLocalEdges_, true);
  eOffsets_ = arcp(new offset_t[nLocalVertices_+1],
		   0, nLocalVertices_+1, true);
  scalar_t *wgts2 = new scalar_t [nLocalEdges_];

  // Form the vertices
  vGids_[0] = comm->getRank();

  // Form the edges
  size_t ptr = 0;
  eOffsets_[0] = ptr;
  for (std::map<int,double>::iterator it = exportpidmap.begin();
       it != exportpidmap.end(); it++) {
    eGids_[ptr] = it->first;
    wgts2[ptr++] = it->second;
  }
  eOffsets_[nLocalVertices_] = ptr;

  // Edge weights
  nWeightsPerEdge_ = 1;
  input_t *wgts = new input_t [nWeightsPerEdge_];
  eWeights_ = arcp(wgts, 0, nWeightsPerEdge_, true);

  for (int w=0; w < nWeightsPerEdge_; w++){
    int stride=0;
    ArrayRCP<const scalar_t> wgtArray = arcp(wgts2, 0, nLocalEdges_, true);
    eWeights_[w] = input_t(wgtArray, stride);
  }

  // Vertex weights
  nWeightsPerVertex_ = 1;
  input_t *weightInfo = new input_t [nWeightsPerVertex_];

  for (int idx=0; idx < nWeightsPerVertex_; idx++){
    scalar_t *wgt = new scalar_t [nLocalVertices_];
    wgt[0] = inGraph->getLocalNumEntries();
    ArrayRCP<const scalar_t> wgtArray = arcp(wgt, 0, nLocalVertices_, true);
    weightInfo[idx] = input_t(wgtArray, 1);
  }

  vWeights_ = arcp<input_t>(weightInfo, 0, nWeightsPerVertex_, true);

  reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
			 &nLocalVertices_, &nGlobalVertices_);
  reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
			 &nLocalEdges_, &nGlobalEdges_);

  // Build vtxDist_ array starting with vGid on each rank
  vtxDist_ = arcp(new size_t[commSize+1], 0, commSize+1, true);
  vtxDist_[0] = 0;
  Teuchos::gatherAll(*comm_, 1, &nLocalVertices_, commSize, &vtxDist_[1]);
  for (int i = 0; i < commSize; i++)
    vtxDist_[i+1] += vtxDist_[i];

  // Migrate the quotient graph into smaller number of MPI ranks (active ranks)
  migrateGraph();
}

template <typename Adapter>
void CommGraphModel<Adapter>::migrateGraph()
{

  // Set default threshold for migration
  threshold_ = 1024;

  // Check if the user set the threshold value
  const ParameterList &pl = env_->getParameters();
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("quotient_threshold");
  if (pe)
    threshold_ = pe->getValue<int>(&threshold_);

  // Compute the sizes of/in the new distribution
  nActiveRanks_ = std::ceil((double) nGlobalVertices_ / threshold_);
  size_t avgVertexShare = nGlobalVertices_ / nActiveRanks_;
  size_t myVertexShare = 0;

  int me = comm_->getRank();
  int commSize = comm_->getSize();

  // Save the original pointers
  ArrayRCP<offset_t> old_eOffsets_ = eOffsets_;
  ArrayRCP<gno_t> old_eGids_ = eGids_;
  size_t old_nLocalEdges_ = nLocalEdges_;
  ArrayRCP<input_t> old_vWeights_ = vWeights_;
  ArrayRCP<input_t> old_eWeights_ = eWeights_;

  // Compute whom to send to
  destRank_ = me / (int) avgVertexShare;
  if(destRank_ >= nActiveRanks_)
    destRank_ = nActiveRanks_ - 1;

  // Start with sending the size of the edge list
  RCP<CommRequest<int>> *requests;
  if(me < nActiveRanks_) {

    // Determine the range of ranks to receive edges from
    // Needs to be updated when chunks are introduced
    startRank_ = me * static_cast<int>(avgVertexShare);
    endRank_ = (me+1) * static_cast<int>(avgVertexShare);
    if(me == nActiveRanks_ - 1 ) // Last rank gets the surplus
      endRank_ = static_cast<int>(nGlobalVertices_);
    myVertexShare = endRank_ - startRank_;

    eOffsets_ = arcp(new offset_t[myVertexShare+1], 0, myVertexShare+1, true);
    eOffsets_[0] = 0;

    // Receive the sizes of their edge list
    requests = new RCP<CommRequest<int>>[myVertexShare];
    for(int i = startRank_; i < endRank_; i++) {
      requests[i-startRank_] = Teuchos::ireceive<int, offset_t>(*comm_,
							   arcp(&eOffsets_[i-startRank_+1], 0, 1, false),
							   i);
    }

    // Send adjacency size  even though this rank will remain active
    Teuchos::send<int, offset_t>(*comm_, 1, &old_eOffsets_[nLocalVertices_], destRank_);

    // Wait
    Teuchos::waitAll<int>(*comm_, Teuchos::arrayView(requests, myVertexShare));

    // Prefix sum over the offsets
    for(size_t i = 1; i <= myVertexShare; i++)
      eOffsets_[i] += eOffsets_[i-1];

    // Recompute the number of local edges
    nLocalEdges_ = eOffsets_[myVertexShare];

    // Reallocate the adjacency array
    eGids_ = arcp(new gno_t[nLocalEdges_], 0, nLocalEdges_, true);


    // Receive the adjacency lists
    for(int i = startRank_; i < endRank_; i++) {
      offset_t adjStartRank_ = eOffsets_[i-startRank_];
      offset_t adjSize = eOffsets_[i-startRank_+1] - adjStartRank_;
      requests[i-startRank_] = Teuchos::ireceive<int, gno_t>(*comm_,
							arcp(&eGids_[adjStartRank_], 0, adjSize, false),
							i);
    }

    // Send adjacency even though this rank will remain active
    Teuchos::send<int, gno_t>(*comm_, old_nLocalEdges_, &old_eGids_[0], destRank_);
    Teuchos::waitAll<int>(*comm_, Teuchos::arrayView(requests, myVertexShare));


    // Migrate vertex weights arrays
    scalar_t *wgts = new scalar_t [myVertexShare];
    for(int i = startRank_; i < endRank_; i++) {
      requests[i-startRank_] = Teuchos::ireceive<int, scalar_t>(*comm_,
							   arcp(&wgts[i-startRank_], 0, 1, false), // assumes one vertex per rank
							   i);
    }

    const scalar_t *wPtr;
    size_t wLen = 0;
    int stride = 0;
    old_vWeights_[0].getStridedList(wLen, wPtr, stride);
    Teuchos::send<int, scalar_t>(*comm_, nLocalVertices_, wPtr, destRank_);

    Teuchos::waitAll<int>(*comm_, Teuchos::arrayView(requests, myVertexShare));

    input_t *weightInfo = new input_t [nWeightsPerVertex_];
    for (int idx=0; idx < nWeightsPerVertex_; idx++){
      ArrayRCP<const scalar_t> wgtArray = arcp(wgts, 0, myVertexShare, true);
      weightInfo[idx] = input_t(wgtArray, 1);
    }
    vWeights_ = arcp<input_t>(weightInfo, 0, nWeightsPerVertex_, true);

    // Migrate edge weights arrays
    scalar_t *ewgts = new scalar_t [nLocalEdges_];
    for(int i = startRank_; i < endRank_; i++) {
      offset_t adjStartRank_ = eOffsets_[i-startRank_];
      offset_t adjSize = eOffsets_[i-startRank_+1] - adjStartRank_;
      requests[i-startRank_] = Teuchos::ireceive<int, scalar_t>(*comm_,
							   arcp(&ewgts[adjStartRank_], 0, adjSize, false), // assumes one vertex per rank
							   i);
    }

    old_eWeights_[0].getStridedList(wLen, wPtr, stride);
    Teuchos::send<int, scalar_t>(*comm_, old_nLocalEdges_, wPtr, destRank_);

    Teuchos::waitAll<int>(*comm_, Teuchos::arrayView(requests, myVertexShare));

    input_t *eweightInfo = new input_t [nWeightsPerEdge_];
    for (int idx=0; idx < nWeightsPerEdge_; idx++){
      ArrayRCP<const scalar_t> ewgtArray = arcp(ewgts, 0, nLocalEdges_, true);
      eweightInfo[idx] = input_t(ewgtArray, 1);
    }
    eWeights_ = arcp<input_t>(eweightInfo, 0, nWeightsPerEdge_, true);


    // Finalize the migration
    vGids_ = arcp(new gno_t[myVertexShare], 0, myVertexShare, true);
    for(int i = startRank_; i < endRank_; i++)
      vGids_[i-startRank_] = i;

    nLocalVertices_ = myVertexShare;


  }
  else {

    // Send adjacency size
    Teuchos::send<int, offset_t>(*comm_, 1, &eOffsets_[nLocalVertices_], destRank_);

    // Send adjacency list
    Teuchos::send<int, gno_t>(*comm_, nLocalEdges_, &eGids_[0], destRank_);

    // Send vertex weights list
    const scalar_t *wPtr;
    size_t wLen = 0;
    int stride = 0;
    vWeights_[0].getStridedList(wLen, wPtr, stride);
    Teuchos::send<int, scalar_t>(*comm_, nLocalVertices_, wPtr, destRank_);

    // Send edge weights list
    eWeights_[0].getStridedList(wLen, wPtr, stride);
    Teuchos::send<int, scalar_t>(*comm_, nLocalEdges_, wPtr, destRank_);

    nLocalVertices_ = 0;
  }

  for (int i = 0; i <= commSize; i++)
    vtxDist_[i] = 0;

  Teuchos::gatherAll(*comm_, 1, &nLocalVertices_, commSize, &vtxDist_[1]);
  for (int i = 0; i < commSize; i++)
    vtxDist_[i+1] += vtxDist_[i];

}

//////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void CommGraphModel<Adapter>::print()
{
  std::ostream *os = env_->getDebugOStream();

  int me = comm_->getRank();

  *os << me
      << " Nvtx  " << nLocalVertices_
      << " Nedge " << nLocalEdges_
      << " NVWgt " << nWeightsPerVertex_
      << " NEWgt " << nWeightsPerEdge_
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

}

}   // namespace Zoltan2


#endif

