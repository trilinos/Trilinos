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

/*! \file Zoltan2_HyperGraphModel.hpp
    \brief Defines the HyperGraphModel interface.
*/

#ifndef _ZOLTAN2_HYPERGRAPHMODEL_HPP_
#define _ZOLTAN2_HYPERGRAPHMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_ModelHelpers.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_MeshAdapter.hpp>

#include <vector>
#include <unordered_map>
#include <queue>
#include <Teuchos_Hashtable.hpp>

namespace Zoltan2 {

  /*! \brief Enumerate the views for the pins:
   *    HYPEREDGE_CENTRIC: pins are the global ids of the vertices as seen by the hyperedges
   *    VERTEX_CENTRIC: pins are the global ids of the hyperedges as seen by the vertices
   */
enum CentricView {
  HYPEREDGE_CENTRIC,
  VERTEX_CENTRIC
};

//////////////////////////////////////////////////////////////////////////
/*!  \brief HyperGraphModel defines the interface required for hyper graph models.

    The constructor of the HyperGraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.
*/
template <typename Adapter>
class HyperGraphModel : public Model<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::offset_t    offset_t;
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::node_t      node_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef Tpetra::Map<lno_t, gno_t>     map_t;
  typedef StridedData<lno_t, scalar_t>  input_t;
#endif

  //!  Destructor
  ~HyperGraphModel() { }

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
  
  HyperGraphModel(const RCP<const MatrixAdapter<user_t,userCoord_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags, CentricView view)
  {
    throw std::runtime_error("Building HyperGraphModel from MatrixAdapter not implemented yet");
  }

  HyperGraphModel(const RCP<const GraphAdapter<user_t,userCoord_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags, CentricView view)
  {
    throw std::runtime_error("Building HyperGraphModel from GraphAdapter not implemented yet");
  }
  
  HyperGraphModel(const RCP<const MeshAdapter<user_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
                  modelFlag_t &modelflags, CentricView view);
  
  HyperGraphModel(const RCP<const VectorAdapter<userCoord_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags, CentricView view)
  {
    throw std::runtime_error("cannot build HyperGraphModel from VectorAdapter");
  }

  HyperGraphModel(const RCP<const IdentifierAdapter<user_t> > &ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &flags, CentricView view)
  {
    throw std::runtime_error("cannot build HyperGraphModel from IdentifierAdapter");
  }
  

  /*! \brief Returns the centric view of the hypergraph
   */
  CentricView getCentricView() const {return view_;}

  /*! \brief Returns true if the vertices are unique false otherwise
   */
  bool areVertexIDsUnique() const {return unique;}
  
  /*! \brief Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return numLocalVertices_; }

  /*! \brief Returns the number vertices on this process that are owned.
   */
  size_t getLocalNumOwnedVertices() const { return numOwnedVertices_; }

  /*! \brief Returns the global number vertices.
   */
  size_t getGlobalNumVertices() const { return numGlobalVertices_; }

  /*! \brief Returns the number of hyper edges on this process.
   *  These are all hyper edges that have an adjacency to at 
   *  least one on process vertex.
   */
  size_t getLocalNumHyperEdges() const { return numLocalEdges_; }

  /*! \brief Returns the global number hyper edges.
   */
  size_t getGlobalNumHyperEdges() const { return numGlobalEdges_; }

  /*! \brief Returns the local number of pins
   */
  size_t getLocalNumPins() const {return numLocalPins_; }

  /*! \brief Returns the number (0 or greater) of weights per vertex
   */
  int getNumWeightsPerVertex() const { return numWeightsPerVertex_; }

  /*! \brief Returns the number (0 or greater) of weights per edge.
   */
  int getNumWeightsPerHyperEdge() const { return nWeightsPerEdge_; }

  /*! \brief Returns the number (0 or greater) of weights per pins.
   */
  int getNumWeightesPerPin() const {return nWeightsPerPin_;}

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
    size_t nv = gids_.size();
    Ids = gids_(0, nv);
    wgts = vWeights_.view(0, numWeightsPerVertex_);
    return nv;
  }

  /*! \brief Sets pointers to this process' vertex coordinates, if available

      \param xyz If vertex coordinate data is available, \c xyz
         will on return point to a StridedData object of coordinates.
   */
  size_t getVertexCoords(ArrayView<input_t> &xyz) const
  {
    size_t nv = gids_.size();
    xyz = vCoords_.view(0, vCoordDim_);
    return nv;
  }

  /*! \brief Sets pointer to the ownership of this processes vertices.

      \param isOwner will on return point to the list of ownership for
        each vertex on this process, true if this process owns the vertex
        false otherwise.
   */
  size_t getOwnedList(ArrayView<bool> &isOwner) const
  {
    size_t nv = isOwner_.size();
    isOwner = isOwner_(0, nv);
    return nv;
  }

  /*! \brief Sets pointers to the vertex map with copies and the vertex map without copies
   *         Note: the pointers will not exist if the hypergraph has unique vertices
   *               check the areVertexIDsUnique() function before calling this function
   *
   *  \param copiesMap on return points to the map of vertices with copies
   *  \param onetooneMap on return points to the map of vertices without copies
   */
  void getVertexMaps(Teuchos::RCP<const map_t>& copiesMap, 
                     Teuchos::RCP<const map_t>& onetooneMap) const 
  {
    copiesMap = mapWithCopies;
    onetooneMap = oneToOneMap;
  }

  /*! \brief Sets pointers to this process' hyperedge Ids and their weights.

      \param Ids will on return point to the list of the global Ids for
        each hyperedge on this process.
      \param wgts If hyperedge weights is available, \c wgts
         will on return point to a StridedData object of weights.
   */
  size_t getEdgeList(
    ArrayView<const gno_t> &Ids,
    ArrayView<input_t> &wgts) const
  {
    size_t nv = edgeGids_.size();
    Ids = edgeGids_(0, nv);
    wgts = eWeights_.view(0, nWeightsPerEdge_);
    return nv;
  }

  /*! \brief Sets pointers to this process' pins global Ids based on 
    the centric view given by getCentricView()

      \param pinIds This is the list of global neighbor Ids corresponding
        to the vertices or hyperedges listed in getVertexList/getEdgeList.
      \param offsets offsets[i] is the offset into pinIds to the start
        of neighbors for ith neighbor.
      \param wgts If pin weights is available, \c wgts
         will on return point to a StridedData object of weights.

      \return The number of ids in the pinIds list.
   */
  size_t getPinList( ArrayView<const gno_t> &pinIds,
    ArrayView<const offset_t> &offsets,
    ArrayView<input_t> &wgts) const
  {
    pinIds = pinGids_(0, numLocalPins_);
    offsets = offsets_.view(0, offsets_.size());
    wgts = pWeights_.view(0, nWeightsPerPin_);
    return pinGids_.size();
  }


  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const { return numLocalVertices_; }

  size_t getGlobalNumObjects() const { return numGlobalVertices_; }

private:

  struct GhostCell {
    lno_t lid; //Assumes lno_t is signed (-1 corresponds to not on this process)
    gno_t gid;
    unsigned int dist;
    GhostCell(lno_t l,gno_t g, unsigned int d) {lid=l;gid=g;dist=d;}
    bool operator<(const struct GhostCell& other) const {return dist>other.dist;}
  };
  template <typename AdapterWithCoords>
  void shared_GetVertexCoords(const AdapterWithCoords *ia);
  

  const RCP<const Environment > env_;
  const RCP<const Comm<int> > comm_;

  CentricView view_;
  bool unique;
  ArrayRCP<const gno_t> gids_;        // vertices of input graph
  
  ArrayRCP<bool> isOwner_;

  int numWeightsPerVertex_;
  ArrayRCP<input_t> vWeights_;

  int vCoordDim_;
  ArrayRCP<input_t> vCoords_;

  ArrayRCP<const gno_t> edgeGids_;
 
  int nWeightsPerEdge_;
  ArrayRCP<input_t> eWeights_;

  ArrayRCP<const gno_t> pinGids_;
  ArrayRCP<const offset_t> offsets_;

  int nWeightsPerPin_;
  ArrayRCP<input_t> pWeights_;

  // For convenience

  size_t numLocalVertices_;
  size_t numOwnedVertices_;
  size_t numGlobalVertices_;
  size_t numLocalEdges_;
  size_t numGlobalEdges_;
  size_t numLocalPins_;
  
  // For unique mapping
  Teuchos::RCP<const map_t> mapWithCopies;
  Teuchos::RCP<const map_t> oneToOneMap;

  // For debugging
  void print();

};


////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//TODO get the weights hyperedges
//GFD Do we need weights for pins too?
template <typename Adapter>
HyperGraphModel<Adapter>::HyperGraphModel(
  const RCP<const MeshAdapter<user_t> > &ia,
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &comm,
  modelFlag_t &modelFlags,
  CentricView view):
       env_(env),
       comm_(comm),
       view_(view),
       gids_(),
       isOwner_(),
       numWeightsPerVertex_(0),
       vWeights_(),
       vCoordDim_(0),
       vCoords_(),
       edgeGids_(),
       nWeightsPerEdge_(0),
       eWeights_(),
       pinGids_(),
       offsets_(),
       nWeightsPerPin_(0),
       pWeights_(),
       numLocalVertices_(0),
       numGlobalVertices_(0),
       numLocalEdges_(0),
       numGlobalEdges_(0),
       numLocalPins_(0)
{
  env_->timerStart(MACRO_TIMERS, "HyperGraphModel constructed from MeshAdapter");
  //Model Type is either traditional or ghosting
  //  Traditional:
  //    vertices == ia->getPrimaryEntityType()
  //    hyperedges == ia->getAdjacencyEntityType()
  //    pins == first adjacency between primary and adjacency types
  //  Ghosting:
  //    vertices == ia->getPrimaryEntityType()
  //    hyperedges == ia->getPrimaryEntityType()
  //    pins == k layers of second adjacency from primary through second adjacency types
  std::string model_type("traditional");
  const Teuchos::ParameterList &pl = env->getParameters();
  const Teuchos::ParameterEntry *pe2 = pl.getEntryPtr("hypergraph_model_type");
  if (pe2){
    model_type = pe2->getValue<std::string>(&model_type);
  }

  // Get the hypergraph types from adapter
  Zoltan2::MeshEntityType primaryEType = ia->getPrimaryEntityType();
  Zoltan2::MeshEntityType adjacencyEType = ia->getAdjacencyEntityType();

  // Get the IDs of the primary entity type; these are hypergraph vertices
  gno_t const *vtxIds=NULL;
  try {
    numLocalVertices_ = ia->getLocalNumOf(primaryEType);
    ia->getIDsViewOf(primaryEType, vtxIds);
    size_t maxId = *(std::max_element(vtxIds,vtxIds+numLocalVertices_));
    reduceAll(*comm_,Teuchos::REDUCE_MAX,1,&maxId,&numGlobalVertices_);
    // TODO:  KDD 1/17 The above computation of numGlobalVertices_ is
    // TODO:  correct only when the vertices are consecutively numbered
    // TODO:  starting at ID 1.  Github #1024
  }
  Z2_FORWARD_EXCEPTIONS;

  gids_ = arcp<const gno_t>(vtxIds, 0, numLocalVertices_, false);

  //A mapping from gids to lids for efficiency
  std::unordered_map<gno_t,lno_t> lid_mapping;
  for (size_t i=0;i<numLocalVertices_;i++)
    lid_mapping[gids_[i]]=i;

  // Define owners for each hypergraph vertex using Tpetra 
  // one to one map. This defines each hypergraph vertex to
  // one process in the case that the adapter has copied 
  // primary entity types
  //If the mesh adapter knows the entities are unique we can optimize out the ownership
  unique = ia->areEntityIDsUnique(ia->getPrimaryEntityType());
  numOwnedVertices_=numLocalVertices_;
  isOwner_ = ArrayRCP<bool>(numLocalVertices_,true);
  if (!unique) {
    
    Tpetra::global_size_t numGlobalCoords = 
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    mapWithCopies = rcp(new map_t(numGlobalCoords, gids_(), 0, comm));
    // TODO KDD 1/17 It would be better to use minimum GID rather than
    // TODO zero in the above Tpetra::Map constructor.  Github #1024
    oneToOneMap = Tpetra::createOneToOne<lno_t, gno_t>(mapWithCopies);

    numOwnedVertices_=oneToOneMap->getNodeNumElements();
    for (size_t i=0;i<numLocalVertices_;i++) {
      isOwner_[i] = oneToOneMap->isNodeGlobalElement(gids_[i]);
    }
  }


  if (model_type=="traditional") {
    // Traditional: Get the IDs of the adjacency entity type; 
    //              these are hypergraph hyperedges
  
    gno_t const *edgeIds=NULL;
    try {
      numLocalEdges_ = ia->getLocalNumOf(adjacencyEType);
      ia->getIDsViewOf(adjacencyEType, edgeIds);
      size_t maxId = *(std::max_element(edgeIds,edgeIds+numLocalEdges_));
      reduceAll(*comm_,Teuchos::REDUCE_MAX,1,&maxId,&numGlobalEdges_);
    }
    Z2_FORWARD_EXCEPTIONS;
    
    edgeGids_ = arcp<const gno_t>(edgeIds, 0, numLocalEdges_, false);
  }
  else if (model_type=="ghosting") {
    // Ghosting: Use the vertices as the hyperedges as well
    numLocalEdges_ = numLocalVertices_;
    edgeGids_ = arcp<const gno_t>(vtxIds, 0, numLocalVertices_, false);
    numGlobalEdges_ = numGlobalVertices_;
  }
 
  //Define the entity types to use for the pins based on the centric view
  Zoltan2::MeshEntityType primaryPinType = primaryEType;
  Zoltan2::MeshEntityType adjacencyPinType = adjacencyEType;
  size_t numPrimaryPins = numLocalVertices_;
  if (view_==HYPEREDGE_CENTRIC) {
    primaryPinType = adjacencyEType;
    adjacencyPinType = primaryEType;
    numPrimaryPins = numLocalEdges_;
  }
  if (model_type=="traditional") {
    //Get the pins from using the traditional method of first adjacency
    gno_t const *nborIds=NULL;
    offset_t const *offsets=NULL;
    
    try {
      ia->getAdjsView(primaryPinType,adjacencyPinType,offsets,nborIds);
    }
    Z2_FORWARD_EXCEPTIONS;
    
    numLocalPins_ = offsets[numPrimaryPins];

    pinGids_ = arcp<const gno_t>(nborIds, 0, numLocalPins_, false);
    offsets_ = arcp<const offset_t>(offsets, 0, numPrimaryPins + 1, false);
  }
  else if (model_type=="ghosting") { 
    // set the view to either since it doesn't matter 
    // vertices==hyperedges
    view_ = VERTEX_CENTRIC;
    // unique set of global ids for the ghosts
    typedef std::set<gno_t> ghost_t;

    // mapping from global id to the set of ghosts 
    typedef std::unordered_map<gno_t,ghost_t> ghost_map_t;
    
    primaryPinType=primaryEType;
    adjacencyPinType =ia->getSecondAdjacencyEntityType();

    // number of layers of ghosting to do
    unsigned int layers=2;
    const Teuchos::ParameterEntry *pe3 = pl.getEntryPtr("ghost_layers");
    if (pe3){
      int l;
      l = pe3->getValue<int>(&l);
      layers = static_cast<unsigned int>(l);
    }
 
    typedef int nonzero_t;  // adjacency matrix doesn't need scalar_t
    typedef Tpetra::CrsMatrix<nonzero_t,lno_t,gno_t,node_t>   sparse_matrix_type;
    
    // Get an adjacency matrix representing the graph on the mesh 
    // using second adjacencies. If second adjacencies are not 
    // provided build the matrix from first adjacencies.
    RCP<sparse_matrix_type> secondAdj;
    if (!ia->avail2ndAdjs(primaryPinType,adjacencyPinType)) {
      secondAdj=Zoltan2::get2ndAdjsMatFromAdjs<user_t>(ia,comm_,primaryPinType, adjacencyPinType);
    }
    else {
      const offset_t* offsets;
      const gno_t* adjacencyIds;
      ia->get2ndAdjsView(primaryPinType,adjacencyPinType,offsets,adjacencyIds);
      if (unique) {
        Tpetra::global_size_t numGlobalCoords = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
        oneToOneMap = rcp(new map_t(numGlobalCoords, gids_(), 0, comm));
        // TODO KDD 1/17 It would be better to use minimum GID rather than
        // TODO zero in the above Tpetra::Map constructor.  Github #1024
      }
      Teuchos::Array<size_t> nPerRow(numLocalVertices_);
      size_t rowcnt = 0;
      for (size_t i=0; i<numLocalVertices_;i++) {
        if (!isOwner_[i])
          continue;
        nPerRow[rowcnt++] = offsets[i+1]-offsets[i];
      }
      secondAdj = rcp(new sparse_matrix_type(oneToOneMap,nPerRow(0,rowcnt)));
      for (size_t i=0; i<numLocalVertices_;i++) {
        if (!isOwner_[i])
          continue;
        gno_t row = gids_[i];
        offset_t num_adjs = offsets[i+1]-offsets[i];
        ArrayRCP<nonzero_t> ones(num_adjs,1);
        ArrayRCP<const gno_t> cols(adjacencyIds,offsets[i],num_adjs,false);
        secondAdj->insertGlobalValues(row,cols(),ones());
      }
      secondAdj->fillComplete();
    }

    //The mapping of the ghosts per hypergraph vertex
    ghost_map_t ghosts; 

    //Read the 1 layer ghosts from the second adjacency matrix
    Array<gno_t> Indices;
    Array<nonzero_t> Values;
    for (unsigned int i=0;i<numLocalEdges_;i++) {
      if (!isOwner_[i])
        continue;
      gno_t gid = edgeGids_[i];
      size_t NumEntries = secondAdj->getNumEntriesInGlobalRow (gid);
      Indices.resize (NumEntries);
      Values.resize (NumEntries);
      secondAdj->getGlobalRowCopy(gid,Indices(),Values(),NumEntries);
      for (size_t j = 0; j < NumEntries; ++j) {
	if(gid != Indices[j]) {
          ghosts[gid].insert(Indices[j]);
	}
      }
    }
    
    // The ith power of the second adjacency matrix is the ith layer of ghosts.
    // Here we compute the ith power of the matrix and add the ith layer ghosts
    // from the new matrix.
    RCP<sparse_matrix_type> mat_old = secondAdj;
    for (unsigned int i=1;i<layers;i++) {
      RCP<sparse_matrix_type> mat_new = 
        rcp (new sparse_matrix_type(secondAdj->getRowMap(),0));
      Tpetra::MatrixMatrix::Multiply(*mat_old,false,*secondAdj,false,*mat_new);
      for (unsigned int j=0;j<numLocalEdges_;j++) {
        if (!isOwner_[j])
          continue;
        gno_t gid = edgeGids_[j];
        size_t NumEntries = mat_new->getNumEntriesInGlobalRow (gid);
        Indices.resize(NumEntries);
        Values.resize(NumEntries);
        mat_new->getGlobalRowCopy(gid,Indices(),Values(),NumEntries);
        for (size_t k = 0; k < NumEntries; ++k) 
          if(gid != Indices[k]) 
            ghosts[gid].insert(Indices[k]);
        
      }
      mat_old = mat_new;
    }

    //Make the pins from the ghosts
    for (size_t i=0;i<numLocalVertices_;i++) {//for each local entity
      numLocalPins_+=ghosts[gids_[i]].size();
    }
    gno_t* temp_pins = new gno_t[numLocalPins_];
    offset_t* temp_offsets = new offset_t[numLocalVertices_+1];
    gno_t j=0;
    for (size_t i=0;i<numLocalVertices_;i++) {//for each local entity
      temp_offsets[i]=j;
      if (!isOwner_[i])
        continue;
      typename ghost_t::iterator itr;
      for (itr=ghosts[gids_[i]].begin();itr!=ghosts[gids_[i]].end();itr++) { //for each ghost of this entity
        temp_pins[j]=*itr;
        j++;
        
      }
    }
    temp_offsets[numLocalVertices_]=numLocalPins_;
    pinGids_ = arcp<const gno_t>(temp_pins,0,numLocalPins_,true);
    offsets_ = arcp<const offset_t>(temp_offsets,0,numLocalVertices_+1,true);
    
    //==============================Ghosting complete=================================
  }


  //Get the vertex weights
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

  //TODO get the weights for edges, and pins(?)

  //Get the vertex coordinates from the primary types
  typedef MeshAdapter<user_t> adapterWithCoords_t;
  shared_GetVertexCoords<adapterWithCoords_t>(&(*ia));

  env_->timerStop(MACRO_TIMERS, "HyperGraphModel constructed from MeshAdapter");
  print();
}

//////////////////////////////////////////////////////////////////////////

template <typename Adapter>
template <typename AdapterWithCoords>
void HyperGraphModel<Adapter>::shared_GetVertexCoords(const AdapterWithCoords *ia)
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
void HyperGraphModel<Adapter>::print()
{
  //only prints the model if debug status is verbose
  if (env_->getDebugLevel() < VERBOSE_DETAILED_STATUS)
    return;

  std::ostream *os = env_->getDebugOStream();
  
  int me = comm_->getRank();
  std::string fn(" ");

  *os << me << fn
      << " Nvtx  " << gids_.size()
      << " Nedge " << edgeGids_.size()
      << " NPins " << numLocalPins_
      << " NVWgt " << numWeightsPerVertex_
      << " NEWgt " << nWeightsPerEdge_
      << " NPWgt " << nWeightsPerPin_
      << " CDim  " << vCoordDim_
      << std::endl;

  for (lno_t i = 0; i < gids_.size(); i++) {
    *os << me << fn << i << " VTXGID " << gids_[i]<<" isOwner: "<<isOwner_[i];
    if (numWeightsPerVertex_==1)
      *os << " weight: " << vWeights_[0][i]; 
    if (view_==VERTEX_CENTRIC) {
      *os <<" pins:";
      for (offset_t j = offsets_[i]; j< offsets_[i+1];j++)
        *os <<" "<<pinGids_[j];
    }
    *os << std::endl;
  }
  for (lno_t i = 0; i<edgeGids_.size();i++) {
    *os << me << fn << i << " EDGEGID " << edgeGids_[i];
    if (view_==HYPEREDGE_CENTRIC) {
      *os <<":";
      for (offset_t j = offsets_[i]; j< offsets_[i+1];j++)
        *os <<" "<<pinGids_[j];
    }
    *os << std::endl;
  }
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

