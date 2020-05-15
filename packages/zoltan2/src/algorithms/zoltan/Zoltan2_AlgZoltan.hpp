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
#ifndef _ZOLTAN2_ALGZOLTAN_HPP_
#define _ZOLTAN2_ALGZOLTAN_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#include <Zoltan2_Model.hpp>

#include <Zoltan2_AlgZoltanCallbacks.hpp>
#include <zoltan_cpp.h>
#include <zoltan_partition_tree.h>

//////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgZoltan.hpp
//! \brief interface to the Zoltan package
//  
//  This first design templates Zoltan's callback functions on the 
//  input adapter.  This approach has the advantage of simplicity and
//  is most similar to current usage of Zoltan (where the callbacks define
//  the model).
//  A better approach might template them on a model, 
//  allowing Zoltan2 greater flexibility in creating models from the input.
//  Alternatively, different callback implementations could be provided to
//  represent different models to Zoltan.
//////////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

template <typename Adapter>
class AlgZoltan : public Algorithm<Adapter>
{

private:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const typename Adapter::base_adapter_t> adapter;
  RCP<const Model<Adapter> > model;
  RCP<Zoltan> zz;

  MPI_Comm mpicomm;
  
  void setMPIComm(const RCP<const Comm<int> > &problemComm__) {
#   ifdef HAVE_ZOLTAN2_MPI
      mpicomm = Teuchos::getRawMpiComm(*problemComm__);
#   else
      mpicomm = MPI_COMM_WORLD;  // taken from siMPI
#   endif
  }

  void zoltanInit() {
    // call Zoltan_Initialize to make sure MPI_Init is called (in MPI or siMPI).
    int argc = 0;
    char **argv = NULL;
    float ver;
    Zoltan_Initialize(argc, argv, &ver);
  }

  void setCallbacksIDs()
  {
    zz->Set_Num_Obj_Fn(zoltanNumObj<Adapter>, (void *) &(*adapter));
    zz->Set_Obj_List_Fn(zoltanObjList<Adapter>, (void *) &(*adapter));

    const part_t *myparts;
    adapter->getPartsView(myparts);
    if (myparts != NULL)
      zz->Set_Part_Multi_Fn(zoltanParts<Adapter>, (void *) &(*adapter));

    char tmp[4];
    sprintf(tmp, "%d", TPL_Traits<ZOLTAN_ID_PTR, gno_t>::NUM_ID);
    zz->Set_Param("NUM_GID_ENTRIES", tmp);
    sprintf(tmp, "%d", TPL_Traits<ZOLTAN_ID_PTR, lno_t>::NUM_ID);
    zz->Set_Param("NUM_LID_ENTRIES", tmp);
  }

  template <typename AdapterWithCoords>
  void setCallbacksGeom(const AdapterWithCoords *ia)
  {
    // Coordinates may be provided by the MeshAdapter or VectorAdapter.
    // VectorAdapter may be provided directly by user or indirectly through
    // GraphAdapter or MatrixAdapter.  So separate template type is needed.
    zz->Set_Num_Geom_Fn(zoltanNumGeom<AdapterWithCoords>, (void *) ia);
    zz->Set_Geom_Multi_Fn(zoltanGeom<AdapterWithCoords>, (void *) ia);
  }

  void setCallbacksGraph(
    const RCP<const GraphAdapter<user_t,userCoord_t> > &/* adp */)
  {
    // std::cout << "NotReadyForGraphCallbacksYet" << std::endl;
    // TODO
  }

  void setCallbacksGraph(
    const RCP<const MatrixAdapter<user_t,userCoord_t> > &adp)
  {
    // std::cout << "NotReadyForGraphCallbacksYet" << std::endl;
    // TODO
  }

  void setCallbacksGraph(
    const RCP<const MeshAdapter<user_t> > &adp)
  {
    // std::cout << "NotReadyForGraphCallbacksYet" << std::endl;
    // TODO
  }

  void setCallbacksHypergraph(
    const RCP<const MatrixAdapter<user_t,userCoord_t> > &adp)
  {
    // TODO:  If add parameter list to this function, can register 
    // TODO:  different callbacks depending on the hypergraph model to use

    zz->Set_HG_Size_CS_Fn(zoltanHGSizeCS_withMatrixAdapter<Adapter>,
                          (void *) &(*adp));
    zz->Set_HG_CS_Fn(zoltanHGCS_withMatrixAdapter<Adapter>,
                     (void *) &(*adp));

    // zz->Set_HG_Size_Edge_Wts_Fn(zoltanHGSizeEdgeWts_withMatrixAdapter<Adapter>,
    //                             (void *) &(*adapter));
    // zz->Set_HG_Edge_Wts_Fn(zoltanHGSizeEdgeWts_withMatrixAdapter<Adapter>,
    //                             (void *) &(*adapter));
  }

  void setCallbacksHypergraph(
    const RCP<const GraphAdapter<user_t,userCoord_t> > &adp)
  {
    zz->Set_HG_Size_CS_Fn(zoltanHGSizeCS_withGraphAdapter<Adapter>,
                          (void *) &(*adp));
    zz->Set_HG_CS_Fn(zoltanHGCS_withGraphAdapter<Adapter>,
                     (void *) &(*adp));

    if (adp->getNumWeightsPerEdge() != 0) {
      if (adp->getNumWeightsPerEdge() > 1) {
        std::cout << "Zoltan2 warning:  getNumWeightsPerEdge() returned "
                  << adp->getNumWeightsPerEdge() << " but PHG supports only "
                  << " one weight per edge; only first weight will be used."
                  << std::endl;
      }
      zz->Set_HG_Size_Edge_Wts_Fn(zoltanHGSizeEdgeWts_withGraphAdapter<Adapter>,
                                  (void *) &(*adapter));
      zz->Set_HG_Edge_Wts_Fn(zoltanHGEdgeWts_withGraphAdapter<Adapter>,
                             (void *) &(*adapter));
    }
  }

  void setCallbacksHypergraph(const RCP<const MeshAdapter<user_t> > &adp)
  {
    
    const Teuchos::ParameterList &pl = env->getParameters();

    const Teuchos::ParameterEntry *pe = pl.getEntryPtr("hypergraph_model_type");
    std::string model_type("traditional");
    if (pe){
      model_type = pe->getValue<std::string>(&model_type);
    }

    if (model_type=="ghosting" || 
        !adp->areEntityIDsUnique(adp->getPrimaryEntityType())) {
      Zoltan2::modelFlag_t flags;
      HyperGraphModel<Adapter>* mdl = new HyperGraphModel<Adapter>(adp, env,
                                                          problemComm, flags,
                                                          HYPEREDGE_CENTRIC);
      model = rcp(static_cast<const Model<Adapter>* >(mdl),true);
      
      zz->Set_Num_Obj_Fn(zoltanHGNumObj_withModel<Adapter>, (void *) &(*mdl));
      zz->Set_Obj_List_Fn(zoltanHGObjList_withModel<Adapter>, (void *) &(*mdl));
      
      zz->Set_HG_Size_CS_Fn(zoltanHGSizeCS_withModel<Adapter>, (void *) &(*mdl));
      zz->Set_HG_CS_Fn(zoltanHGCS_withModel<Adapter>, (void *) &(*mdl));
    }
    else {
      //If entities are unique we dont need the extra cost of the model
      zz->Set_HG_Size_CS_Fn(zoltanHGSizeCS_withMeshAdapter<Adapter>,
                            (void *) &(*adp));
      zz->Set_HG_CS_Fn(zoltanHGCS_withMeshAdapter<Adapter>,
                       (void *) &(*adp));
    }
    // zz->Set_HG_Size_Edge_Wts_Fn(zoltanHGSizeEdgeWts_withMeshAdapter<Adapter>,
    //                               (void *) &(*adp));
    // zz->Set_HG_Edge_Wts_Fn(zoltanHGSizeEdgeWts_withMeshAdapter<Adapter>,
    //                         (void *) &(*adp));
  }
  
  //! \brief  rcb is always binary
  virtual bool isPartitioningTreeBinary() const
  {
    return true;
  }

  //! \brief  handles the building of the splitRangeBeg and splitRangeEnd arrays
  void rcb_recursive_partitionTree_calculations(
                        part_t arrayIndex,
                        part_t numParts,
                        std::vector<part_t> & splitRangeBeg,
                        std::vector<part_t> & splitRangeEnd) const
  {
    // Note the purpose of the recursive method is make sure the children of a
    // node have updated their values for splitRangeBeg and splitRangeEnd
    // Then we can set our own values simply based on the union
    // first load the rcb data for the node
    int parent = -1;
    int left_leaf = -1;
    int right_leaf = -1;
    int err = Zoltan_RCB_Partition_Tree(zz->Get_C_Handle(),
      arrayIndex - numParts + 1, // rcb starts as 1 but does not include terminals
      &parent, &left_leaf, &right_leaf);
    if(err != 0) {
      throw std::logic_error( "Zoltan_RCB_Partition_Tree returned in error." );
    }
    // check that children both have their ranges set and if not, do those
    // range first so we can build them to make our range
    if(left_leaf > 0) { // neg is terminal and always already built
      rcb_recursive_partitionTree_calculations(left_leaf+numParts-1, numParts,
        splitRangeBeg, splitRangeEnd);
    }
    if(right_leaf > 0) { // neg is terminal and always already built
      rcb_recursive_partitionTree_calculations(right_leaf+numParts-1, numParts,
        splitRangeBeg, splitRangeEnd);
    }
    // now we can build our ranges from the children
    // note this exploits the rcb conventions for right and left so we know
    // that left_leaf will be our smaller indices
    int leftIndex = (left_leaf > 0) ? (left_leaf-1+numParts) : (-left_leaf);
    int rightIndex = (right_leaf > 0) ? (right_leaf-1+numParts) : (-right_leaf);
    splitRangeBeg[arrayIndex] = splitRangeBeg[leftIndex];
    splitRangeEnd[arrayIndex] = splitRangeEnd[rightIndex];
    // for debugging sanity check verify left_leaf is a set of indices which
    // goes continuously into the right_leaf
    if(splitRangeBeg[rightIndex] != splitRangeEnd[leftIndex]) { // end is non-inclusive
      throw std::logic_error( "RCB expected left_leaf indices and right leaf"
        " indices to be continuous but it was not so." );
    }
  }

  //! \brief  fill arrays with rcb partition tree info
  void rcb_getPartitionTree(part_t numParts,
                        part_t & numTreeVerts,
                        std::vector<part_t> & permPartNums,
                        std::vector<part_t> & splitRangeBeg,
                        std::vector<part_t> & splitRangeEnd,
                        std::vector<part_t> & treeVertParents) const
  {
    // CALCULATE: numTreeVerts
    // For rcb a tree node always takes 2 nodes and turns them into 1 node
    // That means it takes numParts - 1 nodes to reduce a tree of numParts to
    // a single root node - but we do 2 * numParts - 1 because we are currently
    // treating all of the 'trivial' terminals as tree nodes - something we
    // discussed we may change later
    part_t numTreeNodes = 2 * numParts - 1;
    numTreeVerts = numTreeNodes - 1; // by design convention root not included
    // CALCULATE: permPartNums
    permPartNums.resize(numParts);
    for(part_t n = 0; n < numParts; ++n) {
      permPartNums[n] = n; // for rcb we can assume this is trivial and in order
    }
    // CALCULATE: treeVertParents
    treeVertParents.resize(numTreeNodes); // allocate space for numTreeNodes array
    // scan all the non terminal nodes and set all children to have us as parent
    // that will set all parents except for the root node which we will set to -1
    // track the children of the root and final node for swapping later. Couple
    // ways to do this - all seem a bit awkward but at least this is efficient.
    part_t rootNode = 0; // track index of the root node for swapping
    // a bit awkward but efficient - save the children of root and final node
    // for swap at end to satisfy convention that root is highest index node
    part_t saveRootNodeChildrenA = -1;
    part_t saveRootNodeChildrenB = -1;
    part_t saveFinalNodeChildrenA = -1;
    part_t saveFinalNodeChildrenB = -1;
    for(part_t n = numParts; n < numTreeNodes; ++n) { // scan and set all parents
      int parent = -1;
      int left_leaf = -1;
      int right_leaf = -1;
      int err = Zoltan_RCB_Partition_Tree(zz->Get_C_Handle(),
        static_cast<int>(n - numParts) + 1, // rcb starts as 1 but does not include terminals
        &parent, &left_leaf, &right_leaf);
      if(err != 0) {
        throw std::logic_error("Zoltan_RCB_Partition_Tree returned in error.");
      }
      part_t leftIndex = (left_leaf > 0) ? (left_leaf-1+numParts) : (-left_leaf);
      part_t rightIndex = (right_leaf > 0) ? (right_leaf-1+numParts) : (-right_leaf);
      treeVertParents[leftIndex] = n;
      treeVertParents[rightIndex] = n;
      // save root node for final swap
      if(parent == 1 || parent == -1) { // is it the root?
        rootNode = n; // remember I am the root
        saveRootNodeChildrenA = leftIndex;
        saveRootNodeChildrenB = rightIndex;
      }
      if(n == numTreeNodes-1) {
        saveFinalNodeChildrenA = leftIndex;
        saveFinalNodeChildrenB = rightIndex;
      }
    }
    treeVertParents[rootNode] = -1; // convention parent is root -1
    // splitRangeBeg and splitRangeEnd
    splitRangeBeg.resize(numTreeNodes);
    splitRangeEnd.resize(numTreeNodes);
    // for terminal nodes this is trivial
    for(part_t n = 0; n < numParts; ++n) {
      splitRangeBeg[n] = n;
      splitRangeEnd[n] = n + 1;
    }
    if(numParts > 1) { // not relevant for 1 part
      // build the splitRangeBeg and splitRangeEnd recursively which forces the
      // children of each node to be calculated before the parent so parent can
      // just take the union of the two children
      rcb_recursive_partitionTree_calculations(rootNode, numParts, splitRangeBeg,
        splitRangeEnd);
      // now as a final step handle the swap to root is the highest index node
      // swap the parent of the two nodes
      std::swap(treeVertParents[rootNode], treeVertParents[numTreeNodes-1]);
      // get the children of the swapped nodes to have updated parents
      treeVertParents[saveFinalNodeChildrenA] = rootNode;
      treeVertParents[saveFinalNodeChildrenB] = rootNode;
      // handle case where final node is child of the root
      if(saveRootNodeChildrenA == numTreeNodes - 1) {
        saveRootNodeChildrenA = rootNode;
      }
      if(saveRootNodeChildrenB == numTreeNodes - 1) {
        saveRootNodeChildrenB = rootNode;
      }
      treeVertParents[saveRootNodeChildrenA] = numTreeNodes - 1;
      treeVertParents[saveRootNodeChildrenB] = numTreeNodes - 1;
      // update the beg and end indices - simply swap them
      std::swap(splitRangeBeg[rootNode], splitRangeBeg[numTreeNodes-1]);
      std::swap(splitRangeEnd[rootNode], splitRangeEnd[numTreeNodes-1]);
    }
  }

  //! \brief  fill arrays with rcb partition tree info
  void phg_getPartitionTree(part_t numParts,
                        part_t & numTreeVerts,
                        std::vector<part_t> & permPartNums,
                        std::vector<part_t> & splitRangeBeg,
                        std::vector<part_t> & splitRangeEnd,
                        std::vector<part_t> & treeVertParents) const
  {
    // First thing is to get the length of the tree from zoltan.
    // The tree is a list of pairs (lo,hi) for each node.
    // Here tree_array_size is the number of pairs.
    // In phg indexing the first pairt (i=0) is empty garbage.
    // The second pair (index 1) will be the root.
    // Some nodes will be empty nodes, determined by hi = -1.
    int tree_array_size = -1; // will be set by Zoltan_PHG_Partition_Tree_Size
    int err = Zoltan_PHG_Partition_Tree_Size(
      zz->Get_C_Handle(), &tree_array_size);
    if(err != 0) {
      throw std::logic_error("Zoltan_PHG_Partition_Tree_Size returned error.");
    }
    // Determine the number of valid nodes (PHG will have empty slots)
    // We scan the list of pairs and count each node which does not have hi = -1
    // During the loop we will also construct mapIndex which maps initial n
    // to final n due to some conversions we apply to meet the design specs.
    // The conversions implemented by mapIndex are:
    //    Move all terminals to the beginning (terminals have hi = lo)
    //    Resort the terminals in order (simply map to index lo works)
    //    Move non-terminals after the terminals (they start at index numParts)
    //    Map the first pair (root) to the be last to meet the design spec
    part_t numTreeNodes = 0;
    std::vector<part_t> mapIndex(tree_array_size); // maps n to final index
    part_t trackNonTerminalIndex = numParts; // starts after terminals
    for(part_t n = 0; n < static_cast<part_t>(tree_array_size); ++n) {
      part_t phgIndex = n + 1; // phg indexing starts at 1
      int lo_index = -1;
      int hi_index = -1;
      err = Zoltan_PHG_Partition_Tree(
        zz->Get_C_Handle(), phgIndex, &lo_index, &hi_index);
      if(hi_index != -1) { // hi -1 means it's an unused node
        ++numTreeNodes; // increase the total count because this is a real node
        if(n != 0) { // the root is mapped last but we don't know the length yet
          mapIndex[n] = (lo_index == hi_index) ? // is it a terminal?
            lo_index : // terminals map in sequence - lo_index is correct
            (trackNonTerminalIndex++); // set then bump trackNonTerminalIndex +1
        }
      }
    }
    // now complete the map by setting root to the length-1 for the design specs
    mapIndex[0] = numTreeNodes - 1;
    // CALCULATE: numTreeVerts
    numTreeVerts = numTreeNodes - 1; // this is the design - root not included
    // CALCULATE: permPartNums
    permPartNums.resize(numParts);
    for(part_t n = 0; n < numParts; ++n) {
      permPartNums[n] = n; // for phg we can assume this is trivial and in order
    }
    // CALCULATE: treeVertParents, splitRangeBeg, splitRangeEnd
    // we will determine all of these in this second loop using mapIndex
    // First set the arrays to have the proper length
    treeVertParents.resize(numTreeNodes);
    splitRangeBeg.resize(numTreeNodes);
    splitRangeEnd.resize(numTreeNodes);
    // Now loop a second time
    for(part_t n = 0; n < tree_array_size; ++n) {
      part_t phgIndex = n + 1; // phg indexing starts at 1
      int lo_index = -1; // zoltan Zoltan_PHG_Partition_Tree will set this
      int hi_index = -1; // zoltan Zoltan_PHG_Partition_Tree will set this
      err = Zoltan_PHG_Partition_Tree( // access zoltan phg tree data
        zz->Get_C_Handle(), phgIndex, &lo_index, &hi_index);
      if(err != 0) {
        throw std::logic_error("Zoltan_PHG_Partition_Tree returned in error.");
      }
      if(hi_index != -1) { // hi -1 means it's an unused node (a gap)
        // get final index using mapIndex - so convert from phg to design plan
        part_t finalIndex = mapIndex[n]; // get the index for the final output
        // now determine the parent
        // In the original phg indexing, the parent can be directly calculated
        // from the pair index using the following rules:
        // if phgIndex even, then parent is phgIndex/2
        //   here we determine even by ((phgIndex%2) == 0)
        // if phgIndex odd, then parent is (phgIndex-1)/2
        // but after getting parentPHGIndex we convert back to this array
        // indexing by subtracting 1
        part_t parentPHGIndex =
          ((phgIndex%2) == 0) ? (phgIndex/2) : ((phgIndex-1)/2);
        // map the parent phg index to the final parent index
        // however, for the special case of the root (n=1), set the parent to -1
        treeVertParents[finalIndex] = (n==0) ? -1 : mapIndex[parentPHGIndex-1];
        // set begin (inclusive) and end (non-inclusive), so end is hi+1
        splitRangeBeg[finalIndex] = static_cast<part_t>(lo_index);
        splitRangeEnd[finalIndex] = static_cast<part_t>(hi_index+1);
      }
    }
  }

  //! \brief  fill arrays with partition tree info
  void getPartitionTree(part_t numParts,
                        part_t & numTreeVerts,
                        std::vector<part_t> & permPartNums,
                        std::vector<part_t> & splitRangeBeg,
                        std::vector<part_t> & splitRangeEnd,
                        std::vector<part_t> & treeVertParents) const
  {
    // first check that our parameters requested we keep the tree
    const Teuchos::ParameterList &pl = env->getParameters();
    bool keep_partition_tree = false;
    const Teuchos::ParameterEntry * pe = pl.getEntryPtr("keep_partition_tree");
    if(pe) {
      keep_partition_tree = pe->getValue(&keep_partition_tree);
      if(!keep_partition_tree) {
        throw std::logic_error(
          "Requested tree when param keep_partition_tree is off.");
      }
    }

    // now call the appropriate method based on LB_METHOD in the zoltan
    // parameters list.
    const Teuchos::ParameterList & zoltan_pl = pl.sublist("zoltan_parameters");
    std::string lb_method;
    pe = zoltan_pl.getEntryPtr("LB_METHOD");
    if(pe) {
      lb_method = pe->getValue(&lb_method);
    }
    if(lb_method == "phg") {
      phg_getPartitionTree(numParts, numTreeVerts, permPartNums,
          splitRangeBeg, splitRangeEnd, treeVertParents);
    }
    else if(lb_method == "rcb") {
      rcb_getPartitionTree(numParts, numTreeVerts, permPartNums,
          splitRangeBeg, splitRangeEnd, treeVertParents);
    }
    else {
      throw std::logic_error("Did not recognize LB_METHOD: " + lb_method);
    }
  }

public:

  /*! Zoltan constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter  the user's input adapter
   */
  AlgZoltan(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const IdentifierAdapter<user_t> > &adapter__):
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
    zoltanInit();
    zz = rcp(new Zoltan(mpicomm)); 
    setCallbacksIDs();
  }

  AlgZoltan(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const VectorAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
    zoltanInit();
    zz = rcp(new Zoltan(mpicomm)); 
    setCallbacksIDs();
    setCallbacksGeom(&(*adapter));
  }

  AlgZoltan(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
    zoltanInit();
    zz = rcp(new Zoltan(mpicomm)); 
    setCallbacksIDs();
    setCallbacksGraph(adapter);
    setCallbacksHypergraph(adapter);
    if (adapter->coordinatesAvailable()) {
      setCallbacksGeom(adapter->getCoordinateInput());
    }
  }

  AlgZoltan(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const MatrixAdapter<user_t,userCoord_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
    zoltanInit();
    zz = rcp(new Zoltan(mpicomm)); 
    setCallbacksIDs();
    setCallbacksGraph(adapter);
    setCallbacksHypergraph(adapter);
    if (adapter->coordinatesAvailable()) {
      setCallbacksGeom(adapter->getCoordinateInput());
    }
  }

  AlgZoltan(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const MeshAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
    zoltanInit();
    zz = rcp(new Zoltan(mpicomm)); 
    setCallbacksIDs();
    setCallbacksGraph(adapter);
    //TODO:: check parameter list to see if hypergraph is needed. We dont want to build the model
    //       if we don't have to and we shouldn't as it can take a decent amount of time if the
    //       primary entity is copied
    setCallbacksHypergraph(adapter);
    setCallbacksGeom(&(*adapter));
  }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);
  // void color(const RCP<ColoringSolution<Adapter> > &solution);
};

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgZoltan<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;
  char paramstr[128];

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  sprintf(paramstr, "%lu", numGlobalParts);
  zz->Set_Param("NUM_GLOBAL_PARTS", paramstr);

  int wdim = adapter->getNumWeightsPerID();
  sprintf(paramstr, "%d", wdim);
  zz->Set_Param("OBJ_WEIGHT_DIM", paramstr);

  const Teuchos::ParameterList &pl = env->getParameters();

  double tolerance;
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("imbalance_tolerance");
  if (pe){
    char str[30];
    tolerance = pe->getValue<double>(&tolerance);
    sprintf(str, "%f", tolerance);
    zz->Set_Param("IMBALANCE_TOL", str);
  }
  
  pe = pl.getEntryPtr("partitioning_approach");
  if (pe){
    std::string approach;
    approach = pe->getValue<std::string>(&approach);
    if (approach == "partition")
      zz->Set_Param("LB_APPROACH", "PARTITION");
    else
      zz->Set_Param("LB_APPROACH", "REPARTITION");
  }

  pe = pl.getEntryPtr("partitioning_objective");
  if (pe){
    std::string strChoice = pe->getValue<std::string>(&strChoice);
    if (strChoice == std::string("multicriteria_minimize_total_weight"))
      zz->Set_Param("RCB_MULTICRITERIA_NORM", "1");
    else if (strChoice == std::string("multicriteria_balance_total_maximum"))
      zz->Set_Param("RCB_MULTICRITERIA_NORM", "2");
    else if (strChoice == std::string("multicriteria_minimize_maximum_weight"))
      zz->Set_Param("RCB_MULTICRITERIA_NORM", "3");
  }

  // perhaps make this a bool stored in the AlgZoltan if we want to follow
  // the pattern of multijagged mj_keep_part_boxes for example. However we can
  // collect the error straight from Zoltan if we attempt to access the tree
  // when we never stored it so that may not be necessary
  bool keep_partition_tree = false;
  pe = pl.getEntryPtr("keep_partition_tree");
  if (pe) {
    keep_partition_tree = pe->getValue(&keep_partition_tree);
    if (keep_partition_tree) {
      // need to resolve the organization of this
      // when do we want to use the zoltan parameters directly versus
      // using the zoltan2 parameters like this
      zz->Set_Param("KEEP_CUTS", "1");      // rcb zoltan setting
      zz->Set_Param("PHG_KEEP_TREE", "1");  // phg zoltan setting
    }
  }

  pe = pl.getEntryPtr("rectilinear");
  if (pe) {
    bool val = pe->getValue(&val);
    if (val)
      zz->Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
  }

  // Look for zoltan_parameters sublist; pass all zoltan parameters to Zoltan
  try {
    const Teuchos::ParameterList &zpl = pl.sublist("zoltan_parameters");
    for (ParameterList::ConstIterator iter =  zpl.begin();
                                      iter != zpl.end(); iter++) {
      const std::string &zname = pl.name(iter);
      // Convert the value to a string to pass to Zoltan
      std::string zval = pl.entry(iter).getValue(&zval);
      zz->Set_Param(zname.c_str(), zval.c_str());      
    }
  }
  catch (std::exception &e) {
    // No zoltan_parameters sublist found; no Zoltan parameters to register
  }

  // Get target part sizes
  int pdim = (wdim > 1 ? wdim : 1);
  for (int d = 0; d < pdim; d++) {
    if (!solution->criteriaHasUniformPartSizes(d)) {
      float *partsizes = new float[numGlobalParts];
      int *partidx = new int[numGlobalParts];
      int *wgtidx = new int[numGlobalParts];
      for (size_t i=0; i<numGlobalParts; i++) partidx[i] = i;
      for (size_t i=0; i<numGlobalParts; i++) wgtidx[i] = d;
      for (size_t i=0; i<numGlobalParts; i++)
        partsizes[i] = solution->getCriteriaPartSize(0, i);
      zz->LB_Set_Part_Sizes(1, numGlobalParts, partidx, wgtidx, partsizes);
      delete [] partsizes;
      delete [] partidx;
      delete [] wgtidx;
    }
  }

  // Make the call to LB_Partition
  int changed = 0;
  int nGidEnt = TPL_Traits<ZOLTAN_ID_PTR,gno_t>::NUM_ID;
  int nLidEnt = TPL_Traits<ZOLTAN_ID_PTR,lno_t>::NUM_ID;

  int nDummy = -1;                         // Dummy vars to satisfy arglist
  ZOLTAN_ID_PTR dGids = NULL, dLids = NULL;
  int *dProcs = NULL, *dParts = NULL;
  int nObj = -1;                           // Output vars with new part info
  ZOLTAN_ID_PTR oGids = NULL, oLids = NULL;
  int *oProcs = NULL, *oParts = NULL;

  zz->Set_Param("RETURN_LISTS", "PARTS");  // required format for Zoltan2;
                                           // results in last five arguments

  int ierr = zz->LB_Partition(changed, nGidEnt, nLidEnt,
                              nDummy, dGids, dLids, dProcs, dParts,
                              nObj,   oGids, oLids, oProcs, oParts);

  env->globalInputAssertion(__FILE__, __LINE__, "Zoltan LB_Partition", 
    (ierr==ZOLTAN_OK || ierr==ZOLTAN_WARN), BASIC_ASSERTION, problemComm);

  int numObjects=nObj;
  // The number of objects may be larger than zoltan knows due to copies that 
  // were removed by the hypergraph model
  if (model!=RCP<const Model<Adapter> >() &&
      dynamic_cast<const HyperGraphModel<Adapter>* >(&(*model)) &&
      !(dynamic_cast<const HyperGraphModel<Adapter>* >(&(*model))->areVertexIDsUnique())) {
    numObjects=model->getLocalNumObjects();
  }

  // Load answer into the solution.
  ArrayRCP<part_t> partList(new part_t[numObjects], 0, numObjects, true);
  for (int i = 0; i < nObj; i++) {
    lno_t tmp;
    TPL_Traits<lno_t, ZOLTAN_ID_PTR>::ASSIGN(tmp, &(oLids[i*nLidEnt]));
    partList[tmp] = oParts[i];
  }
  
  if (model!=RCP<const Model<Adapter> >() &&
      dynamic_cast<const HyperGraphModel<Adapter>* >(&(*model)) &&
      !(dynamic_cast<const HyperGraphModel<Adapter>* >(&(*model))->areVertexIDsUnique())) {
    // Setup the part ids for copied entities removed by ownership in 
    // hypergraph model.
    const HyperGraphModel<Adapter>* mdl = 
                    static_cast<const HyperGraphModel<Adapter>* >(&(*model));
    
    typedef typename HyperGraphModel<Adapter>::map_t map_t;
    Teuchos::RCP<const map_t> mapWithCopies;
    Teuchos::RCP<const map_t> oneToOneMap;
    mdl->getVertexMaps(mapWithCopies,oneToOneMap);
    
    typedef Tpetra::Vector<scalar_t, lno_t, gno_t> vector_t;
    vector_t vecWithCopies(mapWithCopies);
    vector_t oneToOneVec(oneToOneMap);

    // Set values in oneToOneVec:  each entry == rank
    assert(nObj == lno_t(oneToOneMap->getNodeNumElements()));
    for (lno_t i = 0; i < nObj; i++)
      oneToOneVec.replaceLocalValue(i, oParts[i]);
    
    // Now import oneToOneVec's values back to vecWithCopies
    Teuchos::RCP<const Tpetra::Import<lno_t, gno_t> > importer = 
      Tpetra::createImport<lno_t, gno_t>(oneToOneMap, mapWithCopies);
    vecWithCopies.doImport(oneToOneVec, *importer, Tpetra::REPLACE);

    // Should see copied vector values when print VEC WITH COPIES
    lno_t nlocal = lno_t(mapWithCopies->getNodeNumElements());
    for (lno_t i = 0; i < nlocal; i++)
      partList[i] = vecWithCopies.getData()[i];
  }
  
  solution->setParts(partList);

  // Clean up
  zz->LB_Free_Part(&oGids, &oLids, &oProcs, &oParts);
}

} // namespace Zoltan2

#endif
