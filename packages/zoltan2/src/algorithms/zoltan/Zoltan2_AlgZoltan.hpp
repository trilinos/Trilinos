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
    const RCP<const GraphAdapter<user_t,userCoord_t> > &adp)
  {
    // std::cout << "NotReadyForGraphYet" << std::endl;
    // TODO
  }

  void setCallbacksGraph(
    const RCP<const MatrixAdapter<user_t,userCoord_t> > &adp)
  {
    // std::cout << "NotReadyForGraphYet" << std::endl;
    // TODO
  }

  void setCallbacksGraph(
    const RCP<const MeshAdapter<user_t> > &adp)
  {
    // std::cout << "NotReadyForGraphYet" << std::endl;
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
  
  void loadRCBPartitionTree(const RCP<PartitioningSolution<Adapter> > &solution)
  {
    // This method is responsible for accessing the rcb formatted partition
    // tree and converting it to a general format all algorithms can use.
    // This would only be called if param "keep_partition_tree" was turned on.

    // For clarity - summarize the conventions - we convert from rcb to general:
    // From RCB:
    //   Each node has two children (binary) called left_leaf or right_leaf
    //   The 0 index node is not part of the tree - allows rcb sign conventions.
    //   A leaf <= 0 means it points to a terminal with part = -leaf
    //   Parent < 0 meant node was the left leaf of the parent. So use -parent.
    //   Parent by convention each node index + 1 so it starts at index 1
    //   However rcb also adds a second + 1 to parent
    //
    // To General:
    //   Each node has 2 or more children nodes - can be non-binary (MJ for ex.)
    //   The children are the indices into the array.
    //   Each node has a parent - also an index into the array.

    // design decision - original plan includes the trivial nodes which add
    // a tree node which points simply to the part ID. In rcb conventions this
    // tree node is not included - the terminal branch points to 2 parts.
    // This flag controls whether to insert them here. They are added at this
    // stage so that all follow up methods on the tree can act naturally
    // without worrying about the convention but it makes the loading a bit
    // convoluted. We have to insert the extra tree nodes at the beginning,
    // offset all node parent and leaf indices to account for the change, and
    // also make sure parent mapping is correct for these trivial nodes.
    bool bLoadTrivialNodes = true;

    // get the total number of global parts
    part_t num_parts = static_cast<part_t>(solution->getTargetGlobalNumberOfParts());

    // Determine how many nodes are in the rcb tree
    part_t numTreeNodes = num_parts - 1; // always true for rcb binary?
    if(bLoadTrivialNodes) { // if using trivial nodes we insert the extra
      numTreeNodes += num_parts; // 1 trivial node for each part is inserted
    }

    // Allocate partitionTreeNode vector and we will set all values below
    std::vector<partitionTreeNode<part_t>> partitionTreeNodes(numTreeNodes);

    // Add the rcb nodes to the general tree, converting conventions as we go
    for(part_t n = 0; n < numTreeNodes; ++n) {
      // is this a trivial node
      // rcb doesn't include these
      // we won't know the parent until we load the parent
      bool bTrivialNode = (bLoadTrivialNodes && n < num_parts);

      if(bTrivialNode) {
        // node parent will be determined when we load it below
        // the trivial terminal node will have it's own node index
        // and also be point to a single child - the part
        // We may need to set a bool flag instead of the neg convnetion
        ArrayRCP<part_t> children(1);
        children[0] = -n;
        partitionTreeNodes[n].children = children;
      }
      else {
        // call zoltan wrapper to get the node information - starts at +1 by conv
        // these are int types because zoltan is defined as int
        int parent = -1;
        int left_leaf = -1;
        int right_leaf = -1;

        int rcbIndex = static_cast<int>(n) + 1; // rcb starts at 1
        if(bLoadTrivialNodes) {
          rcbIndex -= num_parts; // shift back to account for extra trivial nodes
        }

        zz->RCB_Partition_Tree(rcbIndex, parent, left_leaf, right_leaf);

        // convert parent from rcb conventions to our generalized convention
        parent = abs(parent); // remove negative it exists
        parent -= 1; // remove +1 shift from rcb leaving parent=arrayIndex+1

        // if using trivial nodes we need to offset node indices += num_parts
        // also we need to offset parts ids += num_parts
        if(bLoadTrivialNodes) {
          // parent just shifts up by num_parts - though conv is parent
          // for leaf we need to shift in the sign direction because currently
          // we are keeping the convention that a negative sign indicates a
          // part id (terminal), otherwise it's a node index - however in both
          // cases we shift by num_parts
          // part id is determined by leaf <= 0 -> so 0 means shift negative
          if(parent != 0) { // do not shift the root - preserve convention
            parent += num_parts; // offset for parent is easier
          }
          left_leaf += (left_leaf > 0) ? num_parts : -num_parts;
          right_leaf += (right_leaf > 0) ? num_parts : -num_parts;
        }

        // set node parent
        partitionTreeNodes[n].parent = parent;

        // set node children
        ArrayRCP<part_t> children(2);
        children[0] = static_cast<part_t>(left_leaf);
        children[1] = static_cast<part_t>(right_leaf);
        partitionTreeNodes[n].children = children;

        // special case - if we are the parent of a trivial node tell them
        // terminal means leaf <= 0 in which case partID = -leaf
        // but trivials are also inserted first so node index = partID
        // and parent convention is node index + 1, so parent = n + 1
        if(bLoadTrivialNodes) {
          if(left_leaf <= 0) { // is left terminal?
            partitionTreeNodes[-left_leaf - num_parts].parent = n + 1; // [CONVENTION]
            partitionTreeNodes[n].children[0] = -left_leaf - num_parts + 1;
          }
          if(right_leaf <= 0) { // is right terminal?
            partitionTreeNodes[-right_leaf - num_parts].parent = n + 1; // [CONVENTION]
            partitionTreeNodes[n].children[1] = -right_leaf - num_parts + 1;
          }
        }
      }
    }

    // For debugging output the tree with the children for each
    /*
    for(part n = 0; n < static_cast<part_t>(partitionTreeNodes.size()); ++n) {
      const partitionTreeNode<part_t> & node = partitionTreeNodes[n];
      std::cout << "General index: " << n
        << " parent: " << partitionTreeNodes[n].parent << " ";
      for(int c = 0; c < static_cast<int>(node.children.size()); ++c) {
        std::cout << "Child" << c+1 << ": " << node.children[c] << " ";
      }
      std::cout << std::endl; // end the child lines
      std::cout << std::endl; // space between lines
    }
    */

    // might want to avoid this copy - not sure if this is important
    this->setPartitionTreeNodes(partitionTreeNodes);
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

  bool bKeepPartitionTree = false;
  pe = pl.getEntryPtr("keep_partition_tree");
  if (pe) {
    bKeepPartitionTree = pe->getValue(&bKeepPartitionTree);
    if (bKeepPartitionTree)
      zz->Set_Param("KEEP_CUTS", "1");
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
  
  // if the param keep_partition_tree is not set true then the tree
  // will not be available and not loaded - attempts to obtain it will throw
  if(bKeepPartitionTree) {
    loadRCBPartitionTree(solution);
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
