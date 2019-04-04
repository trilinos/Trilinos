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
#ifndef _ZOLTAN2_ALGSCOTCH_HPP_
#define _ZOLTAN2_ALGSCOTCH_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_OrderingSolution.hpp> // BDD: needed by ordering method
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgScotch.hpp
//! \brief interface to the Scotch third-party library

////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

// this function called by the two scotch types below
static inline void getScotchParameters(Teuchos::ParameterList & pl)
{
  // bool parameter
  pl.set("scotch_verbose", false, "verbose output",
    Environment::getBoolValidator());

  RCP<Teuchos::EnhancedNumberValidator<int>> scotch_level_Validator =
    Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 1000, 1, 0) );
  pl.set("scotch_level", 0, "scotch ordering - Level of the subgraph in the "
    "separators tree for the initial graph at the root of the tree",
    scotch_level_Validator);

  pl.set("scotch_imbalance_ratio", 0.2, "scotch ordering - Dissection "
    "imbalance ratio", Environment::getAnyDoubleValidator());

  // bool parameter
  pl.set("scotch_ordering_default", true, "use default scotch ordering "
    "strategy", Environment::getBoolValidator());

  pl.set("scotch_ordering_strategy", "", "scotch ordering - Dissection "
    "imbalance ratio");
}

} // Zoltan2 namespace

#ifndef HAVE_ZOLTAN2_SCOTCH

namespace Zoltan2 {

// Error handling for when Scotch is requested
// but Zoltan2 not built with Scotch.

template <typename Adapter>
class AlgPTScotch : public Algorithm<Adapter>
{
public:
  typedef typename Adapter::base_adapter_t base_adapter_t;
  //AlgPTScotch(const RCP<const Environment> &env,
  //            const RCP<const Comm<int> > &problemComm,
  //            const RCP<GraphModel<typename Adapter::base_adapter_t> > &model
  //) BDD: old inteface for reference
  AlgPTScotch(const RCP<const Environment> &/* env */,
              const RCP<const Comm<int> > &/* problemComm */,
              const RCP<const base_adapter_t> &/* adapter */
  )
  {
    throw std::runtime_error(
          "BUILD ERROR:  Scotch requested but not compiled into Zoltan2.\n"
          "Please set CMake flag Zoltan2_ENABLE_Scotch:BOOL=ON.");
  }

  /*! \brief Set up validators specific to this algorithm
  */
  static void getValidParameters(ParameterList & pl)
  {
    getScotchParameters(pl);
  }
};
}
#endif

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ZOLTAN2_SCOTCH

namespace Zoltan2 {

// stdint.h for int64_t in scotch header
#include <stdint.h>
extern "C" {
#ifndef HAVE_ZOLTAN2_MPI
#include "scotch.h"
#else
#include "ptscotch.h"
#endif
}

#ifdef SHOW_ZOLTAN2_SCOTCH_MEMORY
//
// Scotch keeps track of memory high water mark, but doesn't
// provide a way to get that number.  So add this function:
//   "size_t SCOTCH_getMemoryMax() { return memorymax;}"
// to src/libscotch/common_memory.c
//
// and this macro:
//   "#define HAVE_SCOTCH_GETMEMORYMAX
// to include/ptscotch.h
//
// and compile scotch with -DCOMMON_MEMORY_TRACE
//
#ifdef HAVE_SCOTCH_GETMEMORYMAX
  extern "C"{
    extern size_t SCOTCH_getMemoryMax();
  }
#else
#error "Either turn off ZOLTAN2_ENABLE_SCOTCH_MEMORY_REPORT in cmake configure, or see SHOW_ZOLTAN2_SCOTCH_MEMORY comment in Zoltan2_AlgScotch.hpp"
#endif  // HAVE_SCOTCH_GETMEMORYMAX
#endif  // SHOW_ZOLTAN2_SCOTCH_MEMORY

template <typename Adapter>
class AlgPTScotch : public Algorithm<Adapter>
{
public:

  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef GraphModel<typename Adapter::base_adapter_t> graphModel_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

//  /*! Scotch constructor
//   *  \param env  parameters for the problem and library configuration
//   *  \param problemComm  the communicator for the problem
//   *  \param model a graph
//   *
//   *  Preconditions: The parameters in the environment have been processed.
//   *  TODO:  THIS IS A MINIMAL CONSTRUCTOR FOR NOW.
//   *  TODO:  WHEN ADD SCOTCH ORDERING OR MAPPING, MOVE SCOTCH GRAPH CONSTRUCTION
//   *  TODO:  TO THE CONSTRUCTOR SO THAT CODE MAY BE SHARED.
//   */
//  AlgPTScotch(const RCP<const Environment> &env__,
//              const RCP<const Comm<int> > &problemComm__,
//              const RCP<graphModel_t> &model__) :
//    env(env__), problemComm(problemComm__), 
//#ifdef HAVE_ZOLTAN2_MPI
//    mpicomm(Teuchos::getRawMpiComm(*problemComm__)),
//#endif
//    model(model__) BDD: olde interface for reference

  /*! Scotch constructor
   *  \param env          parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter      the user's input adapter
   *
   *  We're building a graph model, so throw an error if we can't  
   *    build the model from the input adapter passed to constructor
   *  For matrix and mesh adapters, additionally determine which 
   *    objects we wish to partition
   */
  AlgPTScotch(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<const IdentifierAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), 
#ifdef HAVE_ZOLTAN2_MPI
    mpicomm(Teuchos::getRawMpiComm(*problemComm__)),
#endif
    adapter(adapter__)
  {
    std::string errStr = "cannot build GraphModel from IdentifierAdapter, ";
    errStr            += "Scotch requires Graph, Matrix, or Mesh Adapter";
    throw std::runtime_error(errStr);
  }
  
  AlgPTScotch(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<const VectorAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), 
#ifdef HAVE_ZOLTAN2_MPI
    mpicomm(Teuchos::getRawMpiComm(*problemComm__)),
#endif
    adapter(adapter__)
  {
    std::string errStr = "cannot build GraphModel from IdentifierAdapter, ";
    errStr            += "Scotch requires Graph, Matrix, or Mesh Adapter";
    throw std::runtime_error(errStr);
  }
  
  AlgPTScotch(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<const GraphAdapter<user_t, userCoord_t> > &adapter__) :
    env(env__), problemComm(problemComm__), 
#ifdef HAVE_ZOLTAN2_MPI
    mpicomm(Teuchos::getRawMpiComm(*problemComm__)),
#endif
    adapter(adapter__), model_flags()
  {
    this->model_flags.reset();
  }
  
  AlgPTScotch(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<const MatrixAdapter<user_t, userCoord_t> > &adapter__) :
    env(env__), problemComm(problemComm__), 
#ifdef HAVE_ZOLTAN2_MPI
    mpicomm(Teuchos::getRawMpiComm(*problemComm__)),
#endif
    adapter(adapter__), model_flags()
  {
    this->model_flags.reset();
    
    const ParameterList &pl = env->getParameters();
    const Teuchos::ParameterEntry *pe;

    std::string defString("default");
    std::string objectOfInterest(defString);
    pe = pl.getEntryPtr("objects_to_partition");
    if (pe)
      objectOfInterest = pe->getValue<std::string>(&objectOfInterest);

    if (objectOfInterest == defString ||
        objectOfInterest == std::string("matrix_rows") )
      this->model_flags.set(VERTICES_ARE_MATRIX_ROWS);
    else if (objectOfInterest == std::string("matrix_columns"))
      this->model_flags.set(VERTICES_ARE_MATRIX_COLUMNS);
    else if (objectOfInterest == std::string("matrix_nonzeros"))
      this->model_flags.set(VERTICES_ARE_MATRIX_NONZEROS);
  }
  
  AlgPTScotch(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<const MeshAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), 
#ifdef HAVE_ZOLTAN2_MPI
    mpicomm(Teuchos::getRawMpiComm(*problemComm__)),
#endif
    adapter(adapter__), model_flags()
  {
    this->model_flags.reset();
    
    const ParameterList &pl = env->getParameters();
    const Teuchos::ParameterEntry *pe;

    std::string defString("default");
    std::string objectOfInterest(defString);
    pe = pl.getEntryPtr("objects_to_partition");
    if (pe)
      objectOfInterest = pe->getValue<std::string>(&objectOfInterest);

    if (objectOfInterest == defString ||
        objectOfInterest == std::string("mesh_nodes") )
      this->model_flags.set(VERTICES_ARE_MESH_NODES);
    else if (objectOfInterest == std::string("mesh_elements"))
      this->model_flags.set(VERTICES_ARE_MESH_ELEMENTS);
  }

  /*! \brief Set up validators specific to this algorithm
  */
  static void getValidParameters(ParameterList & pl)
  {
    getScotchParameters(pl);
  }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);
  
  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution);
  int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &solution);

private:

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
#ifdef HAVE_ZOLTAN2_MPI
  const MPI_Comm mpicomm;
#endif
  const RCP<const base_adapter_t> adapter;
  modelFlag_t model_flags;
  RCP<graphModel_t > model; // BDD:to be constructed

  void buildModel(bool isLocal);
  void scale_weights(size_t n, StridedData<lno_t, scalar_t> &fwgts,
                     SCOTCH_Num *iwgts);
  static int setStrategy(SCOTCH_Strat * c_strat_ptr, const ParameterList &pList, int procRank);
};

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgPTScotch<Adapter>::buildModel(bool isLocal) {
  HELLO;  
  const ParameterList &pl = env->getParameters();
  const Teuchos::ParameterEntry *pe;

  std::string defString("default");
  std::string symParameter(defString);
  pe = pl.getEntryPtr("symmetrize_graph");
  if (pe){
    symParameter = pe->getValue<std::string>(&symParameter);
    if (symParameter == std::string("transpose"))
      this->model_flags.set(SYMMETRIZE_INPUT_TRANSPOSE);
    else if (symParameter == std::string("bipartite"))
      this->model_flags.set(SYMMETRIZE_INPUT_BIPARTITE);  } 

  bool sgParameter = false;
  pe = pl.getEntryPtr("subset_graph");
  if (pe)
    sgParameter = pe->getValue(&sgParameter);
  if (sgParameter)
      this->model_flags.set(BUILD_SUBSET_GRAPH);

  this->model_flags.set(REMOVE_SELF_EDGES);
  this->model_flags.set(GENERATE_CONSECUTIVE_IDS);
  if (isLocal) {
    HELLO; // only for ordering!
    this->model_flags.set(BUILD_LOCAL_GRAPH);
  }

  this->env->debug(DETAILED_STATUS, "    building graph model");
  this->model = rcp(new graphModel_t(this->adapter, 
        this->env, 
        this->problemComm, 
        this->model_flags));
  this->env->debug(DETAILED_STATUS, "    graph model built");
}

template <typename Adapter>
void AlgPTScotch<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;
  this->buildModel(false);

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  SCOTCH_Num partnbr=0;
  TPL_Traits<SCOTCH_Num, size_t>::ASSIGN(partnbr, numGlobalParts);

#ifdef HAVE_ZOLTAN2_MPI
  int ierr = 0;
  int me = problemComm->getRank();

  const SCOTCH_Num  baseval = 0;  // Base value for array indexing.
                                  // GraphModel returns GNOs from base 0.

  SCOTCH_Strat stratstr;          // Strategy string
                                  // TODO:  Set from parameters
  SCOTCH_stratInit(&stratstr);

  // Allocate and initialize PTScotch Graph data structure.
  SCOTCH_Dgraph *gr = SCOTCH_dgraphAlloc();  // Scotch distributed graph
  ierr = SCOTCH_dgraphInit(gr, mpicomm);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphInit", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Get vertex info
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  size_t nVtx = model->getVertexList(vtxID, vwgts);
  SCOTCH_Num vertlocnbr=0;
  TPL_Traits<SCOTCH_Num, size_t>::ASSIGN(vertlocnbr, nVtx);
  SCOTCH_Num vertlocmax = vertlocnbr; // Assumes no holes in global nums.

  // Get edge info
  ArrayView<const gno_t> edgeIds;
  ArrayView<const offset_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;

  size_t nEdge = model->getEdgeList(edgeIds, offsets, ewgts);

  SCOTCH_Num edgelocnbr=0;
  TPL_Traits<SCOTCH_Num, size_t>::ASSIGN(edgelocnbr, nEdge);
  const SCOTCH_Num edgelocsize = edgelocnbr;  // Assumes adj array is compact.

  SCOTCH_Num *vertloctab;  // starting adj/vtx
  TPL_Traits<SCOTCH_Num, const offset_t>::ASSIGN_ARRAY(&vertloctab, offsets);

  SCOTCH_Num *edgeloctab;  // adjacencies
  TPL_Traits<SCOTCH_Num, const gno_t>::ASSIGN_ARRAY(&edgeloctab, edgeIds);

  // We don't use these arrays, but we need them as arguments to Scotch.
  SCOTCH_Num *vendloctab = NULL;  // Assume consecutive storage for adj
  SCOTCH_Num *vlblloctab = NULL;  // Vertex label array
  SCOTCH_Num *edgegsttab = NULL;  // Array for ghost vertices

  // Get weight info.
  SCOTCH_Num *velotab = NULL;  // Vertex weights
  SCOTCH_Num *edlotab = NULL;  // Edge weights

  int nVwgts = model->getNumWeightsPerVertex();
  int nEwgts = model->getNumWeightsPerEdge();
  if (nVwgts > 1 && me == 0) {
    std::cerr << "Warning:  NumWeightsPerVertex is " << nVwgts 
              << " but Scotch allows only one weight. "
              << " Zoltan2 will use only the first weight per vertex."
              << std::endl;
  }
  if (nEwgts > 1 && me == 0) {
    std::cerr << "Warning:  NumWeightsPerEdge is " << nEwgts 
              << " but Scotch allows only one weight. "
              << " Zoltan2 will use only the first weight per edge."
              << std::endl;
  }

  if (nVwgts) {
    velotab = new SCOTCH_Num[nVtx+1];  // +1 since Scotch wants all procs 
                                       // to have non-NULL arrays
    scale_weights(nVtx, vwgts[0], velotab);
  }

  if (nEwgts) {
    edlotab = new SCOTCH_Num[nEdge+1];  // +1 since Scotch wants all procs 
                                         // to have non-NULL arrays
    scale_weights(nEdge, ewgts[0], edlotab);
  }

  // Build PTScotch distributed data structure
  ierr = SCOTCH_dgraphBuild(gr, baseval, vertlocnbr, vertlocmax,
                            vertloctab, vendloctab, velotab, vlblloctab,
                            edgelocnbr, edgelocsize,
                            edgeloctab, edgegsttab, edlotab);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphBuild", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Create array for Scotch to return results in.
  SCOTCH_Num *partloctab = new SCOTCH_Num[nVtx+1];
    // Note: Scotch does not like NULL arrays, so add 1 to always have non-null.
    //       ParMETIS has this same "feature."  See Zoltan bug 4299.

  // Get target part sizes
  float *partsizes = new float[numGlobalParts];
  if (!solution->criteriaHasUniformPartSizes(0))
    for (size_t i=0; i<numGlobalParts; i++)
      partsizes[i] = solution->getCriteriaPartSize(0, i);
  else
    for (size_t i=0; i<numGlobalParts; i++)
      partsizes[i] = 1.0 / float(numGlobalParts);

  // Allocate and initialize PTScotch target architecture data structure
  SCOTCH_Arch archdat;
  SCOTCH_archInit(&archdat);

  SCOTCH_Num velosum = 0;
  SCOTCH_dgraphSize (gr, &velosum, NULL, NULL, NULL);
  SCOTCH_Num *goalsizes = new SCOTCH_Num[partnbr];
  // TODO: The goalsizes are set as in Zoltan; not sure it is correct there 
  // or here.
  // It appears velosum is global NUMBER of vertices, not global total 
  // vertex weight.  I think we should use the latter.
  // Fix this when we add vertex weights.
  for (SCOTCH_Num i = 0; i < partnbr; i++)
    goalsizes[i] = SCOTCH_Num(ceil(velosum * partsizes[i]));
  delete [] partsizes;

  SCOTCH_archCmpltw(&archdat, partnbr, goalsizes);

  // Call partitioning; result returned in partloctab.
  ierr = SCOTCH_dgraphMap(gr, &archdat, &stratstr, partloctab);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphMap", 
    !ierr, BASIC_ASSERTION, problemComm);

  SCOTCH_archExit(&archdat);
  delete [] goalsizes;

  // TODO - metrics

#ifdef SHOW_ZOLTAN2_SCOTCH_MEMORY
  int me = env->comm_->getRank();
#endif

#ifdef HAVE_SCOTCH_ZOLTAN2_GETMEMORYMAX
  if (me == 0){
    size_t scotchBytes = SCOTCH_getMemoryMax();
    std::cout << "Rank " << me << ": Maximum bytes used by Scotch: ";
    std::cout << scotchBytes << std::endl;
  }
#endif

  // Clean up PTScotch
  SCOTCH_dgraphExit(gr);
  free(gr);
  SCOTCH_stratExit(&stratstr);

  // Load answer into the solution.

  ArrayRCP<part_t> partList;
  if (nVtx > 0)
    TPL_Traits<part_t, SCOTCH_Num>::SAVE_ARRAYRCP(&partList, partloctab, nVtx);
  TPL_Traits<SCOTCH_Num, part_t>::DELETE_ARRAY(&partloctab);

  solution->setParts(partList);

  env->memory("Zoltan2-Scotch: After creating solution");

  // Clean up copies made due to differing data sizes.
  TPL_Traits<SCOTCH_Num, lno_t>::DELETE_ARRAY(&vertloctab);
  TPL_Traits<SCOTCH_Num, gno_t>::DELETE_ARRAY(&edgeloctab);

  if (nVwgts) delete [] velotab;
  if (nEwgts) delete [] edlotab;

#else // DO NOT HAVE MPI

  // TODO:  Handle serial case with calls to Scotch.
  // TODO:  For now, assign everything to rank 0 and assume only one part.
  // TODO:  Can probably use the code above for loading solution,
  // TODO:  instead of duplicating it here.
  // TODO
  // TODO:  Actual logic should call Scotch when number of processes == 1.
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  size_t nVtx = model->getVertexList(vtxID, vwgts);

  ArrayRCP<part_t> partList(new part_t[nVtx], 0, nVtx, true);
  for (size_t i = 0; i < nVtx; i++) partList[i] = 0;

  solution->setParts(partList);

#endif // DO NOT HAVE MPI
}

/////////////////////////////////////////////////////////////////////////////
// Scale and round scalar_t weights (typically float or double) to 
// SCOTCH_Num (typically int or long).
// subject to sum(weights) <= max_wgt_sum.
// Only scale if deemed necessary.
//
// Note that we use ceil() instead of round() to avoid
// rounding to zero weights.
// Based on Zoltan's scale_round_weights, mode 1.

template <typename Adapter>
void AlgPTScotch<Adapter>::scale_weights(
  size_t n,
  StridedData<typename Adapter::lno_t, typename Adapter::scalar_t> &fwgts,
  SCOTCH_Num *iwgts
)
{
  const double INT_EPSILON = 1e-5;

  SCOTCH_Num nonint, nonint_local = 0;
  double sum_wgt, sum_wgt_local = 0.;
  double max_wgt, max_wgt_local = 0.;

  // Compute local sums of the weights 
  // Check whether all weights are integers
  for (size_t i = 0; i < n; i++) {
    double fw = double(fwgts[i]);
    if (!nonint_local){
      SCOTCH_Num tmp = (SCOTCH_Num) floor(fw + .5); /* Nearest int */
      if (fabs((double)tmp-fw) > INT_EPSILON) {
        nonint_local = 1;
      }
    }
    sum_wgt_local += fw;
    if (fw > max_wgt_local) max_wgt_local = fw;
  }

  Teuchos::reduceAll<int,SCOTCH_Num>(*problemComm, Teuchos::REDUCE_MAX, 1,
                              &nonint_local,  &nonint);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_SUM, 1,
                                 &sum_wgt_local, &sum_wgt);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_MAX, 1,
                                 &max_wgt_local, &max_wgt);

  double scale = 1.;
  const double max_wgt_sum = double(SCOTCH_NUMMAX/8);

  // Scaling needed if weights are not integers or weights' 
  // range is not sufficient
  if (nonint || (max_wgt <= INT_EPSILON) || (sum_wgt > max_wgt_sum)) {
    /* Calculate scale factor */
    if (sum_wgt != 0.) scale = max_wgt_sum/sum_wgt;
  }

  /* Convert weights to positive integers using the computed scale factor */
  for (size_t i = 0; i < n; i++)
    iwgts[i] = (SCOTCH_Num) ceil(double(fwgts[i])*scale);

}

template <typename Adapter>
int AlgPTScotch<Adapter>::setStrategy(SCOTCH_Strat * c_strat_ptr, const ParameterList &pList, int procRank) {
  // get ordering parameters from parameter list
  bool bPrintOutput = false; // will be true if rank 0 and verbose is true
  const Teuchos::ParameterEntry *pe;

  if (procRank == 0) {
    pe = pList.getEntryPtr("scotch_verbose");
    if (pe) {
      bPrintOutput = pe->getValue(&bPrintOutput);
    }
  }

  // get parameter setting for ordering default true/false
  bool scotch_ordering_default = true;
  pe = pList.getEntryPtr("scotch_ordering_default");
  if (pe) {
    scotch_ordering_default = pe->getValue(&scotch_ordering_default);
  }
  if (bPrintOutput) {
    std::cout <<
      "Scotch: Ordering default setting (parameter: scotch_ordering_default): "
      << (scotch_ordering_default ? "True" : "False" ) << std::endl;
  }

  // set up error code for return
  int ierr = 1; // will be set 0 if successful

  if (scotch_ordering_default) {
    // get parameter scotch_level
    int scotch_level = 0;
    pe = pList.getEntryPtr("scotch_level");
    if (pe) {
      scotch_level = pe->getValue(&scotch_level);
    }
    if (bPrintOutput) {
      std::cout << "Scotch: Ordering level (parameter: scotch_level): " <<
        scotch_level << std::endl;
    }

    // get parameter scotch_imbalance_ratio
    double scotch_imbalance_ratio = 0.2;
    pe = pList.getEntryPtr("scotch_imbalance_ratio");
    if (pe) {
      scotch_imbalance_ratio = pe->getValue(&scotch_imbalance_ratio);
    }
    if (bPrintOutput) {
      std::cout << "Scotch: Ordering dissection imbalance ratio "
        "(parameter: scotch_imbalance_ratio): "
        << scotch_imbalance_ratio << std::endl;
    }

    SCOTCH_stratInit(c_strat_ptr);

    ierr = SCOTCH_stratGraphOrderBuild( c_strat_ptr,
        SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN |
        SCOTCH_STRATLEAFSIMPLE | SCOTCH_STRATSEPASIMPLE,
        scotch_level, scotch_imbalance_ratio); // based on Siva's example
  }
  else {
    // get parameter scotch_ordering_strategy
    std::string scotch_ordering_strategy_string("");
    pe = pList.getEntryPtr("scotch_ordering_strategy");
    if (pe) {
      scotch_ordering_strategy_string =
        pe->getValue(&scotch_ordering_strategy_string);
    }
    if (bPrintOutput) {
      std::cout << "Scotch ordering strategy"
        " (parameter: scotch_ordering_strategy): " <<
        scotch_ordering_strategy_string << std::endl;
    }

    SCOTCH_stratInit(c_strat_ptr);
    ierr = SCOTCH_stratGraphOrder( c_strat_ptr,
      scotch_ordering_strategy_string.c_str());
  }
  
  return ierr;
} 

template <typename Adapter>
int AlgPTScotch<Adapter>::globalOrder(
    const RCP<GlobalOrderingSolution<gno_t> > &solution) {
  throw std::logic_error("AlgPTScotch does not yet support global ordering.");
}

template <typename Adapter>
int AlgPTScotch<Adapter>::localOrder(
    const RCP<LocalOrderingSolution<lno_t> > &solution) {
  // TODO: translate teuchos sublist parameters to scotch strategy string
  // TODO: construct local graph model
  // TODO: solve with scotch
  // TODO: set solution fields

  HELLO; // say hi so that we know we have called this method
  int me = problemComm->getRank();
  const ParameterList &pl = env->getParameters();
  const Teuchos::ParameterEntry *pe;

  bool isVerbose = false;
  pe = pl.getEntryPtr("scotch_verbose");
  if (pe)
    isVerbose = pe->getValue(&isVerbose);
  
  // build a local graph model
  this->buildModel(true);
  if (isVerbose && me==0) {
    std::cout << "Built local graph model." << std::endl;
  }

  // based off of Siva's example
  SCOTCH_Strat c_strat_ptr; // strategy
  SCOTCH_Graph c_graph_ptr; // pointer to scotch graph
  int ierr;

  ierr = SCOTCH_graphInit( &c_graph_ptr);
  if ( ierr != 0) {
    throw std::runtime_error("Failed to initialize Scotch graph!");
  } else if (isVerbose && me == 0) {
    std::cout << "Initialized Scotch graph." << std::endl;
  }

  // Get vertex info
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  size_t nVtx = model->getVertexList(vtxID, vwgts);
  SCOTCH_Num vertnbr=0;
  TPL_Traits<SCOTCH_Num, size_t>::ASSIGN(vertnbr, nVtx);

  // Get edge info
  ArrayView<const gno_t> edgeIds;
  ArrayView<const offset_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;

  size_t nEdge = model->getEdgeList(edgeIds, offsets, ewgts);
  SCOTCH_Num edgenbr=0;
  TPL_Traits<SCOTCH_Num, size_t>::ASSIGN(edgenbr, nEdge);
  
  SCOTCH_Num *verttab;  // starting adj/vtx
  TPL_Traits<SCOTCH_Num, const offset_t>::ASSIGN_ARRAY(&verttab, offsets);

  SCOTCH_Num *edgetab;  // adjacencies
  TPL_Traits<SCOTCH_Num, const gno_t>::ASSIGN_ARRAY(&edgetab, edgeIds);
  
  // We don't use these arrays, but we need them as arguments to Scotch.
  SCOTCH_Num *vendtab = NULL;  // Assume consecutive storage for adj
  //SCOTCH_Num *vendtab = verttab+1;  // Assume consecutive storage for adj
  // Get weight info.
  SCOTCH_Num *velotab = NULL;  // Vertex weights
  SCOTCH_Num *vlbltab = NULL;  // vertes labels
  SCOTCH_Num *edlotab = NULL;  // Edge weights
  
  int nVwgts = model->getNumWeightsPerVertex();
  int nEwgts = model->getNumWeightsPerEdge();
  if (nVwgts > 1 && me == 0) {
    std::cerr << "Warning:  NumWeightsPerVertex is " << nVwgts 
              << " but Scotch allows only one weight. "
              << " Zoltan2 will use only the first weight per vertex."
              << std::endl;
  }
  if (nEwgts > 1 && me == 0) {
    std::cerr << "Warning:  NumWeightsPerEdge is " << nEwgts 
              << " but Scotch allows only one weight. "
              << " Zoltan2 will use only the first weight per edge."
              << std::endl;
  }
  
  if (nVwgts) {
    velotab = new SCOTCH_Num[nVtx+1];  // +1 since Scotch wants all procs 
                                       // to have non-NULL arrays
    scale_weights(nVtx, vwgts[0], velotab);
  }

  if (nEwgts) {
    edlotab = new SCOTCH_Num[nEdge+1];  // +1 since Scotch wants all procs 
                                         // to have non-NULL arrays
    scale_weights(nEdge, ewgts[0], edlotab);
  }

  // build scotch graph
  int baseval = 0;
  ierr = 1;
  ierr = SCOTCH_graphBuild( &c_graph_ptr, baseval,
                            vertnbr, verttab, vendtab, velotab, vlbltab,
                            edgenbr, edgetab, edlotab);
  if (ierr != 0) {
    throw std::runtime_error("Failed to build Scotch graph!");
  } else if (isVerbose && me == 0) {
    std::cout << "Built Scotch graph." << std::endl;
  }

  ierr = SCOTCH_graphCheck(&c_graph_ptr);
  if (ierr != 0) {
    throw std::runtime_error("Graph is inconsistent!");
  } else if (isVerbose && me == 0) {
    std::cout << "Graph is consistent." << std::endl;
  }

  // set the strategy 
  ierr = AlgPTScotch<Adapter>::setStrategy(&c_strat_ptr, pl, me); 

  if (ierr != 0) {
    throw std::runtime_error("Can't build strategy!");
  }else if (isVerbose && me == 0) {
    std::cout << "Graphing strategy built." << std::endl;
  }

  // Allocate results
  SCOTCH_Num cblk = 0;
  SCOTCH_Num *permtab;  // permutation array
  SCOTCH_Num *peritab;  // inverse permutation array
  SCOTCH_Num *rangetab; // separator range array
  SCOTCH_Num *treetab;  // separator tree

  if (TPL_Traits<lno_t, SCOTCH_Num>::OK_TO_CAST()) {
    permtab = reinterpret_cast<SCOTCH_Num*>(solution->getPermutationView(false));
    peritab = reinterpret_cast<SCOTCH_Num*>(solution->getPermutationView(true));
    rangetab = reinterpret_cast<SCOTCH_Num*>(solution->getSeparatorRangeView());
    treetab = reinterpret_cast<SCOTCH_Num*>(solution->getSeparatorTreeView());
  }
  else {
    permtab = new SCOTCH_Num[nVtx];
    peritab = new SCOTCH_Num[nVtx];
    rangetab = new SCOTCH_Num[nVtx+1];
    treetab = new SCOTCH_Num[nVtx];
  }

  ierr = SCOTCH_graphOrder(&c_graph_ptr, &c_strat_ptr, permtab, peritab, 
                           &cblk, rangetab, treetab);
  if (ierr != 0) {
    throw std::runtime_error("Could not compute ordering!!");
  } else if(isVerbose && me == 0) {
    std::cout << "Ordering computed." << std::endl;
  }
  
  lno_t nSepBlocks;
  TPL_Traits<lno_t, SCOTCH_Num>::ASSIGN(nSepBlocks, cblk);
  solution->setNumSeparatorBlocks(nSepBlocks);

  if (!TPL_Traits<lno_t, SCOTCH_Num>::OK_TO_CAST()) {

    const ArrayRCP<lno_t> arv_perm = solution->getPermutationRCP(false);
    for (size_t i = 0; i < nVtx; i++)
      TPL_Traits<lno_t, SCOTCH_Num>::ASSIGN(arv_perm[i], permtab[i]);
    delete [] permtab;

    const ArrayRCP<lno_t> arv_peri = solution->getPermutationRCP(true);
    for (size_t i = 0; i < nVtx; i++)
      TPL_Traits<lno_t, SCOTCH_Num>::ASSIGN(arv_peri[i], peritab[i]);
    delete [] peritab;

    const ArrayRCP<lno_t> arv_range = solution->getSeparatorRangeRCP();
    for (lno_t i = 0; i <= nSepBlocks; i++)
      TPL_Traits<lno_t, SCOTCH_Num>::ASSIGN(arv_range[i], rangetab[i]);
    delete [] rangetab;

    const ArrayRCP<lno_t> arv_tree = solution->getSeparatorTreeRCP();
    for (lno_t i = 0; i < nSepBlocks; i++)
      TPL_Traits<lno_t, SCOTCH_Num>::ASSIGN(arv_tree[i], treetab[i]);
    delete [] treetab;
  }

  solution->setHaveSeparator(true); 

  // reclaim memory
  // Clean up copies made due to differing data sizes.
  TPL_Traits<SCOTCH_Num, lno_t>::DELETE_ARRAY(&verttab);
  TPL_Traits<SCOTCH_Num, gno_t>::DELETE_ARRAY(&edgetab);

  if (nVwgts) delete [] velotab;
  if (nEwgts) delete [] edlotab;

  SCOTCH_stratExit(&c_strat_ptr);
  SCOTCH_graphFree(&c_graph_ptr);

  if (isVerbose && problemComm->getRank() == 0) {
    std::cout << "Freed Scotch graph!" << std::endl;
  }
  return 0;
}

} // namespace Zoltan2

#endif // HAVE_ZOLTAN2_SCOTCH

////////////////////////////////////////////////////////////////////////

#endif
