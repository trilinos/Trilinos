// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGPULP_HPP_
#define _ZOLTAN2_ALGPULP_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgPuLP.hpp
//! \brief interface to the PuLP third-party library

////////////////////////////////////////////////////////////////////////
#ifndef HAVE_ZOLTAN2_PULP

namespace Zoltan2 {
// Error handling for when PuLP is requested
// but Zoltan2 not built with PuLP.

template <typename Adapter>
class AlgPuLP : public Algorithm<Adapter>
{
public:
  typedef typename Adapter::base_adapter_t base_adapter_t;
  AlgPuLP(const RCP<const Environment> &/* env */,
          const RCP<const Comm<int> > &/* problemComm */,
          const RCP<const base_adapter_t> &/* adapter */
  )
  {
    throw std::runtime_error(
          "BUILD ERROR:  PuLP requested but not compiled into Zoltan2.\n"
          "Please set CMake flag TPL_ENABLE_PuLP:BOOL=ON.");
  }

  /*! \brief Set up validators specific to this algorithm
  */
  static void getValidParameters(ParameterList & /* pl */)
  {
    // No parameters needed in this error-handling version of AlgPuLP
  }
};

}
#endif

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ZOLTAN2_PULP

namespace Zoltan2 {

extern "C" {
// TODO: XtraPuLP
#ifndef HAVE_ZOLTAN2_MPI
#include "pulp.h"
#else
#include "xtrapulp.h"
#endif
}


template <typename Adapter>
class AlgPuLP : public Algorithm<Adapter>
{
public:
  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

  /*! PuLP constructors
   *  \param env          parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter      the user's input adapter
   * 
   *  We're building a graph model, so throw an error if we can't  
   *    build the model from the input adapter passed to constructor
   *  For matrix and mesh adapters, additionally determine which 
   *    objects we wish to partition
   */
  AlgPuLP(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const IdentifierAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    std::string errStr = "cannot build GraphModel from IdentifierAdapter, ";
    errStr            += "PuLP requires Graph, Matrix, or Mesh Adapter";
    throw std::runtime_error(errStr);
  }  

  AlgPuLP(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const VectorAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    std::string errStr = "cannot build GraphModel from VectorAdapter, ";
    errStr            += "PuLP requires Graph, Matrix, or Mesh Adapter";
    throw std::runtime_error(errStr);
  }   

  AlgPuLP(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    modelFlag_t flags;
    flags.reset();

    buildModel(flags);
  }  

  AlgPuLP(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const MatrixAdapter<user_t,userCoord_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  {   
    modelFlag_t flags;
    flags.reset();

    const ParameterList &pl = env->getParameters();
    const Teuchos::ParameterEntry *pe;

    std::string defString("default");
    std::string objectOfInterest(defString);
    pe = pl.getEntryPtr("objects_to_partition");
    if (pe)
      objectOfInterest = pe->getValue<std::string>(&objectOfInterest);

    if (objectOfInterest == defString ||
        objectOfInterest == std::string("matrix_rows") )
      flags.set(VERTICES_ARE_MATRIX_ROWS);
    else if (objectOfInterest == std::string("matrix_columns"))
      flags.set(VERTICES_ARE_MATRIX_COLUMNS);
    else if (objectOfInterest == std::string("matrix_nonzeros"))
      flags.set(VERTICES_ARE_MATRIX_NONZEROS);

    buildModel(flags);
  }

  AlgPuLP(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const MeshAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    modelFlag_t flags;
    flags.reset();

    const ParameterList &pl = env->getParameters();
    const Teuchos::ParameterEntry *pe;

    std::string defString("default");
    std::string objectOfInterest(defString);
    pe = pl.getEntryPtr("objects_to_partition");
    if (pe)
      objectOfInterest = pe->getValue<std::string>(&objectOfInterest);

    if (objectOfInterest == defString ||
        objectOfInterest == std::string("mesh_nodes") )
      flags.set(VERTICES_ARE_MESH_NODES);
    else if (objectOfInterest == std::string("mesh_elements"))
      flags.set(VERTICES_ARE_MESH_ELEMENTS);

    buildModel(flags);
  }

  /*! \brief Set up validators specific to this algorithm
  */
  static void getValidParameters(ParameterList & pl)
  {
    pl.set("pulp_vert_imbalance", 1.1, "vertex imbalance tolerance, ratio of "
      "maximum load over average load",
      Environment::getAnyDoubleValidator());

    pl.set("pulp_edge_imbalance", 1.1, "edge imbalance tolerance, ratio of "
      "maximum load over average load",
      Environment::getAnyDoubleValidator());

    pl.set("pulp_imbalance", 1.1, "multiweight imbalance tolerance, ratio of "
      "maximum load over average load",
      Environment::getAnyDoubleValidator());

    // bool parameter
    pl.set("pulp_lp_init", false, "perform label propagation-based "
      "initialization", Environment::getBoolValidator() );

    // bool parameter
    pl.set("pulp_minimize_maxcut", false, "perform per-part max cut "
      "minimization", Environment::getBoolValidator() );

    // bool parameter
    pl.set("pulp_verbose", false, "verbose output",
      Environment::getBoolValidator() );

    // bool parameter
    pl.set("pulp_do_repart", false, "perform repartitioning",
      Environment::getBoolValidator() );

    pl.set("pulp_seed", 0, "set pulp seed", Environment::getAnyIntValidator());
  }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  void buildModel(modelFlag_t &flags);

  int* scale_weights( size_t n, int nWgt,
                      ArrayView<StridedData<lno_t, scalar_t> > fwgts);

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const base_adapter_t> adapter;
  RCP<const GraphModel<base_adapter_t> > model;
};


/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgPuLP<Adapter>::buildModel(modelFlag_t &flags)
{   
  const ParameterList &pl = env->getParameters();
  const Teuchos::ParameterEntry *pe;

  std::string defString("default");
  std::string symParameter(defString);
  pe = pl.getEntryPtr("symmetrize_graph");
  if (pe){
    symParameter = pe->getValue<std::string>(&symParameter);
    if (symParameter == std::string("transpose"))
      flags.set(SYMMETRIZE_INPUT_TRANSPOSE);
    else if (symParameter == std::string("bipartite"))
      flags.set(SYMMETRIZE_INPUT_BIPARTITE);  } 

  bool sgParameter = false;
  pe = pl.getEntryPtr("subset_graph");
  if (pe)
    sgParameter = pe->getValue(&sgParameter);
  if (sgParameter)
      flags.set(BUILD_SUBSET_GRAPH);

  flags.set(REMOVE_SELF_EDGES);
  flags.set(GENERATE_CONSECUTIVE_IDS);
#ifndef HAVE_ZOLTAN2_MPI
  flags.set(BUILD_LOCAL_GRAPH);
#endif
  this->env->debug(DETAILED_STATUS, "    building graph model");
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env, 
                                            this->problemComm, flags));
  this->env->debug(DETAILED_STATUS, "    graph model built");
}

/* 
NOTE:
  Assumes installed PuLP library is version pulp-0.2
  Assumes installed XtraPuLP library is version xtrapulp-0.3
*/
template <typename Adapter>
void AlgPuLP<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  int num_parts = (int)numGlobalParts;
  //TPL_Traits<int, size_t>::ASSIGN(num_parts, numGlobalParts, env);

  //#ifdef HAVE_ZOLTAN2_MPI
  // TODO: XtraPuLP

  int ierr = 0;
  int np = problemComm->getSize();

  // Get number of vertices and edges
  const size_t modelVerts = model->getLocalNumVertices();
  const size_t modelEdges = model->getLocalNumEdges();
  int num_verts = (int)modelVerts;
  long num_edges = (long)modelEdges;
  //TPL_Traits<int, size_t>::ASSIGN(num_verts, modelVerts, env);
  //TPL_Traits<long, size_t>::ASSIGN(num_edges, modelEdges, env);

  // Get vertex info
  ArrayView<const gno_t> vtxIDs;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  size_t nVtx = model->getVertexList(vtxIDs, vwgts);
  int nVwgts = model->getNumWeightsPerVertex();


#ifndef HAVE_ZOLTAN2_MPI
  // PuLP current only supports a single vertex weight
  if (nVwgts > 1) {
    std::cerr << "Warning:  NumWeightsPerVertex is " << nVwgts 
              << " but PuLP allows only one weight. "
              << " Zoltan2 will use only the first weight per vertex."
              << std::endl;
  }

  std::unique_ptr<int[]> vertex_weights;
  long vertex_weights_sum = 0;
  if (nVwgts) {
    nVwgts = 1;
    vertex_weights = std::unique_ptr<int[]>(scale_weights(nVtx, nVwgts, vwgts));

    for (int i = 0; i < num_verts; ++i)
      vertex_weights_sum += (long)vertex_weights[i];
  } else {
    vertex_weights = std::unique_ptr<int[]>(nullptr);
  }
#else
  // XtraPuLP supports an arbitrary number of vertex weights
  std::unique_ptr<int[]> vertex_weights;
  if (nVwgts) {
    vertex_weights = std::unique_ptr<int[]>(scale_weights(nVtx, nVwgts, vwgts));
  } else {
    vertex_weights = std::unique_ptr<int[]>(nullptr);
  }
#endif

  // Get edge info
  ArrayView<const gno_t> adjs;
  ArrayView<const offset_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;
  size_t nEdge = model->getEdgeList(adjs, offsets, ewgts);
  int nEwgts = model->getNumWeightsPerEdge();
  if (nEwgts > 1) {
    std::cerr << "Warning:  NumWeightsPerEdge is " << nEwgts 
              << " but PuLP/XtraPuLP allows only one weight. "
              << " Zoltan2 will use only the first weight per edge."
              << std::endl;
  }

  std::unique_ptr<int[]> edge_weights;
  if (nEwgts) {
    nEwgts = 1;
    edge_weights = std::unique_ptr<int[]>(scale_weights(nEdge, nEwgts, ewgts));
    if (!nVwgts) {
      // For XtraPulp, we need to fill vertex_weights array if we have
      // any edge weights.
      vertex_weights = std::unique_ptr<int[]>(new int[nVtx]);
      nVwgts = 1;
      for (size_t i = 0; i < nVtx; ++i) {
        vertex_weights[i] = 1;
      }
    }
  } else if (nVwgts) {
    // For XtraPulp, we need to fill edge_weights array if we have
    // any vertex weights.
    edge_weights = std::unique_ptr<int[]>(new int[nEdge]);
    for (size_t i = 0; i < nEdge; ++i) {
      edge_weights[i] = 1;
    }
  } else {
    edge_weights = std::unique_ptr<int[]>(nullptr);
  }

#ifndef HAVE_ZOLTAN2_MPI
  // Create PuLP's graph structure
  int* out_edges = nullptr;
  long* out_offsets = nullptr;
  TPL_Traits<int, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
  TPL_Traits<long, const offset_t>::ASSIGN_ARRAY(&out_offsets, offsets);

  pulp_graph_t g = {num_verts, num_edges, 
                    out_edges, out_offsets,
                    vertex_weights.get(), edge_weights.get(), 
                    vertex_weights_sum};

#else
  // Create XtraPuLP's graph structure
  unsigned long* out_edges = nullptr;
  unsigned long* out_offsets = nullptr;
  TPL_Traits<unsigned long, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
  TPL_Traits<unsigned long, const offset_t>::ASSIGN_ARRAY(&out_offsets, offsets);

  const size_t modelVertsGlobal = model->getGlobalNumVertices();
  const size_t modelEdgesGlobal = model->getGlobalNumEdges();
  unsigned long num_verts_global = (unsigned long)modelVertsGlobal;
  unsigned long num_edges_global = (unsigned long)modelEdgesGlobal;

  unsigned long* global_ids = nullptr;
  TPL_Traits<unsigned long, const gno_t>::ASSIGN_ARRAY(&global_ids, vtxIDs);

  ArrayView<size_t> vtxDist;
  model->getVertexDist(vtxDist);  
  unsigned long* verts_per_rank = new unsigned long[np+1];
  for (int i = 0; i < np+1; ++i)
    verts_per_rank[i] = vtxDist[i];

  dist_graph_t g;
  create_xtrapulp_dist_graph(&g, num_verts_global, num_edges_global, 
                          (unsigned long)num_verts, (unsigned long)num_edges, 
                          out_edges, out_offsets, global_ids, verts_per_rank,
                          nVwgts, vertex_weights.get(), edge_weights.get());
#endif


  // Create array for PuLP to return results in.
  // Or write directly into solution parts array
  ArrayRCP<part_t> partList(new part_t[num_verts], 0, num_verts, true);
  int* parts = nullptr;
  if (num_verts && (sizeof(int) == sizeof(part_t))) {
    // Can write directly into the solution's memory
    parts = (int *)partList.getRawPtr();
  }
  else {
    // Can't use solution memory directly; will have to copy later.
    parts = new int[num_verts];
  }

  // TODO
  // Implement target part sizes

  // Grab options from parameter list
  const Teuchos::ParameterList &pl = env->getParameters();
  const Teuchos::ParameterEntry *pe;

  // figure out which parts of the algorithm we're going to run
  // Default to PuLP with BFS init
  // PuLP - do_edge_min = false, do_maxcut_min = false
  // PuLP-M - do_edge_bal = true, do_maxcut_min = false
  // PuLP-MM - do_edge_bal = true/false, do_maxcut_min = true
  // PuLP-MC - do_edge_bal = false, do_maxcut_min = true/false
  bool do_lp_init = false;
  bool do_bfs_init = true;
  bool do_edge_bal = false;
  bool do_repart = false;
  bool do_maxcut_min = false;
  bool verbose_output = false;

  // Do label propagation initialization instead of bfs?
  pe = pl.getEntryPtr("pulp_lp_init");
  if (pe) do_lp_init = pe->getValue(&do_lp_init);
  if (do_lp_init) do_bfs_init = false;

  // Now look at additional objective
  pe = pl.getEntryPtr("pulp_minimize_maxcut");
  if (pe) {
    do_maxcut_min = pe->getValue(&do_maxcut_min);
    // If we're doing the secondary objective, 
    //   set the additional constraint as well
    if (do_maxcut_min) do_edge_bal = true;
  }

  pe = pl.getEntryPtr("pulp_do_repart");
  if (pe) {
    do_repart = pe->getValue(&do_repart);
    if (do_repart) {
      // Do repartitioning with input parts
      // TODO: read in current parts
      // for (int i = 0; i < num_verts; ++i)
      //   parts[i] = something;
      do_bfs_init = false;
      do_lp_init = false;
      std::string errStr = "repartitioning within (Xtra)PuLP is ";
      errStr            += "currently unsupported.";
      throw std::runtime_error(errStr);
    }
  }

  // Now grab vertex and edge imbalances, defaults at 10%
  double vert_imbalance = 1.1;
  double edge_imbalance = 1.1;
  double imbalance = 1.1;

  pe = pl.getEntryPtr("pulp_vert_imbalance");
  if (pe) vert_imbalance = pe->getValue<double>(&vert_imbalance);
  pe = pl.getEntryPtr("pulp_edge_imbalance");
  if (pe) {
    edge_imbalance = pe->getValue<double>(&edge_imbalance);
    // if manually set edge imbalance, add do_edge_bal flag to true
    do_edge_bal = true;
  }
  pe = pl.getEntryPtr("pulp_imbalance");
  if (pe) imbalance = pe->getValue<double>(&imbalance);

  if (vert_imbalance < 1.0)
    throw std::runtime_error("pulp_vert_imbalance must be '1.0' or greater.");
  if (edge_imbalance < 1.0)
    throw std::runtime_error("pulp_edge_imbalance must be '1.0' or greater.");
  if (imbalance < 1.0)
    throw std::runtime_error("pulp_imbalance must be '1.0' or greater.");

  // verbose output?  
  // TODO: fully implement verbose flag throughout PuLP
  pe = pl.getEntryPtr("pulp_verbose");
  if (pe) verbose_output = pe->getValue(&verbose_output);

  // using pulp seed? 
  int pulp_seed = rand();
  pe = pl.getEntryPtr("pulp_seed");
  if (pe) pulp_seed = pe->getValue(&pulp_seed);


  if (verbose_output) {
    printf("procid: %d, n: %i, m: %li, vb: %f, eb: %f, imb: %f, p: %i\n",
      problemComm->getRank(), 
      num_verts, num_edges, 
      vert_imbalance, edge_imbalance, imbalance, num_parts);
  }

  // Call partitioning; result returned in parts array
#ifndef HAVE_ZOLTAN2_MPI
  // Create PuLP's partitioning data structure
  pulp_part_control_t ppc = {vert_imbalance, edge_imbalance,
    do_lp_init, do_bfs_init, do_repart,
    do_edge_bal, do_maxcut_min,
    verbose_output, pulp_seed};

  ierr = pulp_run(&g, &ppc, parts, num_parts);  

  env->globalInputAssertion(__FILE__, __LINE__, "pulp_run", 
    !ierr, BASIC_ASSERTION, problemComm);
#else
  // Create XtraPuLP's partitioning data structure
  std::unique_ptr<double[]> constraints;
  if (nVwgts > 0) {
    constraints = std::unique_ptr<double[]>(new double[nVwgts]);
    for (int i = 0; i < nVwgts; ++i) {
      constraints[i] = imbalance;
    }
  } else { 
    constraints = std::unique_ptr<double[]>(nullptr);
  }


  pulp_part_control_t ppc = {
      vert_imbalance, edge_imbalance, 
      constraints.get(), nVwgts, 
      do_lp_init, do_bfs_init, do_repart, 
      do_edge_bal, do_maxcut_min,
      verbose_output, pulp_seed};

  ierr = xtrapulp_run(&g, &ppc, parts, num_parts);

  env->globalInputAssertion(__FILE__, __LINE__, "xtrapulp_run", 
    !ierr, BASIC_ASSERTION, problemComm);
#endif


  // Load answer into the solution if necessary
  if ( (sizeof(int) != sizeof(part_t)) || (num_verts == 0) ) {
    for (int i = 0; i < num_verts; i++) 
      partList[i] = parts[i];
    delete [] parts;
  }

  solution->setParts(partList);

  env->memory("Zoltan2-(Xtra)PuLP: After creating solution");

  // Clean up copies made due to differing data sizes.
#ifndef HAVE_ZOLTAN2_MPI
  TPL_Traits<int, const gno_t>::DELETE_ARRAY(&out_edges);
  TPL_Traits<long, const lno_t>::DELETE_ARRAY(&out_offsets);
#else
  TPL_Traits<unsigned long, const gno_t>::DELETE_ARRAY(&out_edges);
  TPL_Traits<unsigned long, const lno_t>::DELETE_ARRAY(&out_offsets);
  TPL_Traits<unsigned long, const gno_t>::DELETE_ARRAY(&global_ids);
#endif


//#endif // DO NOT HAVE_MPI
}

/////////////////////////////////////////////////////////////////////////////
// Scale and round scalar_t weights (typically float or double) to 
// PuLP int
// subject to sum(weights) <= max_wgt_sum.
// Scale only if deemed necessary.
//
// Note that we use ceil() instead of round() to avoid
// rounding to zero weights.
// Based on Zoltan's scale_round_weights, mode 1
//
// Modified from scale_weights() in Zoltan2_AlgParMETIS.hpp
template <typename Adapter>
int* AlgPuLP<Adapter>::scale_weights(
  size_t nVtx, int nWgt,
  ArrayView<StridedData<typename Adapter::lno_t, 
                        typename Adapter::scalar_t> > fwgts
)
{
  const double INT_EPSILON = 1e-5;

  int *iwgts = new int[nVtx*nWgt];
  int *nonint_local = new int[nWgt+nWgt]; 
  int *nonint = nonint_local + nWgt;

  double *sum_wgt_local = new double[nWgt*4];
  double *max_wgt_local = sum_wgt_local + nWgt;
  double *sum_wgt = max_wgt_local + nWgt;
  double *max_wgt = sum_wgt + nWgt;

  for (int i = 0; i < nWgt; i++) {
    nonint_local[i] = 0;
    sum_wgt_local[i] = 0.0;
    max_wgt_local[i] = 0.0; 
  }

  // Compute local sums of the weights 
  // Check whether all weights are integers
  for (int j = 0; j < nWgt; j++) {
    for (size_t i = 0; i < nVtx; i++) {
      double fw = double(fwgts[j][i]);
      if (!nonint_local[j]) {
        int tmp = (int) floor(fw + .5); /* Nearest int */
        if (fabs((double)tmp-fw) > INT_EPSILON) {
          nonint_local[j] = 1;
        }
      }
      sum_wgt_local[j] += fw;
      if (fw > max_wgt_local[j]) max_wgt_local[j] = fw;
    }
  }

  Teuchos::reduceAll<int,int> (*problemComm, Teuchos::REDUCE_MAX, nWgt,
                              nonint_local,  nonint);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_SUM, nWgt,
                                 sum_wgt_local, sum_wgt);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_MAX, nWgt,
                                 max_wgt_local, max_wgt);

  const double max_wgt_sum = double(std::numeric_limits<int>::max()/8);
  for (int j = 0; j < nWgt; j++) {
    double scale = 1.;

    // Scaling needed if weights are not integers or weights' 
    // range is not sufficient
    if (nonint[j] || (max_wgt[j]<=INT_EPSILON) || (sum_wgt[j]>max_wgt_sum)) {
      /* Calculate scale factor */
      if (sum_wgt[j] != 0.) scale = max_wgt_sum/sum_wgt[j];
    }

    /* Convert weights to positive integers using the computed scale factor */
    for (size_t i = 0; i < nVtx; i++)
      iwgts[i*nWgt+j] = (int) ceil(double(fwgts[j][i])*scale);
  }

  delete [] nonint_local;
  delete [] sum_wgt_local;

  return iwgts;
}


} // namespace Zoltan2

#endif // HAVE_ZOLTAN2_PULP

////////////////////////////////////////////////////////////////////////


#endif

