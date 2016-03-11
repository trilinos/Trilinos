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
  AlgPuLP(const RCP<const Environment> &env,
          const RCP<const Comm<int> > &problemComm,
          const RCP<const base_adapter_t> &adapter
  )
  {
    throw std::runtime_error(
          "BUILD ERROR:  PuLP requested but not compiled into Zoltan2.\n"
          "Please set CMake flag Zoltan2_ENABLE_PuLP:BOOL=ON.");
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
//#ifndef HAVE_ZOLTAN2_MPI
#include "pulp.h"
//#else
//endif
}




template <typename Adapter>
class AlgPuLP : public Algorithm<Adapter>
{
public:
  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
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

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  void buildModel(modelFlag_t &flags);

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

  int sgParameter = 0;
  pe = pl.getEntryPtr("subset_graph");
  if (pe)
    sgParameter = pe->getValue<int>(&sgParameter);
  if (sgParameter == 1)
      flags.set(BUILD_SUBSET_GRAPH);

  flags.set(REMOVE_SELF_EDGES);
  flags.set(GENERATE_CONSECUTIVE_IDS);
  flags.set(BUILD_LOCAL_GRAPH);
  this->env->debug(DETAILED_STATUS, "    building graph model");
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env, 
                                            this->problemComm, flags));
  this->env->debug(DETAILED_STATUS, "    graph model built");
}


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
  //int np = problemComm->getSize();

  // Get number of vertices and edges
  const size_t modelVerts = model->getLocalNumVertices();
  const size_t modelEdges = model->getLocalNumEdges();
  int num_verts = (int)modelVerts;
  long num_edges = (long)modelEdges;
  //TPL_Traits<int, size_t>::ASSIGN(num_verts, modelVerts, env);
  //TPL_Traits<long, size_t>::ASSIGN(num_edges, modelEdges, env);

  // Get adjacency list and offset pointers 
  ArrayView<const gno_t> adjs;
  ArrayView<const lno_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts; // ignore this for now
  model->getEdgeList(adjs, offsets, ewgts);

  // Create PuLP's graph structure
  int *out_edges;
  long *out_offsets;
  TPL_Traits<int, const gno_t>::ASSIGN_ARRAY(&out_edges, adjs);
  TPL_Traits<long, const lno_t>::ASSIGN_ARRAY(&out_offsets, offsets);

  pulp_graph_t g = {num_verts, num_edges, out_edges, out_offsets};

  // Create array for PuLP to return results in.
  // Or write directly into solution parts array
  ArrayRCP<part_t> partList(new part_t[num_verts], 0, num_verts, true);
  int* parts = NULL;
  if (num_verts && (sizeof(int) == sizeof(part_t))) {
    // Can write directly into the solution's memory
    parts = (int *) partList.getRawPtr();
  }
  else {
    // Can't use solution memory directly; will have to copy later.
    parts = new int[num_verts];
  }

  // TODO
  // Implement vertex and edge weights

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
  int do_lp_init = 0;
  int do_bfs_init = 1;
  int do_edge_bal = 0;
  int do_maxcut_min = 0;
  int verbose_output = 0;

  // Do label propagation initialization instead of bfs?
  pe = pl.getEntryPtr("pulp_lp_init");
  if (pe) do_lp_init = pe->getValue<int>(&do_lp_init);
  if (do_lp_init) do_bfs_init = 0;

  // Now look at additional objective
  pe = pl.getEntryPtr("pulp_minimize_maxcut");
  if (pe) {
    do_maxcut_min = pe->getValue<int>(&do_maxcut_min);
    // If we're doing the secondary objective, 
    //   set the additional constraint as well
    if (do_maxcut_min) do_edge_bal = 1;
  }

  // Now grab vertex and edge imbalances, defaults at 10%
  double vert_imbalance = 1.1;
  double edge_imbalance = 1.1;

  pe = pl.getEntryPtr("pulp_vert_imbalance");
  if (pe) vert_imbalance = pe->getValue<double>(&vert_imbalance);
  pe = pl.getEntryPtr("pulp_edge_imbalance");
  if (pe) {
    edge_imbalance = pe->getValue<double>(&edge_imbalance);
    // if manually set edge imbalance, add do_edge_bal flag to true
    do_edge_bal = 1;
  }

  if (vert_imbalance < 1.0)
    throw std::runtime_error("pulp_vert_imbalance must be '1.0' or greater.");
  if (edge_imbalance < 1.0)
    throw std::runtime_error("pulp_edge_imbalance must be '1.0' or greater.");

  // verbose output?  
  // TODO: fully implement verbose flag throughout PuLP
  pe = pl.getEntryPtr("pulp_verbose");
  if (pe) verbose_output = pe->getValue<int>(&verbose_output);

  // Create PuLP's partitioning data structure
  pulp_part_control_t ppc = {vert_imbalance, edge_imbalance,
    (bool)do_lp_init, (bool)do_bfs_init, 
    (bool)do_edge_bal, (bool)do_maxcut_min,
    (bool)verbose_output};

  if (verbose_output) {
    printf("n: %i, m: %li, vb: %lf, eb: %lf, p: %i\n",
      num_verts, num_edges, vert_imbalance, edge_imbalance, num_parts);
  }

  // Call partitioning; result returned in parts array
  ierr = pulp_run(&g, &ppc, parts, num_parts);

  env->globalInputAssertion(__FILE__, __LINE__, "pulp_run", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Load answer into the solution if necessary
  if ((sizeof(int) != sizeof(part_t)) || (num_verts == 0)) {
    for (int i = 0; i < num_verts; i++) partList[i] = parts[i];
    delete [] parts;
  }

  solution->setParts(partList);

  env->memory("Zoltan2-PuLP: After creating solution");

  // Clean up copies made due to differing data sizes.
  TPL_Traits<int, const lno_t>::DELETE_ARRAY(&out_edges);
  TPL_Traits<long, const gno_t>::DELETE_ARRAY(&out_offsets);

//#endif // DO NOT HAVE_MPI
}

} // namespace Zoltan2

#endif // HAVE_ZOLTAN2_PULP

////////////////////////////////////////////////////////////////////////


#endif






