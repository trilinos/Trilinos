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
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>

typedef zoltan2_partId_t partId_t;

#ifndef HAVE_ZOLTAN2_SCOTCH

// Error handling for when Scotch is requested
// but Zoltan2 not built with Scotch.

namespace Zoltan2 {

/*! Scotch partitioning method.
 *
 *  \param env  parameters for the problem and library configuration
 *  \param problemComm  the communicator for the problem
 *  \param model a graph
 *  \param solution  a Solution object
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).
 */

template <typename Adapter>
void AlgPTScotch(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpicomm,
#endif
  const RCP<GraphModel<typename Adapter::base_adapter_t> > &model,
  RCP<PartitioningSolution<Adapter> > &solution

)
{
  throw std::runtime_error(
        "BUILD ERROR:  Scotch requested but not compiled into Zoltan2.\n"
        "Please set CMake flag Zoltan2_ENABLE_Scotch:BOOL=ON.");
}
}

#else  //HAVE_ZOLTAN2_SCOTCH


// stdint.h for int64_t in scotch header

#include <stdint.h>
#ifndef HAVE_ZOLTAN2_MPI
#include "scotch.h"
#else
#include "ptscotch.h"
#endif

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

#endif

#endif


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Scotch.hpp
//! \brief Parallel graph partitioning using Scotch.

namespace Zoltan2{

/////////////////////////////////////////////////////////////////////////////
//  Traits struct to handle conversions between gno_t/lno_t and SCOTCH_Num.
/////////////////////////////////////////////////////////////////////////////

// General case:  SCOTCH_Num and gno_t/lno_t (called zno_t here) differ.
template <typename zno_t>
struct SCOTCH_Num_Traits {

  static inline SCOTCH_Num ASSIGN_TO_SCOTCH_NUM(
    SCOTCH_Num &a,
    zno_t b,
    const RCP<const Environment> &env)
  {
    // Assign a = b; make sure SCOTCH_Num is large enough to accept zno_t.
    if (b <= SCOTCH_NUMMAX) a = b;
    else 
      env->localInputAssertion(__FILE__, __LINE__, 
       "Value too large for SCOTCH_Num, Rebuild Scotch with larger SCOTCH_Num",
       false, BASIC_ASSERTION);
    return a;
  }

  static inline void ASSIGN_SCOTCH_NUM_ARRAY(
    SCOTCH_Num **a,
    ArrayView<const zno_t> &b,
    const RCP<const Environment> &env)
  {
    // Allocate array a; copy b values into a.
    size_t size = b.size();
    *a = new SCOTCH_Num[size];
    for (size_t i = 0; i < size; i++) ASSIGN_TO_SCOTCH_NUM((*a)[i], b[i], env);
  }

  static inline void DELETE_SCOTCH_NUM_ARRAY(SCOTCH_Num *a)
  {
    // Delete the copy made in ASSIGN_SCOTCH_NUM_ARRAY.
    delete [] a;
  }

  //TODO static inline ASSIGN_FROM_SCOTCH_NUM(zno_t b, SCOTCH_Num a)
  //TODO {
  //TODO   Make sure zno_t is large enough to accept SCOTCH_Num.
  //TODO   NOT NEEDED AS LONG AS PARTS ARE size_t.
  //TODO }
};


// Special case:  zno_t == SCOTCH_Num. No error checking or copies needed.
template <>
struct SCOTCH_Num_Traits<SCOTCH_Num> {
  static inline SCOTCH_Num ASSIGN_TO_SCOTCH_NUM(
    SCOTCH_Num &a,
    SCOTCH_Num b,
    const RCP<const Environment> &env)
  {
    a = b;
    return a;
  }
  static inline void ASSIGN_SCOTCH_NUM_ARRAY(
    SCOTCH_Num **a,
    ArrayView<const SCOTCH_Num> &b,
    const RCP<const Environment> &env)
  {
    *a = const_cast<SCOTCH_Num *> (b.getRawPtr());
  }
  static inline void DELETE_SCOTCH_NUM_ARRAY(SCOTCH_Num *a) { }
};


///////////////////////////////////////////////////////////////////////
// Now, the actual Scotch algorithm.
///////////////////////////////////////////////////////////////////////

/*! Scotch partitioning method.
 *
 *  \param env  parameters for the problem and library configuration
 *  \param problemComm  the communicator for the problem
 *  \param model a graph
 *  \param solution  a Solution object
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).
 */

template <typename Adapter>
void AlgPTScotch(
  const RCP<const Environment> &env,        // parameters & app comm
  const RCP<const Comm<int> > &problemComm, // problem comm
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpicomm,
#endif
  const RCP<GraphModel<typename Adapter::base_adapter_t> > &model, // the graph
  RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr = 0;

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  SCOTCH_Num partnbr;
  SCOTCH_Num_Traits<size_t>::ASSIGN_TO_SCOTCH_NUM(partnbr, numGlobalParts, env);

#ifdef HAVE_ZOLTAN2_MPI

  const SCOTCH_Num  baseval = 0;  // Base value for array indexing.
                                  // GraphModel returns GNOs from base 0.

  SCOTCH_Strat stratstr;          // Strategy string
                                  // TODO:  Set from parameters
  SCOTCH_stratInit(&stratstr);

  // Allocate & initialize PTScotch data structure.
  SCOTCH_Dgraph *gr = SCOTCH_dgraphAlloc();  // Scotch distributed graph
  ierr = SCOTCH_dgraphInit(gr, mpicomm);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphInit", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Get vertex info
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > xyz;
  ArrayView<StridedData<lno_t, scalar_t> > vtxWt;
  size_t nVtx = model->getVertexList(vtxID, xyz, vtxWt);
  SCOTCH_Num vertlocnbr;
  SCOTCH_Num_Traits<size_t>::ASSIGN_TO_SCOTCH_NUM(vertlocnbr, nVtx, env);
  SCOTCH_Num vertlocmax = vertlocnbr; // Assumes no holes in global nums.

  // Get edge info
  ArrayView<const gno_t> edgeIds;
  ArrayView<const int>   procIds;
  ArrayView<const lno_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;

  size_t nEdges = model->getEdgeList(edgeIds, procIds, offsets, ewgts);

  SCOTCH_Num edgelocnbr;
  SCOTCH_Num_Traits<size_t>::ASSIGN_TO_SCOTCH_NUM(edgelocnbr, nEdges, env);
  const SCOTCH_Num edgelocsize = edgelocnbr;  // Assumes adj array is compact.

  SCOTCH_Num *vertloctab;  // starting adj/vtx
  SCOTCH_Num_Traits<lno_t>::ASSIGN_SCOTCH_NUM_ARRAY(&vertloctab, offsets, env);

  SCOTCH_Num *edgeloctab;  // adjacencies
  SCOTCH_Num_Traits<gno_t>::ASSIGN_SCOTCH_NUM_ARRAY(&edgeloctab, edgeIds, env);

  // We don't use these arrays, but we need them as arguments to Scotch.
  SCOTCH_Num *vendloctab = NULL;  // Assume consecutive storage for adj
  SCOTCH_Num *vlblloctab = NULL;  // Vertex label array
  SCOTCH_Num *edgegsttab = NULL;  // Array for ghost vertices

  // Get weight info.
  // TODO:  Actually get the weights; for now, not using weights.
  SCOTCH_Num *veloloctab = NULL;  // Vertex weights
  SCOTCH_Num *edloloctab = NULL;  // Edge weights
  //TODO int vwtdim = model->getVertexWeightDim();
  //TODO int ewtdim = model->getEdgeWeightDim();
  //TODO if (vwtdim) veloloctab = new SCOTCH_Num[nVtx];
  //TODO if (ewtdim) edloloctab = new SCOTCH_Num[nEdges];
  //TODO scale weights to SCOTCH_Nums.

  // Build PTScotch distributed data structure
  ierr = SCOTCH_dgraphBuild(gr, baseval, vertlocnbr, vertlocmax,
                            vertloctab, vendloctab, veloloctab, vlblloctab,
                            edgelocnbr, edgelocsize,
                            edgeloctab, edgegsttab, edloloctab);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphBuild", 
    !ierr, BASIC_ASSERTION, problemComm);

  // Create array for Scotch to return results in.
  ArrayRCP<partId_t> partList(new partId_t [nVtx], 0, nVtx,true);
  SCOTCH_Num *partloctab;
  if (sizeof(SCOTCH_Num) == sizeof(partId_t)) {
    // Can write directly into the solution's memory
    partloctab = (SCOTCH_Num *) partList.getRawPtr();
  }
  else {
    // Can't use solution memory directly; will have to copy later.
    partloctab = new SCOTCH_Num[nVtx];
  }

  // Call partitioning; result returned in partloctab.
  // TODO:  Use SCOTCH_dgraphMap so can include a machine model in partitioning
  ierr = SCOTCH_dgraphPart(gr, partnbr, &stratstr, partloctab);

  env->globalInputAssertion(__FILE__, __LINE__, "SCOTCH_dgraphPart", 
    !ierr, BASIC_ASSERTION, problemComm);

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
  SCOTCH_stratExit(&stratstr);

  // Load answer into the solution.

  if (sizeof(SCOTCH_Num) != sizeof(partId_t)) {
    for (size_t i = 0; i < nVtx; i++) partList[i] = partloctab[i];
  }

  ArrayRCP<const gno_t> gnos = arcpFromArrayView(vtxID);

  solution->setParts(gnos, partList, true);

  env->memory("Zoltan2-Scotch: After creating solution");

  //if (me == 0) cout << " done." << endl;
  // Clean up Zoltan2
  //TODO if (vwtdim) delete [] velotab;
  //TODO if (ewtdim) delete [] edlotab;

  // Clean up copies made due to differing data sizes.
  if (sizeof(lno_t) != sizeof(SCOTCH_Num)) delete [] vertloctab;
  if (sizeof(gno_t) != sizeof(SCOTCH_Num)) delete [] edgeloctab;

#else // DO NOT HAVE_MPI

  // TODO:  Handle serial case with calls to Scotch.
  // TODO:  For now, assign everything to rank 0 and assume only one part.
  // TODO:  Can probably use the code above for loading solution,
  // TODO:  instead of duplicating it here.
  // TODO
  // TODO:  Actual logic should call Scotch when number of processes == 1.
  ArrayView<const gno_t> vtxID;
  ArrayView<StridedData<lno_t, scalar_t> > xyz;
  ArrayView<StridedData<lno_t, scalar_t> > vtxWt;
  size_t nVtx = model->getVertexList(vtxID, xyz, vtxWt);

  for (size_t i = 0; i < nVtx; i++) partList[i] = 0;


#endif // DO NOT HAVE_MPI
}

}
#endif // HAVE_ZOLTAN2_SCOTCH
#endif
