#ifndef _ZOLTAN2_ALGSCOTCH_HPP_
#define _ZOLTAN2_ALGSCOTCH_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#ifdef HAVE_SCOTCH
#ifndef HAVE_MPI
#include "scotch.h"
#else
#include "ptscotch.h"
#endif
#endif

#ifdef SHOW_LINUX_MEMINFO
extern "C"{
static char *z2_meminfo=NULL;
extern void Zoltan_get_linux_meminfo(char *msg, char **result);
}
#endif

#ifdef SHOW_SCOTCH_HIGH_WATER_MARK
extern "C"{
//
// Scotch keeps track of memory high water mark, but doesn't
// provide a way to get that number.  So add this function:  
//   "size_t SCOTCH_getMemoryMax() { return memorymax;}"
// to src/libscotch/common_memory.c
// and compile scotch with -DCOMMON_MEMORY_TRACE
//
extern size_t SCOTCH_getMemoryMax();
}
#endif


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Scotch.hpp
//! \brief Parallel graph partitioning using Scotch.


// Placeholder for real error handling.
#define KDD_HANDLE_ERROR {\
    cout << __func__ << ":" << __LINE__ << " KDDERROR" << endl;\
    }

namespace Zoltan2{

template <typename Adapter>
void AlgPTScotch(
  const size_t nParts,
  const RCP<GraphModel<Adapter> > &model, 
  RCP<PartitioningSolution<typename Adapter::gid_t, 
                           typename Adapter::lid_t,
                           typename Adapter::lno_t> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<const Teuchos::Comm<int> > &comm,
  const RCP<const Environment> &env
) 
{
#ifndef HAVE_SCOTCH
  throw std::runtime_error(
        "BUILD ERROR:  Scotch requested but not compiled into Zoltan2.\n" 
        "Please set CMake flag Zoltan2_ENABLE_Scotch:BOOL=ON.");

#else

  HELLO;

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr = 0;
  const SCOTCH_Num partnbr = nParts;

#ifdef HAVE_MPI

  if (sizeof(lno_t) != sizeof(SCOTCH_Num) || 
      sizeof(gno_t) != sizeof(SCOTCH_Num)) {
    throw std::runtime_error(
          "Incompatible Scotch build; "
          "sizeof SCOTCH_Num must match sizeof(lno_t) and sizeof(gno_t)\n");
    return;
  }

  const SCOTCH_Num  baseval = 0;  // Base value for array indexing.
                                  // GraphModel returns GNOs from base 0.

  SCOTCH_Strat stratstr;          // Strategy string
                                  // TODO:  Set from parameters
  SCOTCH_stratInit(&stratstr);

  MPI_Comm mpicomm = *(rcp_dynamic_cast<const Teuchos::MpiComm<int> > 
                                       (comm)->getRawMpiComm());

  // Allocate & initialize PTScotch data structure.
  SCOTCH_Dgraph *gr = SCOTCH_dgraphAlloc();  // Scotch distributed graph 
  ierr = SCOTCH_dgraphInit(gr, mpicomm);  
  if (ierr) {
    KDD_HANDLE_ERROR;
  }

  // Get vertex info
  ArrayView<const gno_t> vtxID;
  ArrayView<const scalar_t> xyz;
  ArrayView<const scalar_t> vtxWt;
  size_t nVtx = model->getVertexList(vtxID, xyz, vtxWt);
  const SCOTCH_Num vertlocnbr = nVtx;
  const SCOTCH_Num vertlocmax = vertlocnbr; // Assumes no holes in global nums.

  // Get edge info
  ArrayView<const gno_t> edgeIds;
  ArrayView<const int>   procIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> ewgts;
  size_t nEdges = model->getEdgeList(edgeIds, procIds, offsets, ewgts);
  const SCOTCH_Num edgelocnbr = nEdges;
  const SCOTCH_Num edgelocsize = edgelocnbr;  // Assumes adj array is compact.


  SCOTCH_Num *vertloctab = 
              const_cast<SCOTCH_Num *>(offsets.getRawPtr()); // starting adj/vtx

  SCOTCH_Num *edgeloctab = 
              const_cast<SCOTCH_Num *>(edgeIds.getRawPtr()); // adjacencies

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
  if (ierr) {
    KDD_HANDLE_ERROR;
  }

  // Create array for Scotch to return results in.
  SCOTCH_Num *partloctab;
  if (sizeof(SCOTCH_Num) == sizeof(size_t)) {
    // Can write directly into the solution's memory
    partloctab = (SCOTCH_Num *) (solution->getPartsRCP().getRawPtr());
  }
  else {
    // Can't use solution memory directly; will have to copy later.
    partloctab = new SCOTCH_Num[nVtx];
  }

  // Call partitioning; result returned in partloctab.
  // TODO:  Use SCOTCH_dgraphMap so can include a machine model in partitioning
  ierr = SCOTCH_dgraphPart(gr, partnbr, &stratstr, partloctab);
  if (ierr) KDD_HANDLE_ERROR;

#ifdef SHOW_SCOTCH_HIGH_WATER_MARK
  int me = comm->getRank();
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

  if (sizeof(SCOTCH_Num) != sizeof(size_t)) {
    // Need to copy parts if didn't use parts array directly from Solution.
    ArrayRCP<size_t> solnParts = solution->getPartsRCP();
    for (size_t i = 0; i < nVtx; i++) solnParts[i] = partloctab[i];
    delete [] partloctab;
  }

  if (!(model->getIdentifierMap().is_null())) {
    // TODO Need to translate vtxIDs to GIDs for solution
    Z2_LOCAL_BUG_ASSERTION(*env, "GID translation not yet implemented", 1, 0);
  }
  else {
    // Copy vtxIDs to Solution, assuming vtxIDs == GIDs.
    // TODO:  Handling of LIDs is not yet clear.
    ArrayRCP<gid_t> solnGIDs = solution->getGidsRCP();
    ArrayRCP<lid_t> solnLIDs = solution->getLidsRCP();
    for (size_t i = 0; i < nVtx; i++) solnGIDs[i] = vtxID[i];
    for (size_t i = 0; i < solnLIDs.size(); i++) solnLIDs[i] = i; // TODO NOT CORRECT YET. KDD
  }

#ifdef SHOW_LINUX_MEMINFO
  if (me==0){
    Zoltan_get_linux_meminfo("After creating solution", &z2_meminfo);
    if (z2_meminfo){
      std::cout << "Rank " << me << ": " << z2_meminfo << std::endl;
      free(z2_meminfo);
      z2_meminfo=NULL;
    }
  }
#endif

  //if (me == 0) cout << " done." << endl;
  // Clean up Zoltan2
  //TODO if (vwtdim) delete [] velotab;
  //TODO if (ewtdim) delete [] edlotab;

#else // DO NOT HAVE_MPI

  // TODO:  Handle serial case with calls to Scotch.
  // TODO:  For now, assign everything to rank 0 and assume only one part.
  // TODO:  Can probably use the code above for loading solution,
  // TODO:  instead of duplicating it here.
  // TODO
  // TODO:  Actual logic should call Scotch when number of processes == 1.
  ArrayView<const gno_t> vtxID;
  ArrayView<const scalar_t> xyz;
  ArrayView<const scalar_t> vtxWt;
  size_t nVtx = model->getVertexList(vtxID, xyz, vtxWt);

  ArrayRCP<size_t> solnParts = solution->getParts();
  for (size_t i = 0; i < nVtx; i++) solnParts[i] = 0;

  if (!(model->getIdentifierMap().isNull())) {
    // TODO Need to translate vtxIDs to GIDs for solution
    Z2_LOCAL_BUG_ASSERTION(*env, "GID translation not yet implemented", 1, 0);
  }
  else {
    // Copy vtxIDs to Solution, assuming vtxIDs == GIDs.
    // TODO:  Handling of LIDs is not yet clear.
    RCP<gid_t> solnGIDs = solution->getGidRCP();
    RCP<lid_t> solnLIDs = solution->getLidRCP();
    for (size_t i = 0; i < nVtx; i++) solnGIDs[i] = vtxID[i];
    for (size_t i = 0; i < solnLIDs->size(); i++) solnLIDs[i] = i; // TODO NOT CORRECT YET. KDD
  }

#endif // DO NOT HAVE_MPI
#endif // HAVE_SCOTCH
}

}
#endif
