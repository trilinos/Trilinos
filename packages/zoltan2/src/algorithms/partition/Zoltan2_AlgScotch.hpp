#ifndef _ZOLTAN2_ALGSCOTCH_HPP_
#define _ZOLTAN2_ALGSCOTCH_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#ifdef HAVE_SCOTCH
#include "ptscotch.h"
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
int AlgPTScotch(
  const RCP<GraphModel<Adapter> > &model, 
  const RCP<PartitioningSolution<Adapter> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<const Teuchos::Comm<int> > &comm
) 
{
#ifndef HAVE_SCOTCH
  cout << "Scotch requested but not compiled into Zoltan2." << endl
       << "Please set CMake flag Zoltan2_ENABLE_Scotch:BOOL=ON." << endl;
  return 1;
#else
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::scalar_t scalar_t;

  HELLO;
  int ierr = 0;
  int me = comm->getRank();

  if (sizeof(lno_t) != sizeof(SCOTCH_Num) || 
      sizeof(gno_t) != sizeof(SCOTCH_Num)) {
    cout << "Incompatible Scotch build." << endl;
    KDD_HANDLE_ERROR;
  }

  const SCOTCH_Num  baseval = 0;  // Base value for array indexing.
                                  // GraphModel returns GNOs from base 0.

  SCOTCH_Strat stratstr;          // Strategy string
                                  // TODO:  Set from parameters
  SCOTCH_stratInit(&stratstr);

  MPI_Comm mpicomm = *(rcp_dynamic_cast<const Teuchos::MpiComm<int> > 
                                       (comm)->getRawMpiComm());

  // Allocate & initialize PTScotch data structure.
  if (me == 0) cout << __func__ << ": SCOTCH_dgraphAlloc ";
  SCOTCH_Dgraph *gr = SCOTCH_dgraphAlloc();  // Scotch distributed graph 
  if (me == 0) cout << " done." << endl;

  if (me == 0) cout << __func__ << ": SCOTCH_dgraphInit ";
  ierr = SCOTCH_dgraphInit(gr, mpicomm);  // TODO Handle non-MPI builds.
  if (ierr) KDD_HANDLE_ERROR;
  if (me == 0) cout << " done." << endl;
  
  // Get vertex info
  const SCOTCH_Num vertlocnbr = model->getLocalNumVertices();
  const SCOTCH_Num vertlocmax = vertlocnbr; // Assumes no holes in
                                            // global numbering

  // Get edge info
  const SCOTCH_Num edgelocnbr = model->getLocalNumEdges();
  const SCOTCH_Num edgelocsize = edgelocnbr;// Assumes adj array is compact.

  // Get edges
  ArrayView<const gno_t> edgeIds;
  ArrayView<const int>   procIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> ewgts;
  size_t nEdges = model->getEdgeList(edgeIds, procIds, offsets, ewgts);
  assert(nEdges == edgelocnbr);

  SCOTCH_Num *vertloctab = 
              const_cast<SCOTCH_Num *>(offsets.getRawPtr()); // starting adj/vtx

  SCOTCH_Num *edgeloctab = 
              const_cast<SCOTCH_Num *>(edgeIds.getRawPtr()); // adjacencies

  // We don't need these arrays, I think.
  SCOTCH_Num *vendloctab = NULL;  // Unneeded; assume consecutive 
                                        // storage for adj
  SCOTCH_Num *vlblloctab = NULL;  // Vertex label array; not needed?
  SCOTCH_Num *edgegsttab = NULL;  // Array for ghost vertices; not used


  // Get weight info.
  // TODO:  Actually get the weights; for now, not using weights.
  SCOTCH_Num *veloloctab = NULL;  // Vertex weights
  SCOTCH_Num *edloloctab = NULL;  // Edge weights
  
  //TODO int vwtdim = model->getVertexWeightDim();
  //TODO int ewtdim = model->getEdgeWeightDim();
  //TODO if (vwtdim) veloloctab = new SCOTCH_Num[vertlocnbr];
  //TODO if (ewtdim) edloloctab = new SCOTCH_Num[edgelocnbr];
  //TODO scale weights to SCOTCH_Nums.

  // Build PTScotch distributed data structure
  if (me == 0) cout << __func__ << ": SCOTCH_dgraphBuild ";
  ierr = SCOTCH_dgraphBuild(gr, baseval, vertlocnbr, vertlocmax, 
                            vertloctab, vendloctab, veloloctab, vlblloctab, 
                            edgelocnbr, edgelocsize, 
                            edgeloctab, edgegsttab, edloloctab);
  if (ierr) KDD_HANDLE_ERROR;
  if (me == 0) cout << " done." << endl;

  // Call partitioning; result returned in partloctab.
  // TODO:  Use SCOTCH_dgraphMap so can include a machine model in partitioning
  const SCOTCH_Num partnbr = comm->getSize();  // TODO read from params later.
  SCOTCH_Num *partloctab = new SCOTCH_Num[vertlocnbr];

  if (me == 0) cout << __func__ << ": SCOTCH_dgraphPart ";
  ierr = SCOTCH_dgraphPart(gr, partnbr, &stratstr, partloctab);
  if (ierr) KDD_HANDLE_ERROR;
  if (me == 0) cout << " done." << endl;

  // Clean up PTScotch
  if (me == 0) cout << __func__ << ": SCOTCH_dgraphExit ";
  SCOTCH_dgraphExit(gr);
  SCOTCH_stratExit(&stratstr);
  if (ierr) KDD_HANDLE_ERROR;
  if (me == 0) cout << " done." << endl;

  // Load answer into the solution.
  // TODO May move getVertexList call above when need weights.
  if (me == 0) cout << __func__ << ": Load solution ";
  ArrayView<const gno_t> vtxID;
  ArrayView<const scalar_t> xyz;
  ArrayView<const scalar_t> vtxWt;
  size_t nVtx = model->getVertexList(vtxID, xyz, vtxWt);

  size_t *parts;
  if (sizeof(SCOTCH_Num) == sizeof(size_t))
    parts = (size_t *) partloctab;
  else {  // Need a copy.
    parts = new size_t[vertlocnbr];
    for (size_t i = 0; i < vertlocnbr; i++) parts[i] = partloctab[i];
    delete [] partloctab;
  }

  solution->setPartition((size_t) partnbr, nVtx,
               (gid_t *) (vtxID.getRawPtr()), // TODO Use IdentifierMap
                                              //      instead of cast.
               (lid_t *) NULL,                // TODO Use User's LIDs
               parts);

  if (me == 0) cout << " done." << endl;

  // Clean up Zoltan2
  //TODO if (vwtdim) delete [] velotab;
  //TODO if (ewtdim) delete [] edlotab;

  return ierr;
#endif // HAVE_SCOTCH
}

}
#endif
