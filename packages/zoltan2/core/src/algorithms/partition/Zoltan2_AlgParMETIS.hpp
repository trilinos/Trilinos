// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGPARMETIS_HPP_
#define _ZOLTAN2_ALGPARMETIS_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_CommGraphModel.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>

/////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgParMETIS.hpp
//! \brief Interface to the third-party library ParMETIS
/////////////////////////////////////////////////////////////////////////////

#ifndef HAVE_ZOLTAN2_PARMETIS

// Error handling for when ParMETIS is requested
// but Zoltan2 not built with ParMETIS.

namespace Zoltan2 {
template <typename Adapter, typename Model=GraphModel<typename Adapter::base_adapter_t>>
class AlgParMETIS : public Algorithm<Adapter>
{
public:
    AlgParMETIS(const RCP<const Environment> &,
                const RCP<const Comm<int> > &,
                const RCP<const typename Adapter::base_adapter_t> &,
                const modelFlag_t& graphFlags_ = modelFlag_t())
  {
    throw std::runtime_error(
          "BUILD ERROR:  ParMETIS requested but not compiled into Zoltan2.\n"
          "Please set CMake flag Zoltan2_ENABLE_ParMETIS:BOOL=ON.");
  }
};
}

#endif

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ZOLTAN2_PARMETIS

#ifndef HAVE_ZOLTAN2_MPI

// ParMETIS requires compilation with MPI.
// If MPI is not available, make compilation fail.
#error "TPL ParMETIS requires compilation with MPI.  Configure with -DTPL_ENABLE_MPI:BOOL=ON or -DZoltan2_ENABLE_ParMETIS:BOOL=OFF"

#else

extern "C" {
#include "parmetis.h"
}

#if (PARMETIS_MAJOR_VERSION < 4)

// Zoltan2 requires ParMETIS v4.x.
// Make compilation fail for earlier versions of ParMETIS.
#error "Specified version of ParMETIS is not compatible with Zoltan2; upgrade to ParMETIS v4 or later, or build Zoltan2 without ParMETIS."

#else

// MPI and ParMETIS version requirements are met.  Proceed.

namespace Zoltan2 {

template <typename Adapter, typename Model=GraphModel<typename Adapter::base_adapter_t>>
class AlgParMETIS : public Algorithm<Adapter>
{
public:

  typedef GraphModel<typename Adapter::base_adapter_t> graphModel_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;

  typedef idx_t  pm_idx_t;
  typedef real_t pm_real_t;

  /*! ParMETIS constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter an adapter used to create the model
   *
   *  Preconditions: The parameters in the environment have been processed.
   *  TODO:  THIS IS A MINIMAL CONSTRUCTOR FOR NOW.
   *  TODO:  WHEN ADD PARMETIS ORDERING, MOVE PARMETIS GRAPH CONSTRUCTION
   *  TODO:  TO THE CONSTRUCTOR SO THAT CODE MAY BE SHARED.
   */
  AlgParMETIS(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<const typename Adapter::base_adapter_t> &adapter__,
              const modelFlag_t& graphFlags_ = modelFlag_t()) :
    env(env__), problemComm(problemComm__),
    adapter(adapter__), graphFlags(graphFlags_)
  { }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const typename Adapter::base_adapter_t> adapter;
  modelFlag_t graphFlags;

  void scale_weights(size_t n, ArrayView<StridedData<lno_t, scalar_t> > &fwgts,
                     pm_idx_t *iwgts);
};


/////////////////////////////////////////////////////////////////////////////
  template <typename Adapter, typename Model>
  void AlgParMETIS<Adapter, Model>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  int me = problemComm->getRank();
  int np = problemComm->getSize();

  // Get vertex info
  ArrayView<const gno_t> vtxgnos;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;

  const auto model = rcp(new Model(adapter, env, problemComm, graphFlags));

  int nVwgt = model->getNumWeightsPerVertex();
  size_t nVtx = model->getVertexList(vtxgnos, vwgts);
  pm_idx_t pm_nVtx;
  TPL_Traits<pm_idx_t,size_t>::ASSIGN(pm_nVtx, nVtx);

  pm_idx_t *pm_vwgts = NULL;
  if (nVwgt) {
    pm_vwgts = new pm_idx_t[nVtx*nVwgt];
    scale_weights(nVtx, vwgts, pm_vwgts);
  }

  // Get edge info
  ArrayView<const gno_t> adjgnos;
  ArrayView<const offset_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;
  int nEwgt = model->getNumWeightsPerEdge();
  size_t nEdge = model->getEdgeList(adjgnos, offsets, ewgts);

  pm_idx_t *pm_ewgts = NULL;
  if (nEwgt) {
    pm_ewgts = new pm_idx_t[nEdge*nEwgt];
    scale_weights(nEdge, ewgts, pm_ewgts);
  }

  // Convert index types for edges, if needed
  pm_idx_t *pm_offsets;
  TPL_Traits<pm_idx_t,const offset_t>::ASSIGN_ARRAY(&pm_offsets, offsets);
  pm_idx_t *pm_adjs;
  pm_idx_t pm_dummy_adj;
  if (nEdge)
    TPL_Traits<pm_idx_t,const gno_t>::ASSIGN_ARRAY(&pm_adjs, adjgnos);
  else
    pm_adjs = &pm_dummy_adj;  // ParMETIS does not like NULL pm_adjs;


  // Build vtxdist
  pm_idx_t *pm_vtxdist;
  ArrayView<size_t> vtxdist;
  model->getVertexDist(vtxdist);
  TPL_Traits<pm_idx_t,size_t>::ASSIGN_ARRAY(&pm_vtxdist, vtxdist);

  // ParMETIS does not like processors having no vertices.
  // Inspect vtxdist and remove from communicator procs that have no vertices
  RCP<Comm<int> > subcomm;
  MPI_Comm mpicomm;  // Note:  mpicomm is valid only while subcomm is in scope

  int nKeep = 0;
  if (np > 1) {
    Array<int> keepRanks(np);
    for (int i = 0; i < np; i++) {
      if ((pm_vtxdist[i+1] - pm_vtxdist[i]) > 0) {
        keepRanks[nKeep] = i;
        pm_vtxdist[nKeep] = pm_vtxdist[i];
        nKeep++;
      }
    }
    pm_vtxdist[nKeep] = pm_vtxdist[np];
    if (nKeep < np) {
      subcomm = problemComm->createSubcommunicator(keepRanks.view(0,nKeep));
      if (subcomm != Teuchos::null)
        mpicomm = Teuchos::getRawMpiComm(*subcomm);
      else
        mpicomm = MPI_COMM_NULL;
    }
    else {
      mpicomm = Teuchos::getRawMpiComm(*problemComm);
    }
  }
  else {
    mpicomm = Teuchos::getRawMpiComm(*problemComm);
  }

  // Create array for ParMETIS to return results in.
  pm_idx_t *pm_partList = NULL;
  if (nVtx) pm_partList = new pm_idx_t[nVtx];
  for (size_t i = 0; i < nVtx; i++) pm_partList[i] = 0;
  int pm_return = METIS_OK;

  if (mpicomm != MPI_COMM_NULL) {
    // If in ParMETIS' communicator (i.e., have vertices), call ParMETIS

    // Get target part sizes
    pm_idx_t pm_nCon = (nVwgt == 0 ? 1 : pm_idx_t(nVwgt));
    pm_real_t *pm_partsizes = new pm_real_t[numGlobalParts*pm_nCon];
    for (pm_idx_t dim = 0; dim < pm_nCon; dim++) {
      if (!solution->criteriaHasUniformPartSizes(dim))
        for (size_t i=0; i<numGlobalParts; i++)
          pm_partsizes[i*pm_nCon+dim] =
                       pm_real_t(solution->getCriteriaPartSize(dim,i));
      else
        for (size_t i=0; i<numGlobalParts; i++)
          pm_partsizes[i*pm_nCon+dim] = pm_real_t(1.)/pm_real_t(numGlobalParts);
    }

    // Get imbalance tolerances
    double tolerance = 1.1;
    const Teuchos::ParameterList &pl = env->getParameters();
    const Teuchos::ParameterEntry *pe = pl.getEntryPtr("imbalance_tolerance");
    if (pe) tolerance = pe->getValue<double>(&tolerance);

    // ParMETIS requires tolerance to be greater than 1.0;
    // fudge it if condition is not met
    if (tolerance <= 1.0) {
      if (me == 0)
        std::cerr << "Warning:  ParMETIS requires imbalance_tolerance > 1.0; "
                  << "to comply, Zoltan2 reset imbalance_tolerance to 1.01."
                  << std::endl;
      tolerance = 1.01;
    }

    pm_real_t *pm_imbTols = new pm_real_t[pm_nCon];
    for (pm_idx_t dim = 0; dim < pm_nCon; dim++)
      pm_imbTols[dim] = pm_real_t(tolerance);

    std::string parmetis_method("PARTKWAY");
    pe = pl.getEntryPtr("partitioning_approach");
    if (pe){
      std::string approach;
      approach = pe->getValue<std::string>(&approach);
      if ((approach == "repartition") || (approach == "maximize_overlap")) {
        if (nKeep > 1)
          // ParMETIS_V3_AdaptiveRepart requires two or more processors
          parmetis_method = "ADAPTIVE_REPART";
        else
          // Probably best to do PartKway if nKeep == 1;
          // I think REFINE_KWAY won't give a good answer in most use cases
          // parmetis_method = "REFINE_KWAY";
          parmetis_method = "PARTKWAY";
      }
    }

    // Other ParMETIS parameters?
    pm_idx_t pm_wgtflag = 2*(nVwgt > 0) + (nEwgt > 0);
    pm_idx_t pm_numflag = 0;
    pm_idx_t pm_edgecut = -1;
    pm_idx_t pm_options[METIS_NOPTIONS];
    pm_options[0] = 1;   // Use non-default options for some ParMETIS options
    for (int i = 1; i < METIS_NOPTIONS; i++)
      pm_options[i] = 0; // Default options
    pm_options[2] = 15;  // Matches default value used in Zoltan

    pm_idx_t pm_nPart;
    TPL_Traits<pm_idx_t,size_t>::ASSIGN(pm_nPart, numGlobalParts);

    if (parmetis_method == "PARTKWAY") {
      pm_return = ParMETIS_V3_PartKway(pm_vtxdist, pm_offsets, pm_adjs,
                                       pm_vwgts, pm_ewgts, &pm_wgtflag,
                                       &pm_numflag, &pm_nCon, &pm_nPart,
                                       pm_partsizes, pm_imbTols, pm_options,
                                       &pm_edgecut, pm_partList, &mpicomm);
    }
    else if (parmetis_method == "ADAPTIVE_REPART") {
      // Get object sizes:  pm_vsize
      // TODO:  get pm_vsize info from input adapter or graph model
      // TODO:  This is just a placeholder
      pm_idx_t *pm_vsize = new pm_idx_t[nVtx];
      for (size_t i = 0; i < nVtx; i++) pm_vsize[i] = 1;

      pm_real_t itr = 100.;  // Same default as in Zoltan
      pm_return = ParMETIS_V3_AdaptiveRepart(pm_vtxdist, pm_offsets, pm_adjs,
                                             pm_vwgts,
                                             pm_vsize, pm_ewgts, &pm_wgtflag,
                                             &pm_numflag, &pm_nCon, &pm_nPart,
                                             pm_partsizes, pm_imbTols,
                                             &itr, pm_options, &pm_edgecut,
                                             pm_partList, &mpicomm);
      delete [] pm_vsize;
    }
    // else if (parmetis_method == "REFINE_KWAY") {
    //   We do not currently have an execution path that calls REFINE_KWAY.
    //   pm_return = ParMETIS_V3_RefineKway(pm_vtxdist, pm_offsets, pm_adjs,
    //                                      pm_vwgts, pm_ewgts, &pm_wgtflag,
    //                                     &pm_numflag, &pm_nCon, &pm_nPart,
    //                                    pm_partsizes, pm_imbTols, pm_options,
    //                                      &pm_edgecut, pm_partList, &mpicomm);
    // }
    else {
      // We should not reach this condition.
      throw std::logic_error("\nInvalid ParMETIS method requested.\n");
    }

    // Clean up
    delete [] pm_partsizes;
    delete [] pm_imbTols;
  }

  // Load answer into the solution.

  ArrayRCP<part_t> partList;
  if (nVtx)
    TPL_Traits<part_t, pm_idx_t>::SAVE_ARRAYRCP(&partList, pm_partList, nVtx);
  TPL_Traits<pm_idx_t, part_t>::DELETE_ARRAY(&pm_partList);

  solution->setParts(partList);

  env->memory("Zoltan2-ParMETIS: After creating solution");

  // Clean up copies made due to differing data sizes.
  TPL_Traits<pm_idx_t,size_t>::DELETE_ARRAY(&pm_vtxdist);
  TPL_Traits<pm_idx_t,const lno_t>::DELETE_ARRAY(&pm_offsets);
  if (nEdge)
    TPL_Traits<pm_idx_t,const gno_t>::DELETE_ARRAY(&pm_adjs);

  if (nVwgt) delete [] pm_vwgts;
  if (nEwgt) delete [] pm_ewgts;

  if (pm_return != METIS_OK) {
    throw std::runtime_error(
          "\nParMETIS returned an error; no valid partition generated.\n"
          "Look for 'PARMETIS ERROR' in your output for more details.\n");
  }
}

/////////////////////////////////////////////////////////////////////////////
// Scale and round scalar_t weights (typically float or double) to
// ParMETIS' idx_t (typically int or long).
// subject to sum(weights) <= max_wgt_sum.
// Scale only if deemed necessary.
//
// Note that we use ceil() instead of round() to avoid
// rounding to zero weights.
// Based on Zoltan's scale_round_weights, mode 1

  template <typename Adapter, typename Model>
  void AlgParMETIS<Adapter, Model>::scale_weights(
  size_t n,
  ArrayView<StridedData<typename Adapter::lno_t,
                        typename Adapter::scalar_t> > &fwgts,
  pm_idx_t *iwgts
)
{
  const double INT_EPSILON = 1e-5;
  const int nWgt = fwgts.size();

  int *nonint_local = new int[nWgt+nWgt];
  int *nonint = nonint_local + nWgt;

  double *sum_wgt_local = new double[nWgt*4];
  double *max_wgt_local = sum_wgt_local + nWgt;
  double *sum_wgt = max_wgt_local + nWgt;
  double *max_wgt = sum_wgt + nWgt;

  for (int i = 0; i < nWgt; i++) {
    nonint_local[i] = 0;
    sum_wgt_local[i] = 0.;
    max_wgt_local[i] = 0;
  }

  // Compute local sums of the weights
  // Check whether all weights are integers
  for (int j = 0; j < nWgt; j++) {
    for (size_t i = 0; i < n; i++) {
      double fw = double(fwgts[j][i]);
      if (!nonint_local[j]) {
        pm_idx_t tmp = (pm_idx_t) floor(fw + .5); /* Nearest int */
        if (fabs((double)tmp-fw) > INT_EPSILON) {
          nonint_local[j] = 1;
        }
      }
      sum_wgt_local[j] += fw;
      if (fw > max_wgt_local[j]) max_wgt_local[j] = fw;
    }
  }

  Teuchos::reduceAll<int,int>(*problemComm, Teuchos::REDUCE_MAX, nWgt,
                              nonint_local,  nonint);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_SUM, nWgt,
                                 sum_wgt_local, sum_wgt);
  Teuchos::reduceAll<int,double>(*problemComm, Teuchos::REDUCE_MAX, nWgt,
                                 max_wgt_local, max_wgt);

  const double max_wgt_sum = double(std::numeric_limits<pm_idx_t>::max()/8);
  for (int j = 0; j < nWgt; j++) {
    double scale = 1.;

    // Scaling needed if weights are not integers or weights'
    // range is not sufficient
    if (nonint[j] || (max_wgt[j]<=INT_EPSILON) || (sum_wgt[j]>max_wgt_sum)) {
      /* Calculate scale factor */
      if (sum_wgt[j] != 0.) scale = max_wgt_sum/sum_wgt[j];
    }

    /* Convert weights to positive integers using the computed scale factor */
    for (size_t i = 0; i < n; i++)
      iwgts[i*nWgt+j] = (pm_idx_t) ceil(double(fwgts[j][i])*scale);
  }
  delete [] nonint_local;
  delete [] sum_wgt_local;
}

} // namespace Zoltan2

#endif // PARMETIS VERSION 4 OR HIGHER CHECK

#endif // HAVE_ZOLTAN2_MPI

#endif // HAVE_ZOLTAN2_PARMETIS

#endif
