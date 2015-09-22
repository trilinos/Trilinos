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
#ifndef _ZOLTAN2_ALGPARMETIS_HPP_
#define _ZOLTAN2_ALGPARMETIS_HPP_

#include <Zoltan2_GraphModel.hpp>
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
template <typename Adapter>
class AlgParMETIS : public Algorithm<Adapter>
{
public:
  AlgParMETIS(const RCP<const Environment> &env,
              const RCP<const Comm<int> > &problemComm,
              const RCP<GraphModel<typename Adapter::base_adapter_t> > &model
  )
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

template <typename Adapter>
class AlgParMETIS : public Algorithm<Adapter>
{
public:

  typedef GraphModel<typename Adapter::base_adapter_t> graphModel_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;

  typedef idx_t  pm_idx_t;
  typedef real_t pm_real_t;

  /*! ParMETIS constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param model a graph
   *
   *  Preconditions: The parameters in the environment have been processed.
   *  TODO:  THIS IS A MINIMAL CONSTRUCTOR FOR NOW.
   *  TODO:  WHEN ADD PARMETIS ORDERING, MOVE PARMETIS GRAPH CONSTRUCTION
   *  TODO:  TO THE CONSTRUCTOR SO THAT CODE MAY BE SHARED.
   */
  AlgParMETIS(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<graphModel_t> &model__) :
    env(env__), problemComm(problemComm__), 
    mpicomm(TeuchosConst2MPI(problemComm__)),
    model(model__)
  { }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  MPI_Comm mpicomm;
  const RCP<GraphModel<typename Adapter::base_adapter_t> > model;

  void scale_weights(size_t n, ArrayView<StridedData<lno_t, scalar_t> > &fwgts,
                     pm_idx_t *iwgts);
};


/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgParMETIS<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

  int np = problemComm->getSize();

  // Get vertex info
  ArrayView<const gno_t> vtxgnos;
  ArrayView<StridedData<lno_t, scalar_t> > xyz;
  ArrayView<StridedData<lno_t, scalar_t> > vwgts;
  int nVwgt = model->getNumWeightsPerVertex();
  size_t nVtx = model->getVertexList(vtxgnos, xyz, vwgts);
  pm_idx_t pm_nVtx;
  TPL_Traits<pm_idx_t,size_t>::ASSIGN_TPL_T(pm_nVtx, nVtx, env);

  pm_idx_t *pm_vwgts = NULL;
  if (nVwgt) {
    pm_vwgts = new pm_idx_t[nVtx*nVwgt];
    scale_weights(nVtx, vwgts, pm_vwgts);
  }

  // Get edge info
  ArrayView<const gno_t> adjgnos;
  ArrayView<const int>   procs;
  ArrayView<const lno_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > ewgts;
  int nEwgt = model->getNumWeightsPerEdge();
  size_t nEdge = model->getEdgeList(adjgnos, procs, offsets, ewgts);

  pm_idx_t *pm_ewgts = NULL;
  if (nEwgt) {
    pm_ewgts = new pm_idx_t[nEdge*nEwgt]; 
    scale_weights(nEdge, ewgts, pm_ewgts);
  }

  // Convert index types for edges, if needed
  pm_idx_t *pm_offsets;  
  TPL_Traits<pm_idx_t,lno_t>::ASSIGN_TPL_T_ARRAY(&pm_offsets, offsets, env);
  pm_idx_t *pm_adjs;  
  TPL_Traits<pm_idx_t,gno_t>::ASSIGN_TPL_T_ARRAY(&pm_adjs, adjgnos, env);

  // Build vtxdist
  pm_idx_t *pm_vtxdist = new pm_idx_t[np+1];
  pm_vtxdist[0] = 0;
  Teuchos::gatherAll(*problemComm, 1, &pm_nVtx, np, &(pm_vtxdist[1]));
  for (int i = 2; i <= np; i++)
    pm_vtxdist[i] += pm_vtxdist[i-1];

  // Create array for ParMETIS to return results in.
  // Note:  ParMETIS does not like NULL arrays,
  //        so add 1 to always have non-null.
  //        See Zoltan bug 4299.
  pm_idx_t *pm_partList = new pm_idx_t[nVtx+1];

  // Get target part sizes and imbalance tolerances

  pm_idx_t pm_nCon = (nVwgt == 0 ? 1 : pm_idx_t(nVwgt));
  pm_real_t *pm_partsizes = new pm_real_t[numGlobalParts*pm_nCon];
  for (pm_idx_t dim = 0; dim < pm_nCon; dim++) {
    if (!solution->criteriaHasUniformPartSizes(dim))
      for (size_t i=0; i<numGlobalParts; i++)
        pm_partsizes[i*pm_nCon+dim] = 
                     pm_real_t(solution->getCriteriaPartSize(dim,i));
    else
      for (size_t i=0; i<numGlobalParts; i++)
        pm_partsizes[i*pm_nCon+dim] = pm_real_t(1.) / pm_real_t(numGlobalParts);
  }
  pm_real_t *pm_imbTols = new pm_real_t[pm_nCon];
  for (pm_idx_t dim = 0; dim < pm_nCon; dim++)
    pm_imbTols[dim] = 1.05;  // TODO:  GET THE PARAMETER

  std::string parmetis_method("PARTKWAY");
  pm_idx_t pm_wgtflag = 2*(nVwgt > 0) + (nEwgt > 0);
  pm_idx_t pm_numflag = 0;

  pm_idx_t pm_nPart;
  TPL_Traits<pm_idx_t,size_t>::ASSIGN_TPL_T(pm_nPart, numGlobalParts, env);

  if (parmetis_method == "PARTKWAY") {

    pm_idx_t pm_edgecut = -1;
    pm_idx_t pm_options[3];
    pm_options[0] = 0;   // Use default options
    pm_options[1] = 0;   // Debug level (ignored if pm_options[0] == 0)
    pm_options[2] = 0;   // Seed (ignored if pm_options[0] == 0)

    ParMETIS_V3_PartKway(pm_vtxdist, pm_offsets, pm_adjs, pm_vwgts, pm_ewgts,
                         &pm_wgtflag, &pm_numflag, &pm_nCon, &pm_nPart,
                         pm_partsizes, pm_imbTols, pm_options,
                         &pm_edgecut, pm_partList, &mpicomm);
  }
  else if (parmetis_method == "ADAPTIVE_REPART") {
    // Get object sizes
    std::cout << "NOT READY FOR ADAPTIVE_REPART YET" << std::endl;
    exit(-1);
  }
  else if (parmetis_method == "PART_GEOM") {
    // Get coordinate info, too.
    std::cout << "NOT READY FOR PART_GEOM YET" << std::endl;
    exit(-1);
  }

  // Clean up 
  delete [] pm_vtxdist;
  delete [] pm_partsizes;
  delete [] pm_imbTols;

  // Load answer into the solution.

  ArrayRCP<part_t> partList;
  if (TPL_Traits<pm_idx_t, part_t>::OK_TO_CAST_TPL_T()) {
    partList = ArrayRCP<part_t>((part_t *)pm_partList, 0, nVtx, true);
  }
  else {
    // TODO Probably should have a TPL_Traits function to do the following
    partList = ArrayRCP<part_t>(new part_t[nVtx], 0, nVtx, true);
    for (size_t i = 0; i < nVtx; i++) {
      partList[i] = part_t(pm_partList[i]);
    }
    delete [] pm_partList;
  }

  solution->setParts(partList);

  env->memory("Zoltan2-ParMETIS: After creating solution");

  // Clean up copies made due to differing data sizes.
  TPL_Traits<pm_idx_t,lno_t>::DELETE_TPL_T_ARRAY(&pm_offsets);
  TPL_Traits<pm_idx_t,gno_t>::DELETE_TPL_T_ARRAY(&pm_adjs);

  if (nVwgt) delete [] pm_vwgts;
  if (nEwgt) delete [] pm_ewgts;
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

template <typename Adapter>
void AlgParMETIS<Adapter>::scale_weights(
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
