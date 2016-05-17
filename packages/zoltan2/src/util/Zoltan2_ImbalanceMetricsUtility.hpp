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

/*! \file Zoltan2_ImbalanceMetricsUtility.hpp
 */

#ifndef ZOLTAN2_IMBALANCEMETRICSUTILITY_HPP
#define ZOLTAN2_IMBALANCEMETRICSUTILITY_HPP

#include <Zoltan2_ImbalanceMetrics.hpp>
#include <Zoltan2_MetricUtility.hpp>

namespace Zoltan2{

/*! \brief Given the local partitioning, compute the global sums in each part.
 *
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param part   \c part[i] is the part ID for local object \c i
 *   \param vwgts  \c vwgts[w] is the StridedData object
 *       representing weight index \c w. The number of weights
 *       (which must be at least one  TODO  WHY?) is taken to be \c vwgts.size().
 *   \param mcNorm the multiCriteria norm, to be used if the number of weights is
 *             greater than one.
 *   \param numParts  on return this is the global number of parts.
 *   \param numNonemptyParts  on return this is the number of those
 *          parts that are non-empty.
 *   \param metrics on return points to a list of named MetricValues objects
 *     that each contains the global min, max and avg over parts of
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 0", "weight 1" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one set of non-uniform weights were given, then
 *     "object count" and "weight 0" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 0", "weight 1", and so on.
 *   \param globalSums If weights are uniform, the globalSums is the
 *      \c numParts totals of global number of objects in each part.
 *     Suppose the number of weights is \c W.  If
 *     W is 1, then on return this is an array of length \c 2*numParts .
 *     The first \c numParts entries are the count of objects in each part,
 *     and the second is the total weight in each part.
 *     If \c W is greater than one, then the length of this array is
 *     \c (2+W)*numParts .
 *     The first \c numParts entries are the count of objects in each part.
 *     The next \c numParts entries are the sum of the normed weights in
 *     each part.
 *     The final entries are the sum of the individual weights in each part,
 *     by weight index by part number.  The array is allocated here.
 *
 *
() must be called by all processes in \c comm.
 * The imbalance metrics are not yet set in the MetricValues objects,
 * because they require part size information.
 */

template <typename scalar_t, typename lno_t, typename part_t>
  void globalSumsByPart(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const ArrayView<const part_t> &part,
    int vwgtDim,
    const ArrayView<StridedData<lno_t, scalar_t> > &vwgts,
    multiCriteriaNorm mcNorm,
    part_t &numParts,
    part_t &numNonemptyParts,
    ArrayRCP<RCP<BaseClassMetrics<scalar_t> > > &metrics,
    ArrayRCP<scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalSumsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  numParts = numNonemptyParts = 0;

  int numMetrics = 1;                       // "object count"
  if (vwgtDim) numMetrics++;                // "normed weight" or "weight 0"
  if (vwgtDim > 1) numMetrics += vwgtDim;   // "weight n"

  // add some more metrics to the array
  typedef ImbalanceMetrics<scalar_t> mv_t;
  typedef typename ArrayRCP<RCP<BaseClassMetrics<scalar_t> > >::size_type array_size_type;
  metrics.resize( metrics.size() + numMetrics );
  for( array_size_type n = metrics.size() - numMetrics; n < metrics.size(); ++n )
  {
	  mv_t * newMetric = new mv_t;									// allocate the new memory

	  // moved this here because we now allocate the polymorphic classes
	  // we should probably reorganize these functions so all data setup is done on the derived classes
	  // then as a last step we can insert them into the general array of MetricBase types
	  if (vwgtDim > 1)
	    newMetric->setNorm(multiCriteriaNorm(mcNorm));

	  env->localMemoryAssertion(__FILE__,__LINE__,1,newMetric);		// check errors
	  metrics[n] = rcp( newMetric ); 				// create the new members
  }
  array_size_type next = metrics.size() - numMetrics; // MDM - this is most likely temporary to preserve the format here - we are now filling a larger array so we may not have started at 0

  //////////////////////////////////////////////////////////
  // Figure out the global number of parts in use.
  // Verify number of vertex weights is the same everywhere.

  lno_t localNumObj = part.size();
  part_t localNum[2], globalNum[2];
  localNum[0] = static_cast<part_t>(vwgtDim);
  localNum[1] = 0;

  for (lno_t i=0; i < localNumObj; i++)
    if (part[i] > localNum[1]) localNum[1] = part[i];

  try{
    reduceAll<int, part_t>(*comm, Teuchos::REDUCE_MAX, 2,
      localNum, globalNum);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->globalBugAssertion(__FILE__,__LINE__,
    "inconsistent number of vertex weights",
    globalNum[0] == localNum[0], DEBUG_MODE_ASSERTION, comm);

  part_t nparts = globalNum[1] + 1;

  part_t globalSumSize = nparts * numMetrics;
  scalar_t * sumBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, sumBuf);
  globalSums = arcp(sumBuf, 0, globalSumSize);

  //////////////////////////////////////////////////////////
  // Calculate the local totals by part.

  scalar_t *localBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, localBuf);
  memset(localBuf, 0, sizeof(scalar_t) * globalSumSize);

  scalar_t *obj = localBuf;              // # of objects

  for (lno_t i=0; i < localNumObj; i++)
    obj[part[i]]++;

  if (numMetrics > 1){

    scalar_t *wgt = localBuf + nparts; // single normed weight or weight 0
    try{
      normedPartWeights<scalar_t, lno_t, part_t>(env, nparts,
        part, vwgts, mcNorm, wgt);
    }
    Z2_FORWARD_EXCEPTIONS

    // This code assumes the solution has the part ordered the
    // same way as the user input.  (Bug 5891 is resolved.)
    if (vwgtDim > 1){
      wgt += nparts;         // individual weights
      for (int vdim = 0; vdim < vwgtDim; vdim++){
        for (lno_t i=0; i < localNumObj; i++)
          wgt[part[i]] += vwgts[vdim][i];
        wgt += nparts;
      }
    }
  }

  // Metric: local sums on process
  metrics[next]->setName("object count");
  metrics[next]->setMetricValue("local sum", localNumObj);

  next++;

  if (numMetrics > 1){
    scalar_t *wgt = localBuf + nparts; // single normed weight or weight 0
    scalar_t total = 0.0;

    for (int p=0; p < nparts; p++){
      total += wgt[p];
    }

    if (vwgtDim == 1)
      metrics[next]->setName("weight 0");
    else
      metrics[next]->setName("normed weight");

    metrics[next]->setMetricValue("local sum", total);

    next++;

    if (vwgtDim > 1){
      for (int vdim = 0; vdim < vwgtDim; vdim++){
        wgt += nparts;
        total = 0.0;
        for (int p=0; p < nparts; p++){
          total += wgt[p];
        }

        std::ostringstream oss;
        oss << "weight " << vdim;

        metrics[next]->setName(oss.str());
        metrics[next]->setMetricValue("local sum", total);

        std::cout << "***** Logging ImbalanceMetric Data Index: " << next << "  Name:" << oss.str() << "   local sum: " << total << std::endl;

        next++;
      }
    }

  }

  //////////////////////////////////////////////////////////
  // Obtain global totals by part.

  try{
    reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, globalSumSize,
      localBuf, sumBuf);
  }
  Z2_THROW_OUTSIDE_ERROR(*env);

  delete [] localBuf;

  //////////////////////////////////////////////////////////
  // Global sum, min, max, and average over all parts

  obj = sumBuf;                     // # of objects
  scalar_t min=0, max=0, sum=0;
  next = metrics.size() - numMetrics; // MDM - this is most likely temporary to preserve the format here - we are now filling a larger array so we may not have started at 0


  ArrayView<scalar_t> objVec(obj, nparts);
  getStridedStats<scalar_t>(objVec, 1, 0, min, max, sum);

  metrics[next]->setMetricValue("global minimum", min);
  metrics[next]->setMetricValue("global maximum", max);
  metrics[next]->setMetricValue("global sum", sum);
  metrics[next]->setMetricValue("global average", sum / nparts);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + nparts;        // single normed weight or weight 0

    ArrayView<scalar_t> normedWVec(wgt, nparts);
    getStridedStats<scalar_t>(normedWVec, 1, 0, min, max, sum);

    metrics[next]->setMetricValue("global minimum", min);
    metrics[next]->setMetricValue("global maximum", max);
    metrics[next]->setMetricValue("global sum", sum);
    metrics[next]->setMetricValue("global average", sum / nparts);
    next++;

    if (vwgtDim > 1){
      for (int vdim=0; vdim < vwgtDim; vdim++){
        wgt += nparts;       // individual weights
        ArrayView<scalar_t> fromVec(wgt, nparts);
        getStridedStats<scalar_t>(fromVec, 1, 0, min, max, sum);

        metrics[next]->setMetricValue("global minimum", min);
        metrics[next]->setMetricValue("global maximum", max);
        metrics[next]->setMetricValue("global sum", sum);
        metrics[next]->setMetricValue("global average", sum / nparts);
        next++;
      }
    }
  }

  //////////////////////////////////////////////////////////
  // How many parts do we actually have.

  numParts = nparts;
  obj = sumBuf;               // # of objects

  /*for (part_t p=nparts-1; p > 0; p--){
    if (obj[p] > 0) break;
    numParts--;
    }*/

  numNonemptyParts = numParts;

  for (part_t p=0; p < numParts; p++)
    if (obj[p] == 0) numNonemptyParts--;

  env->debug(DETAILED_STATUS, "Exiting globalSumsByPart");
}

/*! \brief Compute imbalance metrics for a distribution.
 *
 *   \param env   The problem environment.
 *   \param comm  The problem communicator.
 *   \param ia the InputAdapter object which corresponds to the Solution.
 *   \param solution the PartitioningSolution to be evaluated.
 *   \param mcNorm  is the multicriteria norm to use if the number of weights
 *           is greater than one.  See the multiCriteriaNorm enumerator for
 *           \c mcNorm values.
 *   \param graphModel the graph model.
 *   \param numParts on return is the global number of parts in the solution
 *   \param numNonemptyParts on return is the global number of parts to which
 *                                objects are assigned.
 *   \param metrics on return points to a list of named MetricValues objects
 *     that each contains the global min, max and avg over parts and
 *     also imbalance measures of
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 0", "weight 1" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one set of non-uniform weights were given, then
 *     "object count" and "weight 0" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 0", "weight 1", and so on.
 *
 *  objectMetrics() must be called by all processes in \c comm.
 *  See the metricOffset enumerator in the MetricValues class for the
 *  interpretation of the metric quantities.
 *   \todo check that part sizes sum to one if we're doing COMPLEX_ASSERTION
 */

template <typename Adapter>
  void objectMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    multiCriteriaNorm mcNorm,
    const Adapter *ia,
    const PartitioningSolution<Adapter> *solution,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel,
    typename Adapter::part_t &numParts,
    typename Adapter::part_t &numNonemptyParts,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics)
{
  env->debug(DETAILED_STATUS, "Entering objectMetrics");

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef StridedData<lno_t, scalar_t> sdata_t;

  // Local number of objects.

  size_t numLocalObjects = ia->getLocalNumIDs();

  // Parts to which objects are assigned.

  const part_t *parts;
  if (solution) {
    // User provided a partitioning solution; use it.
    parts = solution->getPartListView();
    env->localInputAssertion(__FILE__, __LINE__, "parts not set",
      ((numLocalObjects == 0) || parts), BASIC_ASSERTION);
  } else {
    // User did not provide a partitioning solution;
    // Use input adapter partition.

    parts = NULL;
    ia->getPartsView(parts);
    if (parts == NULL) {
      // User has not provided input parts in input adapter
      part_t *procs = new part_t [numLocalObjects];
      for (size_t i = 0; i < numLocalObjects; i++) procs[i] = comm->getRank();
      parts = procs;
    }
  }
  ArrayView<const part_t> partArray(parts, numLocalObjects);

  // Weights, if any, for each object.

  int nWeights = ia->getNumWeightsPerID();
  int numCriteria = (nWeights > 0 ? nWeights : 1);
  Array<sdata_t> weights(numCriteria);

  if (nWeights == 0){
    // One set of uniform weights is implied.
    // StridedData default constructor creates length 0 strided array.
    weights[0] = sdata_t();
  }
  else{
    // whether vertex degree is ever used as vertex weight.
    enum BaseAdapterType adapterType = ia->adapterType();
    bool useDegreeAsWeight = false;
    if (adapterType == GraphAdapterType) {
      useDegreeAsWeight = reinterpret_cast<const GraphAdapter
	<typename Adapter::user_t, typename Adapter::userCoord_t> *>(ia)->
	useDegreeAsWeight(0);
    } else if (adapterType == MatrixAdapterType) {
      useDegreeAsWeight = reinterpret_cast<const MatrixAdapter
	<typename Adapter::user_t, typename Adapter::userCoord_t> *>(ia)->
	useDegreeAsWeight(0);
    } else if (adapterType == MeshAdapterType) {
      useDegreeAsWeight =
	reinterpret_cast<const MeshAdapter<typename Adapter::user_t> *>(ia)->
	useDegreeAsWeight(0);
    }
    if (useDegreeAsWeight) {
      ArrayView<const gno_t> Ids;
      ArrayView<sdata_t> vwgts;
      if (graphModel == Teuchos::null) {
	std::bitset<NUM_MODEL_FLAGS> modelFlags;
	RCP<GraphModel<base_adapter_t> > graph;
	const RCP<const base_adapter_t> bia =
	  rcp(dynamic_cast<const base_adapter_t *>(ia), false);
	graph = rcp(new GraphModel<base_adapter_t>(bia,env,comm,modelFlags));
	graph->getVertexList(Ids, vwgts);
      } else {
	graphModel->getVertexList(Ids, vwgts);
      }
      scalar_t *wgt = new scalar_t[numLocalObjects];
      for (int i=0; i < nWeights; i++){
	for (size_t j=0; j < numLocalObjects; j++) {
	  wgt[j] = vwgts[i][j];
	}
	ArrayRCP<const scalar_t> wgtArray(wgt,0,numLocalObjects,false);
	weights[i] = sdata_t(wgtArray, 1);
      }
    } else {
      for (int i=0; i < nWeights; i++){
	const scalar_t *wgt;
	int stride;
	ia->getWeightsView(wgt, stride, i);
	ArrayRCP<const scalar_t> wgtArray(wgt,0,stride*numLocalObjects,false);
	weights[i] = sdata_t(wgtArray, stride);
      }
    }
  }

  // Relative part sizes, if any, assigned to the parts.

  part_t targetNumParts = comm->getSize();

  if (solution)
    targetNumParts = solution->getTargetGlobalNumberOfParts();

  scalar_t *psizes = NULL;

  ArrayRCP<ArrayRCP<scalar_t> > partSizes(numCriteria);
  for (int dim=0; dim < numCriteria; dim++){
    if (solution)
    if (solution->criteriaHasUniformPartSizes(dim) != true){
      psizes = new scalar_t [targetNumParts];
      env->localMemoryAssertion(__FILE__, __LINE__, numParts, psizes);
      for (part_t i=0; i < targetNumParts; i++){
        psizes[i] = solution->getCriteriaPartSize(dim, i);
      }
      partSizes[dim] = arcp(psizes, 0, targetNumParts, true);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Get number of parts, and the number that are non-empty.
  // Get sums per part of objects, individual weights, and normed weight sums.

  ArrayRCP<scalar_t> globalSums;

  try{
    globalSumsByPart<scalar_t, lno_t, part_t>(env, comm,
      partArray, nWeights, weights.view(0, numCriteria), mcNorm,
      numParts, numNonemptyParts, metrics, globalSums);
  }
  Z2_FORWARD_EXCEPTIONS

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the object count.
  // (Use first index of part sizes.)

  scalar_t *objCount  = globalSums.getRawPtr();
  scalar_t min, max, avg;
  psizes=NULL;

  if (partSizes[0].size() > 0)
    psizes = partSizes[0].getRawPtr();

  computeImbalances<scalar_t, part_t>(numParts, targetNumParts, psizes,
      metrics[0]->getMetricValue("global sum"), objCount,
      min, max, avg);

  // MDM - note that this indexing works because we have Metrics first - but we should generalize this - perhaps we will not start at index 0 in the new scheme
  metrics[0]->setMetricValue("maximum imbalance", 1.0 + max);
  metrics[0]->setMetricValue("average imbalance", avg);

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the normed weight sum.

  scalar_t *wgts = globalSums.getRawPtr() + numParts;

  if (metrics.size() > 1){

    computeImbalances<scalar_t, part_t>(numParts, targetNumParts,
      numCriteria, partSizes.view(0, numCriteria),
      metrics[1]->getMetricValue("global sum"), wgts,
      min, max, avg);

    metrics[1]->setMetricValue("maximum imbalance", 1.0 + max);
    metrics[1]->setMetricValue("average imbalance", avg);

    if (metrics.size() > 2){

    ///////////////////////////////////////////////////////////////////////////
    // Compute imbalances for each individual weight.

      int next = 2;	// MDM - also generalize here for the new scheme where we build the list - currently still fine because ImbalanceMetrics happen to be loaded first

      for (int vdim=0; vdim < numCriteria; vdim++){
        wgts += numParts;
        psizes = NULL;

        if (partSizes[vdim].size() > 0)
           psizes = partSizes[vdim].getRawPtr();

        computeImbalances<scalar_t, part_t>(numParts, targetNumParts, psizes,
          metrics[next]->getMetricValue("global sum"), wgts, min, max, avg);

        metrics[next]->setMetricValue("maximum imbalance", 1.0 + max);
        metrics[next]->setMetricValue("average imbalance", avg);
        next++;
      }
    }

  }
  env->debug(DETAILED_STATUS, "Exiting objectMetrics");
}

/*! \brief Print out header info for imbalance metrics.
 */
template <typename scalar_t, typename part_t>
void printImbalanceMetricsHeader(std::ostream &os, part_t targetNumParts, part_t numParts, part_t numNonemptyParts)
{
  os << "Imbalance Metrics: Number of parts is " << numParts;
  if (numNonemptyParts < numParts) {
    os << " (" << numNonemptyParts << " of which are non-empty)";
  }
  os << std::endl;
  if (targetNumParts != numParts) {
    os << "Target number of parts is " << targetNumParts << std::endl;
  }
  ImbalanceMetrics<scalar_t>::printHeader(os);
}

/*! \brief Print out list of imbalance metrics.
 */
template <typename scalar_t, typename part_t>
void printImbalanceMetrics(std::ostream &os, part_t targetNumParts, part_t numParts, part_t numNonemptyParts, const ArrayView<RCP<BaseClassMetrics<scalar_t>>> &infoList)
{
  printImbalanceMetricsHeader<scalar_t, part_t>(os, targetNumParts, numParts, numNonemptyParts);
  for (int i=0; i < infoList.size(); i++) {
    if (infoList[i]->getName() != METRICS_UNSET_STRING) {
      infoList[i]->printLine(os);
    }
  }
  os << std::endl;
}

/*! \brief Print out header and a single imbalance metric.
 */
template <typename scalar_t, typename part_t>
void printImbalanceMetrics(std::ostream &os, part_t targetNumParts, part_t numParts, part_t numNonemptyParts, RCP<BaseClassMetrics<scalar_t>> metricValue)
{
  printImbalanceMetricsHeader<scalar_t, part_t>(os, targetNumParts, numParts, numNonemptyParts);
  metricValue->printLine(os);
}

} //namespace Zoltan2


#endif
