// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
 *       (which must be at least one  TODO WHY?) is taken to be \c vwgts.size().
 *   \param mcNorm the multiCriteria norm to be used if the number of weights is
 *             greater than one.
 *   \param targetNumParts  input:  number of requested parts
 *   \param numExistingParts  on return this is the maximum part ID + 1.
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
 *      \c numExistingParts totals of global number of objects in each part.
 *     Suppose the number of weights is \c W.  If
 *     W is 1, then on return this is an array of length \c 2*numExistingParts .
 *     The first \c numExistingParts entries are the count of objects in each 
 *     part and the second is the total weight in each part.
 *     If \c W is greater than one, then the length of this array is
 *     \c (2+W)*numExistingParts .
 *     The first \c numExistingParts entries are the count of objects in each 
 *     part.
 *     The next \c numExistingParts entries are the sum of the normed weights in
 *     each part.
 *     The final entries are the sum of the individual weights in each part,
 *     by weight index by part number.  The array is allocated here.
 *
 * () must be called by all processes in \c comm.
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
    part_t targetNumParts,
    part_t &numExistingParts,
    part_t &numNonemptyParts,
    ArrayRCP<RCP<BaseClassMetrics<scalar_t> > > &metrics,
    ArrayRCP<scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalSumsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  numExistingParts = numNonemptyParts = 0;

  int numMetrics = 1;                       // "object count"
  if (vwgtDim) numMetrics++;                // "normed weight" or "weight 0"
  if (vwgtDim > 1) numMetrics += vwgtDim;   // "weight n"

  auto next = metrics.size(); // where we will start filling
  typedef ImbalanceMetrics<scalar_t> im_t;
  for(int n = 0; n < numMetrics; ++n) {
    RCP<im_t> newMetric = addNewMetric<im_t, scalar_t>(env, metrics);
    if (vwgtDim > 1) {
      newMetric->setNorm(multiCriteriaNorm(mcNorm));
    }
  }

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

  part_t maxPartPlusOne = globalNum[1] + 1;  // Range of possible part IDs:
                                             // [0,maxPartPlusOne)

  part_t globalSumSize = maxPartPlusOne * numMetrics;

  scalar_t * sumBuf = new scalar_t[globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, sumBuf);
  globalSums = arcp(sumBuf, 0, globalSumSize);

  //////////////////////////////////////////////////////////
  // Calculate the local totals by part.

  scalar_t *localBuf = new scalar_t[globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, localBuf);
  memset(localBuf, 0, sizeof(scalar_t) * globalSumSize);

  scalar_t *obj = localBuf;              // # of objects

  for (lno_t i=0; i < localNumObj; i++)
    obj[part[i]]++;

  if (numMetrics > 1){

    scalar_t *wgt = localBuf+maxPartPlusOne; // single normed weight or weight 0
    try{
      normedPartWeights<scalar_t, lno_t, part_t>(env, maxPartPlusOne,
        part, vwgts, mcNorm, wgt);
    }
    Z2_FORWARD_EXCEPTIONS

    // This code assumes the solution has the part ordered the
    // same way as the user input.  (Bug 5891 is resolved.)
    if (vwgtDim > 1){
      wgt += maxPartPlusOne;         // individual weights
      for (int vdim = 0; vdim < vwgtDim; vdim++){
        for (lno_t i=0; i < localNumObj; i++)
          wgt[part[i]] += vwgts[vdim][i];
        wgt += maxPartPlusOne;
      }
    }
  }

  // Metric: local sums on process
  metrics[next]->setName("object count");
  metrics[next]->setMetricValue("local sum", localNumObj);

  next++;

  if (numMetrics > 1){
    scalar_t *wgt = localBuf+maxPartPlusOne; // single normed weight or weight 0
    scalar_t total = 0.0;

    for (int p=0; p < maxPartPlusOne; p++){
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
        wgt += maxPartPlusOne;
        total = 0.0;
        for (int p=0; p < maxPartPlusOne; p++){
          total += wgt[p];
        }

        std::ostringstream oss;
        oss << "weight " << vdim;

        metrics[next]->setName(oss.str());
        metrics[next]->setMetricValue("local sum", total);

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
  next = metrics.size() - numMetrics; // MDM - this is most likely temporary 
                                      // to preserve the format here - we are 
                                      // now filling a larger array so we may 
                                      // not have started at 0


  ArrayView<scalar_t> objVec(obj, maxPartPlusOne);
  getStridedStats<scalar_t>(objVec, 1, 0, min, max, sum);
  if (maxPartPlusOne < targetNumParts)
    min = scalar_t(0);  // Some of the target parts are empty

  metrics[next]->setMetricValue("global minimum", min);
  metrics[next]->setMetricValue("global maximum", max);
  metrics[next]->setMetricValue("global sum", sum);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + maxPartPlusOne; // single normed weight or weight 0

    ArrayView<scalar_t> normedWVec(wgt, maxPartPlusOne);
    getStridedStats<scalar_t>(normedWVec, 1, 0, min, max, sum);
    if (maxPartPlusOne < targetNumParts)
      min = scalar_t(0);  // Some of the target parts are empty

    metrics[next]->setMetricValue("global minimum", min);
    metrics[next]->setMetricValue("global maximum", max);
    metrics[next]->setMetricValue("global sum", sum);
    next++;

    if (vwgtDim > 1){
      for (int vdim=0; vdim < vwgtDim; vdim++){
        wgt += maxPartPlusOne;       // individual weights
        ArrayView<scalar_t> fromVec(wgt, maxPartPlusOne);
        getStridedStats<scalar_t>(fromVec, 1, 0, min, max, sum);
        if (maxPartPlusOne < targetNumParts)
          min = scalar_t(0);  // Some of the target parts are empty

        metrics[next]->setMetricValue("global minimum", min);
        metrics[next]->setMetricValue("global maximum", max);
        metrics[next]->setMetricValue("global sum", sum);
        next++;
      }
    }
  }

  //////////////////////////////////////////////////////////
  // How many parts do we actually have.

  numExistingParts = maxPartPlusOne;
  obj = sumBuf;               // # of objects

  /*for (part_t p=nparts-1; p > 0; p--){
    if (obj[p] > 0) break;
    numExistingParts--;
    }*/

  numNonemptyParts = numExistingParts;

  for (part_t p=0; p < numExistingParts; p++)
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
 *   \param numExistingParts on return is the max Part ID + 1.
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
  void imbalanceMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    multiCriteriaNorm mcNorm,
    const Adapter *ia,
    const PartitioningSolution<Adapter> *solution,
    const ArrayView<const typename Adapter::part_t> &partArray,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel,
    typename Adapter::part_t &numExistingParts,
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
      RCP<GraphModel<base_adapter_t> > graph;
      if (graphModel == Teuchos::null) {
	std::bitset<NUM_MODEL_FLAGS> modelFlags;
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
  scalar_t *psizes = NULL;

  ArrayRCP<ArrayRCP<scalar_t> > partSizes(numCriteria);
  if (solution) {
    targetNumParts = solution->getTargetGlobalNumberOfParts();
    for (int dim=0; dim < numCriteria; dim++){
      if (solution->criteriaHasUniformPartSizes(dim) != true){
        psizes = new scalar_t [targetNumParts];
        env->localMemoryAssertion(__FILE__, __LINE__, targetNumParts, psizes);
        for (part_t i=0; i < targetNumParts; i++){
          psizes[i] = solution->getCriteriaPartSize(dim, i);
        }
        partSizes[dim] = arcp(psizes, 0, targetNumParts, true);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Get number of parts, and the number that are non-empty.
  // Get sums per part of objects, individual weights, and normed weight sums.

  ArrayRCP<scalar_t> globalSums;

  int initialMetricCount = metrics.size();
  try{
    globalSumsByPart<scalar_t, lno_t, part_t>(env, comm,
      partArray, nWeights, weights.view(0, numCriteria), mcNorm,
      targetNumParts, numExistingParts, numNonemptyParts, metrics, globalSums);
  }
  Z2_FORWARD_EXCEPTIONS

  int addedMetricsCount = metrics.size() - initialMetricCount;

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the object count.
  // (Use first index of part sizes.)

  int index = initialMetricCount;

  scalar_t *objCount  = globalSums.getRawPtr();
  scalar_t min, max, avg;
  psizes=NULL;

  if (partSizes[0].size() > 0)
    psizes = partSizes[0].getRawPtr();

  scalar_t gsum = metrics[index]->getMetricValue("global sum");
  computeImbalances<scalar_t, part_t>(numExistingParts, targetNumParts, psizes,
      gsum, objCount, min, max, avg);

  metrics[index]->setMetricValue("global average", gsum / targetNumParts);

  metrics[index]->setMetricValue("maximum imbalance", 1.0 + max);
  metrics[index]->setMetricValue("average imbalance", avg);

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the normed weight sum.

  scalar_t *wgts = globalSums.getRawPtr() + numExistingParts;

  if (addedMetricsCount > 1){
    ++index;
    gsum = metrics[index]->getMetricValue("global sum");
    computeImbalances<scalar_t, part_t>(numExistingParts, targetNumParts,
      numCriteria, partSizes.view(0, numCriteria), gsum, wgts, min, max, avg);

    metrics[index]->setMetricValue("global average", gsum / targetNumParts);

    metrics[index]->setMetricValue("maximum imbalance", 1.0 + max);
    metrics[index]->setMetricValue("average imbalance", avg);

    if (addedMetricsCount > 2){

    ///////////////////////////////////////////////////////////////////////////
    // Compute imbalances for each individual weight.

      ++index;

      for (int vdim=0; vdim < numCriteria; vdim++){
        wgts += numExistingParts;
        psizes = NULL;

        if (partSizes[vdim].size() > 0)
           psizes = partSizes[vdim].getRawPtr();

        gsum = metrics[index]->getMetricValue("global sum");
        computeImbalances<scalar_t, part_t>(numExistingParts, targetNumParts, 
                                            psizes, gsum, wgts, min, max, avg);

        metrics[index]->setMetricValue("global average", gsum / targetNumParts);

        metrics[index]->setMetricValue("maximum imbalance", 1.0 + max);
        metrics[index]->setMetricValue("average imbalance", avg);
        index++;
      }
    }

  }
  env->debug(DETAILED_STATUS, "Exiting objectMetrics");
}

/*! \brief Print out header info for imbalance metrics.
 */
template <typename scalar_t, typename part_t>
void printImbalanceMetricsHeader(
  std::ostream &os, 
  part_t targetNumParts, 
  part_t numExistingParts, 
  part_t numNonemptyParts)
{
  os << "Imbalance Metrics: (" << numExistingParts << " existing parts)";
  if (numNonemptyParts < numExistingParts) {
    os << " (" << numNonemptyParts << " of which are non-empty)";
  }
  os << std::endl;
  if (targetNumParts != numExistingParts) {
    os << "Target number of parts is " << targetNumParts << std::endl;
  }
  ImbalanceMetrics<scalar_t>::printHeader(os);
}

/*! \brief Print out list of imbalance metrics.
 */
template <typename scalar_t, typename part_t>
void printImbalanceMetrics(
  std::ostream &os, 
  part_t targetNumParts, 
  part_t numExistingParts, 
  part_t numNonemptyParts, 
  const ArrayView<RCP<BaseClassMetrics<scalar_t> > > &infoList)
{
  printImbalanceMetricsHeader<scalar_t, part_t>(os, targetNumParts, 
                                                numExistingParts,
                                                numNonemptyParts);
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
void printImbalanceMetrics(
  std::ostream &os, 
  part_t targetNumParts, 
  part_t numExistingParts, 
  part_t numNonemptyParts, 
  RCP<BaseClassMetrics<scalar_t>> metricValue)
{
  printImbalanceMetricsHeader<scalar_t, part_t>(os, targetNumParts,
                                                numExistingParts,
                                                numNonemptyParts);
  metricValue->printLine(os);
}

} //namespace Zoltan2


#endif
