// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_MetricFunctions.hpp
 *  \functions which were with the metric classes but do not explicitly depend on them
 */

#ifndef ZOLTAN2_METRICFUNCTIONS_HPP
#define ZOLTAN2_METRICFUNCTIONS_HPP

#include <Zoltan2_StridedData.hpp>

namespace Zoltan2{

// this utlity method is used to allocate more metrics in the metric array
// the array is composed of an array of ptrs to BaseClassMetric
// but each ptr is allocated to the specific derived metric type
// So the array can access the polymorphic hierarchy
//
// Note this is currently only relevant to EvaluatePartition and the
// GraphMetrics and ImbalanceMetrics calculations
template <typename metric_t, typename scalar_t>
RCP<metric_t> addNewMetric(const RCP<const Environment> &env,
  ArrayRCP<RCP<BaseClassMetrics<scalar_t> > > &metrics)
{
  metrics.resize(metrics.size() + 1); // increase array size by 1
  metric_t * newMetric = new metric_t;  // allocate
  env->localMemoryAssertion(__FILE__,__LINE__,1,newMetric); // check errors
  RCP<metric_t> newRCP = rcp(newMetric);       // rcp of derived class
  metrics[metrics.size()-1] = newRCP; 				 // create the new rcp
  return newRCP;
}

///////////////////////////////////////////////////////////////////
// Namespace methods to compute metric values
///////////////////////////////////////////////////////////////////

/*! \brief Find min, max and sum of metric values.
 *   \param v  a list of values
 *   \param stride  the value such that \c v[offset + stride*i]
 *             will be included in the calculation for all possible i.
 *   \param offset  the offset at which calculation will begin.
 *   \param min  on return, min will hold the minimum of the values.
 *   \param max  on return, max will hold the maximum of the values.
 *   \param sum on return, sum will hold the sum of the values.
 */
template <typename scalar_t>
 void getStridedStats(const ArrayView<scalar_t> &v, int stride,
   int offset, scalar_t &min, scalar_t &max, scalar_t &sum)
{
  if (v.size() < 1) return;
  min = max = sum = v[offset];

  for (int i=offset+stride; i < v.size(); i += stride){
    if (v[i] < min) min = v[i];
    else if (v[i] > max) max = v[i];
    sum += v[i];
  }
}

/*! \brief Find max and sum of graph metric values.
 *   \param v  a list of values
 *   \param stride  the value such that \c v[offset + stride*i]
 *             will be included in the calculation for all possible i.
 *   \param offset  the offset at which calculation will begin.
 *   \param max  on return, max will hold the maximum of the values.
 *   \param sum on return, sum will hold the sum of the values.
 */
template <typename scalar_t>
 void getStridedStats(const ArrayView<scalar_t> &v, int stride,
   int offset, scalar_t &max, scalar_t &sum)
{
  if (v.size() < 1) return;
  max = sum = v[offset];

  for (int i=offset+stride; i < v.size(); i += stride){
    if (v[i] > max) max = v[i];
    sum += v[i];
  }
}

/*! \brief Compute the total weight in each part on this process.
 *
 * \param env the Environment for error messages
 * \param numberOfParts the number of Parts with respect to
 *          which weights should be computed.
 * \param parts the part Id for each object, which may range
 *         from 0 to one less than \c numberOfParts
 * \param vwgts \c vwgts[w] is the StridedData object
 *    representing weight index
 *    \c w.  vwgts[w][i] is the \c w'th weight for object \c i.
 * \param mcNorm the multiCriteria norm, to be used if the number of weights
 *           is greater than one.
 * \param weights on return, \c weights[p] is the total weight for part
      \c p. \c weights is allocated by the caller
 *
 * \todo - Zoltan_norm() in Zoltan may scale the weight. Do we ever need this?
 */

template <typename scalar_t, typename lno_t, typename part_t>
  void normedPartWeights(
    const RCP<const Environment> &env,
    part_t numberOfParts,
    const ArrayView<const part_t> &parts,
    const ArrayView<StridedData<lno_t, scalar_t> > &vwgts,
    multiCriteriaNorm mcNorm,
    scalar_t *weights)
{
  env->localInputAssertion(__FILE__, __LINE__, "parts or weights",
    numberOfParts > 0 && vwgts.size() > 0, BASIC_ASSERTION);

  int numObjects = parts.size();
  int vwgtDim = vwgts.size();

  memset(weights, 0, sizeof(scalar_t) * numberOfParts);

  if (numObjects == 0)
    return;

  if (vwgtDim == 0) {
    for (lno_t i=0; i < numObjects; i++){
      weights[parts[i]]++;
    }
  }
  else if (vwgtDim == 1){
    for (lno_t i=0; i < numObjects; i++){
      weights[parts[i]] += vwgts[0][i];
    }
  }
  else{
    switch (mcNorm){
      case normMinimizeTotalWeight:   /*!< 1-norm = Manhattan norm */
        for (int wdim=0; wdim < vwgtDim; wdim++){
          for (lno_t i=0; i < numObjects; i++){
            weights[parts[i]] += vwgts[wdim][i];
          }
        }  // next weight index
        break;

      case normBalanceTotalMaximum:   /*!< 2-norm = sqrt of sum of squares */
        for (lno_t i=0; i < numObjects; i++){
          scalar_t ssum = 0;
          for (int wdim=0; wdim < vwgtDim; wdim++)
            ssum += (vwgts[wdim][i] * vwgts[wdim][i]);
          weights[parts[i]] += sqrt(ssum);
        }
        break;

      case normMinimizeMaximumWeight: /*!< inf-norm = maximum norm */
        for (lno_t i=0; i < numObjects; i++){
          scalar_t max = 0;
          for (int wdim=0; wdim < vwgtDim; wdim++)
            if (vwgts[wdim][i] > max)
              max = vwgts[wdim][i];
          weights[parts[i]] += max;
        }
        break;

      default:
        env->localBugAssertion(__FILE__, __LINE__, "invalid norm", false,
          BASIC_ASSERTION);
        break;
    }
  }
}

/*! \brief Compute the imbalance
 *  \param numExistingParts the max Part ID + 1, which is the
 *             length of the \c vals array.
 *  \param targetNumParts the number of parts desired, which is the
 *             length of the \c psizes array if it is defined.
 *  \param psizes  if part sizes are not uniform then <tt> psizes[p]</tt>
 *        is the part size for part \c p, for \c p ranging from zero
 *        to one less than \c targetNumParts.  Part sizes must sum to one.
 *        If part sizes are uniform, \c psizes should be NULL.
 *  \param sumVals is the sum of the values in the \c vals list.
 *  \param vals  <tt> vals[p] </tt> is the amount in part \c p, for \c p
 *          ranging from zero to one less than \c numExistingParts.
 *  \param min  on return, min will be the minimum (best)
 *           imbalance of all the parts.
 *  \param max  on return, max will be the maximum imbalance of all the parts.
 *  \param avg  on return avg will be the average imbalance across the parts.
 *
 * Imbalance is a value between zero and one.  If \c target is the desired
 * amount in part \c p and \c actual is the actual amount in part \c p, then
 * the imbalance is:
    \code
           abs(target - actual) / target
    \endcode
 *
 * If the part is supposed to be empty (\c target is zero), then no
 * imbalance is computed for that part.  If \c actual for that part
 * is non-zero, then other parts are too small and the imbalance will
 * be found in those other parts.
 */

template <typename scalar_t, typename part_t>
void computeImbalances(
  part_t numExistingParts,  // Max Part ID + 1
  part_t targetNumParts,    // comm.size() or requested global # parts from soln
  const scalar_t *psizes,
  scalar_t sumVals,
  const scalar_t *vals,
  scalar_t &min,
  scalar_t &max,
  scalar_t &avg)
{
  min = sumVals;
  max = avg = 0;

  if (sumVals <= 0 || targetNumParts < 1 || numExistingParts < 1)
    return;

  if (targetNumParts==1) {
    min = max = avg = 0;  // 0 imbalance
    return;
  }

  if (!psizes){
    scalar_t target = sumVals / targetNumParts;
    for (part_t p=0; p < numExistingParts; p++){
      scalar_t diff = vals[p] - target;
      scalar_t adiff = (diff >= 0 ? diff : -diff);
      scalar_t tmp = diff / target;
      scalar_t atmp = adiff / target;
      avg += atmp;
      if (tmp > max) max = tmp;
      if (tmp < min) min = tmp;
    }
    part_t emptyParts = targetNumParts - numExistingParts;
    if (emptyParts > 0){
      if (max < 1.0)
        max = 1.0;       // target divided by target
      avg += emptyParts;
    }
  }
  else{
    for (part_t p=0; p < targetNumParts; p++){
      if (psizes[p] > 0){
        if (p < numExistingParts){
          scalar_t target = sumVals * psizes[p];
          scalar_t diff = vals[p] - target;
          scalar_t adiff = (diff >= 0 ? diff : -diff);
          scalar_t tmp = diff / target;
          scalar_t atmp = adiff / target;
          avg += atmp;
          if (tmp > max) max = tmp;
          if (tmp < min) min = tmp;
        }
        else{
          if (max < 1.0)
            max = 1.0;  // target divided by target
          avg += 1.0;
        }
      }
    }
  }

  avg /= targetNumParts;
}

/*! \brief Compute the imbalance in the case of multiple part sizes.
 *
 *  \param numExistingParts the max Part ID + 1, which is the
 *             length of the \c vals array.
 *  \param targetNumParts the number of parts desired, which is the
 *             length of the \c psizes array if it is defined.
 *  \param numSizes the number of part size arrays
 *  \param psizes  is an array of \c numSizes pointers to part size arrays.
 *         If the part sizes for index \c w are uniform, then
 *         <tt>psizes[w]</tt> should be NULL.  Otherwise it should
 *         point to \c targetNumParts sizes, and the sizes for each
 *         index should sum to one.
 *  \param sumVals is the sum of the values in the \c vals list.
 *  \param vals  <tt> vals[p] </tt> is the amount in part \c p, for \c p
 *          ranging from zero to one less than \c numExistingParts.
 *  \param min  on return, min will be the minimum (best) imbalance
 *          of all the parts.
 *  \param max  on return, max will be the maximum imbalance of all the parts.
 *  \param avg  on return avg will be the average imbalance across the parts.
 *
 * Imbalance is a value between zero and one.  If \c target is the desired
 * amount in part \c p and \c actual is the actual amount in part \c p, then
 * the imbalance is:
    \code
           abs(target - actual) / target
    \endcode
 *
 * If the part is supposed to be empty (\c target is zero), then no
 * imbalance is computed for that part.  If \c actual for that part
 * is non-zero, then other parts are too small and the imbalance will
 * be found in those other parts.
 */

template <typename scalar_t, typename part_t>
void computeImbalances(
  part_t numExistingParts,
  part_t targetNumParts,
  int numSizes,
  ArrayView<ArrayRCP<scalar_t> > psizes,
  scalar_t sumVals,
  const scalar_t *vals,
  scalar_t &min,
  scalar_t &max,
  scalar_t &avg)
{
  min = sumVals;
  max = avg = 0;

  if (sumVals <= 0 || targetNumParts < 1 || numExistingParts < 1)
    return;

  if (targetNumParts==1) {
    min = max = avg = 0;  // 0 imbalance
    return;
  }

  bool allUniformParts = true;
  for (int i=0; i < numSizes; i++){
    if (psizes[i].size() != 0){
      allUniformParts = false;
      break;
    }
  }

  if (allUniformParts){
    computeImbalances<scalar_t, part_t>(numExistingParts, targetNumParts, NULL,
      sumVals, vals, min, max, avg);
    return;
  }

  double uniformSize = 1.0 / targetNumParts;
  std::vector<double> sizeVec(numSizes);
  for (int i=0; i < numSizes; i++){
    sizeVec[i] = uniformSize;
  }

  for (part_t p=0; p < numExistingParts; p++){

    // If we have objects in parts that should have 0 objects,
    // we don't compute an imbalance.  It means that other
    // parts have too few objects, and the imbalance will be
    // picked up there.

    if (p >= targetNumParts)
      break;

    // Vector of target amounts: T

    double targetNorm = 0;
    for (int i=0; i < numSizes; i++) {
      if (psizes[i].size() > 0)
        sizeVec[i] = psizes[i][p];
      sizeVec[i] *= sumVals;
      targetNorm += (sizeVec[i] * sizeVec[i]);
    }
    targetNorm = sqrt(targetNorm);

    // If part is supposed to be empty, we don't compute an
    // imbalance.  Same argument as above.

    if (targetNorm > 0){

      // Vector of actual amounts: A

      std::vector<double> actual(numSizes);
      double actualNorm = 0.;
      for (int i=0; i < numSizes; i++) {
        actual[i] = vals[p] * -1.0;
        actual[i] += sizeVec[i];
        actualNorm += (actual[i] * actual[i]);
      }
      actualNorm = sqrt(actualNorm);

      //  |A - T| / |T|

      scalar_t imbalance = actualNorm / targetNorm;

      if (imbalance < min)
        min = imbalance;
      else if (imbalance > max)
        max = imbalance;
      avg += imbalance;
    }
  }

  part_t numEmptyParts = 0;

  for (part_t p=numExistingParts; p < targetNumParts; p++){
    bool nonEmptyPart = false;
    for (int i=0; !nonEmptyPart && i < numSizes; i++)
      if (psizes[i].size() > 0 && psizes[i][p] > 0.0)
        nonEmptyPart = true;

    if (nonEmptyPart){
      // The partition has no objects for this part, which
      // is supposed to be non-empty.
      numEmptyParts++;
    }
  }

  if (numEmptyParts > 0){
    avg += numEmptyParts;
    if (max < 1.0)
      max = 1.0;       // target divided by target
  }

  avg /= targetNumParts;
}

/*! \brief Compute the norm of the vector of weights.
 */
template <typename scalar_t>
  scalar_t normedWeight(ArrayView <scalar_t> weights,
    multiCriteriaNorm norm)
{
  size_t dim = weights.size();
  if (dim == 0)
    return 0.0;
  else if (dim == 1)
    return weights[0];

  scalar_t nweight = 0;

  switch (norm){
    case normMinimizeTotalWeight:   /*!< 1-norm = Manhattan norm */

      for (size_t i=0; i <dim; i++)
        nweight += weights[i];

      break;

    case normBalanceTotalMaximum:   /*!< 2-norm = sqrt of sum of squares */
      for (size_t i=0; i <dim; i++)
        nweight += (weights[i] * weights[i]);

      nweight = sqrt(nweight);

      break;

    case normMinimizeMaximumWeight: /*!< inf-norm = maximum norm */

      nweight = weights[0];
      for (size_t i=1; i <dim; i++)
        if (weights[i] > nweight)
          nweight = weights[i];

      break;

    default:
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__ << std::endl;
      emsg << "bug: " << "invalid norm" << std::endl;
      throw std::logic_error(emsg.str());
  }

  return nweight;
}

/*! \brief Compute the norm of the vector of weights stored as StridedData.
 */
template<typename lno_t, typename scalar_t>
  scalar_t normedWeight(ArrayView<StridedData<lno_t, scalar_t> > weights,
     lno_t idx, multiCriteriaNorm norm)
{
  size_t dim = weights.size();
  if (dim < 1)
    return 0;

  Array<scalar_t> vec(dim, 1.0);
  for (size_t i=0; i < dim; i++)
    if (weights[i].size() > 0)
       vec[dim] = weights[i][idx];

  return normedWeight(vec, norm);
}

} //namespace Zoltan2

#endif
