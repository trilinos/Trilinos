// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Metric.hpp
 *  \brief Metric class and namespace methods to compute quality metrics.
 *  \todo Add graph and hypergraph metrics.
 */

#ifndef ZOLTAN2_METRIC_HPP
#define ZOLTAN2_METRIC_HPP

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_StridedData.hpp>

#include <Epetra_SerialDenseVector.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace Zoltan2{

///////////////////////////////////////////////////////////////////
// Enumerators
///////////////////////////////////////////////////////////////////

/*!\brief  Enumerator for handling of weight dimension > 1
 */
enum multiCriteriaNorm{
  normMinimizeTotalWeight,   /*!< 1-norm = Manhattan norm */
  normBalanceTotalMaximum,   /*!< 2-norm = sqrt of sum of squares */
  normMinimizeMaximumWeight, /*!< inf-norm = maximum norm */
  normNumNorms
};

///////////////////////////////////////////////////////////////////
// Classes
///////////////////////////////////////////////////////////////////

/*! \brief A class containing the metrics for one measurable item.
 */

template <typename scalar_t>
  class MetricValues{

private:
  void resetValues(){
    scalar_t *tmp = new scalar_t [evalNumMetrics];
    memset(tmp, 0, sizeof(scalar_t) * evalNumMetrics);
    values = arcp(tmp, 0, evalNumMetrics);
  }
  ArrayRCP<scalar_t> values;
  string metricName;

public:

/*! \brief  Enumerator for offsets into metric data.
 *
 *  When part sizes are all uniform, it is sufficient to 
 *  look at totals per part.  For non-uniform part sizes, the
 *  total is not really significant, but rather the min, max and
 *  average part imbalances.  We provide both types of metrics.
 */
enum metricOffset{
  evalLocalSum,    /*!< the total on this process */
  evalGlobalSum,   /*!< the global total on all parts */
  evalGlobalMin,   /*!< the minimum across all parts */
  evalGlobalMax,   /*!< the maximum across all parts */
  evalGlobalAvg,   /*!< the global sum divided by the number of parts */
  evalMinImbalance, /*!< the imbalance of best balanced part */
  evalAvgImbalance, /*!< the average of the part imbalances */
  evalMaxImbalance, /*!< the worst, which is the overall imbalance */
  evalNumMetrics    /*!< the number of metric values */
};

/*! \brief Print a standard header 
     \todo This method is static and the header is the same for
  every type of MetricValues object.  But when we add metrics for
  graph and hypergraph, we may need to look at which type of metric
  are set in order to decide what header to print.  So it would
  be a member method instead of static.
 */
static void printHeader(ostream &os);

/*! \brief Print a standard line of data that fits under the header. */
void printLine(ostream &os) const;

/*! \brief Constructor */
MetricValues(string mname) : metricName(mname) { resetValues();}

/*! \brief Constructor */
MetricValues() : metricName("unset") { resetValues();}

/*! \brief Set or reset the name.  */
void setName(string name) { metricName = name;}

/*! \brief Set the sum on the local process.  */
void setLocalSum(scalar_t x) { values[evalLocalSum] = x;}

/*! \brief Set the global sum.  */
void setGlobalSum(scalar_t x) { values[evalGlobalSum] = x;}

/*! \brief Set the global minimum across parts.  */
void setGlobalMin(scalar_t x) { values[evalGlobalMin] = x;}

/*! \brief Set the global maximum across parts.  */
void setGlobalMax(scalar_t x) { values[evalGlobalMax] = x;}

/*! \brief Set the global average (sum / numParts).  */
void setGlobalAvg(scalar_t x) { values[evalGlobalAvg] = x;}

/*! \brief Set the imbalance of the least imbalanced part. */
void setMinImbalance(scalar_t x) { values[evalMinImbalance] = x;}

/*! \brief Set the imbalance of the worst imbalanced part. 
     This is what we normally call the imbalance of a partition.
*/
void setMaxImbalance(scalar_t x) { values[evalMaxImbalance] = x;}

/*! \brief Set the average imbalance of all parts. */
void setAvgImbalance(scalar_t x) { values[evalAvgImbalance] = x;}

/*! \brief Get the name of the item measured. */
const string &getName() const { return metricName; }

/*! \brief Get the sum on the local process. */
scalar_t getLocalSum() const { return values[evalLocalSum];}

/*! \brief Get the global sum for all parts. */
scalar_t getGlobalSum() const { return values[evalGlobalSum];}

/*! \brief Get the global minimum across all parts. */
scalar_t getGlobalMin() const { return values[evalGlobalMin];}

/*! \brief Get the global maximum across all parts. */
scalar_t getGlobalMax() const { return values[evalGlobalMax];}

/*! \brief Get the average of the sum over all parts. */
scalar_t getGlobalAvg() const { return values[evalGlobalAvg];}

/*! \brief Get the imbalance of the least imbalanced part. */
scalar_t getMinImbalance() const { return values[evalMinImbalance];}

/*! \brief Get the imbalance of the most imbalanced part. 
     This is what we normally call the imbalance of a partition.
*/
scalar_t getMaxImbalance() const { return values[evalMaxImbalance];}

/*! \brief Get the average of the part imbalances. */
scalar_t getAvgImbalance() const { return values[evalAvgImbalance];}

};  // end class


template <typename scalar_t>
  void MetricValues<scalar_t>::printLine(ostream &os) const
{
  os << setw(20) << metricName;
  os << setw(12) << setprecision(4) << scientific << values[evalGlobalMin];
  os << setw(12) << setprecision(4) << scientific << values[evalGlobalAvg];
  os << setw(12) << setprecision(4) << scientific << values[evalGlobalMax];
  os << setw(6) << setprecision(4) << values[evalMinImbalance];
  os << setw(6) << setprecision(4) << values[evalAvgImbalance];
  os << setw(6) << setprecision(4) << values[evalMaxImbalance];
  os << endl;
}

template <typename scalar_t>
  void MetricValues<scalar_t>::printHeader(ostream &os)
{
  os << setw(20) << " ";
  os << setw(36) << "SUM PER PART";
  os << setw(18) << "IMBALANCE PER PART";
  os << endl;

  os << setw(20) << " ";
  os << setw(12) << "min" << "max" << "avg";
  os << setw(6) << "min" << "max" << "avg";
  os << endl;
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
 *   \param sum on return, max will hold the sum of the values.
 */
template <typename scalar_t>
 void getStridedStats(const ArrayView<scalar_t> &v, int stride, int offset,
                      scalar_t &min, scalar_t &max, scalar_t &sum)
{
  if (v.size() < 1) return;
  min = max = sum = v[offset];

  for (int i=offset+stride; i < v.size(); i += stride){
    if (v[i] < min) min = v[i];
    else if (v[i] > max) max = v[i];
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
 *   \param partNumMin   If only a range of part IDs are to be counted,
 *          then \c partNumMin is the lowest part ID to count.  
 *          If \c partNumMin
 *          and \c partNumMax are the same, then all part IDs present in the
 *          \c part list are counted.
 *   \param partNumMax   If only a range of part IDs are to be counted,
 *          then \c partNumMax is the highest part ID to count.
 *          If \c partNumMin
 *          and \c partNumMax are the same, then all part IDs present in the
 *          \c part list are counted.
 * \param vwgts \c vwgts[w] is the StridedData object 
 *    representing weight dimension
 *    \c w.  vwgts[w][i] is the \c w'th weight for object \c i.
 * \param mcNorm the multiCriteria norm, to be used if weight dimension is 
 *             greater than one.
 * \param weights on return, \c weights[p] is the total weight for part \c p. 
 *                    \c weights is allocated by the caller
 *
 * \todo - Zoltan_norm() in Zoltan may scale the weight. Do we ever need this?
 */

template <typename scalar_t, typename pnum_t, typename lno_t>
  void normedPartWeights(
    const RCP<const Environment> &env,
    int numberOfParts,
    const ArrayView<pnum_t> &parts,
    pnum_t partNumMin,
    pnum_t partNumMax,
    const ArrayView<StridedData<lno_t, scalar_t> > &vwgts,
    multiCriteriaNorm mcNorm,
    scalar_t *weights)
{
  env->localInputAssertion(__FILE__, __LINE__, "parts or weights", 
    numberOfParts > 0 && vwgts.size() > 0, BASIC_ASSERTION);

  bool checkNum = (partNumMin < partNumMax);
  int numObjects = parts.size();
  int vwgtDim = vwgts.size();
  scalar_t *partWeights = new scalar_t [numberOfParts];
  env->localMemoryAssertion(__FILE__, __LINE__, numberOfParts, partWeights);
  memset(partWeights, 0, sizeof(scalar_t) * numberOfParts);

  memset(weights, 0, sizeof(scalar_t) * numberOfParts);

  if (numObjects == 0)
    return;

  bool haveUniform = false;
  bool haveNonUniform = false;

  for (int wdim=0; wdim < vwgtDim; wdim++){
    if (vwgts[wdim].size() > 0)
      haveNonUniform = true;
    if (vwgts[wdim].size() == 0)
      haveUniform = true;
  }

  if (!haveNonUniform){
    for (lno_t i=0; i < numObjects; i++){
      if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
        continue;
      weights[parts[i]]++;
    }
  }
  else if (vwgtDim == 1){
    for (lno_t i=0; i < numObjects; i++){
      if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
        continue;
      weights[parts[i]] += vwgts[0][i];
    }
  }
  else{
    switch (mcNorm){
      case normMinimizeTotalWeight:   /*!< 1-norm = Manhattan norm */

        for (int wdim=0; wdim < vwgtDim; wdim++){
          if (vwgts[wdim].size() == 0){
            for (lno_t i=0; i < numObjects; i++){
              if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
                continue;
              weights[parts[i]]++;
            }
          }
          else{
            for (lno_t i=0; i < numObjects; i++){
              if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
                continue;
                weights[parts[i]] += vwgts[wdim][i];
            }
          }
        }  // next weight dimension
       
      case normBalanceTotalMaximum:   /*!< 2-norm = sqrt of sum of squares */
        if (!haveUniform){
          for (lno_t i=0; i < numObjects; i++){
            if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
              continue;
            scalar_t ssum = 0;
            for (int wdim=0; wdim < vwgtDim; wdim++)
              ssum += (vwgts[wdim][i] * vwgts[wdim][i]);
            weights[parts[i]] += sqrt(ssum);
          }
        }
        else{
          scalar_t uniformFactor = 0;
          for (int wdim=0; wdim < vwgtDim; wdim++)
            if (vwgts[wdim].size() == 0)
              uniformFactor++;
              
          for (lno_t i=0; i < numObjects; i++){
            if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
              continue;
            scalar_t ssum = 0;
            for (int wdim=0; wdim < vwgtDim; wdim++){
              if (vwgts[wdim].size() > 0)
                ssum += (vwgts[wdim][i] * vwgts[wdim][i]);
            }
            ssum += uniformFactor;
            weights[parts[i]] += sqrt(ssum);
          }
        }

      case normMinimizeMaximumWeight: /*!< inf-norm = maximum norm */

        if (!haveUniform){

          for (lno_t i=0; i < numObjects; i++){
            if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
              continue;
            scalar_t max = 0;
            for (int wdim=0; wdim < vwgtDim; wdim++)
              if (vwgts[wdim][i] > max)
                max = vwgts[wdim][i];
            weights[parts[i]] += max;
          }
        }
        else{

          for (lno_t i=0; i < numObjects; i++){
            if (checkNum && (parts[i] < partNumMin || parts[i] > partNumMax))
              continue;
            scalar_t max = 1.0;
            for (int wdim=0; wdim < vwgtDim; wdim++)
              if (vwgts[wdim].size() > 0 && vwgts[wdim][i] > max)
                max = vwgts[wdim][i];
            weights[parts[i]] += max;
          }
        }

      default:
        env->localBugAssertion(__FILE__, __LINE__, "invalid norm", false,
          BASIC_ASSERTION);
    } 
  } 
}

/*! \brief Given the local partitioning, compute the global sums in each part.
 *
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param part   \c part[i] is the part ID for local object \c i
 *   \param partNumMin   If only a range of part IDs are to be counted,
 *          then \c partNumMin is the lowest part ID to count.  
 *          If \c partNumMin
 *          and \c partNumMax are the same, then all part IDs present in the
 *          \c part list are counted.
 *   \param partNumMax   If only a range of part IDs are to be counted,
 *          then \c partNumMax is the highest part ID to count.
 *          If \c partNumMin
 *          and \c partNumMax are the same, then all part IDs present in the
 *          \c part list are counted.
 *   \param vwgts  \c vwgts[w] is the StridedData object 
 *       representing weight dimension \c w. The weight dimension
 *       (which must be at least one) is taken to be \c vwgts.size().  
 *       If <tt>vwgts[wdim].size()</tt> is
 *       zero, we assume uniform weights for weight dimension \c wdim.
 *   \param mcNorm the multiCriteria norm, to be used if weight dimension is
 *             greater than one.
 *   \param numParts  on return this is the total number of parts globally.
 *   \param numNonemptyParts  on return this is the number of 
 *          parts that are non-empty.
 *   \param metrics on return points to a list of named MetricValues objects 
 *     that each contains the global min, max and avg over parts of 
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 1", "weight 2" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one dimension of non-uniform weights were given, then
 *     "object count" and "weight 1" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 1", "weight 2", and so on.
 *   \param globalSums If weights are uniform, the globalSums is the
 *      \c numParts totals of global number of objects in each part.
 *     Suppose the weight dimension is \c W.  If
 *     W is 1, then on return this is an array of length \c 2*numParts .
 *     The first \c numParts entries are the count of objects in each part,  
 *     and the second is the total weight in each part.
 *     If \c W is greater than one, then the length of this array is 
 *     \c (2+W)*numParts .
 *     The first \c numParts entries are the count of objects in each part.  
 *     The next \c numParts entries are the sum of the normed weights in 
 *     each part.  
 *     The final entries are the sum of the individual weights in each part,
 *     by weight dimension by part number.  The array is allocated here.
 *
 * globalSumsByPart() must be called by all processes in \c comm.
 * The imbalance metrics are not yet set in the MetricValues objects, 
 * because they require part size information.
 */

template <typename scalar_t, typename pnum_t, typename lno_t>
  void globalSumsByPart( 
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm, 
    const ArrayView<pnum_t> &part, 
    pnum_t partNumMin,
    pnum_t partNumMax,
    const ArrayView<StridedData<lno_t, scalar_t> > &vwgts,
    multiCriteriaNorm mcNorm,
    partId_t &numParts, 
    partId_t &numNonemptyParts,
    ArrayRCP<MetricValues<scalar_t> > &metrics,
    ArrayRCP<scalar_t> &globalSums)
{
  //////////////////////////////////////////////////////////
  // Initialize return values

  numParts = numNonemptyParts = 0;

  int vwgtDim = vwgts.size();

  int numMetrics = 1;        // "object count"
  if (vwgts[0].size() > 0)
    numMetrics++;            // "normed weight" or "weight 1"
  if (vwgtDim > 1)
    numMetrics += vwgtDim;   // "weight n"

  metrics = arcp(new MetricValues<scalar_t> [numMetrics], 0, numMetrics);

  bool checkNum = (partNumMin < partNumMax);

  //////////////////////////////////////////////////////////
  // Figure out the global number of parts in use.
  // Verify vertex weight dim is the same everywhere.

  lno_t localNumObj = part.size();
  partId_t localNum[2], globalNum[2];
  localNum[0] = static_cast<partId_t>(vwgtDim);  
  localNum[1] = 0;

  if (!checkNum){
    for (lno_t i=0; i < localNumObj; i++)
      if (part[i] > localNum[1]) localNum[1] = part[i];
  }
  else{
    for (lno_t i=0; i < localNumObj; i++){
      if (part[i] < partNumMin || part[i] > partNumMax)
        continue;
      if (part[i] > localNum[1]) 
        localNum[1] = part[i];
    }
  }

  try{
    reduceAll<int, partId_t>(*comm, Teuchos::REDUCE_MAX, 2, 
      localNum, globalNum);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->globalBugAssertion(__FILE__, __LINE__, "inconsistent vertex dimension",
    globalNum[0] > 0 && globalNum[0] == localNum[0], 
    DEBUG_MODE_ASSERTION, comm);

  partId_t nparts = globalNum[1] + 1;

  int globalSumSize = nparts * numMetrics;
  scalar_t * sumBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, sumBuf);
  globalSums = arcp(sumBuf, 0, globalSumSize);

  //////////////////////////////////////////////////////////
  // Calculate the local totals by part.

  scalar_t *localBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, localBuf);
  memset(localBuf, 0, sizeof(scalar_t) * globalSumSize);

  scalar_t *obj = localBuf;              // # of objects

  if (!checkNum){
    for (lno_t i=0; i < localNumObj; i++)
      obj[part[i]]++;
  }
  else{
    for (lno_t i=0; i < localNumObj; i++)
      if (part[i] >= partNumMin && part[i] <= partNumMax)
        obj[part[i]]++;
  }

  if (numMetrics > 1){

    scalar_t *wgt = localBuf + nparts; // single normed weight
    try{
      normedPartWeights<scalar_t, pnum_t, lno_t>(env, nparts, 
        part, partNumMin, partNumMax, vwgts, mcNorm, wgt);
    }
    Z2_FORWARD_EXCEPTIONS
  
    if (vwgtDim > 1){
      wgt += nparts;         // individual weights
      for (int vdim = 0; vdim < vwgtDim; vdim++){
        if (vwgts[vdim].size()){
          if (!checkNum){
            for (lno_t i=0; i < localNumObj; i++)
              wgt[part[i]] += vwgts[vdim][i];
          }
          else{
            for (lno_t i=0; i < localNumObj; i++)
              if (part[i] >= partNumMin && part[i] <= partNumMax)
                wgt[part[i]] += vwgts[vdim][i];
          }
        }
        else{  // uniform weights
          for (int p=0; p < nparts; p++)
            wgt[p] = obj[p];
        }
        wgt += nparts;
      }
    }
  }

  // Metric: local sums on process

  int next = 0;

  metrics[next].setName("object count");
  metrics[next].setLocalSum(localNumObj);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = localBuf + nparts; // single normed weight
    scalar_t total = 0.0;
  
    for (int p=0; p < nparts; p++){
      total += wgt[p];
    }

    if (vwgtDim == 1)
      metrics[next].setName("weight 1");
    else
      metrics[next].setName("normed weight");

    metrics[next].setLocalSum(total);
    next++;
  
    if (vwgtDim > 1){
      for (int vdim = 0; vdim < vwgtDim; vdim++){
        wgt += nparts;
        total = 0.0;
        for (int p=0; p < nparts; p++){
          total += wgt[p];
        }

        ostringstream oss;
        oss << "weight " << vdim+1;

        metrics[next].setName(oss.str());
        metrics[next].setLocalSum(total);
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
  scalar_t min, max, sum;
  next = 0;

  ArrayView<scalar_t> objVec(obj, nparts);
  getStridedStats<scalar_t>(objVec, 1, 0, min, max, sum);

  metrics[next].setGlobalMin(min);
  metrics[next].setGlobalMax(max);
  metrics[next].setGlobalSum(sum);
  metrics[next].setGlobalAvg(sum / nparts);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + nparts;        // single normed weight
  
    ArrayView<scalar_t> normedWVec(wgt, nparts);
    getStridedStats<scalar_t>(normedWVec, 1, 0, min, max, sum);

    metrics[next].setGlobalMin(min);
    metrics[next].setGlobalMax(max);
    metrics[next].setGlobalSum(sum);
    metrics[next].setGlobalAvg(sum / nparts);
    next++;
  
    if (vwgtDim > 1){
      for (int vdim=0; vdim < vwgtDim; vdim++){
        wgt += nparts;       // individual weights
        ArrayView<scalar_t> fromVec(wgt, nparts);
        getStridedStats<scalar_t>(fromVec, 1, 0, min, max, sum);

        metrics[next].setGlobalMin(min);
        metrics[next].setGlobalMax(max);
        metrics[next].setGlobalSum(sum);
        metrics[next].setGlobalAvg(sum / nparts);
        next++;
      }
    }
  }

  //////////////////////////////////////////////////////////
  // How many parts do we actually have.

  numParts = nparts;
  obj = sumBuf;               // # of objects

  for (partId_t p=nparts-1; p > 0; p--){
    if (obj[p] > 0) break;
    numParts--;
  }

  numNonemptyParts = numParts; 

  for (partId_t p=0; p < numParts; p++)
    if (obj[p] == 0) numNonemptyParts--;

  delete [] sumBuf;
}

/*! \brief Compute the imbalance
 *  \param numParts the number of parts supplied, which is the
 *             length of the \c vals array.
 *  \param targetNumParts the number of parts desired, which is the
 *             length of the \c psizes array if it is defined.
 *  \param psizes  if part sizes are not uniform then <tt> psizes[p]</tt>
 *        is the part size for part \c p, for \c p ranging from zero
 *        to one less than \c targetNumParts.  Part sizes must sum to one.
 *        If part sizes are uniform, \c psizes should be NULL.
 *  \param sumVals is the sum of the values in the \c vals list.
 *  \param vals  <tt> vals[p] </tt> is the amount in part \c p, for \c p
 *          ranging from zero to one less than \c numParts.
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

template <typename scalar_t>
  void computeImbalances(partId_t numParts, partId_t targetNumParts,
    const scalar_t *psizes, scalar_t sumVals , const scalar_t *vals, 
    scalar_t &min, scalar_t &max, scalar_t &avg)
{
  min = sumVals;
  max = avg = 0;

  if (sumVals <= 0 || targetNumParts < 1 || numParts < 1)
    return;

  if (!psizes){
    scalar_t target = sumVals / targetNumParts;
    for (partId_t p=0; p < numParts; p++){
      scalar_t diff = abs(vals[p] - target);
      scalar_t tmp = diff / target;
      avg += tmp;
      if (tmp > max) max = tmp;
      if (tmp < min) min = tmp;
    }
    partId_t emptyParts = targetNumParts - numParts;  
    if (emptyParts > 0){
      if (max < 1.0)
        max = 1.0;       // target divided by target
      avg += emptyParts;
    }
  }
  else{
    for (partId_t p=0; p < targetNumParts; p++){
      if (psizes[p] > 0){
        if (p < numParts){
          scalar_t target = sumVals * psizes[p];
          scalar_t diff = abs(vals[p] - target);
          scalar_t tmp = diff / target;
          avg += tmp;
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
 *  \param numParts the number of parts supplied, which is the
 *             length of the \c vals array.
 *  \param targetNumParts the number of parts desired, which is the
 *             length of the \c psizes array if it is defined.
 *  \param numSizes the number of part size arrays
 *  \param psizes  is an array of \c numSizes pointers to part size arrays.
 *         If the part sizes for dimension \c w are uniform, then
 *         <tt>psizes[w]</tt> should be NULL.  Otherwise it should
 *         point to \c targetNumParts sizes, and the sizes for each
 *         dimension should sum to one.
 *  \param sumVals is the sum of the values in the \c vals list.
 *  \param vals  <tt> vals[p] </tt> is the amount in part \c p, for \c p
 *          ranging from zero to one less than \c numParts.
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

template <typename scalar_t>
 void computeImbalances(partId_t numParts, partId_t targetNumParts,
   int numSizes, const scalar_t * const *psizes, 
   scalar_t sumVals , const scalar_t *vals, 
   scalar_t &min, scalar_t &max, scalar_t &avg)
{
  min = sumVals;
  max = avg = 0;

  if (sumVals <= 0 || targetNumParts < 1 || numParts < 1)
    return;

  bool allUniformParts = true;
  for (int i=0; i < numSizes; i++){
    if (psizes[i] != NULL){
      allUniformParts = false;
      break;
    }
  }

  if (allUniformParts){
    computeImbalances<scalar_t>(numParts, targetNumParts, NULL,
      sumVals, vals, min, max, avg);
    return;
  }

  scalar_t uniformSize = 1.0 / targetNumParts;
  Array<double> sizeVec(numSizes, uniformSize);

  for (partId_t p=0; p < numParts; p++){

    // If we have objects in parts that should have 0 objects,
    // we don't compute an imbalance.  It means that other
    // parts have too few objects, and the imbalance will be
    // picked up there.

    if (p >= targetNumParts)
      break;

    // Vector of target amounts: T

    for (int i=0; i < numSizes; i++)
      if (psizes[i])
        sizeVec[i] = psizes[i][p];

    Epetra_SerialDenseVector target(View, sizeVec.getRawPtr(), numSizes);
    target.Scale(sumVals);
    double targetNorm = target.Norm2();

    // If part is supposed to be empty, we don't compute an
    // imbalance.  Same argument as above.

    if (targetNorm > 0){

      // Vector of actual amounts: A

      Epetra_SerialDenseVector actual(numSizes);
      for (int i=0; i < numSizes; i++)
        actual[i] = vals[p] * -1.0;
      
      actual += target;

      //  |A - T| / |T|

      double imbalance = actual.Norm2() / targetNorm;

      if (imbalance < min)
        min = imbalance;
      else if (imbalance > max)
        max = imbalance;
      avg += imbalance; 
    }
  }

  partId_t numEmptyParts = 0;

  for (partId_t p=numParts; p < targetNumParts; p++){
    bool nonEmptyPart = false;
    for (int i=0; !nonEmptyPart && i < numSizes; i++)
      if (psizes[i] && psizes[i][p] > 0.0)
        nonEmptyPart = true;

    if (nonEmptyPart){
      // The partitioning has no objects for this part, which
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

/*! \brief Compute imbalance metrics for a distribution.
 *
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param targetNumParts is usually the same as \c numParts, which is the
 *         global number of different parts to which objects have been 
 *         assigned in the \c parts array.  
 *         However if we are calling objectMetrics() before
 *         and after partitioning, it is possible that before partitioning
 *         the number of parts which have objects is not the same as
 *         targetNumParts.  For example, all objects may be on process zero.
 *         We calculate the imbalance with respect to \c targetNumParts 
 *         and not \c numParts.
 *   \param partSizes \c partSizes[w][p] is the desired part size for 
 *         weight dimension \w for part \c p. For each weight there should
 *         be \c targetNumParts part sizes. However if \c partSizes[w].size() 
 *         is zero we assume uniform part sizes. Weight dimension must be at
 *         least one.  If part sizes are supplied, they must sum
 *         to 1.0 for any weight dimension.
 *   \param parts  \c parts[i] is the part assigned to object \c i.
 *   \param vwgts  \c vwgts[w][i] is the weight for weight dimension \c w for
 *                      object \c i. If \c vwgts[w].size() is zero we
 *                     assume uniform weights. Weight dimension must be at
 *                     least one.
 *   \param mcNorm  is the multicriteria norm to use if the weight dimension
 *           is greater than one.  See the multiCriteriaNorm enumerator for
 *           \c mcNorm values.
 *   \param numParts on return is the global number of parts
 *   \param numNonemptyParts on return is the global number of parts to which 
 *                                objects have been assigned.
 *   \param metrics on return points to a list of named MetricValues objects 
 *     that each contains the global min, max and avg over parts and
 *     also imbalance measures of 
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 1", "weight 2" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one dimension of non-uniform weights were given, then
 *     "object count" and "weight 1" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 1", "weight 2", and so on.
 *
 *  objectMetrics() must be called by all processes in \c comm.  
 *  See the metricOffset enumerator in the MetricValues class for the 
 *  interpretation of the metric quantities.
 *   \todo check that part sizes sum to one if we're doing COMPLEX_ASSERTION
 */

template <typename scalar_t, typename lno_t>
  void objectMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    partId_t targetNumParts, 
    const ArrayView<ArrayView<scalar_t> > &partSizes,
    const ArrayView<partId_t> &parts,
    const ArrayView<StridedData<lno_t, scalar_t> > &vwgts,
    multiCriteriaNorm mcNorm,
    partId_t &numParts,
    partId_t &numNonemptyParts,
    ArrayRCP<MetricValues<scalar_t> > &metrics)
{
  ///////////////////////////////////////////////////////////////////////////
  // Get number of parts, and the number that are non-empty.
  // Get sums per part of objects, individual weights, and normed weight sums.

  ArrayRCP<scalar_t> globalSums;

  try{
    globalSumsByPart<scalar_t, partId_t, lno_t>(env, comm, 
      parts, 0, 0, vwgts, mcNorm,
      numParts, numNonemptyParts, metrics, globalSums);
  }
  Z2_FORWARD_EXCEPTIONS

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the object count.  
  // (Use first dimension of part sizes.) 

  scalar_t *objCount  = globalSums.getRawPtr();
  scalar_t *psizes = NULL;
  if (partSizes[0].size() > 0)
    psizes = partSizes[0].getRawPtr();

  scalar_t min, max, avg;

  computeImbalances<scalar_t>(numParts, targetNumParts, psizes, 
      metrics[0].getGlobalSum(), objCount, 
      min, max, avg);

  metrics[0].setMinImbalance(min);
  metrics[0].setMaxImbalance(max);
  metrics[0].setAvgImbalance(avg);

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the normed weight sum.

  int vwgtDim = vwgts.size();
  scalar_t *wgts = globalSums.getRawPtr() + numParts;

  if (metrics.size() > 1){
    Array<scalar_t *>psizesMultiple(vwgtDim);
  
    for (int vdim=0; vdim < vwgtDim; vdim++){
      if (partSizes[vdim].size() > 0)
        psizesMultiple[vdim] = partSizes[vdim].getRawPtr();
      else
        psizesMultiple[vdim] = NULL;
    }
  
    computeImbalances<scalar_t>(numParts, targetNumParts, 
      vwgtDim, psizesMultiple.getRawPtr(),
      metrics[1].getGlobalSum(), wgts,
      min, max, avg);
  
    metrics[1].setMinImbalance(min);
    metrics[1].setMaxImbalance(max);
    metrics[1].setAvgImbalance(avg);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for each individual weight dimension.

  if (metrics.size() > 2){
    int next = 2;

    for (int vdim=0; vdim < vwgtDim; vdim++){
      wgts += numParts;
      psizes = NULL;
      if (partSizes[vdim].size() > 0)
        psizes = partSizes[vdim].getRawPtr();

      computeImbalances<scalar_t>(numParts, targetNumParts, psizes,
        metrics[next].getGlobalSum(), wgts,
        min, max, avg);

      metrics[next].setMinImbalance(min);
      metrics[next].setMaxImbalance(max);
      metrics[next].setAvgImbalance(avg);
      next++;
    }
  }
}


/*! \brief Compute object metrics when weight dimension is one and 
 *           there are no part sizes.
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param targetNumParts is usually the same as \c numParts, which is the
 *         global number of different parts to which objects have been
 *         assigned in the \c parts array.
 *         However if we are calling objectMetrics() before
 *         and after partitioning, it is possible that before partitioning
 *         the number of parts which have objects is not the same as
 *         targetNumParts.  For example, all objects may be on process zero.
 *         We calculate the imbalance with respect to \c targetNumParts
 *         and not \c numParts.
 *   \param parts  \c parts[i] is the part assigned to object \c i.
 *   \param vwgts  is the StridedData object containing the weights.  
 *            If \c vwgts.size()
 *                    is zero, we assume uniform weights.
 *   \param numParts on return is the global number of parts
 *   \param numNonemptyParts on return is the global number of parts to which
 *                                objects have been assigned.
 *   \param metrics on return points to a list of named MetricValues objects 
 *     that each contains the global min, max and avg over parts and
 *     also imbalance measures of 
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 1", "weight 2" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one dimension of non-uniform weights were given, then
 *     "object count" and "weight 1" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 1", "weight 2", and so on.
 *
 */

template <typename scalar_t, typename lno_t>
  void objectMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    partId_t targetNumParts, 
    const ArrayView<partId_t> &parts,
    const StridedData<lno_t, scalar_t>  &vwgts,
    partId_t &numParts,
    partId_t &numNonemptyParts,
    ArrayRCP<MetricValues<scalar_t> > &metrics)  
{
  Array<ArrayView<scalar_t> > psArray(1);
  Array<StridedData<lno_t, scalar_t> > wgts(1, vwgts);

  return objectMetrics(env, comm, targetNumParts,
    psArray.view(0,1), parts, wgts.view(0, 1), normMinimizeTotalWeight,
    numParts, numNonemptyParts, metrics);
}

/*! \brief Compute object metrics when weight dimension is one and we
 *                 have non-uniform part sizes.
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param targetNumParts is usually the same as \c numParts, which is the
 *         global number of different parts to which objects have been
 *         assigned in the \c parts array.
 *         However if we are calling objectMetrics() before
 *         and after partitioning, it is possible that before partitioning
 *         the number of parts which have objects is not the same as
 *         targetNumParts.  For example, all objects may be on process zero.
 *         We calculate the imbalance with respect to \c targetNumParts
 *         and not \c numParts.
 *   \param partSizes \c partSizes[p] is the desired part size for
 *         part \c p.  However if \c partSizes.size()
 *         is zero we assume uniform part sizes.
 *   \param parts  \c parts[i] is the part assigned to object \c i.
 *   \param vwgts  is the StridedData object containing the object weights. 
 *                 If \c vwgts.size() is zero we assume uniform weights.
 *   \param mcNorm  is the multicriteria norm to use if the weight dimension
 *           is greater than one.  See the multiCriteriaNorm enumerator for
 *           \c mcNorm values.
 *   \param numParts on return is the global number of parts
 *   \param numNonemptyParts on return is the global number of parts to which
 *                                objects have been assigned.
 *   \param metrics on return points to a list of named MetricValues objects 
 *     that each contains the global min, max and avg over parts and
 *     also imbalance measures of 
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 1", "weight 2" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one dimension of non-uniform weights were given, then
 *     "object count" and "weight 1" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 1", "weight 2", and so on.
 *
 */

template <typename scalar_t, typename lno_t>
  void objectMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    partId_t targetNumParts, 
    const ArrayView<scalar_t> &partSizes,
    const ArrayView<partId_t> &parts,
    const StridedData<lno_t, scalar_t>  &vwgts,
    partId_t &numParts,
    partId_t &numNonemptyParts,
    ArrayRCP<MetricValues<scalar_t> > &metrics)
{
  Array<ArrayView<scalar_t> > psArray(1, partSizes);
  Array<StridedData<lno_t, scalar_t> > wgts(1, vwgts);

  return objectMetrics(env, comm, targetNumParts,
    psArray.view(0,1), parts, wgts.view(0, 1), normMinimizeTotalWeight,
    numParts, numNonemptyParts, metrics);
  
}

/*! \brief Compute object metrics when there are no part sizes or weights.
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param targetNumParts is usually the same as \c numParts, which is the
 *         global number of different parts to which objects have been
 *         assigned in the \c parts array.
 *         However if we are calling objectMetrics() before
 *         and after partitioning, it is possible that before partitioning
 *         the number of parts which have objects is not the same as
 *         targetNumParts.  For example, all objects may be on process zero.
 *         We calculate the imbalance with respect to \c targetNumParts
 *         and not \c numParts.
 *   \param parts  \c parts[i] is the part assigned to object \c i.
 *   \param numParts on return is the global number of parts
 *   \param numNonemptyParts on return is the global number of parts to which
 *                                objects have been assigned.
 *   \param metrics on return points to a list of named MetricValues objects 
 *     that each contains the global min, max and avg over parts and
 *     also imbalance measures of 
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 1", "weight 2" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one dimension of non-uniform weights were given, then
 *     "object count" and "weight 1" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 1", "weight 2", and so on.
 *
 */
template <typename scalar_t, typename lno_t>
  void objectMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    partId_t targetNumParts, 
    const ArrayView<partId_t> &parts,
    partId_t &numParts,
    partId_t &numNonemptyParts,
    ArrayRCP<MetricValues<scalar_t> > &metrics)
{
  Array<ArrayView<scalar_t> > noPartSizes(1);
  Array<StridedData<lno_t, scalar_t> > noWeights(1);

  return objectMetrics(env, comm, targetNumParts,
    noPartSizes.view(0,1), parts, noWeights.view(0, 1), 
    normMinimizeTotalWeight,
    numParts, numNonemptyParts, metrics);
}

/*! \brief Compute object metrics when there are no part sizes but there are
 *                 multiple weights.
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param targetNumParts is usually the same as \c numParts, which is the
 *         global number of different parts to which objects have been
 *         assigned in the \c parts array.
 *         However if we are calling objectMetrics() before
 *         and after partitioning, it is possible that before partitioning
 *         the number of parts which have objects is not the same as
 *         targetNumParts.  For example, all objects may be on process zero.
 *         We calculate the imbalance with respect to \c targetNumParts
 *         and not \c numParts.
 *   \param parts  \c parts[i] is the part assigned to object \c i.
 *   \param vwgts  \c vwgts[w][i] is the weight for weight dimension \c w for
 *                      object \c i. If \c vwgts[w].size() is zero we
 *                     assume uniform weights. Weight dimension must be at
 *                     least one.
 *   \param mcNorm  is the multicriteria norm to use if the weight dimension
 *           is greater than one.  See the multiCriteriaNorm enumerator for
 *           \c mcNorm values.
 *   \param numParts on return is the global number of parts
 *   \param numNonemptyParts on return is the global number of parts to which
 *                                objects have been assigned.
 *   \param metrics on return points to a list of named MetricValues objects 
 *     that each contains the global min, max and avg over parts and
 *     also imbalance measures of 
 *     the item being measured. The list may contain "object count",
 *     "normed weight", "weight 1", "weight 2" and so on in that order.
 *     If uniform weights were given, then only "object count" appears.
 *     If one dimension of non-uniform weights were given, then
 *     "object count" and "weight 1" appear.  Finally, if multiple
 *     weights were given, we have "object count", then "normed weight",
 *     then the individual weights "weight 1", "weight 2", and so on.
 *
 *
 *  See the metricOffset enumerator for the interpretation of the
 *  evalNumMetrics quantities.
 */
template <typename scalar_t, typename lno_t>
    void objectMetrics(
      const RCP<const Environment> &env,
      const RCP<const Comm<int> > &comm,
      partId_t targetNumParts, 
      const ArrayView<partId_t> &parts,
      const ArrayView<StridedData<lno_t, scalar_t> > &vwgts,
      multiCriteriaNorm mcNorm,
      partId_t &numParts,
      partId_t &numNonemptyParts,
      ArrayRCP<MetricValues<scalar_t> > &metrics)
{
  int vwgtDim = vwgts.size();
  Array<ArrayView<scalar_t> > noPartSizes(vwgtDim);

  return objectMetrics(env, comm, targetNumParts,
    noPartSizes.view(0,vwgtDim), parts, vwgts, mcNorm,
    numParts, numNonemptyParts, metrics);
}

/*! \brief Print out a header and the list metric values.
 */

template <typename scalar_t>
  void printMetrics( ostream &os,
    partId_t targetNumParts, partId_t numParts, partId_t numNonemptyParts, 
    const ArrayView<MetricValues<scalar_t> > &infoList)
{
  os << "NUMBER OF PARTS IS " << numParts;
  if (numNonemptyParts < numParts)
    os << " (" << numNonemptyParts << " of which are non-empty)";
  os << endl;
  if (targetNumParts != numParts)
    os << "TARGET NUMBER OF PARTS IS " << targetNumParts << std::endl;

  string unset("unset");

  MetricValues<scalar_t>::printHeader(os);

  for (int i=0; i < infoList.size(); i++)
    if (infoList[i].getName() != unset)
      infoList[i].printLine(os);

  os << endl;
}

/*! \brief Print out a header and the metric values.
 */

template <typename scalar_t>
  void printMetrics( ostream &os,
    partId_t targetNumParts, partId_t numParts, partId_t numNonemptyParts, 
    const MetricValues<scalar_t> &info)
{
  ArrayView<MetricValues<scalar_t> > infoList(&info, 1);
  printMetrics( os, targetNumParts, numParts, numNonemptyParts, infoList);
}





} //namespace Zoltan2
#endif
