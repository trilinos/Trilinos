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

/*! \file Zoltan2_Metric.hpp
 *  \brief Metric class and namespace methods to compute quality metrics.
 *  \todo Add graph and hypergraph metrics.
 */

#ifndef ZOLTAN2_METRIC_HPP
#define ZOLTAN2_METRIC_HPP

#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Zoltan2_GraphModel.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <set>

namespace Zoltan2{

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
    values_ = arcp(tmp, 0, evalNumMetrics, true);
  }
  ArrayRCP<scalar_t> values_;
  std::string metricName_;
  multiCriteriaNorm mcnorm_;   // store "actualNorm + 1"
  static std::set<std::string> metrics_;

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
  evalNumMetrics    /*!< the number of metric values_ */
};

/*! \brief Print a standard header 
     \todo This method is static and the header is the same for
  every type of MetricValues object.  But when we add metrics for
  graph and hypergraph, we may need to look at which type of metric
  are set in order to decide what header to print.  So it would
  be a member method instead of static.
 */
static void printHeader(std::ostream &os);

/*! \brief Print a standard line of data that fits under the header. */
void printLine(std::ostream &os) const;

/*! \brief Constructor */
MetricValues(std::string mname) : 
  values_(), metricName_(mname), mcnorm_(multiCriteriaNorm(0)) { 
    resetValues();}

/*! \brief Constructor */
MetricValues() : 
  values_(), metricName_("unset"), mcnorm_(multiCriteriaNorm(0)) { 
    resetValues();}

/*! \brief Set or reset the name.  */
void setName(std::string name) { metricName_ = name;}

/*! \brief Set or reset the norm.  */
void setNorm(multiCriteriaNorm normVal) { 
  mcnorm_ = multiCriteriaNorm(normVal+1);}

/*! \brief Set the sum on the local process.  */
void setLocalSum(scalar_t x) { values_[evalLocalSum] = x;}

/*! \brief Set the global sum.  */
void setGlobalSum(scalar_t x) { values_[evalGlobalSum] = x;}

/*! \brief Set the global minimum across parts.  */
void setGlobalMin(scalar_t x) { values_[evalGlobalMin] = x;}

/*! \brief Set the global maximum across parts.  */
void setGlobalMax(scalar_t x) { values_[evalGlobalMax] = x;}

/*! \brief Set the global average (sum / numParts).  */
void setGlobalAvg(scalar_t x) { values_[evalGlobalAvg] = x;}

/*! \brief Set the imbalance of the least imbalanced part. */
void setMinImbalance(scalar_t x) { values_[evalMinImbalance] = x;}

/*! \brief Set the imbalance of the worst imbalanced part. 
     This is what we normally call the imbalance of a partition.
*/
void setMaxImbalance(scalar_t x) { values_[evalMaxImbalance] = x;}

/*! \brief Set the average imbalance of all parts. */
void setAvgImbalance(scalar_t x) { values_[evalAvgImbalance] = x;}

/*! \brief Get the name of the item measured. */
const std::string &getName() const { return metricName_; }

/*! \brief Get the norm.  */
multiCriteriaNorm getNorm() { return multiCriteriaNorm(mcnorm_-1);}

/*! \brief Get the sum on the local process. */
scalar_t getLocalSum() const { return values_[evalLocalSum];}

/*! \brief Get the global sum for all parts. */
scalar_t getGlobalSum() const { return values_[evalGlobalSum];}

/*! \brief Get the global minimum across all parts. */
scalar_t getGlobalMin() const { return values_[evalGlobalMin];}

/*! \brief Get the global maximum across all parts. */
scalar_t getGlobalMax() const { return values_[evalGlobalMax];}

/*! \brief Get the average of the sum over all parts. */
scalar_t getGlobalAvg() const { return values_[evalGlobalAvg];}

/*! \brief Get the imbalance of the least imbalanced part. */
scalar_t getMinImbalance() const { return values_[evalMinImbalance];}

/*! \brief Get the imbalance of the most imbalanced part. 
     This is what we normally call the imbalance of a partition.
*/
scalar_t getMaxImbalance() const { return values_[evalMaxImbalance];}

/*! \brief Get the average of the part imbalances. */
scalar_t getAvgImbalance() const { return values_[evalAvgImbalance];}

/// \brief Return a metric value specified by name
///
/// @param metric_name Name of metric to return
/// @param[out] value metric value returned by reference
///
/// @return Returns a boolean indicated whether or not the metric was returned 
scalar_t getMetricValue(const std::string & metric_name) const {
  if (metric_name == "local sum") {
    return this->getLocalSum();
  } else if (metric_name == "global sum") {
    return this->getGlobalSum();
  } else if (metric_name == "global maximum") {
    return this->getGlobalMax();
  } else if (metric_name == "global minimum") {
    return this->getGlobalMin();
  } else if (metric_name == "global average") {
    return this->getGlobalAvg();
  } else if (metric_name == "minimum imbalance") {
    return this->getMinImbalance();
  } else if (metric_name == "maximum imbalance") {
    return this->getMaxImbalance();
  } else if (metric_name == "average imbalance") {
    return this->getAvgImbalance();
  } else {
    return 0.0; // throw error
  }
}

bool hasMetricValue(const std::string & metric_name) const {
  return MetricValues<scalar_t>::metrics_.find(metric_name) !=
         MetricValues<scalar_t>::metrics_.end();
}
};  // end class

template <typename scalar_t>
std::set<std::string> MetricValues<scalar_t>::metrics_ = {
  "local sum",
  "global sum",
  "global maximum",
  "global minimum",
  "global average",
  "minimum imbalance",
  "maximum imbalance",
  "average imbalance",
};

/*! \brief A class containing the metrics for one measurable item.
 */

template <typename scalar_t>
  class GraphMetricValues{

private:
  void resetValues(){
    scalar_t *tmp = new scalar_t [evalNumMetrics];
    memset(tmp, 0, sizeof(scalar_t) * evalNumMetrics);
    values_ = arcp(tmp, 0, evalNumMetrics, true);
  }
  ArrayRCP<scalar_t> values_;
  std::string metricName_;
  static std::set<std::string> metrics_;

public:

/*! \brief  Enumerator for offsets into metric data.
 */
enum metricOffset{
  evalGlobalSum,   /*!< the global total on all parts */
  evalGlobalMax,   /*!< the maximum across all parts */
  evalNumMetrics    /*!< the number of metric values_ */
};

/*! \brief Print a standard header
 */
static void printHeader(std::ostream &os);

/*! \brief Print a standard line of data that fits under the header. */
void printLine(std::ostream &os) const;

/*! \brief Constructor */
GraphMetricValues(std::string mname) :
  values_(), metricName_(mname) {
  resetValues();}

/*! \brief Constructor */
GraphMetricValues() : 
  values_(), metricName_("unset") { 
    resetValues();}

/*! \brief Set or reset the name.  */
void setName(std::string name) { metricName_ = name;}

/*! \brief Set the global sum.  */
void setGlobalSum(scalar_t x) { values_[evalGlobalSum] = x;}

/*! \brief Set the global maximum across parts.  */
void setGlobalMax(scalar_t x) { values_[evalGlobalMax] = x;}

/*! \brief Get the name of the item measured. */
const std::string &getName() const { return metricName_; }

/*! \brief Get the global sum of edge cuts for all parts. */
scalar_t getGlobalSum() const { return values_[evalGlobalSum];}

/*! \brief Get the global maximum of edge cuts per part across all parts. */
scalar_t getGlobalMax() const { return values_[evalGlobalMax];}

/// \brief Return a metric value specified by name
///
/// @param metric_name Name of metric to return
/// @param[out] value metric value returned by reference
///
/// @return Returns a boolean indicated whether or not the metric was returned 
scalar_t getMetricValue(const std::string & metric_name) const {
  if (metric_name == "global maximum") {
    return this->getGlobalMax();
  } else if (metric_name == "global sum") {
    return this->getGlobalSum();
  } else {
    return 0.0; // throw error
  }
}

bool hasMetricValue(const std::string & metric_name) const {
  return GraphMetricValues<scalar_t>::metrics_.find(metric_name) !=
         GraphMetricValues<scalar_t>::metrics_.end();
}
};  // end class

template <typename scalar_t>
std::set<std::string> GraphMetricValues<scalar_t>::metrics_ = {
  "global sum",
  "global maximum",
};

template <typename scalar_t>
  void MetricValues<scalar_t>::printLine(std::ostream &os) const
{
  std::string label(metricName_);
  if (mcnorm_ > 0){
    multiCriteriaNorm realNorm = multiCriteriaNorm(mcnorm_ - 1);
    std::ostringstream oss;
    switch (realNorm){
      case normMinimizeTotalWeight:   // 1-norm = Manhattan norm 
        oss << metricName_ << " (1)";
        break;
      case normBalanceTotalMaximum:   // 2-norm = sqrt of sum of squares
        oss << metricName_ << " (2)";
        break;
      case normMinimizeMaximumWeight: // inf-norm = maximum norm 
        oss << metricName_ << " (inf)";
        break;
      default:
        oss << metricName_ << " (?)";
        break;
    }

    label = oss.str();
  }

  os << std::setw(20) << label;
  os << std::setw(11) << std::setprecision(4) << values_[evalGlobalMin];
  os << std::setw(11) << std::setprecision(4) << values_[evalGlobalMax];
  os << std::setw(11) << std::setprecision(4) << values_[evalGlobalAvg];
  os << std::setw(2) << " ";
  os << std::setw(10) << std::setprecision(4) << values_[evalMinImbalance];
  os << std::setw(10) << std::setprecision(4) << values_[evalMaxImbalance];
  os << std::setw(2) << char(177);
  os << std::setprecision(3) << values_[evalAvgImbalance];
  os << std::endl;
}

template <typename scalar_t>
  void MetricValues<scalar_t>::printHeader(std::ostream &os)
{
  os << std::setw(20) << " ";
  os << std::setw(36) << "------------SUM PER PART-----------";
  os << std::setw(2) << " ";
  os << std::setw(24) << "---IMBALANCE PER PART---";
  os << std::endl;

  os << std::setw(20) << " ";
  os << std::setw(11) << "min" << std::setw(11) << "max" << std::setw(11) << "avg";
  os << std::setw(2) << " ";
  os << std::setw(10) << "lightest" << std::setw(10) << "heaviest" << std::setw(6) << "avg";
  os << std::endl;
}

template <typename scalar_t>
  void GraphMetricValues<scalar_t>::printLine(std::ostream &os) const
{
  std::string label(metricName_);

  os << std::setw(20) << label;
  os << std::setw(12) << std::setprecision(4) << values_[evalGlobalSum];
  os << std::setw(12) << std::setprecision(4) << values_[evalGlobalMax];
  os << std::endl;
}

template <typename scalar_t>
  void GraphMetricValues<scalar_t>::printHeader(std::ostream &os)
{
  os << std::setw(20) << " ";
  os << std::setw(24) << "----------SUM----------";
  os << std::endl;

  os << std::setw(20) << " ";
  os << std::setw(12) << "sum" << std::setw(12) << "max";
  os << std::endl;
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
 * globalSumsByPart() must be called by all processes in \c comm.
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
    ArrayRCP<MetricValues<scalar_t> > &metrics,
    ArrayRCP<scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalSumsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  numParts = numNonemptyParts = 0;

  int numMetrics = 1;                       // "object count"
  if (vwgtDim) numMetrics++;                // "normed weight" or "weight 0"
  if (vwgtDim > 1) numMetrics += vwgtDim;   // "weight n"

  typedef MetricValues<scalar_t> mv_t;
  mv_t *newMetrics = new mv_t [numMetrics];
  env->localMemoryAssertion(__FILE__, __LINE__, numMetrics, newMetrics); 
  ArrayRCP<mv_t> metricArray(newMetrics, 0, numMetrics, true);

  metrics = metricArray;

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

  int next = 0;

  metrics[next].setName("object count");
  metrics[next].setLocalSum(localNumObj);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = localBuf + nparts; // single normed weight or weight 0
    scalar_t total = 0.0;
  
    for (int p=0; p < nparts; p++){
      total += wgt[p];
    }

    if (vwgtDim == 1)
      metrics[next].setName("weight 0");
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

        std::ostringstream oss;
        oss << "weight " << vdim;

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
  scalar_t min=0, max=0, sum=0;
  next = 0;

  ArrayView<scalar_t> objVec(obj, nparts);
  getStridedStats<scalar_t>(objVec, 1, 0, min, max, sum);

  metrics[next].setGlobalMin(min);
  metrics[next].setGlobalMax(max);
  metrics[next].setGlobalSum(sum);
  metrics[next].setGlobalAvg(sum / nparts);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + nparts;        // single normed weight or weight 0
  
    ArrayView<scalar_t> normedWVec(wgt, nparts);
    getStridedStats<scalar_t>(normedWVec, 1, 0, min, max, sum);

    if (vwgtDim > 1)
      metrics[next].setNorm(multiCriteriaNorm(mcNorm));

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

  /*for (part_t p=nparts-1; p > 0; p--){
    if (obj[p] > 0) break;
    numParts--;
    }*/

  numNonemptyParts = numParts; 

  for (part_t p=0; p < numParts; p++)
    if (obj[p] == 0) numNonemptyParts--;

  env->debug(DETAILED_STATUS, "Exiting globalSumsByPart");
}

/*! \brief Given the local partitioning, compute the global weighted cuts in each part.
 *
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param graph Graph model
 *   \param part   \c part[i] is the part ID for local object \c i
 *   \param numParts  on return this is the global number of parts.
 *   \param metrics on return points to a list of named GraphMetricValues cuts 
 *     that each contains the global max and sum over parts of 
 *     the item being measured. The list may contain "cut count", or
 *     "weight 0", "weight 1" and so on in that order.
 *     If uniform weights were given, then only "cut count" appears.
 *     If one set of non-uniform weights were given, then
 *     "weight 0" appear.  Finally, if multiple
 *     weights were given, we have
 *     the individual weights "weight 0", "weight 1", and so on.
 *   \param globalSums If weights are uniform, the globalSums is the
 *      \c numParts totals of global number of cuts in each part.
 *     Suppose the number of weights is \c W.  If
 *     W is 1, then on return this is an array of length \c numParts .
 *     The \c numParts entries are the total weight in each part.
 *     If \c W is greater than one, then the length of this array is 
 *     \c W*numParts .
 *     The entries are the sum of the individual weights in each part,
 *     by weight index by part number.  The array is allocated here.
 *
 * globalWeightedCutsByPart() must be called by all processes in \c comm.
 */

template <typename Adapter>
  void globalWeightedCutsByPart( 
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm, 
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graph,
    const ArrayView<const typename Adapter::part_t> &part, 
    typename Adapter::part_t &numParts, 
    ArrayRCP<GraphMetricValues<typename Adapter::scalar_t> > &metrics,
    ArrayRCP<typename Adapter::scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalWeightedCutsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  numParts = 0;

  int ewgtDim = graph->getNumWeightsPerEdge();

  int numMetrics = 1;                   // "cut count"
  if (ewgtDim) numMetrics += ewgtDim;   // "weight n"

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::node_t node_t;
  typedef typename Adapter::part_t part_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  typedef GraphMetricValues<scalar_t> mv_t;
  typedef Tpetra::CrsMatrix<part_t,lno_t,gno_t,node_t>  sparse_matrix_type;
  typedef Tpetra::Vector<part_t,lno_t,gno_t,node_t>     vector_t;
  typedef Tpetra::Map<lno_t, gno_t, node_t>                map_type;
  typedef Tpetra::global_size_t GST;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  using Teuchos::as;

  mv_t *newMetrics = new mv_t [numMetrics];
  env->localMemoryAssertion(__FILE__,__LINE__,numMetrics,newMetrics); 
  ArrayRCP<mv_t> metricArray(newMetrics, 0, numMetrics, true);

  metrics = metricArray;

  //////////////////////////////////////////////////////////
  // Figure out the global number of parts in use.
  // Verify number of vertex weights is the same everywhere.

  lno_t localNumObj = part.size();
  part_t localNum[2], globalNum[2];
  localNum[0] = static_cast<part_t>(ewgtDim);  
  localNum[1] = 0;

  for (lno_t i=0; i < localNumObj; i++)
    if (part[i] > localNum[1]) localNum[1] = part[i];

  try{
    reduceAll<int, part_t>(*comm, Teuchos::REDUCE_MAX, 2, 
      localNum, globalNum);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->globalBugAssertion(__FILE__,__LINE__,
    "inconsistent number of edge weights",
    globalNum[0] == localNum[0], DEBUG_MODE_ASSERTION, comm);

  part_t nparts = globalNum[1] + 1;

  part_t globalSumSize = nparts * numMetrics;
  scalar_t * sumBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, sumBuf);
  globalSums = arcp(sumBuf, 0, globalSumSize);

  //////////////////////////////////////////////////////////
  // Calculate the local totals by part.

  scalar_t *localBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__,__LINE__,globalSumSize,localBuf);
  memset(localBuf, 0, sizeof(scalar_t) * globalSumSize);

  scalar_t *cut = localBuf;              // # of cuts

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> vwgts;
  //size_t nv =
  graph->getVertexList(Ids, vwgts);

  ArrayView<const gno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<input_t> wgts;
  //size_t numLocalEdges =
  graph->getEdgeList(edgeIds, offsets, wgts);
  // **************************************************************************
  // *************************** BUILD MAP FOR ADJS ***************************
  // **************************************************************************

  RCP<const map_type> vertexMapG;

  // Build a list of the global vertex ids...
  gno_t min = std::numeric_limits<gno_t>::max();
  size_t maxcols = 0;
  for (lno_t i = 0; i < localNumObj; ++i) {
    if (Ids[i] < min) min = Ids[i];
    size_t ncols = offsets[i+1] - offsets[i];
    if (ncols > maxcols) maxcols = ncols;
  }

  gno_t gmin;
  Teuchos::reduceAll<int, gno_t>(*comm,Teuchos::REDUCE_MIN,1,&min,&gmin);

  //Generate Map for vertex
  vertexMapG = rcp(new map_type(INVALID, Ids, gmin, comm));

  // **************************************************************************
  // ************************** BUILD GRAPH FOR ADJS **************************
  // **************************************************************************

  RCP<sparse_matrix_type> adjsMatrix;

  // Construct Tpetra::CrsGraph objects.
  adjsMatrix = rcp (new sparse_matrix_type (vertexMapG, 0));

  Array<part_t> justOneA(maxcols, 1);

  for (lno_t localElement=0; localElement<localNumObj; ++localElement){
    // Insert all columns for global row Ids[localElement] 
    size_t ncols = offsets[localElement+1] - offsets[localElement];
    adjsMatrix->insertGlobalValues(Ids[localElement],
                                   edgeIds(offsets[localElement], ncols),
                                   justOneA(0, ncols));
  }

  //Fill-complete adjs Graph
  adjsMatrix->fillComplete ();

  // Compute part
  RCP<vector_t> scaleVec = Teuchos::rcp( new vector_t(vertexMapG,false) );
  for (lno_t localElement=0; localElement<localNumObj; ++localElement) {
    scaleVec->replaceLocalValue(localElement,part[localElement]);
  }

  // Postmultiply adjsMatrix by part
  adjsMatrix->rightScale(*scaleVec);
  Array<gno_t> Indices;
  Array<part_t> Values;

  for (lno_t i=0; i < localNumObj; i++) {
    const gno_t globalRow = Ids[i];
    size_t NumEntries = adjsMatrix->getNumEntriesInGlobalRow (globalRow);
    Indices.resize (NumEntries);
    Values.resize (NumEntries);
    adjsMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

    for (size_t j=0; j < NumEntries; j++)
      if (part[i] != Values[j])
	cut[part[i]]++;
  }

  if (numMetrics > 1) {

    scalar_t *wgt = localBuf + nparts; // weight 0

    // This code assumes the solution has the part ordered the
    // same way as the user input.  (Bug 5891 is resolved.)
    for (int edim = 0; edim < ewgtDim; edim++){
      for (lno_t i=0; i < localNumObj; i++) {
	const gno_t globalRow = Ids[i];
	size_t NumEntries = adjsMatrix->getNumEntriesInGlobalRow (globalRow);
	Indices.resize (NumEntries);
	Values.resize (NumEntries);
	adjsMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

	for (size_t j=0; j < NumEntries; j++)
	  if (part[i] != Values[j])
	    wgt[part[i]] += wgts[edim][offsets[i] + j];
      }
      wgt += nparts;         // individual weights
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
  // Global max and sum over all parts

  cut = sumBuf;                     // # of cuts
  scalar_t max=0, sum=0;
  int next = 0;

  ArrayView<scalar_t> cutVec(cut, nparts);
  getStridedStats<scalar_t>(cutVec, 1, 0, max, sum);

  metrics[next].setName("cut count");
  metrics[next].setGlobalMax(max);
  metrics[next].setGlobalSum(sum);
  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + nparts;        // weight 0
  
    for (int edim=0; edim < ewgtDim; edim++){
      ArrayView<scalar_t> fromVec(wgt, nparts);
      getStridedStats<scalar_t>(fromVec, 1, 0, max, sum);

      std::ostringstream oss;
      oss << "weight " << edim;

      metrics[next].setName(oss.str());
      metrics[next].setGlobalMax(max);
      metrics[next].setGlobalSum(sum);
      next++;
      wgt += nparts;       // individual weights
    }
  }

  numParts = nparts;

  env->debug(DETAILED_STATUS, "Exiting globalWeightedCutsByPart");
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

template <typename scalar_t, typename part_t>
  void computeImbalances(part_t numParts, part_t targetNumParts,
    const scalar_t *psizes, scalar_t sumVals , const scalar_t *vals, 
    scalar_t &min, scalar_t &max, scalar_t &avg)
{
  min = sumVals;
  max = avg = 0;

  if (sumVals <= 0 || targetNumParts < 1 || numParts < 1)
    return;

  if (targetNumParts==1 || numParts==1){
    min = max = avg = 0;  // 0 imbalance
    return;
  }

  if (!psizes){
    scalar_t target = sumVals / targetNumParts;
    for (part_t p=0; p < numParts; p++){
      scalar_t diff = vals[p] - target;
      scalar_t adiff = (diff >= 0 ? diff : -diff);
      scalar_t tmp = diff / target;
      scalar_t atmp = adiff / target;
      avg += atmp;
      if (tmp > max) max = tmp;
      if (tmp < min) min = tmp;
    }
    part_t emptyParts = targetNumParts - numParts;  
    if (emptyParts > 0){
      if (max < 1.0)
        max = 1.0;       // target divided by target
      avg += emptyParts;
    }
  }
  else{
    for (part_t p=0; p < targetNumParts; p++){
      if (psizes[p] > 0){
        if (p < numParts){
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
 *  \param numParts the number of parts supplied, which is the
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

template <typename scalar_t, typename part_t>
 void computeImbalances(part_t numParts, part_t targetNumParts,
   int numSizes, ArrayView<ArrayRCP<scalar_t> > psizes,
   scalar_t sumVals , const scalar_t *vals, 
   scalar_t &min, scalar_t &max, scalar_t &avg)
{
  min = sumVals;
  max = avg = 0;

  if (sumVals <= 0 || targetNumParts < 1 || numParts < 1)
    return;

  if (targetNumParts==1 && numParts==1){
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
    computeImbalances<scalar_t, part_t>(numParts, targetNumParts, NULL,
      sumVals, vals, min, max, avg);
    return;
  }

  double uniformSize = 1.0 / targetNumParts;
  std::vector<double> sizeVec(numSizes);
  for (int i=0; i < numSizes; i++){
    sizeVec[i] = uniformSize;
  }

  for (part_t p=0; p < numParts; p++){

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

  for (part_t p=numParts; p < targetNumParts; p++){
    bool nonEmptyPart = false;
    for (int i=0; !nonEmptyPart && i < numSizes; i++)
      if (psizes[i].size() > 0 && psizes[i][p] > 0.0)
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
    ArrayRCP<MetricValues<typename Adapter::scalar_t> > &metrics)
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
      metrics[0].getGlobalSum(), objCount, 
      min, max, avg);

  metrics[0].setMinImbalance(1.0 + min);
  metrics[0].setMaxImbalance(1.0 + max);
  metrics[0].setAvgImbalance(avg);

  ///////////////////////////////////////////////////////////////////////////
  // Compute imbalances for the normed weight sum.

  scalar_t *wgts = globalSums.getRawPtr() + numParts;

  if (metrics.size() > 1){
  
    computeImbalances<scalar_t, part_t>(numParts, targetNumParts, 
      numCriteria, partSizes.view(0, numCriteria),
      metrics[1].getGlobalSum(), wgts,
      min, max, avg);

    metrics[1].setMinImbalance(1.0 + min);
    metrics[1].setMaxImbalance(1.0 + max);
    metrics[1].setAvgImbalance(avg);

    if (metrics.size() > 2){

    ///////////////////////////////////////////////////////////////////////////
    // Compute imbalances for each individual weight.

      int next = 2;
  
      for (int vdim=0; vdim < numCriteria; vdim++){
        wgts += numParts;
        psizes = NULL;
  
        if (partSizes[vdim].size() > 0)
           psizes = partSizes[vdim].getRawPtr();
         
        computeImbalances<scalar_t, part_t>(numParts, targetNumParts, psizes,
          metrics[next].getGlobalSum(), wgts, min, max, avg);
  
        metrics[next].setMinImbalance(1.0 + min);
        metrics[next].setMaxImbalance(1.0 + max);
        metrics[next].setAvgImbalance(avg);
        next++;
      }
    }

  }
  env->debug(DETAILED_STATUS, "Exiting objectMetrics");
}

/*! \brief Print out a header and the values for a list of metrics.
 */

template <typename scalar_t, typename part_t>
  void printMetrics( std::ostream &os,
    part_t targetNumParts, part_t numParts, part_t numNonemptyParts, 
    const ArrayView<MetricValues<scalar_t> > &infoList)
{
  os << "NUMBER OF PARTS IS " << numParts;
  if (numNonemptyParts < numParts)
    os << " (" << numNonemptyParts << " of which are non-empty)";
  os << std::endl;
  if (targetNumParts != numParts)
    os << "TARGET NUMBER OF PARTS IS " << targetNumParts << std::endl;

  std::string unset("unset");

  MetricValues<scalar_t>::printHeader(os);

  for (int i=0; i < infoList.size(); i++)
    if (infoList[i].getName() != unset)
      infoList[i].printLine(os);

  os << std::endl;
}

/*! \brief Print out a header and the values for a list of graph metrics.
 */

template <typename scalar_t, typename part_t>
  void printMetrics( std::ostream &os,
    part_t targetNumParts, part_t numParts, 
    const ArrayView<GraphMetricValues<scalar_t> > &infoList)
{
  os << "NUMBER OF PARTS IS " << numParts;
  os << std::endl;
  if (targetNumParts != numParts)
    os << "TARGET NUMBER OF PARTS IS " << targetNumParts << std::endl;

  std::string unset("unset");

  GraphMetricValues<scalar_t>::printHeader(os);

  for (int i=0; i < infoList.size(); i++)
    if (infoList[i].getName() != unset)
      infoList[i].printLine(os);

  os << std::endl;
}

/*! \brief Print out a header and the values for a single metric.
 */

template <typename scalar_t, typename part_t>
  void printMetrics( std::ostream &os,
    part_t targetNumParts, part_t numParts, part_t numNonemptyParts, 
    const MetricValues<scalar_t> &info)
{
  ArrayView<MetricValues<scalar_t> > infoList(&info, 1);
  printMetrics( os, targetNumParts, numParts, numNonemptyParts, infoList);
}

/*! \brief Print out a header and the values for a single metric.
 */

template <typename scalar_t, typename part_t>
  void printMetrics( std::ostream &os,
    part_t targetNumParts, part_t numParts, 
    const GraphMetricValues<scalar_t> &info)
{
  ArrayView<GraphMetricValues<scalar_t> > infoList(&info, 1);
  printMetrics( os, targetNumParts, numParts, infoList);
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
      Environment env;  // default environment
      env.localBugAssertion(__FILE__, __LINE__, "invalid norm", false,
        BASIC_ASSERTION);
      break;
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
