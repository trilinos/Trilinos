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

/*! \file Zoltan2_MetricValues.hpp
 */

#ifndef ZOLTAN2_METRICVALUES_HPP
#define ZOLTAN2_METRICVALUES_HPP

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_MetricBase.hpp>
#include <Zoltan2_MetricFunctions.hpp>

namespace Zoltan2{

/*! \MetricValues class containing the metrics for one measurable item. */
template <typename scalar_t>
  class MetricValues : public MetricBase<scalar_t> {

private:
  multiCriteriaNorm mcnorm_;   // store "actualNorm + 1"

public:

/*! \brief Constructor */
MetricValues(std::string mname) : MetricBase<scalar_t>(evalNumMetrics, mname),mcnorm_(multiCriteriaNorm(0)) {}

/*! \brief Constructor */
MetricValues() : MetricBase<scalar_t>(evalNumMetrics), mcnorm_(multiCriteriaNorm(0)) {}

/*! \brief Print a standard header */
static void printHeader(std::ostream &os);	// MDM - note we may want to make this a regular method

/*! \brief Print a standard line of data that fits under the header. */
virtual void printLine(std::ostream &os) const;

/*! \brief Set or reset the norm.  */
void setNorm(multiCriteriaNorm normVal) { mcnorm_ = multiCriteriaNorm(normVal+1);}

/*! \brief Get the norm.  */
multiCriteriaNorm getNorm() { return multiCriteriaNorm(mcnorm_-1);}

/*! \brief Set the sum on the local process. */
void setLocalSum(scalar_t x) { this->values_[evalLocalSum] = x;}

/*! \brief Set the global sum.  */
void setGlobalSum(scalar_t x) { this->values_[evalGlobalSum] = x;}

/*! \brief Set the global minimum across parts.  */
void setGlobalMin(scalar_t x) { this->values_[evalGlobalMin] = x;}

/*! \brief Set the global maximum across parts.  */
void setGlobalMax(scalar_t x) { this->values_[evalGlobalMax] = x;}

/*! \brief Set the global average (sum / numParts).  */
void setGlobalAvg(scalar_t x) { this->values_[evalGlobalAvg] = x;}

/*! \brief Set the imbalance of the least imbalanced part. */
void setMinImbalance(scalar_t x) { this->values_[evalMinImbalance] = x;}

/*! \brief Set the imbalance of the worst imbalanced part.
     This is what we normally call the imbalance of a partition.
*/
void setMaxImbalance(scalar_t x) { this->values_[evalMaxImbalance] = x;}

/*! \brief Set the average imbalance of all parts. */
void setAvgImbalance(scalar_t x) { this->values_[evalAvgImbalance] = x;}

/*! \brief Get the sum on the local process. */
scalar_t getLocalSum() const { return this->values_[evalLocalSum];}

/*! \brief Get the global sum for all parts. */
scalar_t getGlobalSum() const { return this->values_[evalGlobalSum];}

/*! \brief Get the global minimum across all parts. */
scalar_t getGlobalMin() const { return this->values_[evalGlobalMin];}

/*! \brief Get the global maximum across all parts. */
scalar_t getGlobalMax() const { return this->values_[evalGlobalMax];}

/*! \brief Get the average of the sum over all parts. */
scalar_t getGlobalAvg() const { return this->values_[evalGlobalAvg];}

/*! \brief Get the imbalance of the least imbalanced part. */
scalar_t getMinImbalance() const { return this->values_[evalMinImbalance];}

/*! \brief Get the imbalance of the most imbalanced part.
     This is what we normally call the imbalance of a partition.
*/
scalar_t getMaxImbalance() const { return this->values_[evalMaxImbalance];}

/*! \brief Get the average of the part imbalances. */
scalar_t getAvgImbalance() const { return this->values_[evalAvgImbalance];}

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

/// \brief Return a metric value specified by name
///
/// @param metric_name Name of metric to return
/// @param[out] value metric value returned by reference
///
/// @return Returns a boolean indicated whether or not the metric was returned
virtual scalar_t getMetricValue(const std::string & metric_name) const {
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

/*! \setup a static set of strings. */
static std::set<std::string> static_metrics_;

/*! \this method is enforced by the base class. */
virtual const std::set<std::string> & getMetrics() const { return MetricValues<scalar_t>::static_metrics_; }

};  // end class

// MDM - Need to figure out a better way to encapsulate this
template <typename scalar_t>
std::set<std::string> MetricValues<scalar_t>::static_metrics_ = {
  "local sum",
  "global sum",
  "global maximum",
  "global minimum",
  "global average",
  "minimum imbalance",
  "maximum imbalance",
  "average imbalance",
};

// MDM - These are methods I will move to the .cpp
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
  void MetricValues<scalar_t>::printLine(std::ostream &os) const
{
  std::string label( MetricBase<scalar_t>::getName() );
  if (mcnorm_ > 0){
    multiCriteriaNorm realNorm = multiCriteriaNorm(mcnorm_ - 1);
    std::ostringstream oss;
    switch (realNorm){
      case normMinimizeTotalWeight:   // 1-norm = Manhattan norm
        oss << MetricBase<scalar_t>::getName() << " (1)";
        break;
      case normBalanceTotalMaximum:   // 2-norm = sqrt of sum of squares
        oss << MetricBase<scalar_t>::getName() << " (2)";
        break;
      case normMinimizeMaximumWeight: // inf-norm = maximum norm
        oss << MetricBase<scalar_t>::getName() << " (inf)";
        break;
      default:
        oss << MetricBase<scalar_t>::getName() << " (?)";
        break;
    }

    label = oss.str();
  }

  os << std::setw(20) << label;
  os << std::setw(11) << std::setprecision(4) << this->values_[evalGlobalMin];
  os << std::setw(11) << std::setprecision(4) << this->values_[evalGlobalMax];
  os << std::setw(11) << std::setprecision(4) << this->values_[evalGlobalAvg];
  os << std::setw(2) << " ";
  os << std::setw(10) << std::setprecision(4) << this->values_[evalMinImbalance];
  os << std::setw(10) << std::setprecision(4) << this->values_[evalMaxImbalance];
  os << std::setw(2) << char(177);
  os << std::setprecision(3) << this->values_[evalAvgImbalance];
  os << std::endl;
}

} //namespace Zoltan2
#endif
