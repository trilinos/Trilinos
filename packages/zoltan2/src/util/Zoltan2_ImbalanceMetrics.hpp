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

/*! \file Zoltan2_ImbalanceMetrics.hpp
 */

#ifndef ZOLTAN2_IMBALANCEMETRICS_HPP
#define ZOLTAN2_IMBALANCEMETRICS_HPP

#include <Zoltan2_BaseClassMetrics.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_MetricUtility.hpp>

#define IMBALANCE_METRICS_TYPE_NAME "ImbalanceMetrics"

namespace Zoltan2{

/*! \ImbalanceMetrics class containing the metrics for one measurable item. */
template <typename scalar_t>
  class ImbalanceMetrics : public BaseClassMetrics<scalar_t> {

private:
  multiCriteriaNorm mcnorm_;   // store "actualNorm + 1"

public:
/*! \brief Constructor */
ImbalanceMetrics(std::string mname) : BaseClassMetrics<scalar_t>(static_metricNames_.size(), mname),mcnorm_(multiCriteriaNorm(0)) {}

/*! \brief Constructor */
ImbalanceMetrics() : BaseClassMetrics<scalar_t>(static_metricNames_.size()), mcnorm_(multiCriteriaNorm(0)) {}

/*! \brief Get the class type of the metric. */
virtual const std::string & getMetricType() const { return static_metricTypeName_; }

/*! \brief Print a standard header */
static void printHeader(std::ostream &os);

/*! \brief Print a standard line of data that fits under the header. */
virtual void printLine(std::ostream &os) const;

/*! \brief Set or reset the norm.  */
void setNorm(multiCriteriaNorm normVal) { mcnorm_ = multiCriteriaNorm(normVal+1);}

/*! \brief Get the norm.  */
multiCriteriaNorm getNorm() { return multiCriteriaNorm(mcnorm_-1);}

/*! \brief Set the sum on the local process. */
void setLocalSum(scalar_t x) { this->setMetricValue("local sum", x);}

/*! \brief Set the global sum.  */
void setGlobalSum(scalar_t x) { this->setMetricValue("global sum", x );}

/*! \brief Set the global minimum across parts.  */
void setGlobalMin(scalar_t x) { this->setMetricValue("global minimum", x );}

/*! \brief Set the global maximum across parts.  */
void setGlobalMax(scalar_t x) { this->setMetricValue("global maximum", x );}

/*! \brief Set the imbalance of the worst imbalanced part. This is what we normally call the imbalance of a partition. */
void setMaxImbalance(scalar_t x) { this->setMetricValue("maximum imbalance", x);}

/*! \brief Set the average imbalance of all parts. */
void setAvgImbalance(scalar_t x) { this->setMetricValue("average imbalance", x);}

/*! \brief Get the sum on the local process. */
scalar_t getLocalSum() const { return this->getMetricValue("local sum");}

/*! \brief Get the global sum for all parts. */
scalar_t getGlobalSum() const { return this->getMetricValue("global sum");}

/*! \brief Get the global minimum across all parts. */
scalar_t getGlobalMin() const { return this->getMetricValue("global minimum");}

/*! \brief Get the global maximum across all parts. */
scalar_t getGlobalMax() const { return this->getMetricValue("global maximum");}

/*! \brief Get the imbalance of the most imbalanced part.
     This is what we normally call the imbalance of a partition.
*/
scalar_t getMaxImbalance() const { return this->getMetricValue("maximum imbalance");}

/*! \brief Get the average of the part imbalances. */
scalar_t getAvgImbalance() const { return this->getMetricValue("average imbalance");}

/*! \this method is enforced by the base class. */
virtual const std::vector<std::string> & getMetrics() const { return ImbalanceMetrics<scalar_t>::static_metricNames_; }

/*! \setup a static string name indicating my class name. */
static std::string static_metricTypeName_;

/*! \setup a static vector of strings. */
static std::vector<std::string> static_metricNames_;
};  // end class

/*! \static class name for string - used to identify by parameter lists. */
template <typename scalar_t> std::string ImbalanceMetrics<scalar_t>::static_metricTypeName_ = IMBALANCE_METRICS_TYPE_NAME;

/*! \lists all metrics for this class. */
template <typename scalar_t>
std::vector<std::string> ImbalanceMetrics<scalar_t>::static_metricNames_ = {
  "local sum",
  "global sum",
  "global minimum",
  "global maximum",
  "global average",
  "average imbalance",
  "maximum imbalance",
};

template <typename scalar_t>
  void ImbalanceMetrics<scalar_t>::printHeader(std::ostream &os)
{
  os << std::setw(20) << " ";
  os << std::setw(11) << "min" << std::setw(11) << "max" << std::setw(11) << "avg";
  os << std::setw(2) << " ";
  os << std::setw(10) << "imbalance";
  os << std::endl;
}

template <typename scalar_t>
  void ImbalanceMetrics<scalar_t>::printLine(std::ostream &os) const
{
  std::string label( this->getName() );
  if (mcnorm_ > 0){
    multiCriteriaNorm realNorm = multiCriteriaNorm(mcnorm_ - 1);
    std::ostringstream oss;
    switch (realNorm) {
      case normMinimizeTotalWeight:   // 1-norm = Manhattan norm
        oss << this->getName() << " (1)";
        break;
      case normBalanceTotalMaximum:   // 2-norm = sqrt of sum of squares
        oss << this->getName() << " (2)";
        break;
      case normMinimizeMaximumWeight: // inf-norm = maximum norm
        oss << this->getName() << " (inf)";
        break;
      default:
        oss << this->getName() << " (?)";
        break;
    }

    label = oss.str();
  }

  os << std::setw(20) << label;
  os << std::setw(11) << std::setprecision(4) 
     << this->getMetricValue("global minimum");
  os << std::setw(11) << std::setprecision(4) 
     << this->getMetricValue("global maximum");
  os << std::setw(11) << std::setprecision(4) 
     << this->getMetricValue("global average");

  os << std::setw(2) << " ";
  os << std::setw(10) << std::setprecision(4) 
     << this->getMetricValue("maximum imbalance");

  os << std::endl;
}
} // namespace Zoltan2
#endif
