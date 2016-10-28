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

/*! \file Zoltan2_OrderingMetrics.hpp
 */

#ifndef ZOLTAN2_ORDERINGMETRICS_HPP
#define ZOLTAN2_ORDERINGMETRICS_HPP

#include <Zoltan2_BaseClassMetrics.hpp>
#include <Zoltan2_MetricUtility.hpp>

namespace Zoltan2{

/*! \ImbalanceMetrics class containing the metrics for one measurable item. */
template <typename scalar_t>
  class OrderingMetrics : public BaseClassMetrics<scalar_t> {

private:

public:
/*! \brief Constructor */
OrderingMetrics(std::string mname) : OrderingMetrics<scalar_t>
  (static_metricNames_.size(), mname) {}

/*! \brief Constructor */
OrderingMetrics() : BaseClassMetrics<scalar_t>(static_metricNames_.size()) {}

/*! \brief Get the class type of the metric. */
virtual const std::string & getMetricType()
  const { return static_metricTypeName_; }

/*! \brief Print a standard header */
static void printHeader(std::ostream &os);

/*! \brief Print a standard line of data that fits under the header. */
virtual void printLine(std::ostream &os) const;

/*! \this method is enforced by the base class. */
virtual const std::vector<std::string> & getMetrics() const
  { return OrderingMetrics<scalar_t>::static_metricNames_; }

/*! \setup a static string name indicating my class name. */
static std::string static_metricTypeName_;

/*! \setup a static vector of strings. */
static std::vector<std::string> static_metricNames_;
};  // end class

/*! \static class name for string - used to identify by parameter lists. */
template <typename scalar_t> std::string
  OrderingMetrics<scalar_t>::static_metricTypeName_ = ORDERING_METRICS_TYPE_NAME;

/*! \lists all metrics for this class. */
template <typename scalar_t>
std::vector<std::string> OrderingMetrics<scalar_t>::static_metricNames_ = {
  "unsolved",
  "solved"
};

template <typename scalar_t>
  void OrderingMetrics<scalar_t>::printHeader(std::ostream &os)
{
 // TO DO...
  os << std::setw(20) << " ";
  os << std::setw(11) << "unsolved";
  os << std::setw(11) << "solved";
  os << std::endl;
}

template <typename scalar_t>
  void OrderingMetrics<scalar_t>::printLine(std::ostream &os) const
{
  std::string label( this->getName() );

  os << std::setw(20) << label;
  os << std::setw(11) << std::setprecision(4)
     << this->getMetricValue("unsolved");
  os << std::setw(11) << std::setprecision(4)
     << this->getMetricValue("solved");

  os << std::endl;
}
} // namespace Zoltan2
#endif
